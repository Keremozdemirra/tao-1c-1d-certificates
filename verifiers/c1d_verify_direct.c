/* Exact score of a C_1d certificate by direct O(n^2) convolution in 128-bit integers.
 *
 * Written independently of verify_kronecker.py (different algorithm, different language).
 * Input: a text file with n on the first line, then n integers (the cell values v_i).
 * Output: n, S = sum v, M = max_m |sum_{i+j=m} v_i v_j|, the argmax, and the exact ratio
 * 2n*M / S^2 as a reduced fraction printed in decimal.
 *
 * Overflow guard: requires max|v| < 2^44 and n < 2^20, so |c_m| < n * 2^88 < 2^108.
 *
 * Build: cc -O2 -o verify_direct verify_direct.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
typedef unsigned __int128 u128;

static void print_u128(u128 x) {
    char buf[64];
    int k = 0;
    if (x == 0) buf[k++] = '0';
    while (x) { buf[k++] = (char)('0' + (int)(x % 10)); x /= 10; }
    while (k) putchar(buf[--k]);
}

/* arbitrary-size unsigned decimal multiply/divide helpers for the final fraction */
typedef struct { unsigned int d[40]; int len; } big; /* base 1e9 limbs, little endian */

static void big_from_u128(big *b, u128 x) {
    b->len = 0;
    memset(b->d, 0, sizeof b->d);
    if (x == 0) { b->len = 1; return; }
    while (x) { b->d[b->len++] = (unsigned int)(x % 1000000000u); x /= 1000000000u; }
}

static void big_mul(big *r, const big *a, const big *b) {
    unsigned long long t[80] = {0};
    for (int i = 0; i < a->len; i++)
        for (int j = 0; j < b->len; j++) {
            t[i + j] += (unsigned long long)a->d[i] * b->d[j];
            for (int k = i + j; t[k] >= 1000000000ull; k++) { t[k + 1] += t[k] / 1000000000ull; t[k] %= 1000000000ull; }
        }
    int len = a->len + b->len;
    while (len > 1 && t[len - 1] == 0) len--;
    memset(r->d, 0, sizeof r->d);
    for (int i = 0; i < len; i++) r->d[i] = (unsigned int)t[i];
    r->len = len;
}

static void big_print(const big *b) {
    printf("%u", b->d[b->len - 1]);
    for (int i = b->len - 2; i >= 0; i--) printf("%09u", b->d[i]);
}

static u128 gcd128(u128 a, u128 b) { while (b) { u128 t = a % b; a = b; b = t; } return a; }

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s values.txt\n", argv[0]); return 2; }
    FILE *fp = fopen(argv[1], "r");
    if (!fp) { perror("open"); return 2; }
    long n;
    if (fscanf(fp, "%ld", &n) != 1 || n <= 0 || n >= (1L << 20)) { fprintf(stderr, "bad n\n"); return 2; }
    long long *v = malloc(sizeof(long long) * n);
    long long lim = 1LL << 44;
    i128 S = 0;
    for (long i = 0; i < n; i++) {
        if (fscanf(fp, "%lld", &v[i]) != 1) { fprintf(stderr, "short file\n"); return 2; }
        if (v[i] >= lim || v[i] <= -lim) { fprintf(stderr, "value too large\n"); return 2; }
        S += v[i];
    }
    fclose(fp);
    if (S <= 0) { fprintf(stderr, "need positive sum\n"); return 2; }
    u128 M = 0;
    long arg = -1, count = 0;
    i128 total = 0;
    for (long m = 0; m < 2 * n - 1; m++) {
        long lo = m - (n - 1) > 0 ? m - (n - 1) : 0;
        long hi = m < n - 1 ? m : n - 1;
        i128 c = 0;
        for (long i = lo; i <= hi; i++) c += (i128)v[i] * v[m - i];
        total += c;
        u128 a = c < 0 ? (u128)(-c) : (u128)c;
        if (a > M) { M = a; arg = m; count = 1; }
        else if (a == M) count++;
    }
    if (total != S * S) { fprintf(stderr, "check failed: sum c != S^2\n"); return 1; }
    printf("n = %ld\nS = ", n); print_u128((u128)S);
    printf("\nM = "); print_u128(M);
    printf("\nargmax m = %ld (%ld knots attain M)\n", arg, count);
    /* ratio = 2n*M / S^2: reduce by g = gcd(2n*M, S) twice (S^2 = S*S) */
    u128 num = (u128)(2 * n) * M, s = (u128)S;
    u128 g1 = gcd128(num, s); num /= g1; u128 s1 = s / g1;
    u128 g2 = gcd128(num, s); num /= g2; u128 s2 = s / g2;
    big bn, ba, bb, bd;
    big_from_u128(&bn, num); big_from_u128(&ba, s1); big_from_u128(&bb, s2);
    big_mul(&bd, &ba, &bb);
    printf("ratio = "); big_print(&bn); printf(" / "); big_print(&bd); printf("\n");
    printf("ratio ~ %.18Lf\n", (long double)num / ((long double)s1 * (long double)s2));
    free(v);
    return 0;
}
