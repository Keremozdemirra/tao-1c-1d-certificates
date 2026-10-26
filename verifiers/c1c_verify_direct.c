/* Exact score of a C_1c certificate by direct summation of every lag in 128-bit integers.
 *
 * Independent of verify_kronecker.py (different algorithm and language).
 * Input: text file "K n" on the first line, then n non-negative integer weights.
 * Output: S = sum w, min_{0<=k<=K} A_k with A_k = sum_i w_i w_{i+k}, the argmin, and the
 * score K * min A / S^2 as numerator and denominator (S^2 printed via 256-bit-free long multiplication).
 *
 * Overflow guard: w_i < 2^44, n < 2^20, so A_k < 2^108.
 * Build: cc -O2 -o verify_direct verify_direct.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef unsigned __int128 u128;

static void print_u128(u128 x) {
    char buf[64];
    int k = 0;
    if (x == 0) buf[k++] = '0';
    while (x) { buf[k++] = (char)('0' + (int)(x % 10)); x /= 10; }
    while (k) putchar(buf[--k]);
}

/* print a*b for a, b < 2^64 exactly: 128-bit product */
static u128 gcd128(u128 a, u128 b) { while (b) { u128 t = a % b; a = b; b = t; } return a; }

/* decimal big multiply of two u128 values, printed exactly */
static void print_product(u128 a, u128 b) {
    unsigned long long da[8], db[8], t[17] = {0};
    int la = 0, lb = 0;
    const unsigned long long B = 1000000000ull;
    if (a == 0 || b == 0) { printf("0"); return; }
    while (a) { da[la++] = (unsigned long long)(a % B); a /= B; }
    while (b) { db[lb++] = (unsigned long long)(b % B); b /= B; }
    for (int i = 0; i < la; i++)
        for (int j = 0; j < lb; j++) {
            t[i + j] += da[i] * db[j];
            for (int k = i + j; t[k] >= B; k++) { t[k + 1] += t[k] / B; t[k] %= B; }
        }
    int len = la + lb;
    while (len > 1 && t[len - 1] == 0) len--;
    printf("%llu", t[len - 1]);
    for (int i = len - 2; i >= 0; i--) printf("%09llu", t[i]);
}

int main(int argc, char **argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s weights.txt\n", argv[0]); return 2; }
    FILE *fp = fopen(argv[1], "r");
    if (!fp) { perror("open"); return 2; }
    long K, n;
    if (fscanf(fp, "%ld %ld", &K, &n) != 2 || n <= 0 || n >= (1L << 20) || K <= 0) { fprintf(stderr, "bad header\n"); return 2; }
    unsigned long long *w = malloc(sizeof *w * n);
    u128 S = 0;
    for (long i = 0; i < n; i++) {
        long long x;
        if (fscanf(fp, "%lld", &x) != 1) { fprintf(stderr, "short file\n"); return 2; }
        if (x < 0 || x >= (1LL << 44)) { fprintf(stderr, "weight out of range\n"); return 2; }
        w[i] = (unsigned long long)x;
        S += w[i];
    }
    fclose(fp);
    if (S == 0) { fprintf(stderr, "zero mass\n"); return 2; }
    u128 mn = 0, A0 = 0;
    long arg = -1, cnt = 0;
    for (long k = 0; k <= K; k++) {
        u128 A = 0;
        for (long i = 0; i + k < n; i++) A += (u128)w[i] * w[i + k];
        if (k == 0) A0 = A;
        if (arg < 0 || A < mn) { mn = A; arg = k; cnt = 1; }
        else if (A == mn) cnt++;
    }
    (void)A0;
    printf("K = %ld, cells = %ld\nS = ", K, n); print_u128(S);
    printf("\nmin A = "); print_u128(mn);
    printf("\nargmin k = %ld (%ld lags attain it)\n", arg, cnt);
    u128 num = (u128)K * mn;
    u128 g1 = gcd128(num, S); num /= g1; u128 s1 = S / g1;
    u128 g2 = gcd128(num, S); num /= g2; u128 s2 = S / g2;
    printf("score = "); print_u128(num); printf(" / "); print_product(s1, s2); printf("\n");
    printf("score ~ %.18Lf\n", (long double)num / ((long double)s1 * (long double)s2));
    free(w);
    return 0;
}
