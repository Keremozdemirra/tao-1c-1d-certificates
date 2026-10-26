/* Independent exact cross-check helper for the cold verifiers (written from scratch).
 *
 * Reads signed integers (one per line, each |x| < 2^62) and accumulates in __int128.
 * Writes one exact decimal integer per line.
 *
 *   xcheck conv  A.txt OUT.txt            full linear convolution: out[k] = sum_{i+j=k} a_i a_j, k=0..2n-2
 *   xcheck acorr A.txt OUT.txt MAXLAG     out[k] = sum_i a_i a_{i+k}, k=0..MAXLAG
 *   xcheck xcorr A.txt B.txt OUT.txt      out[l] = sum_i a_i b_{i+k}, k = -(n-1) .. (m-1) ascending
 *
 * The caller (python) is responsible for bounding n * max|x|^2 < 2^126.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef __int128 i128;

static void put128(FILE *f, i128 v) {
    char buf[64];
    int neg = 0, p = 63;
    buf[p] = 0;
    unsigned __int128 u;
    if (v < 0) { neg = 1; u = (unsigned __int128)(-(v + 1)) + 1; } else u = (unsigned __int128)v;
    if (u == 0) { buf[--p] = '0'; }
    while (u > 0) { buf[--p] = (char)('0' + (int)(u % 10)); u /= 10; }
    if (neg) buf[--p] = '-';
    fputs(buf + p, f);
    fputc('\n', f);
}

static int64_t *readints(const char *path, long *n) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(2); }
    long cap = 1024, len = 0;
    int64_t *a = malloc(cap * sizeof(int64_t));
    long long x;
    while (fscanf(f, "%lld", &x) == 1) {
        if (len == cap) { cap *= 2; a = realloc(a, cap * sizeof(int64_t)); }
        a[len++] = (int64_t)x;
    }
    fclose(f);
    *n = len;
    return a;
}

int main(int argc, char **argv) {
    if (argc < 4) { fprintf(stderr, "usage\n"); return 2; }
    const char *mode = argv[1];
    long n, m;
    if (!strcmp(mode, "conv")) {
        int64_t *a = readints(argv[2], &n);
        FILE *o = fopen(argv[3], "w");
        for (long k = 0; k <= 2 * n - 2; k++) {
            long lo = k - (n - 1); if (lo < 0) lo = 0;
            long hi = k < n - 1 ? k : n - 1;
            i128 acc = 0;
            for (long i = lo; i <= hi; i++) acc += (i128)a[i] * (i128)a[k - i];
            put128(o, acc);
        }
        fclose(o);
    } else if (!strcmp(mode, "acorr")) {
        if (argc < 5) return 2;
        int64_t *a = readints(argv[2], &n);
        long maxlag = atol(argv[4]);
        FILE *o = fopen(argv[3], "w");
        for (long k = 0; k <= maxlag; k++) {
            i128 acc = 0;
            for (long i = 0; i + k < n; i++) acc += (i128)a[i] * (i128)a[i + k];
            put128(o, acc);
        }
        fclose(o);
    } else if (!strcmp(mode, "xcorr")) {
        if (argc < 5) return 2;
        int64_t *a = readints(argv[2], &n);
        int64_t *b = readints(argv[3], &m);
        FILE *o = fopen(argv[4], "w");
        for (long k = -(n - 1); k <= m - 1; k++) {
            long lo = k < 0 ? -k : 0;
            long hi = (m - 1 - k) < (n - 1) ? (m - 1 - k) : (n - 1);
            i128 acc = 0;
            for (long i = lo; i <= hi; i++) acc += (i128)a[i] * (i128)b[i + k];
            put128(o, acc);
        }
        fclose(o);
    } else { fprintf(stderr, "bad mode\n"); return 2; }
    return 0;
}
