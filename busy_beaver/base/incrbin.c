// incrbin <records.bin> <indices.int64> <start> <end> <log2T> <out.txt>
// For each indexed record: run 2^log2T steps, collect the steps at which the span grows; print
// "index n_incr factor fingerprint" where factor = (t_b/t_a)^(1/q) over the last q increments after
// 2^14 steps, q the largest multiple of 6 available (0 if none), and fingerprint = FNV-1a hash of all
// increment times (equal behaviour gives equal hashes).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"), *fi = fopen(argv[2], "rb"); long a = atol(argv[3]), b = atol(argv[4]);
    long T = 1L << atoi(argv[5]); FILE *o = fopen(argv[6], "w");
    static unsigned char tape[1 << 20]; static long ts[1 << 16]; unsigned char rec[30]; int64_t idx;
    fseek(fi, 8 * a, SEEK_SET);
    for (long j = a; j < b; j++) {
        if (fread(&idx, 8, 1, fi) != 1) break;
        fseek(f, 30 * idx, SEEK_SET); if (fread(rec, 1, 30, f) != 30) break;
        int W[10], M[10], G[10];
        for (int k = 0; k < 10; k++) { W[k] = rec[3*k]; M[k] = rec[3*k+1] ? -1 : 1; G[k] = (int)rec[3*k+2] - 1; }
        long pos = 1 << 19, lo = pos, hi = pos; int st = 0, n = 0;
        for (long t = 1; t <= T; t++) {
            int c = 2 * st + tape[pos]; if (G[c] < 0) break;
            tape[pos] = W[c]; pos += M[c]; st = G[c];
            if (pos < lo || pos > hi) { if (pos < lo) lo = pos; else hi = pos; if (n < (1 << 16)) ts[n++] = t; }
            if (pos < 2 || pos > (1 << 20) - 2) break;
        }
        int first = 0; while (first < n && ts[first] <= (1 << 14)) first++;
        int late = n - first, q = (late - 1) / 6 * 6; double F = 0;
        if (q >= 6) F = pow((double)ts[n - 1] / ts[n - 1 - q], 1.0 / q);
        uint64_t h = 1469598103934665603ULL; for (int k = 0; k < n; k++) { h ^= (uint64_t)ts[k]; h *= 1099511628211ULL; }
        fprintf(o, "%lld %d %.7f %016llx\n", (long long)idx, n, F, (unsigned long long)h);
        memset(tape + lo, 0, hi - lo + 1);
    }
    fclose(o); return 0;
}
