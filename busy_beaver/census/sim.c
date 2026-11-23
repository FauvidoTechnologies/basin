// sim <records.bin> <start> <end> <log2T> <out.bin>
// For each machine: int32 halt_step (-1 if none), then span[K], ones[K] at t_k = round(2^(10+k/4)).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); long a = atol(argv[2]), b = atol(argv[3]);
    int L = atoi(argv[4]); FILE *o = fopen(argv[5], "wb");
    long T = 1L << L; int K = 4 * (L - 10) + 1;
    long *ck = malloc(K * sizeof(long));
    for (int k = 0; k < K; k++) ck[k] = lround(pow(2.0, 10 + k / 4.0));
    unsigned char *tape = calloc(2 * T + 16, 1);
    int32_t *row = malloc((1 + 2 * K) * sizeof(int32_t));
    unsigned char rec[30];
    fseek(f, 30 * a, SEEK_SET);
    for (long i = a; i < b; i++) {
        if (fread(rec, 1, 30, f) != 30) break;
        int W[10], M[10], G[10];
        for (int s = 0; s < 5; s++) for (int r = 0; r < 2; r++) {
            unsigned char *p = rec + 6 * s + 3 * r;
            W[2*s+r] = p[0]; M[2*s+r] = p[1] ? -1 : 1; G[2*s+r] = (int)p[2] - 1;
        }
        long pos = T + 8, lo = pos, hi = pos, ones = 0, t = 0; int st = 0, k = 0;
        row[0] = -1;
        while (t < T) {
            int c = 2 * st + tape[pos];
            if (G[c] < 0) { row[0] = (int32_t)(t + 1); break; }
            ones += W[c] - tape[pos]; tape[pos] = W[c]; pos += M[c]; st = G[c]; t++;
            if (pos < lo) lo = pos; else if (pos > hi) hi = pos;
            if (t == ck[k]) { row[1 + k] = hi - lo + 1; row[1 + K + k] = ones; k++; }
        }
        for (; k < K; k++) { row[1 + k] = -1; row[1 + K + k] = -1; }
        fwrite(row, sizeof(int32_t), 1 + 2 * K, o);
        memset(tape + lo, 0, hi - lo + 1);
    }
    fclose(o); return 0;
}
