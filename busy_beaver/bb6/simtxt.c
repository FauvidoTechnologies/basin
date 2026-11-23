// simtxt <machines.txt> <log2T> <out.bin>: one machine string per line (last whitespace-separated
// field), up to 8 states; '---' or a target outside A..H means halt. Output rows as census/sim.c:
// int32 halt_step (-1 if none), span[K], ones[K] at t_k = round(2^(10+k/4)); -1 after halting.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); int L = atoi(argv[2]); FILE *o = fopen(argv[3], "wb");
    long T = 1L << L; int K = 4 * (L - 10) + 1;
    long *ck = malloc(K * sizeof(long));
    for (int k = 0; k < K; k++) ck[k] = lround(pow(2.0, 10 + k / 4.0));
    unsigned char *tape = calloc(2 * T + 16, 1);
    int32_t *row = malloc((1 + 2 * K) * sizeof(int32_t));
    char line[512];
    while (fgets(line, sizeof line, f)) {
        char *m = NULL, *tok = strtok(line, " \t\r\n");
        while (tok) { m = tok; tok = strtok(NULL, " \t\r\n"); }
        if (!m) continue;
        int W[16], M[16], G[16], n = (int)(strlen(m) + 1) / 7;
        for (int s = 0; s < n; s++) for (int r = 0; r < 2; r++) {
            char *p = m + 7 * s + 3 * r;
            W[2*s+r] = p[0] == '1'; M[2*s+r] = p[1] == 'L' ? -1 : 1;
            G[2*s+r] = (p[2] >= 'A' && p[2] < 'A' + n) ? p[2] - 'A' : -1;
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
