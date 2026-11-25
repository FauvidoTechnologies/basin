// spandense <machines.txt> <log2T> <out.txt>: span at t_j = round(2^(10 + j/64)), j = 0..64*(log2T-10).
// One machine string per line (last field). Line out: "machine halt_step span_0 ... span_J" (-1 after
// halting or after the head leaves the 2^26-cell tape).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); int L = atoi(argv[2]); FILE *o = fopen(argv[3], "w");
    long T = 1L << L, N = 1L << 26; int J = 64 * (L - 10) + 1;
    long *ck = malloc(J * sizeof(long)); for (int j = 0; j < J; j++) ck[j] = lround(pow(2.0, 10 + j / 64.0));
    unsigned char *tape = calloc(N, 1); long *sp = malloc(J * sizeof(long));
    char line[512];
    while (fgets(line, sizeof line, f)) {
        char *m = NULL, *tok = strtok(line, " ,\t\r\n"); while (tok) { m = tok; tok = strtok(NULL, " ,\t\r\n"); }
        if (!m || strlen(m) < 13) continue;
        int W[16], M[16], G[16], n = (int)(strlen(m) + 1) / 7;
        for (int s = 0; s < n; s++) for (int r = 0; r < 2; r++) { char *p = m + 7 * s + 3 * r;
            W[2*s+r] = p[0] == '1'; M[2*s+r] = p[1] == 'L' ? -1 : 1; G[2*s+r] = (p[2] >= 'A' && p[2] < 'A' + n) ? p[2] - 'A' : -1; }
        long pos = N / 2, lo = pos, hi = pos, t = 0, halt = -1; int st = 0, j = 0;
        while (j < J) {
            int c = 2 * st + tape[pos];
            if (G[c] < 0) { halt = t + 1; break; }
            tape[pos] = W[c]; pos += M[c]; st = G[c]; t++;
            if (pos < lo) lo = pos; else if (pos > hi) hi = pos;
            if (pos < 1 || pos > N - 2) break;
            while (j < J && t == ck[j]) sp[j++] = hi - lo + 1;
        }
        for (; j < J; j++) sp[j] = -1;
        fprintf(o, "%s %ld", m, halt); for (int k = 0; k < J; k++) fprintf(o, " %ld", sp[k]); fprintf(o, "\n");
        memset(tape + lo, 0, hi - lo + 1);
    }
    fclose(o); return 0;
}
