// incr <machines.txt> <log2T>: for each machine, the steps at which the span grows (a new cell is
// visited), up to 2^log2T steps. Prints "machine n t_1 ... t_n" (only the last 60 increments).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); long T = 1L << atoi(argv[2]);
    static unsigned char tape[1 << 20]; static long ts[1 << 16];
    char line[256];
    while (fgets(line, sizeof line, f)) {
        char *m = strtok(line, " \t\r\n"); if (!m) continue;
        int W[16], M[16], G[16], n = (int)(strlen(m) + 1) / 7;
        for (int s = 0; s < n; s++) for (int r = 0; r < 2; r++) { char *p = m + 7*s + 3*r;
            W[2*s+r] = p[0] == '1'; M[2*s+r] = p[1] == 'L' ? -1 : 1; G[2*s+r] = (p[2] >= 'A' && p[2] < 'A' + n) ? p[2] - 'A' : -1; }
        memset(tape, 0, sizeof tape); long pos = 1 << 19, lo = pos, hi = pos; int st = 0, k = 0;
        for (long t = 1; t <= T; t++) {
            int c = 2 * st + tape[pos]; if (G[c] < 0) break;
            tape[pos] = W[c]; pos += M[c]; st = G[c];
            if (pos < lo || pos > hi) { if (pos < lo) lo = pos; else hi = pos; if (k < (1 << 16)) ts[k++] = t; }
            if (pos < 2 || pos > (1 << 20) - 2) break;
        }
        printf("%s %d", m, k); for (int j = (k > 60 ? k - 60 : 0); j < k; j++) printf(" %ld", ts[j]); printf("\n");
    }
    return 0;
}
