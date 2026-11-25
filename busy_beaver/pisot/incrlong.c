// incrlong <machine> <log2T>: every step at which the span grows, up to 2^log2T steps.
// Prints "machine n t_1 ... t_n". Tape is 1<<16 cells (enough for log-growth machines).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int argc, char **argv) {
    char *m = argv[1]; long T = 1L << atoi(argv[2]);
    static unsigned char tape[1 << 16]; static long ts[1 << 16];
    int W[16], M[16], G[16], n = (int)(strlen(m) + 1) / 7;
    for (int s = 0; s < n; s++) for (int r = 0; r < 2; r++) { char *p = m + 7*s + 3*r;
        W[2*s+r] = p[0] == '1'; M[2*s+r] = p[1] == 'L' ? -1 : 1; G[2*s+r] = (p[2] >= 'A' && p[2] < 'A' + n) ? p[2] - 'A' : -1; }
    long pos = 1 << 15, lo = pos, hi = pos, t; int st = 0, k = 0;
    for (t = 1; t <= T; t++) {
        int c = 2 * st + tape[pos]; if (G[c] < 0) break;
        tape[pos] = W[c]; pos += M[c]; st = G[c];
        if (pos < lo || pos > hi) { if (pos < lo) lo = pos; else hi = pos; ts[k++] = t;
            if (pos < 2 || pos > (1 << 16) - 3 || k == (1 << 16)) break; }
    }
    printf("%s %d", m, k); for (int j = 0; j < k; j++) printf(" %ld", ts[j]); printf("\n");
    return 0;
}
