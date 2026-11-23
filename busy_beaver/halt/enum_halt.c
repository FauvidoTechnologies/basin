// enum_halt <n> <LIM> <MINH> <P> <p> <out.txt>: TNF tree (first 1RB), as small/enum.c, but simulating
// each node only up to LIM steps. Every halting leaf (an undefined transition first reached at step h)
// is counted in a log2 histogram of h; leaves with h >= MINH are written as "machine h".
// Work is split across P processes by subtree at depth SPLIT (defined transitions).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define SPLIT 4
static int N, P, PI; static long LIM, MINH; static FILE *out;
static unsigned char W[16], M[16], G[16];
static long hist[64], ctr = 0, nonhalt = 0;
static unsigned char *tape; static long TL;
static void emit(long h) {
    char s[128]; int q = 0;
    for (int x = 0; x < N; x++) { for (int r = 0; r < 2; r++) { int k = 2 * x + r;
        if (!G[k]) { s[q++] = '-'; s[q++] = '-'; s[q++] = '-'; }
        else { s[q++] = '0' + W[k]; s[q++] = M[k] ? 'L' : 'R'; s[q++] = 'A' + G[k] - 1; } }
        if (x < N - 1) s[q++] = '_'; }
    s[q] = 0; fprintf(out, "%s %ld\n", s, h);
}
static void rec(void) {
    int defined = 0, maxs = 0;
    for (int k = 0; k < 2 * N; k++) if (G[k]) { defined++; if (G[k] - 1 > maxs) maxs = G[k] - 1; }
    if (defined == SPLIT && (ctr++ % P) != PI) return;
    memset(tape, 0, TL); long pos = TL / 2, lo = pos, hi = pos; int st = 0;
    for (long t = 0; t < LIM; t++) {
        int c = 2 * st + tape[pos];
        if (!G[c]) {
            if (defined >= SPLIT || PI == 0) {  // count each halting leaf once
                int b = 0; while ((1L << (b + 1)) <= t + 1) b++; hist[b]++;
                if (t + 1 >= MINH) emit(t + 1);
            }
            if (defined == 2 * N - 1) return;
            int top = maxs + 1 < N - 1 ? maxs + 1 : N - 1;
            for (int w = 0; w < 2; w++) for (int m = 0; m < 2; m++) for (int g = 0; g <= top; g++) {
                W[c] = w; M[c] = m; G[c] = g + 1; rec(); G[c] = 0;
            }
            return;
        }
        tape[pos] = W[c]; pos += M[c] ? -1 : 1; st = G[c] - 1;
    }
    nonhalt++;
}
int main(int argc, char **argv) {
    N = atoi(argv[1]); LIM = atol(argv[2]); MINH = atol(argv[3]); P = atoi(argv[4]); PI = atoi(argv[5]);
    out = fopen(argv[6], "w"); TL = 2 * LIM + 16; tape = malloc(TL);
    W[0] = 1; M[0] = 0; G[0] = 2;
    rec();
    fprintf(out, "# nonhalting_at_LIM %ld\n# hist", nonhalt);
    for (int b = 0; b < 40; b++) fprintf(out, " %ld", hist[b]);
    fprintf(out, "\n"); fclose(out); return 0;
}
