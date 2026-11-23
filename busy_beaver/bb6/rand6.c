// rand6 <nstates> <count> <seed> <LIM> <out.txt>: random n-state machines with A0 = 1RB, exactly one
// halting slot ('---') and every other transition uniform over {0,1}x{L,R}x{A..}. Kept if the machine
// visits every state and runs LIM steps without halting. Duplicates are skipped.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned long long rs; static unsigned rnd(unsigned n) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs % n; }
int main(int argc, char **argv) {
    int N = atoi(argv[1]), C = atoi(argv[2]); rs = strtoull(argv[3], 0, 10) * 2654435761ULL + 1; long LIM = atol(argv[4]);
    FILE *o = fopen(argv[5], "w");
    static unsigned char tape[1 << 18];
    int got = 0; long tries = 0;
    char (*seen)[128] = malloc((size_t)C * 128);
    while (got < C) {
        tries++;
        int W[16], M[16], G[16];
        for (int k = 0; k < 2 * N; k++) { W[k] = rnd(2); M[k] = rnd(2) ? 1 : -1; G[k] = rnd(N); }
        W[0] = 1; M[0] = 1; G[0] = 1;
        G[1 + rnd(2 * N - 1)] = -1;  // the one halting slot, never A0
        memset(tape, 0, sizeof tape); long pos = 1 << 17, t = 0; int st = 0, ok = 1, vis = 1;
        while (t < LIM) {
            int c = 2 * st + tape[pos];
            if (G[c] < 0) { ok = 0; break; }
            tape[pos] = W[c]; pos += M[c]; st = G[c]; vis |= 1 << st; t++;
            if (pos < 2 || pos > (1 << 18) - 2) { ok = 0; break; }
        }
        if (!ok || vis != (1 << N) - 1) continue;
        char s[128]; int q = 0;
        for (int x = 0; x < N; x++) { for (int r = 0; r < 2; r++) { int k = 2 * x + r;
            if (G[k] < 0) { s[q++] = '-'; s[q++] = '-'; s[q++] = '-'; }
            else { s[q++] = '0' + W[k]; s[q++] = M[k] > 0 ? 'R' : 'L'; s[q++] = 'A' + G[k]; } }
            if (x < N - 1) s[q++] = '_'; }
        s[q] = 0;
        int dup = 0; for (int j = 0; j < got && !dup; j++) dup = !strcmp(seen[j], s);
        if (dup) continue;
        strcpy(seen[got], s); fprintf(o, "%s\n", s); got++;
    }
    fclose(o); fprintf(stderr, "%d kept of %ld random machines\n", got, tries); return 0;
}
