// enum <nstates> <S(n)> <out.bin>: TNF tree enumeration of n-state 2-symbol machines starting 1RB.
// Writes every machine that runs S(n)+1 steps without reaching an undefined transition
// (hence never halts) in bbchallenge's 30-byte format (5 states, unused ones undefined).
#include <stdio.h>
#include <string.h>
static int N, LIM; static FILE *out;
static unsigned char W[10], M[10], G[10];   // G = 0 undefined, else 1..n
static long nonhalt = 0, halting = 0, nodes = 0;
static unsigned char tape[512];
static void rec(void) {
    nodes++;
    memset(tape, 0, sizeof tape);
    int pos = 256, st = 0, maxs = 0;
    for (int k = 0; k < 2 * N; k++) if (G[k] && G[k] - 1 > maxs) maxs = G[k] - 1;
    for (int t = 0; t < LIM; t++) {
        int c = 2 * st + tape[pos];
        if (!G[c]) {
            halting++;  // the leaf where this transition halts
            int defined = 0; for (int k = 0; k < 2 * N; k++) defined += G[k] != 0;
            if (defined == 2 * N - 1) return;  // last free slot: only the halting choice keeps TNF finite
            int top = maxs + 1 < N - 1 ? maxs + 1 : N - 1;
            for (int w = 0; w < 2; w++) for (int m = 0; m < 2; m++) for (int g = 0; g <= top; g++) {
                W[c] = w; M[c] = m; G[c] = g + 1; rec(); G[c] = 0;
            }
            return;
        }
        tape[pos] = W[c]; pos += M[c] ? -1 : 1; st = G[c] - 1;
        if (pos < 1 || pos > 510) break;  // cannot happen within S(4)+1 steps
    }
    unsigned char rec30[30] = {0};
    for (int k = 0; k < 2 * N; k++) if (G[k]) { rec30[3*k] = W[k]; rec30[3*k+1] = M[k]; rec30[3*k+2] = G[k]; }
    fwrite(rec30, 1, 30, out); nonhalt++;
}
int main(int argc, char **argv) {
    sscanf(argv[1], "%d", &N); sscanf(argv[2], "%d", &LIM); LIM += 1; out = fopen(argv[3], "wb");
    W[0] = 1; M[0] = 0; G[0] = 2;  // A0 = 1RB
    rec();
    fclose(out);
    printf("n=%d: %ld non-halting machines, %ld halting leaves, %ld tree nodes\n", N, nonhalt, halting, nodes);
    return 0;
}
