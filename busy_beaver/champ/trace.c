// Run the BB(5) champion; at each predicted configuration time verify the tape is 0^inf <A 1^n 0^inf
// and print (t, n, span); also print span at 2000 log-spaced times between 2^10 and 2^26.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(void) {
    const char *m = "1RB1LC_1RC1RB_1RD0LE_1LA1LD_1RZ0LA";
    int W[10], M[10], G[10];
    for (int s = 0; s < 5; s++) for (int r = 0; r < 2; r++) { const char *p = m + 7*s + 3*r;
        W[2*s+r] = p[0] == '1'; M[2*s+r] = p[1] == 'L' ? -1 : 1; G[2*s+r] = (p[2] >= 'A' && p[2] <= 'E') ? p[2] - 'A' : -1; }
    long steps[] = {15, 73, 277, 907, 2757, 7957, 22777, 64407, 180307, 504027, 1403967, 3906393, 10861903, 30196527};
    long N = 1 << 16; unsigned char *tape = calloc(2 * N, 1); long pos = N, lo = pos, hi = pos; int st = 0;
    long cfg = 0, ci = 0; long ns[] = {6,16,34,64,114,196,334,564,946,1584,2646,4416,7366,12284};
    int K = 2000; double lt0 = 10, lt1 = 26; int k = 0; long next = 1 << 10;
    for (long t = 0; ; t++) {
        if (ci < 14 && t == cfg + steps[ci]) {
            cfg = t; long n = 0; int ok = st == 0 && tape[pos] == 0;
            long q = pos + 1; while (tape[q]) { n++; q++; }
            for (long x = lo; x <= hi; x++) if ((x < pos + 1 || x > pos + n) && tape[x]) ok = 0;
            printf("CFG t=%ld n=%ld expected=%ld form_ok=%d span=%ld\n", t, n, ns[ci], ok, hi - lo + 1);
            ci++;
        }
        while (k < K && t == next) { printf("S %ld %ld\n", t, hi - lo + 1); k++; next = llround(pow(2, lt0 + (lt1 - lt0) * k / (K - 1))); if (next <= t) next = t + 1; }
        int c = 2 * st + tape[pos]; if (G[c] < 0) { printf("HALT %ld\n", t + 1); break; }
        tape[pos] = W[c]; pos += M[c]; st = G[c];
        if (pos < lo) lo = pos; else if (pos > hi) hi = pos;
    }
    return 0;
}
