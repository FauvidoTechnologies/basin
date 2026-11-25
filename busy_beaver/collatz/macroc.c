/* macroc.c: C version of macro2.py's shape-run extraction, for batch use.
   usage: macroc <machines.txt> <steps>   (each line: "<name> <machine>")
   For each machine: pass 1 finds the (state, side) edge visit most common in the second half;
   pass 2 compresses the trimmed tape at every such visit into tokens (a unit of 1-4 cells
   repeated 3+ times -> (u)^k, else literal; same rule as macro2.py's regex) and prints one line
   per maximal run of equal shapes:  R <t_first> <count> <shape> | <exponents at first event>
   Output starts with "M <name> <machine> <key> <status>" and ends with "E". */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define TAPE (1 << 21)
#define MAXS 8
#define WIDE 100000  /* the screen is about sqrt(t) machines (span ~2e4 at 2^26); skip wide tapes */
static int W[MAXS][2], D[MAXS][2], N[MAXS][2];
static unsigned char tape[TAPE];
static int ns;
static int parse(const char *m) {
    ns = 0;
    for (const char *p = m; *p; ) {
        for (int b = 0; b < 2; b++, p += 3) {
            if (p[0] == '-' || p[2] == 'Z' || p[2] == 'H') { N[ns][b] = -1; W[ns][b] = 1; D[ns][b] = 1; }
            else { W[ns][b] = p[0] - '0'; D[ns][b] = p[1] == 'R' ? 1 : -1; N[ns][b] = p[2] - 'A'; }
        }
        ns++; if (*p == '_') p++;
    }
    return ns;
}
static char sbuf[TAPE + 2], shape[1 << 20], prev[1 << 20];
static long expo[1 << 16]; static int nexp;
static void tokenize(int n) {        /* sbuf[0..n) of '0'/'1' -> shape string + exponents */
    int i = 0, lit = -1; char *o = shape; nexp = 0;
    while (i < n) {
        int got = 0;
        for (int u = 1; u <= 4 && i + 3 * u <= n; u++) {
            int k = 1;
            while (i + (k + 1) * u <= n && !memcmp(sbuf + i, sbuf + i + k * u, u)) k++;
            if (k >= 3) {
                if (lit >= 0) { o += sprintf(o, "%.*s ", i - lit, sbuf + lit); lit = -1; }
                o += sprintf(o, "(%.*s)^n ", u, sbuf + i); expo[nexp++] = k; i += k * u; got = 1; break;
            }
        }
        if (!got) { if (lit < 0) lit = i; i++; }
    }
    if (lit >= 0) o += sprintf(o, "%.*s ", n - lit, sbuf + lit);
    if (o > shape) o[-1] = 0; else *o = 0;
}
static long cnt[MAXS][2];
/* one line per shape run; shapes longer than 200 characters are printed as LONG<len> (they are never
   single-block resets, which is all the rule analysis uses) and at most 8 exponents are printed */
static void emit(long t, long c, long *e, int ne) {
    size_t L = strlen(prev);
    if (L <= 200) printf("R %ld %ld %s |", t, c, prev); else printf("R %ld %ld LONG%zu |", t, c, L);
    for (int i = 0; i < ne && i < 8; i++) printf(" %ld", e[i]);
    printf("\n");
}
static int run(long steps, int pass, int key) {   /* returns 0 ok, 1 halted, 2 left tape, 3 wider than WIDE */
    memset(tape, 0, sizeof tape);
    long pos = TAPE / 2, lo = pos, hi = pos; int st = 0;
    long tfirst = 0, c = 0; long e0[64]; int ne0 = 0; prev[0] = 1; prev[1] = 0;
    for (long t = 1; t <= steps; t++) {
        int b = tape[pos], nx = N[st][b];
        if (nx < 0) return 1;
        tape[pos] = W[st][b]; pos += D[st][b]; st = nx;
        if (pos <= 0 || pos >= TAPE - 1) return 2;
        if (pos < lo) lo = pos; if (pos > hi) hi = pos;
        if (hi - lo > WIDE) return 3;
        if (pos != lo && pos != hi) continue;
        int side = pos == lo ? 0 : 1;
        if (pass == 1) { if (t > steps / 2) cnt[st][side]++; continue; }
        if (st * 2 + side != key) continue;
        long a = lo, z = hi; while (a <= z && !tape[a]) a++; while (z >= a && !tape[z]) z--;
        int n = 0; for (long i = a; i <= z; i++) sbuf[n++] = '0' + tape[i];
        tokenize(n);
        if (prev[0] != 1 && !strcmp(shape, prev)) { c++; continue; }
        if (prev[0] != 1) emit(tfirst, c, e0, ne0);
        strcpy(prev, shape); tfirst = t; c = 1; ne0 = nexp < 64 ? nexp : 64; memcpy(e0, expo, ne0 * sizeof(long));
    }
    if (pass == 2 && prev[0] != 1) emit(tfirst, c, e0, ne0);
    return 0;
}
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "r"); long steps = (long)atof(argv[2]); char name[256], m[256];
    while (fscanf(f, "%255s %255s", name, m) == 2) {
        parse(m); memset(cnt, 0, sizeof cnt);
        int r = run(steps, 1, -1), key = 0;
        for (int s = 0; s < ns * 2; s++) if (cnt[s / 2][s % 2] > cnt[key / 2][key % 2]) key = s;
        printf("M %s %s %c%c %s\n", name, m, 'A' + key / 2, key % 2 ? 'R' : 'L', r == 1 ? "halt" : r == 2 ? "offtape" : r == 3 ? "wide" : "ok");
        if (r == 0) run(steps, 2, key);
        printf("E\n"); fflush(stdout);
    }
    return 0;
}
