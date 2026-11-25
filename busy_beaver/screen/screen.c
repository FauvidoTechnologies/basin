// screen: span at t_j = round(2^(10 + j/16)), j = 0..16*(L-10), for many machines; int32 rows, -1 after
// halting or leaving the tape.
//   screen rec <records.bin> <indices.int64> <start> <end> <L> <out.bin>   (bbchallenge 30-byte records)
//   screen txt <machines.txt> <L> <out.bin>                                 (machine strings, last field)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
static long ck[1024]; static int J; static unsigned char *tape; static long N = 1L << 26;
static void run(int n, int *W, int *M, int *G, long T, int32_t *row) {
    long pos = N / 2, lo = pos, hi = pos, t = 0; int st = 0, j = 0;
    while (j < J) {
        int c = 2 * st + tape[pos]; if (G[c] < 0) break;
        tape[pos] = W[c]; pos += M[c]; st = G[c]; t++;
        if (pos < lo) lo = pos; else if (pos > hi) hi = pos;
        if (pos < 1 || pos > N - 2) break;
        while (j < J && t == ck[j]) row[j++] = (int32_t)(hi - lo + 1);
    }
    for (; j < J; j++) row[j] = -1;
    memset(tape + lo, 0, hi - lo + 1);
}
int main(int argc, char **argv) {
    tape = calloc(N, 1); int rec = !strcmp(argv[1], "rec");
    int L = atoi(rec ? argv[6] : argv[3]); J = 16 * (L - 10) + 1;
    for (int j = 0; j < J; j++) ck[j] = lround(pow(2.0, 10 + j / 16.0));
    FILE *o = fopen(rec ? argv[7] : argv[4], "wb"); int32_t row[1024]; int W[16], M[16], G[16];
    if (rec) {
        FILE *f = fopen(argv[2], "rb"), *fi = fopen(argv[3], "rb"); long a = atol(argv[4]), b = atol(argv[5]);
        unsigned char r[30]; int64_t idx; fseek(fi, 8 * a, SEEK_SET);
        for (long i = a; i < b; i++) {
            if (fread(&idx, 8, 1, fi) != 1) break;
            fseek(f, 30 * idx, SEEK_SET); if (fread(r, 1, 30, f) != 30) break;
            for (int k = 0; k < 10; k++) { W[k] = r[3*k]; M[k] = r[3*k+1] ? -1 : 1; G[k] = (int)r[3*k+2] - 1; }
            run(5, W, M, G, 0, row); fwrite(row, 4, J, o);
        }
    } else {
        FILE *f = fopen(argv[2], "r"); char line[512];
        while (fgets(line, sizeof line, f)) {
            char *m = NULL, *tok = strtok(line, " ,\t\r\n"); while (tok) { m = tok; tok = strtok(NULL, " ,\t\r\n"); }
            if (!m || strlen(m) < 13) continue;
            int n = (int)(strlen(m) + 1) / 7;
            for (int s = 0; s < n; s++) for (int q = 0; q < 2; q++) { char *p = m + 7 * s + 3 * q;
                W[2*s+q] = p[0] == '1'; M[2*s+q] = p[1] == 'L' ? -1 : 1; G[2*s+q] = (p[2] >= 'A' && p[2] < 'A' + n) ? p[2] - 'A' : -1; }
            run(n, W, M, G, 0, row); fwrite(row, 4, J, o);
        }
    }
    fclose(o); return 0;
}
