// scan <records.bin> <start> <end> <out.bin>: span (int32) at 2^12, 2^14, 2^16 steps for each 30-byte record.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
int main(int argc, char **argv) {
    FILE *f = fopen(argv[1], "rb"); long a = atol(argv[2]), b = atol(argv[3]); FILE *o = fopen(argv[4], "wb");
    const long T = 1 << 16; static unsigned char tape[2 * (1 << 16) + 16]; unsigned char rec[30];
    fseek(f, 30 * a, SEEK_SET);
    for (long i = a; i < b; i++) {
        if (fread(rec, 1, 30, f) != 30) break;
        int W[10], M[10], G[10];
        for (int k = 0; k < 10; k++) { W[k] = rec[3*k]; M[k] = rec[3*k+1] ? -1 : 1; G[k] = (int)rec[3*k+2] - 1; }
        long pos = T + 8, lo = pos, hi = pos; int st = 0; int32_t out[3] = {-1, -1, -1};
        for (long t = 1; t <= T; t++) {
            int c = 2 * st + tape[pos]; if (G[c] < 0) break;
            tape[pos] = W[c]; pos += M[c]; st = G[c];
            if (pos < lo) lo = pos; else if (pos > hi) hi = pos;
            if (t == 1 << 12) out[0] = hi - lo + 1; else if (t == 1 << 14) out[1] = hi - lo + 1; else if (t == T) out[2] = hi - lo + 1;
        }
        fwrite(out, 4, 3, o); memset(tape + lo, 0, hi - lo + 1);
    }
    fclose(o); return 0;
}
