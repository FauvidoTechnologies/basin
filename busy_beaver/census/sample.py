# Uniformly sample machine records from the bbchallenge seed database (streamed from the zip).
import sys, zipfile, numpy as np
db, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
TOTAL, NTIME = 88664064, 14322029
idx = np.sort(np.random.default_rng(0).choice(TOTAL, n, replace=False))
z = zipfile.ZipFile(db); f = z.open(z.infolist()[0]); f.read(30)
CH = 1 << 20  # machines per chunk
recs, j, base = [], 0, 0
while j < n:
    buf = f.read(30 * CH)
    if not buf: break
    m = len(buf) // 30
    while j < n and idx[j] < base + m:
        k = idx[j] - base; recs.append(buf[30*k:30*k+30]); j += 1
    base += m
open(out + ".bin", "wb").write(b"".join(recs))
np.save(out + "_idx.npy", idx)
print(len(recs), "sampled;", int((idx < NTIME).sum()), "from time-limit group")
