# Same fits as small/fit_small.py, for several sim files; uses only checkpoints before any halt.
import os, sys
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
import basin
K = 49
t = np.round(2.0 ** (10 + np.arange(K) / 4))
JOBS = []  # (tag, sim file, indices)
for spec in sys.argv[2:]:
    tag, path, n = spec.split(":"); r = np.fromfile(path, np.int32).reshape(-1, 1 + 2 * K)
    idx = range(len(r)) if n == "all" else np.random.default_rng(0).choice(len(r), int(n), replace=False)
    JOBS += [(tag, path, int(i)) for i in idx]
CACHE = {}
def rows(path):
    if path not in CACHE: CACHE[path] = np.fromfile(path, np.int32).reshape(-1, 1 + 2 * K)
    return CACHE[path]
def one(job):
    tag, path, i = job; warnings.filterwarnings("ignore")
    r = rows(path)[i]; span = r[1:1 + K].astype(float)
    ok = np.arange(K) >= 8; ok &= span > 0
    y, x = span[ok], t[ok]
    out = dict(set=tag, i=i, span=int(y[-1]), last_t=int(x[-1]), halt=int(r[0]))
    if np.ptp(y) == 0: out.update(cls="bounded"); return out
    slope = np.log(y[-1] / y[-9]) / np.log(x[-1] / x[-9])
    out["raw_slope"] = float(slope)
    for tg, xx, sym in [("t", x, "t"), ("u", np.log2(x), "u")]:
        if tg == "u" and slope > 0.4:
            out[tg] = dict(law=None, cv=float("inf"), err="skipped: raw slope > 0.4"); continue
        try:
            law, d = basin.discover(xx[:, None], y, [sym], details=True)
            out[tg] = dict(law=str(law), cv=float(d["cv_rel_mse"]))
        except Exception as e:
            out[tg] = dict(law=None, cv=float("inf"), err=repr(e)[:80])
    return out
if __name__ == "__main__":
    with Pool(int(sys.argv[1])) as p, open("bb6_fits.jsonl", "w") as f:
        for res in p.imap_unordered(one, JOBS, chunksize=1):
            f.write(json.dumps(res) + "\n"); f.flush()
