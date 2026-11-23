# fit_window <workers> <kmin> <kmax> <out.jsonl> <tag:simfile:K:n|all> ...
# Basin fits of span on checkpoints kmin..kmax (t = 2^(10+k/4)), same rules as bb6/fit_any.py.
import os, sys
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
import basin
W, KMIN, KMAX, OUT = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
JOBS, FILES = [], {}
for spec in sys.argv[5:]:
    tag, path, K, n = spec.split(":"); K = int(K)
    r = np.fromfile(path, np.int32).reshape(-1, 1 + 2 * K); FILES[path] = K
    idx = range(len(r)) if n == "all" else np.random.default_rng(0).choice(len(r), int(n), replace=False)
    JOBS += [(tag, path, int(i)) for i in idx]
CACHE = {}
def one(job):
    tag, path, i = job; warnings.filterwarnings("ignore"); K = FILES[path]
    if path not in CACHE: CACHE[path] = np.fromfile(path, np.int32).reshape(-1, 1 + 2 * K)
    r = CACHE[path][i]; t = np.round(2.0 ** (10 + np.arange(K) / 4))
    sl = slice(KMIN, KMAX + 1); y, x = r[1:1 + K].astype(float)[sl], t[sl]
    out = dict(set=tag, i=i, halt=int(r[0]), span=int(y[-1]))
    if (y <= 0).any(): out.update(cls="halted_in_window"); return out
    if np.ptp(y) == 0: out.update(cls="bounded"); return out
    slope = np.log(y[-1] / y[-9]) / np.log(x[-1] / x[-9]); out["raw_slope"] = float(slope)
    for tg, xx, sym in [("t", x, "t"), ("u", np.log2(x), "u")]:
        if tg == "u" and slope > 0.4: out[tg] = dict(law=None, cv=float("inf")); continue
        try:
            law, d = basin.discover(xx[:, None], y, [sym], details=True)
            out[tg] = dict(law=str(law), cv=float(d["cv_rel_mse"]))
        except Exception as e:
            out[tg] = dict(law=None, cv=float("inf"), err=repr(e)[:80])
    return out
if __name__ == "__main__":
    with Pool(W) as p, open(OUT, "w") as f:
        for res in p.imap_unordered(one, JOBS, chunksize=4):
            f.write(json.dumps(res) + "\n")
