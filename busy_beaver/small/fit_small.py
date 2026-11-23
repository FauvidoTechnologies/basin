# Same fits as census/fit_all.py; the log-t fit is skipped when the raw 2^20..2^22 slope exceeds 0.4.
import os, sys
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import json, warnings, numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
import basin
SIM, OUT, WORKERS = sys.argv[1], sys.argv[2], int(sys.argv[3])
K = 49
r = np.fromfile(SIM, np.int32).reshape(-1, 1 + 2 * K)
t = np.round(2.0 ** (10 + np.arange(K) / 4)); sel = slice(8, K)
def one(i):
    warnings.filterwarnings("ignore")
    y = r[i, 1:1 + K].astype(float)[sel]; ones = r[i, 1 + K:].astype(float)[sel]
    out = dict(i=int(i), span=int(y[-1]), ones=int(ones[-1]))
    if np.ptp(y) == 0:
        out.update(cls="bounded"); return out
    slope = np.log(y[-1] / y[-9]) / np.log(4)
    for tag, x, sym in [("t", t[sel], "t"), ("u", np.log2(t[sel]), "u")]:
        if tag == "u" and slope > 0.4:
            out[tag] = dict(law=None, cv=float("inf"), err="skipped: raw slope > 0.4"); continue
        try:
            law, d = basin.discover(x[:, None], y, [sym], details=True)
            out[tag] = dict(law=str(law), cv=float(d["cv_rel_mse"]))
        except Exception as e:
            out[tag] = dict(law=None, cv=float("inf"), err=repr(e)[:80])
    return out
if __name__ == "__main__":
    with Pool(WORKERS) as p, open(OUT, "w") as f:
        for res in p.imap_unordered(one, range(len(r)), chunksize=8):
            f.write(json.dumps(res) + "\n")
