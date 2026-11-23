import os
os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
import json, warnings, numpy as np, sympy as sp
from multiprocessing import Pool
warnings.filterwarnings("ignore")
import basin
K = 49
r = np.fromfile("sim22.out", np.int32).reshape(-1, 1 + 2 * K)
t = np.round(2.0 ** (10 + np.arange(K) / 4)); sel = slice(8, K)
tt, uu = sp.Symbol("t", positive=True), sp.Symbol("u", positive=True)
def one(i):
    warnings.filterwarnings("ignore")
    y = r[i, 1:1 + K].astype(float)[sel]; ones = r[i, 1 + K:].astype(float)[sel]
    out = dict(i=int(i), span=int(y[-1]), ones=int(ones[-1]))
    if np.ptp(y) == 0:
        out.update(cls="bounded"); return out
    for tag, x, sym in [("t", t[sel], "t"), ("u", np.log2(t[sel]), "u")]:
        try:
            law, d = basin.discover(x[:, None], y, [sym], details=True)
            out[tag] = dict(law=str(law), cv=float(d["cv_rel_mse"]))
        except Exception as e:
            out[tag] = dict(law=None, cv=float("inf"), err=repr(e)[:80])
    return out
if __name__ == "__main__":
    with Pool(20) as p, open("fits.jsonl", "w") as f:
        for k, res in enumerate(p.imap_unordered(one, range(len(r)), chunksize=8)):
            f.write(json.dumps(res) + "\n")
            if k % 2000 == 0: print(k, flush=True)
