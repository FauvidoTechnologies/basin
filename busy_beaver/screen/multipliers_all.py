"""multipliers_all.py: multiplier (teeth ratio) and behaviour fingerprint for every scale-invariant
record of the full screen. Writes all_scaleinv.csv."""
import numpy as np, pandas as pd, hashlib
from multiprocessing import Pool
from multipliers import lam, near
if __name__ == '__main__':
    d = pd.read_csv('all_scores.csv'); S = np.memmap('all.bin', dtype=np.int32, mode='r').reshape(-1, 257)
    k = np.where(d.S < 0.25)[0]; d = d.iloc[k].copy(); rows = np.array(S[k])
    d['fp'] = [hashlib.md5(r.tobytes()).hexdigest()[:12] for r in rows]
    with Pool(22) as p: res = p.map(lam, list(rows), chunksize=200)
    d['lam'] = [r[0] for r in res]; d['teeth'] = [r[1] for r in res]; d['ratios'] = [r[2] for r in res]
    d['frac'] = d.lam.apply(near)
    d.to_csv('all_scaleinv.csv', index=False)
    print('scale-invariant records', len(d), '| distinct span sequences', d.fp.nunique())
