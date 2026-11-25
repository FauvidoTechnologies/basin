"""multipliers.py: for scale-invariant machines in a scored sample, the multiplier lambda from the
span ratio between successive teeth (sharp minima of the residual), and distinct behaviours."""
import sys, numpy as np, pandas as pd, hashlib
from fractions import Fraction
from scipy.signal import find_peaks
from multiprocessing import Pool
def lam(s):
    s = np.asarray(s, float); u = 10 + np.arange(len(s)) / 16; ok = (s > 0) & (u >= 12)
    x, y = u[ok], np.log(s[ok]); r = y - np.polyval(np.polyfit(x, y, 1), x)
    amp = np.percentile(r, 95) - np.percentile(r, 5)
    i, _ = find_peaks(-r, prominence=0.25 * amp, distance=3)
    v = s[ok][i]; rat = v[1:] / v[:-1]
    return (np.median(rat[-4:]) if len(rat) >= 2 else np.nan), len(i), ' '.join(f'{q:.3f}' for q in rat[-6:])
def near(x):
    if not np.isfinite(x): return 'none'
    f, e = min(((Fraction(a, b), abs(a / b / x - 1)) for b in range(1, 9) for a in range(b + 1, 9 * b)), key=lambda z: z[1])
    return str(f) if e < 0.005 else 'none'
if __name__ == '__main__':
    d = pd.read_csv('sampleA_scores.csv'); S = np.fromfile('sampleA.bin', dtype=np.int32).reshape(-1, 257)
    k = np.where(d.S < 0.25)[0]
    d = d.iloc[k].copy(); rows = S[k]
    d['fp'] = [hashlib.md5(r.tobytes()).hexdigest()[:12] for r in rows]
    with Pool(20) as p: res = p.map(lam, list(rows))
    d['lam'] = [r[0] for r in res]; d['teeth'] = [r[1] for r in res]; d['ratios'] = [r[2] for r in res]
    d['frac'] = d.lam.apply(near)
    d.to_csv('sampleA_scaleinv.csv', index=False)
    print('scale-invariant records', len(d), '| distinct span sequences', d.fp.nunique())
    u = d.drop_duplicates('fp')
    print(u.frac.value_counts().head(15).to_string())
    print('records per fraction:'); print(d.frac.value_counts().head(10).to_string())
