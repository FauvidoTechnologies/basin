"""decay.py: does the relative wobble of the span decay? For each machine in ../wobble/dense32.out:
residual of ln span after a straight-line fit in ln t (t >= 2^12), high-passed by subtracting a
running mean over 2 octaves; amplitude = std of the high-passed residual in an early window
[2^13, 2^18] and a late window [2^27, 2^32]. D = log2(late/early). Bouncers: wobble ~ 1/span, so
D ~ -p * 14 (about -7). Multiplicative (Collatz-like) macro steps: D ~ 0."""
import numpy as np, pandas as pd
from multiprocessing import Pool
def stat(line):
    f = line.split(); s = np.array(f[2:], float); m = f[0]
    j = np.arange(len(s)); u = 10 + j / 64; ok = (s > 0) & (u >= 12)
    if ok.sum() < 64 * 19: return dict(machine=m, p=np.nan, D=np.nan, early=np.nan, late=np.nan)
    x, y = u[ok], np.log(s[ok]); p = np.polyfit(x * np.log(2), y, 1)[0]
    r = y - np.polyval(np.polyfit(x, y, 1), x)
    k = 128; hp = r - np.convolve(r, np.ones(k) / k, mode='same')
    core = (x >= 13) & (x <= 31)                      # avoid the running mean's edge effects
    e = hp[core & (x <= 18)].std(); l = hp[core & (x >= 26)].std()
    return dict(machine=m, p=p, D=np.log2(max(l, 1e-12) / max(e, 1e-12)), early=e, late=l)
if __name__ == '__main__':
    with Pool(20) as pool: rows = pool.map(stat, open('../wobble/dense32.out').readlines())
    d = pd.read_csv('../wobble/machines.csv').merge(pd.DataFrame(rows).drop_duplicates('machine'), on='machine')
    d.to_csv('decay.csv', index=False)
    pd.set_option('display.width', 200)
    print(d.groupby('set').D.describe(percentiles=[.1, .5, .9]).round(2))
    print(d[d.set.isin(['cryptid', 'potential', 'champion6'])][['name', 'p', 'D', 'late']].round(3).to_string())
