"""refine.py: for each detection (and the BB(5) controls), refine the period on a fine grid and
measure lambda directly as the ratio of span at successive minima of the detrended wobble."""
import numpy as np, pandas as pd
from multiprocessing import Pool
spans = {l.split()[0]: np.array(l.split()[2:], float) for l in open('dense32.out')}
def fitP(x, r, P, H=3):
    w = 2*np.pi*x/P
    B = np.vstack([np.ones_like(x)] + [f(h*w) for h in range(1, H+1) for f in (np.cos, np.sin)]).T
    c, *_ = np.linalg.lstsq(B, r, rcond=None); return 1 - (r - B@c).var()/r.var()
def work(m):
    s = spans[m]; j = np.arange(len(s)); t = 2.0**(10 + j/64); ok = (s > 0) & (t >= 2**12)
    x, y = np.log(t[ok]), np.log(s[ok]); p, q = np.polyfit(x, y, 1); r = y - (p*x + q)
    d0 = float(pd.read_csv('wobble_results.csv').groupby('machine').P.first()[m])
    Ps = d0*np.exp(np.linspace(-0.03, 0.03, 601)); R = [fitP(x, r, P) for P in Ps]; P = Ps[int(np.argmax(R))]
    # minima of the wobble: one per cycle, searched in windows of length P
    mins = []; x0 = x[0]
    while x0 + P <= x[-1]:
        k = (x >= x0) & (x < x0 + P); i = np.argmin(r[k]); mins.append((x[k][i], s[ok][k][i])); x0 += P
    mins = mins[1:]  # drop the first (start-up) cycle
    sr = [b[1]/a[1] for a, b in zip(mins, mins[1:])]
    return dict(machine=m, P_fine=P, factor_fine=np.exp(P), lam_fit=np.exp(p*P), lam_minima=np.median(sr) if sr else np.nan,
                minima_ratios=' '.join(f'{v:.4f}' for v in sr))
if __name__ == '__main__':
    det = pd.read_csv('detections.csv')
    d = pd.read_csv('wobble_results.csv')
    extra = d[(d.set == 'halter5') & (d.name.isin(['BB(5) halter #1', 'BB(5) halter #3']))]
    bb5 = d[(d.set == 'bb5_sqrt') & (d.R2 > 0.5)]
    todo = pd.concat([det[['set', 'name', 'machine']], extra[['set', 'name', 'machine']], bb5[['set', 'name', 'machine']]])
    with Pool(20) as pool: res = pool.map(work, list(todo.machine))
    out = todo.merge(pd.DataFrame(res), on='machine'); out.to_csv('refined.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 60)
    print(out[out.set != 'holdout'][['name', 'factor_fine', 'lam_fit', 'lam_minima', 'minima_ratios']].round(4).to_string())
    h = out[out.set == 'holdout']
    print(h[['factor_fine', 'lam_fit', 'lam_minima']].round(3).sort_values('lam_minima').to_string())
