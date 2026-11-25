"""teeth.py: for named machines, find the teeth of the wobble (sharp minima of the residual of
ln span after a straight-line fit in ln t) and print the span ratio between successive teeth."""
import numpy as np, pandas as pd
from scipy.signal import find_peaks
spans = {l.split()[0]: np.array(l.split()[2:], float) for l in open('../wobble/dense32.out')}
def teeth(m, lo=12):
    s = spans[m]; u = 10 + np.arange(len(s)) / 64; ok = (s > 0) & (u >= lo)
    x, y = u[ok], np.log(s[ok]); r = y - np.polyval(np.polyfit(x, y, 1), x)
    amp = np.percentile(r, 95) - np.percentile(r, 5)
    idx, _ = find_peaks(-r, prominence=0.25 * amp, distance=8)
    return x[idx], s[ok][idx]
if __name__ == '__main__':
    sp = pd.read_csv('../bb6/special.csv'); sp = sp[sp.set.isin(['cryptid', 'potential'])]
    for _, row in sp.iterrows():
        x, sv = teeth(row.machine); rat = sv[1:] / sv[:-1]
        print(f"{row['name']:<17} teeth={len(x):>3}  ratios: " + ' '.join(f'{v:.3f}' for v in rat[-18:]))
