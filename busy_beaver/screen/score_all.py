"""score_all.py: scale-invariance index for all 2,767,634 sqrt(t)-like records, at 2^26 and, from the
first 225 samples of the same rows, at 2^24. Writes all_scores.csv."""
import numpy as np, pandas as pd
from multiprocessing import Pool
from scalestat import stat
J24 = 16 * (24 - 10) + 1
def chunk(a): return [stat(s) + stat(s[:J24], L=24) for s in a]
if __name__ == '__main__':
    S = np.fromfile('all.bin', dtype=np.int32).reshape(-1, 257); idx = np.fromfile('all.int64', dtype=np.int64)
    with Pool(22) as p: r = sum(p.map(chunk, np.array_split(S, 2000), chunksize=4), [])
    d = pd.DataFrame(r, columns=['p', 'D', 'S', 'p24', 'D24', 'S24']); d.insert(0, 'db_index', idx)
    d['span_2_26'] = S[:, -1]; d['span_2_24'] = S[:, J24 - 1]; d['sampleA'] = np.arange(len(d)) < 300000
    d.to_csv('all_scores.csv', index=False); print(len(d))
