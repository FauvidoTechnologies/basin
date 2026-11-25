"""score.py <spans.bin> <indices.int64> <out.csv>: scale-invariance statistic for every machine."""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
from scalestat import stat
def chunk(a): return [stat(s) for s in a]
if __name__ == '__main__':
    S = np.fromfile(sys.argv[1], dtype=np.int32).reshape(-1, 257); idx = np.fromfile(sys.argv[2], dtype=np.int64)[:len(S)]
    with Pool(20) as p: r = sum(p.map(chunk, np.array_split(S, 200)), [])
    d = pd.DataFrame(r, columns=['p', 'D', 'S']); d.insert(0, 'db_index', idx); d['span_2_26'] = S[:, -1]
    d.to_csv(sys.argv[3], index=False); print(len(d))
