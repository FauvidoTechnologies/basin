import numpy as np, pandas as pd
from multiprocessing import Pool
from wobble import analyse
meta = pd.read_csv('machines.csv')
def work(line):
    f = line.split(); a = analyse(np.array(f[2:], dtype=float)) or {}
    return dict(machine=f[0], halt=int(f[1]), **a)
if __name__ == '__main__':
    with Pool(20) as p: rows = p.map(work, open('dense32.out').readlines())
    r = pd.DataFrame(rows)
    df = meta.merge(r, on='machine', how='left').drop_duplicates(['set', 'machine'])
    df['factor'] = np.exp(df.P); df['lam'] = np.exp(df.p * df.P)
    df.to_csv('wobble_results.csv', index=False); print(len(df))
