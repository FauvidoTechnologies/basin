# batch_weights.py <csv> <steps>: weight fit for every machine in the csv; summary per machine.
import sys, numpy as np, pandas as pd
from weights import run, fit
d=pd.read_csv(sys.argv[1]); N=int(float(sys.argv[2]))
for m in d.machine:
    L,R=run(m,N)
    end,snaps=('L',L) if len(L)>len(R) else ('R',R)
    w,err,k=fit(snaps,end)
    big=np.where(np.abs(w)>1e-6)[0]
    # base estimate: ratio of weights two cells apart, and one cell apart, at the top end (skip edge cells)
    top=[i for i in range(3,19) if abs(w[i])>0.3]
    r2=np.median([w[i+2]/w[i] for i in top if i+2<len(w) and abs(w[i])>1e-6])
    r1=np.median([w[i+1]/w[i] for i in top if i+1<len(w) and abs(w[i])>1e-6])
    print(f'{m} end={end} incr={len(snaps)} maxres={err:.2g} ratio1={r1:.5f} ratio2={r2:.5f}')
