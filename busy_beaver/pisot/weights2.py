# weights2.py: like weights.py, but an "increment" is the head arriving (moving towards the end) at
# cell k from either end, for k = 0..7; reports the trigger with the best (most exact) fit.
import sys, numpy as np
from weights import parse
def run(m,N,K=8):
    T=parse(m); tape=np.zeros(1<<16,np.int8); pos=lo=hi=1<<15; st=0
    snaps={(e,k):[] for e in 'LR' for k in range(K)}; prev=pos
    for t in range(1,N+1):
        a=T[(st,tape[pos])]
        if a is None: return None
        tape[pos]=a[0]; pos+=a[1]; st=a[2]
        lo=min(lo,pos); hi=max(hi,pos)
        if pos<prev and pos-lo<K: snaps[('L',pos-lo)].append((lo,hi,tape[lo:hi+1].copy()))
        if pos>prev and hi-pos<K: snaps[('R',hi-pos)].append((lo,hi,tape[lo:hi+1].copy()))
        prev=pos
    return snaps
def fit(snaps,end,skip=0.25):
    snaps=snaps[int(len(snaps)*skip):]
    if len(snaps)<200: return None
    W=max(len(s[2]) for s in snaps); X=np.zeros((len(snaps),W))
    for k,(lo,hi,c) in enumerate(snaps):
        c=c if end=='L' else c[::-1]; X[k,:len(c)]=c
    n=np.arange(len(snaps),dtype=float); A=np.hstack([np.ones((len(n),1)),X])
    w,*_=np.linalg.lstsq(A,n,rcond=None)
    return w[1:],np.abs(n-A@w).max(),len(snaps)
def best(m,N):
    S=run(m,N); out=[]
    for (e,k),s in S.items():
        r=fit(s,e)
        if r: out.append((r[1],e,k,r[2],r[0]))
    out.sort(key=lambda x:x[0]); return out[0]
if __name__=='__main__':
    for m in sys.argv[2:]:
        err,e,k,n,w=best(m,int(float(sys.argv[1])))
        mid=[w[i+2]/w[i] for i in range(4,20) if i+2<len(w) and abs(w[i])>0.3]
        print(f'{m} trigger={e}{k} increments={n} maxres={err:.2g} ratio2={np.median(mid) if mid else float("nan"):.5f}')
        print('  weights:',' '.join('%.4g'%x for x in w[:26]))
