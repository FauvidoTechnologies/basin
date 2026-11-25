# phasefit.py <steps> <machine>...: find the longest stretch with no span growth in the last half of
# the run and fit the weights only within it (the counter restarts after each overflow).
import sys, numpy as np
from weights import parse
from weights2 import fit
def go(m,N,K=8):
    T=parse(m); tape=np.zeros(1<<16,np.int8); pos=lo=hi=1<<15; st=0; prev=pos; grow=[0]
    hist=[]   # (t, end, k, lo, hi, tape)
    for t in range(1,N+1):
        a=T[(st,tape[pos])]; tape[pos]=a[0]; pos+=a[1]; st=a[2]
        if pos<lo or pos>hi: grow.append(t)
        lo=min(lo,pos); hi=max(hi,pos)
        if t>N//2:
            if pos<prev and pos-lo<K: hist.append((t,'L',pos-lo,lo,hi,tape[lo:hi+1].copy()))
            if pos>prev and hi-pos<K: hist.append((t,'R',hi-pos,lo,hi,tape[lo:hi+1].copy()))
        prev=pos
    g=np.array(grow+[N]); g=g[g>N//2] if (g>N//2).sum()>1 else g
    i=np.argmax(np.diff(g)); a,b=g[i],g[i+1]
    res=[]
    for e in 'LR':
        for k in range(K):
            s=[(lo,hi,c) for (t,ee,kk,lo,hi,c) in hist if a<t<b and ee==e and kk==k]
            r=fit(s,e,skip=0.05)
            if r: res.append((r[1],e,k,r[2],r[0]))
    res.sort(key=lambda x:x[0]); return a,b,res[0]
for m in sys.argv[2:]:
    a,b,(err,e,k,n,w)=go(m,int(float(sys.argv[1])))
    mid=[w[i+1]/w[i] for i in range(2,30) if i+1<len(w) and abs(w[i])>0.3 and abs(w[i+1])>0.3]
    print(f'{m} phase {a}-{b} trigger={e}{k} increments={n} maxres={err:.2g} ratio1(median)={np.median(mid) if mid else float("nan"):.5f}')
    print('  weights:',' '.join('%.4g'%x for x in w[:28]))
