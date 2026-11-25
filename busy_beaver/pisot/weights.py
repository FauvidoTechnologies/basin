# weights.py <machine> <steps>: snapshot the tape each time the head arrives at the more-visited
# end of the visited region (one "increment"), then fit  n = c + sum_i w_i * cell_i  (cells counted
# from that end). A counter in a positional system fits exactly; the weights give the base.
import sys, numpy as np
def parse(m):
    T={}
    for s,blk in enumerate(m.split('_')):
        for r in range(2):
            p=blk[3*r:3*r+3]
            T[(s,r)]=None if p[2] not in 'ABCDEFGH' else (int(p[0]),-1 if p[1]=='L' else 1,ord(p[2])-65)
    return T
def run(m,N):
    T=parse(m); tape=np.zeros(1<<16,np.int8); pos=lo=hi=1<<15; st=0
    snapsL=[]; snapsR=[]; prev=pos
    for t in range(1,N+1):
        a=T[(st,tape[pos])]
        if a is None: return None
        tape[pos]=a[0]; pos+=a[1]; st=a[2]
        if pos<lo: lo=pos
        if pos>hi: hi=pos
        if pos==lo and prev!=lo: snapsL.append((t,lo,hi,tape[lo:hi+1].copy()))
        if pos==hi and prev!=hi: snapsR.append((t,lo,hi,tape[lo:hi+1].copy()))
        prev=pos
    return snapsL,snapsR
def fit(snaps,end,skip=0.25):
    snaps=snaps[int(len(snaps)*skip):]           # drop the start-up phase
    W=max(len(s[3]) for s in snaps)
    X=np.zeros((len(snaps),W))
    for k,(t,lo,hi,c) in enumerate(snaps):
        c=c if end=='L' else c[::-1]             # index cells from the low end
        X[k,:len(c)]=c
    n=np.arange(len(snaps),dtype=float)
    A=np.hstack([np.ones((len(n),1)),X])
    w,*_=np.linalg.lstsq(A,n,rcond=None)
    res=n-A@w
    return w[1:],np.abs(res).max(),len(snaps)
if __name__=='__main__':
    m=sys.argv[1]; N=int(float(sys.argv[2]))
    L,R=run(m,N)
    end,snaps=('L',L) if len(L)>len(R) else ('R',R)
    w,err,k=fit(snaps,end)
    print(m,'low end',end,'increments',len(snaps),'fitted',k,'max residual %.3g'%err)
    print('weights:',' '.join('%.4g'%x for x in w[:24]))
