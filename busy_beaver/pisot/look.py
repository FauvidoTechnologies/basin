# look.py <machine> <steps> [every]: print the tape (lo..hi) each time the span grows.
import sys
def parse(m):
    T={}
    for s,blk in enumerate(m.split('_')):
        for r in range(2):
            p=blk[3*r:3*r+3]
            T[(s,r)]=None if p[2] not in 'ABCDEFGH' else (int(p[0]),-1 if p[1]=='L' else 1,ord(p[2])-65)
    return T
m=sys.argv[1]; N=int(float(sys.argv[2])); T=parse(m)
tape={}; pos=lo=hi=0; st=0
for t in range(1,N+1):
    a=T[(st,tape.get(pos,0))]
    if a is None: print('halt',t); break
    tape[pos]=a[0]; pos+=a[1]; st=a[2]
    if pos<lo or pos>hi:
        lo=min(lo,pos); hi=max(hi,pos)
        s=''.join(str(tape.get(i,0)) for i in range(lo,hi+1))
        print(f'{t:>10} {hi-lo+1:>3} {"L" if pos==lo else "R"} {s}')
