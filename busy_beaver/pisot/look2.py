import sys
def parse(m):
    T={}
    for s,blk in enumerate(m.split('_')):
        for r in range(2):
            p=blk[3*r:3*r+3]
            T[(s,r)]=None if p[2] not in 'ABCDEFGH' else (int(p[0]),-1 if p[1]=='L' else 1,ord(p[2])-65)
    return T
m=sys.argv[1]; t0=int(float(sys.argv[2])); t1=int(float(sys.argv[3])); T=parse(m)
tape={}; pos=lo=hi=0; st=0; last=0
for t in range(1,t1+1):
    a=T[(st,tape.get(pos,0))]
    tape[pos]=a[0]; d=a[1]; 
    if t>=t0 and d!=last:
        s=''.join(('[%d]'%tape.get(i,0)) if i==pos else str(tape.get(i,0)) for i in range(lo,hi+1))
        print(f'{t:>9} {"ABCDE"[st]}{">" if d>0 else "<"} {s}')
    last=d; pos+=d; st=a[2]; lo=min(lo,pos); hi=max(hi,pos)
