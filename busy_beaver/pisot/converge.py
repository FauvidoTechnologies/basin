# converge.py <long36.out>: time factor per new cell, measured (as in base/incrbin.c) over the last
# q increments (q a multiple of 6) that end before 2^L, for L = 24, 28, 32, 36; error against the
# expected base for each group.
import sys, numpy as np
phi=(1+5**.5)/2
target={'r32':1.5**.5,'r3_2':1.5,'ctl_2':2,'ctl_sqrt2':2**.5,'ctl_phi':phi,
        'ctl_plastic':1.324717957244746,'ctl_narayana':1.465571231876768}
rows=[]
for line in open(sys.argv[1]):
    f=line.split(); g,m,n=f[0],f[1],int(f[2]); ts=np.array(f[3:],dtype=float)
    out=[g,m,n]
    for L in (24,28,32,36):
        t=ts[ts<=2**L]; late=t[t>2**14]; q=(len(late)-1)//6*6
        out.append((t[-1]/t[-1-q])**(1/q) if q>=6 else np.nan)
    rows.append(out)
print('%-12s %4s  %-11s %s'%('group','n','base','median |rel err| at 2^24, 2^28, 2^32, 2^36   (worst at 2^36)'))
for g in target:
    R=[r for r in rows if r[0]==g]
    if not R: continue
    E=np.array([[abs(x/target[g]-1) for x in r[3:]] for r in R])
    print('%-12s %4d  %-11.7f '%(g,len(R),target[g])+'  '.join('%.1e'%np.nanmedian(E[:,j]) for j in range(4))+'   (%.1e)'%np.nanmax(E[:,3]))
