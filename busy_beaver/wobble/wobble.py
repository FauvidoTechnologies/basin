"""wobble.py <dense.out> [min_log2t]: look for a wobble that repeats on a log scale.
For each machine: x = ln t, y = ln span on t >= 2^min (default 12), detrend with a straight line
(slope p = growth exponent), then for each trial period P (in ln t) fit 3 harmonics to the residual
and keep the best P. Reports R2 (share of residual variance the periodic fit explains), P, the time
factor e^P per cycle, and the implied multiplier lam = e^(p*P) of a rule n -> lam*n with span ~ n."""
import sys, numpy as np
def analyse(span, lo=12, H=3):
    j = np.arange(len(span)); t = 2.0 ** (10 + j / 64)
    ok = (span > 0) & (t >= 2 ** lo)
    if ok.sum() < 200: return None
    x, y = np.log(t[ok]), np.log(span[ok].astype(float))
    A = np.vstack([np.ones_like(x), x]).T; c, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ c; p = c[1]; v = r.var()
    if v < 1e-12: return dict(p=p, R2=0.0, P=np.nan, amp=0.0, n=ok.sum())
    span_x = x[-1] - x[0]; best = (0, np.nan, 0)
    for P in np.exp(np.linspace(np.log(0.15), np.log(span_x / 2.5), 600)):
        w = 2 * np.pi * x / P
        B = np.vstack([np.ones_like(x)] + [f(h * w) for h in range(1, H + 1) for f in (np.cos, np.sin)]).T
        cc, *_ = np.linalg.lstsq(B, r, rcond=None); R2 = 1 - (r - B @ cc).var() / v
        if R2 > best[0]: best = (R2, P, np.std(B @ cc))
    R2, P, amp = best
    return dict(p=p, R2=R2, P=P, amp=amp, n=ok.sum())
if __name__ == "__main__":
    lo = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    for line in open(sys.argv[1]):
        f = line.split(); m, halt = f[0], int(f[1]); span = np.array(f[2:], dtype=float)
        a = analyse(span, lo)
        if a is None: print(m, halt, "too few points"); continue
        print(f"{m} halt={halt} p={a['p']:.4f} R2={a['R2']:.3f} P={a['P']:.4f} factor={np.exp(a['P']):.4f} "
              f"lam={np.exp(a['p'] * a['P']):.4f} amp={a['amp']:.2e} n={a['n']}")
