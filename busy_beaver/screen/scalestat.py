"""stat.py: scale-invariance statistic for span series sampled at t = 2^(10 + j/16), j = 0..16*(L-10).
Residual of ln span after a straight-line fit in log2 t on [2^12, 2^L]; high-pass by subtracting a
running mean over 2 octaves; D = log2(std late / std early), early = [2^13, 2^18], late =
[2^(L-5.5), 2^(L-1)]. A wobble that fades like 1/span gives D = -p * (gap between window centres);
S = D / that prediction: S ~ 1 for additive growth (bouncers), S ~ 0 for scale-invariant growth."""
import numpy as np
def stat(s, L=26):
    s = np.asarray(s, float); u = 10 + np.arange(len(s)) / 16
    ok = (s > 0) & (u >= 12)
    if ok.sum() < 16 * (L - 12) or s[-1] <= 0: return np.nan, np.nan, np.nan
    x, y = u[ok], np.log(s[ok]); c = np.polyfit(x, y, 1); p = c[0] / np.log(2)
    r = y - np.polyval(c, x); hp = r - np.convolve(r, np.ones(32) / 32, mode='same')
    e = hp[(x >= 13) & (x <= 18)].std(); l = hp[(x >= L - 5.5) & (x <= L - 1)].std()
    D = np.log2(max(l, 1e-12) / max(e, 1e-12)); gap = (L - 3.25) - 15.5
    return p, D, D / (-gap * p) if p > 0.05 else np.nan
