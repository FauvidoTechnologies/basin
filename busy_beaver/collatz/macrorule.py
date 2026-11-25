"""macrorule.py <macroc output> <out.csv>: infer each machine's macro rule from its shape runs.

For every shape with exactly one repeated block, take the block length n at each new occurrence
(first event of each run). Consecutive occurrences give transitions x -> y. For each such shape with
at least MINT transitions and max n >= NMIN:
  * slope a: least-squares slope of y on x, snapped to the nearest p/q (q <= 8) if within 2%;
  * 'single affine' if q*y - p*x is the same for every transition (exact, no rounding);
  * otherwise test residue rules: for d = 2..6, every class of x mod d with >= 2 transitions must
    be exactly affine with the same slope p/q (offsets may differ), and at least 2 classes must
    have >= 2 transitions.
Machine class, from the shape with the largest log-range of n:
  additive        a = 1 (n -> n + c): counting / bouncer-like
  affine xa       a != 1, one exact formula for all n (e.g. n -> 2n + 6)
  residue xa      a != 1, exact only per residue class mod d: Collatz-like
  irregular xa    a != 1 but no exact affine or residue fit with the data available
  no reset        no single-block shape recurs often enough"""
import sys, re, math, collections
from fractions import Fraction
import numpy as np, pandas as pd
MINT, NMIN = 3, 40

def blocks(path):
    cur = None
    for line in open(path):
        if line.startswith('M '):
            _, name, m, key, status = line.split(); cur = dict(name=name, machine=m, key=key, status=status, runs=[])
        elif line.startswith('R ') and cur is not None:
            head, _, ex = line.rstrip('\n').partition(' |')
            _, t, c, *sh = head.split(' '); cur['runs'].append((int(t), int(c), ' '.join(sh), [int(v) for v in ex.split()]))
        elif line.startswith('E') and cur is not None:
            yield cur; cur = None

def snap(a):
    best = None
    for q in range(1, 9):
        p = round(a * q)
        if p > 0 and (best is None or abs(p / q - a) < abs(best - a)): best = Fraction(p, q)
    return best if best is not None and abs(float(best) - a) <= 0.02 * a else None

def exact(x, y, f):
    return len({q * b - p * a for a, b in zip(x, y) for p, q in [(f.numerator, f.denominator)]}) == 1

def rule(x, y):
    x, y = np.array(x), np.array(y)
    a = np.polyfit(x, y, 1)[0] if np.ptp(x) > 0 else y.mean() / x.mean()
    f = snap(a)
    out = dict(slope=a, frac=str(f) if f else '', kind='')
    if f is None: out['kind'] = 'irregular'; return out
    if exact(x, y, f):
        out['kind'] = 'additive' if f == 1 else 'affine'
        out['offset'] = str(Fraction(int(y[0] * f.denominator - x[0] * f.numerator), f.denominator)); return out
    if f == 1: out['kind'] = 'additive-irregular'; return out
    for d in range(2, 7):
        cls = collections.defaultdict(list)
        for a_, b_ in zip(x, y): cls[a_ % d].append((a_, b_))
        big = [v for v in cls.values() if len(v) >= 2]
        if len(big) >= 2 and all(exact([u for u, _ in v], [w for _, w in v], f) for v in big):
            out['kind'] = 'residue'; out['mod'] = d
            out['offsets'] = ' '.join(f'{r}:{Fraction(v[0][1] * f.denominator - v[0][0] * f.numerator, f.denominator)}'
                                      for r, v in sorted(cls.items()))
            return out
    out['kind'] = 'irregular'; return out

def analyse(b):
    # occurrences of each shape with 1-3 blocks: the exponent vector at the first event of each run
    occ = collections.defaultdict(list)
    for t, c, sh, ex in b['runs']:
        if 1 <= len(ex) <= 3 and not sh.startswith('LONG'): occ[sh].append(ex)
    best = None
    for sh, vs in occ.items():
        if len(vs) < MINT + 1: continue
        for j in range(len(vs[0])):
            ns = [v[j] for v in vs]
            if max(ns) < NMIN: continue
            x, y = ns[:-1], ns[1:]
            # skip the early transient: keep transitions whose source is at least NMIN/2
            k = [i for i in range(len(x)) if x[i] >= NMIN // 2]
            if len(k) < MINT: continue
            x, y = [x[i] for i in k], [y[i] for i in k]
            # a macro coordinate must grow along the run of occurrences
            if sum(b_ > a_ for a_, b_ in zip(x, y)) < 0.75 * len(x): continue
            r = rule(x, y); r.update(shape=sh, coord=j, nblocks=len(vs[0]), ntrans=len(x), nmin=min(x),
                                     nmax=max(y), seq=' '.join(map(str, ns[-8:])),
                                     seq_all=' '.join(map(str, ns)))
            r['range'] = math.log(max(y) / max(min(x), 1))
            # prefer exact rules, then multiplicative ones, then the widest range
            pri = (r['kind'] in ('affine', 'residue', 'additive'), r['kind'] not in ('additive', 'additive-irregular'),
                   r['range'])
            if best is None or pri > best[0]: best = (pri, r)
    row = dict(name=b['name'], machine=b['machine'], key=b['key'], status=b['status'], nshapes=len(b['runs']))
    if best is None: row['kind'] = 'no reset'
    else: row.update(best[1])
    if row['kind'] in ('affine', 'residue', 'irregular') and row.get('frac'):
        row['cls'] = f"{row['kind']} x{row['frac']}"
    else: row['cls'] = row['kind']
    return row

if __name__ == '__main__':
    rows = [analyse(b) for b in blocks(sys.argv[1])]
    d = pd.DataFrame(rows); d.to_csv(sys.argv[2], index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 40)
    cols = [c for c in ['name', 'cls', 'shape', 'coord', 'ntrans', 'mod', 'offset', 'offsets', 'seq'] if c in d]
    print(d[cols].to_string(index=False))
