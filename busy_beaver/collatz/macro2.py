"""macro2.py <machine> <steps> [rows]: macro structure with repeated units.
At each visit of the head to an edge (edge zeros trimmed), compress the tape into tokens: a unit u of
1-4 cells repeated 3+ times becomes (u)^k, anything else stays literal. The tape's *shape* is the
token list with each exponent replaced by 'n'. For the (state, side) visited most in the second half,
print the run of events grouped by shape: for each maximal run of equal shapes, its first and last
event with the exponents. Shape changes are where macro steps happen."""
import sys, re, collections
from macro import parse
TOK = re.compile(r'([01]{1,4}?)\1{2,}')
def tokens(s):
    out, i = [], 0
    for mm in TOK.finditer(s):
        if mm.start() > i: out.append((s[i:mm.start()], 1))
        out.append((mm.group(1), len(mm.group(0)) // len(mm.group(1)))); i = mm.end()
    if i < len(s): out.append((s[i:], 1))
    return out
def shape(tk): return ' '.join(f'({u})^n' if k > 1 else u for u, k in tk)
def expo(tk): return [k for u, k in tk if k > 1]
m = sys.argv[1]; N = int(float(sys.argv[2])); rows = int(sys.argv[3]) if len(sys.argv) > 3 else 30
T = parse(m); tape = collections.defaultdict(int); pos = lo = hi = 0; st = 0; ev = []
for t in range(1, N + 1):
    a = T[(st, tape[pos])]
    if a is None: print('halt', t); break
    tape[pos] = a[0]; pos += a[1]; st = a[2]; lo = min(lo, pos); hi = max(hi, pos)
    if pos == lo or pos == hi:
        ev.append((t, 'ABCDEF'[st] + ('L' if pos == lo else 'R'), lo, hi, pos))
        if len(ev) > 1 and ev[-2][0] == t - 1 and ev[-2][1][1] == ev[-1][1][1]: pass
key = collections.Counter(e[1] for e in ev[len(ev)//2:]).most_common(1)[0][0]
# re-run to take snapshots only at that key (cheaper than storing tapes for every event)
T = parse(m); tape = collections.defaultdict(int); pos = lo = hi = 0; st = 0; groups = []
for t in range(1, N + 1):
    a = T[(st, tape[pos])]
    if a is None: break
    tape[pos] = a[0]; pos += a[1]; st = a[2]; lo = min(lo, pos); hi = max(hi, pos)
    if (pos == lo or pos == hi) and ('ABCDEF'[st] + ('L' if pos == lo else 'R')) == key:
        tk = tokens(''.join(str(tape[i]) for i in range(lo, hi + 1)).strip('0')); sh = shape(tk)
        if groups and groups[-1][0] == sh: groups[-1][2] = (t, expo(tk)); groups[-1][3] += 1
        else: groups.append([sh, (t, expo(tk)), (t, expo(tk)), 1])
print(f'key {key}; {len(groups)} shape runs; last {rows}:')
for sh, (t0, e0), (t1, e1), c in groups[-rows:]:
    print(f'  {c:>5}x  {sh:<40.40} first t={t0} {e0}   last t={t1} {e1}')
