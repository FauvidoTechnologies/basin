"""macro.py <machine> <steps> [k]: find a machine's macro configurations.
At every visit of the head to an edge of the visited region (in state S, on side L/R), write the tape
as runs with edge zeros trimmed. For each (S, side) take the smallest run count seen in the last
half of the run; the configurations with that count are candidate macro configurations. Prints the
(S, side) with the most such events and its last k configurations as run-length blocks."""
import sys, itertools, collections
def parse(m):
    T = {}
    for s, blk in enumerate(m.split('_')):
        for r in range(2):
            p = blk[3*r:3*r+3]
            T[(s, r)] = None if p[2] not in 'ABCDEFGH' else (int(p[0]), -1 if p[1] == 'L' else 1, ord(p[2]) - 65)
    return T
def runs(tape, lo, hi):
    a, b = lo, hi
    while a <= b and tape[a] == 0: a += 1
    while b >= a and tape[b] == 0: b -= 1
    return [(k, len(list(g))) for k, g in itertools.groupby(tape[i] for i in range(a, b + 1))], a, b
def collect(m, N, maxb=14):
    T = parse(m); tape = collections.defaultdict(int); pos = lo = hi = 0; st = 0; bnd = 0; ev = []
    for t in range(1, N + 1):
        a = T[(st, tape[pos])]
        if a is None: return ev, t
        w = a[0]
        if w != tape[pos]:
            for j in (pos - 1, pos):
                if lo <= j < hi: bnd += 1 if tape[j] == tape[j+1] else -1
            tape[pos] = w
        pos += a[1]; st = a[2]
        if pos < lo: lo = pos; bnd += (tape[lo] != tape[lo+1])
        if pos > hi: hi = pos; bnd += (tape[hi-1] != tape[hi])
        if (pos == lo or pos == hi) and bnd <= maxb:
            rl, a0, b0 = runs(tape, lo, hi)
            ev.append((t, 'ABCDEFGH'[st] + ('L' if pos == lo else 'R'), len(rl), rl, pos - a0, b0 - pos))
    return ev, None
if __name__ == '__main__':
    m = sys.argv[1]; N = int(float(sys.argv[2])); k = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    ev, halt = collect(m, N)
    late = [e for e in ev if e[0] > N // 2]
    by = collections.defaultdict(list)
    for e in late: by[e[1]].append(e)
    best = {key: min(e[2] for e in v) for key, v in by.items()}
    cand = sorted(by, key=lambda key: (best[key], -sum(e[2] == best[key] for e in by[key])))
    print('halt', halt, '| min runs per (state,side):', {key: best[key] for key in cand[:6]})
    key = cand[0]; full = [e for e in ev if e[1] == key and e[2] == best[key]]
    print(f'{key}: {len(full)} configurations with {best[key]} runs; last {k}:')
    for e in full[-k:]:
        print(f'  t={e[0]:>10}  ' + ' '.join(f'{s}^{n}' for s, n in e[3]))
