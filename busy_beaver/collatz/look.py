# look.py <machine> <steps> <n>: print the last n tapes (edge zeros trimmed) at visits of the head to
# an edge, compressed: any unit u of length 1-6 repeated 3+ times is written (u)^k.
import sys, re, collections
from macro import parse
def compress(s):
    return re.sub(r'([01]{1,6}?)\1{2,}', lambda mm: f'({mm.group(1)})^{len(mm.group(0))//len(mm.group(1))} ', s)
m = sys.argv[1]; N = int(float(sys.argv[2])); n = int(sys.argv[3]); T = parse(m)
tape = collections.defaultdict(int); pos = lo = hi = 0; st = 0; ring = collections.deque(maxlen=n)
for t in range(1, N + 1):
    a = T[(st, tape[pos])]
    if a is None: print('halt', t); break
    tape[pos] = a[0]; pos += a[1]; st = a[2]; lo = min(lo, pos); hi = max(hi, pos)
    if (pos == lo or pos == hi) and t > N - 400000:
        s = ''.join(str(tape[i]) for i in range(lo, hi + 1)).strip('0')
        ring.append(f'{t} {"ABCDEF"[st]}{"L" if pos == lo else "R"} {compress(s)}')
for l in ring: print(l[:230])
