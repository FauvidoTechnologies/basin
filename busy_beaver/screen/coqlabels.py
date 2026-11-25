"""Label screened time-limit records with Coq-BB5 decider IDs.
usage: coqlabels.py <indices.int64> <out.csv>   (streams BB5_verified_enumeration.csv.zst)"""
import numpy as np, subprocess, sys, os
H = os.path.dirname(os.path.abspath(__file__))
def tm_strings(idx):
    rec = np.memmap(H + '/../bbdata/timelimit.bin', dtype=np.uint8, mode='r').reshape(-1, 30)[np.sort(idx)]
    out = []
    for r in rec:
        t = ['---' if r[3*i+2] == 0 else f"{r[3*i]}{'RL'[r[3*i+1]]}{chr(64+r[3*i+2])}" for i in range(10)]
        out.append('_'.join(a + b for a, b in zip(t[::2], t[1::2])))
    return np.sort(idx), out
if __name__ == '__main__':
    idx, s = tm_strings(np.fromfile(sys.argv[1], dtype=np.int64))
    want = dict(zip(s, idx)); print(len(want), 'strings, e.g.', s[0], flush=True)
    p = subprocess.Popen(['zstdcat', H + '/../bbdata/labels/BB5_verified_enumeration.csv.zst'],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    found = {}; n = 0
    for line in p.stdout:
        n += 1
        m = line[:34]
        if m in want: found[m] = line.rstrip().split(',', 2)[1:]
    print(n, 'enumerated machines;', len(found), 'of', len(want), 'matched')
    with open(sys.argv[2], 'w') as f:
        f.write('db_index,machine,status,decider\n')
        for m, i in want.items():
            st, de = found.get(m, ('', ''))
            f.write(f'{i},{m},{st},{de}\n')
