"""Attach bbchallenge decider labels (bb5_decided_indexes) to screened records."""
import numpy as np, pandas as pd, glob, os, sys
L = os.path.join(os.path.dirname(__file__), '../bbdata/labels')
NAMES = {'cyclers':'cyclers','translated':'translated cyclers','backward':'backward reasoning',
         'halting':'halting segment','bouncers':'bouncers','FAR':'FAR','coq':'Coq-BB5 (holdouts)'}
def load():
    out = {}
    for f in glob.glob(L + '/*-run-*'):
        k = next(v for k, v in NAMES.items() if os.path.basename(f).startswith(k))
        a = np.fromfile(f, dtype='>u4').astype(np.int64)
        out[k] = np.union1d(out.get(k, np.empty(0, np.int64)), a)
    return out
def attach(df, sets):
    for k, a in sets.items(): df[k] = np.isin(df.db_index.values, a)
    return df
if __name__ == '__main__':
    sets = load()
    for k, a in sets.items(): print(f'{k:22s}{len(a):>10,}  in time-limit range {np.sum(a<14322029):>10,}')
    d = pd.read_csv(sys.argv[1])
    d = attach(d, sets)
    d['cls'] = np.select([d.S < .25, d.S <= .75], ['scale-inv', 'intermediate'], 'fades 1/span')
    d.loc[d.S.isna(), 'cls'] = 'unscored'
    cols = list(sets)
    d['nlab'] = d[cols].sum(1)
    print('\nrecords with no label:', (d.nlab == 0).sum(), ' with >1 label:', (d.nlab > 1).sum())
    print('\nfraction of each class carrying each label (%):')
    print((d.groupby('cls')[cols].mean() * 100).round(2).T.to_string())
    print('\nclass counts:'); print(d.cls.value_counts().to_string())
    print('\nfraction scale-invariant among records carrying each label (%):')
    for k in cols:
        s = d[d[k]]; print(f'{k:22s} n={len(s):>7,}  scale-inv {100*(s.cls=="scale-inv").mean():6.2f}  interm {100*(s.cls=="intermediate").mean():5.2f}')
    d.to_csv(sys.argv[2], index=False)
