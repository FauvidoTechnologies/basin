"""Compare basin's roughness (cv_rel_mse) with the scale-invariance index S on the BB(5) census's
time-limit sqrt(t) machines. Writes census_sqrt_S.csv and census_flagged.txt (for macroc)."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from scalestat import stat
cen = pd.read_csv('../results/bb5_census.csv')
m = [l.split()[-1] for l in open('census_sqrt.txt')]
sp = np.fromfile('census_sqrt.bin', dtype=np.int32).reshape(len(m), -1)
d = pd.DataFrame([stat(s) for s in sp], columns=['p', 'D', 'S']); d['machine'] = m
d = d.merge(cen[(cen.group == 'time') & (cen.growth_class == '1/2')][['machine', 'cv_rel_mse', 'law']], on='machine')
d['cls'] = np.select([d.S < .25, d.S <= .75], ['scale-inv', 'interm'], 'fades'); d.loc[d.S.isna(), 'cls'] = 'unscored'
d['cv_pct'] = d.cv_rel_mse.rank(pct=True) * 100
print(d.cls.value_counts().to_dict())
print(d.groupby('cls').cv_rel_mse.describe()[['count', '25%', '50%', '75%']].to_string())
a, b = d[d.cls == 'scale-inv'].cv_rel_mse, d[d.cls == 'fades'].cv_rel_mse
if len(a): print('AUC(rough | scale-inv vs fades) = %.3f  p = %.1g' % (mannwhitneyu(a, b).statistic / len(a) / len(b), mannwhitneyu(a, b).pvalue))
print('spearman(S, cv) = %.2f' % spearmanr(d.S, d.cv_rel_mse, nan_policy='omit').statistic)
top = d.sort_values('cv_rel_mse', ascending=False)
print('\nroughest 15:'); print(top.head(15)[['machine', 'S', 'cls', 'cv_rel_mse', 'cv_pct']].to_string(index=False))
print('\nscale-inv machines:'); print(d[d.cls != 'fades'].sort_values('S')[['machine', 'S', 'cls', 'cv_rel_mse', 'cv_pct']].to_string(index=False))
d.to_csv('census_sqrt_S.csv', index=False)
fl = pd.concat([d[d.cls != 'fades'], top.head(15)]).drop_duplicates('machine')
with open('census_flagged.txt', 'w') as f:
    for i, r in fl.iterrows(): f.write(f'{r.cls}_{i} {r.machine}\n')
