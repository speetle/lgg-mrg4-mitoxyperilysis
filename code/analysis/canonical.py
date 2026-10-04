# -*- coding: utf-8 -*-
import json
import pandas as pd, numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
pd.set_option('display.width',220)
m=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
m['idhwt']=(m['idh']=='IDHwt').astype(int); m['g3']=(m['grade']=='G3').astype(int); m['male']=(m['sex'].str[0]=='M').astype(int)
m['risk_z']=(m['risk']-m['risk'].mean())/m['risk'].std()
m['high']=(m['组']=='高风险').astype(int)
def cox(df,cols,lab):
    d=df[['t','e']+cols].dropna(); c=CoxPHFitter().fit(d,'t','e')
    print(f'\n[{lab}] n={len(d)} deaths={int(d["e"].sum())} C={concordance_index(d["t"],-c.predict_partial_hazard(d),d["e"]):.3f}')
    print(c.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4).to_string())
    return c
print('===== FULL COHORT (n=510) =====')
cox(m,['risk_z'],'unadjusted: risk score per SD')
cox(m,['risk_z','age','male'],'adjusted (age+sex only)')
cox(m,['high'],'BINARY risk group (unadjusted)')
r=multivariate_logrank_test(m['t'],m['组'],m['e'])
print('\nmedian-split log-rank P = %.3e'%r.p_value)
k=KaplanMeierFitter()
for g in ['低风险','高风险']:
    s=m[m['组']==g]; k.fit(s['t'],s['e']); print(f'  {g}: n={len(s)} deaths={int(s["e"].sum())} median OS={k.median_survival_time_:.1f} mo')
k.fit(m['t'],1-m['e']); print('median follow-up (reverse KM) = %.1f mo'%k.median_survival_time_)
print('\n===== COMPLETE ANNOTATION SUBSET (n=502) =====')
a=m.dropna(subset=['subtype','grade']).copy()
cox(a,['risk_z','age','male','idhwt','g3'],'adjusted (age+sex+IDH+grade), continuous')
d2=a.copy()
cox(d2,['high','age','male','idhwt','g3'],'adjusted (age+sex+IDH+grade), binary')
cox(a,['idhwt'],'IDH-wildtype alone')
cox(a,['g3'],'WHO grade 3 alone')
out={'full_n':510,'full_deaths':125,'adj_n':502,'adj_deaths':int(a['e'].sum())}
json.dump(out,open('canonical.json','w')); print('\nsaved')
