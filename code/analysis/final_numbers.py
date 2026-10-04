# -*- coding: utf-8 -*-
import json
import pandas as pd, numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
from scipy.stats import spearmanr, mannwhitneyu
out={}

m=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
m['male']=(m['sex'].str[0]=='M').astype(int); m['idhwt']=(m['idh']=='IDHwt').astype(int); m['g3']=(m['grade']=='G3').astype(int)
m['risk_z']=(m['risk']-m['risk'].mean())/m['risk'].std()

def fu(df,ecol='e'):
    kmf=KaplanMeierFitter(); kmf.fit(df['t'],1-df[ecol]); return float(kmf.median_survival_time_)
a=m.dropna(subset=['subtype','grade']).dropna(subset=['risk','t','e','age','male','idhwt','g3'])
out['tcga']={'n_total':int(len(m)),'n_analysis':int(len(a)),'deaths':int(a['e'].sum()),
 'age_median':round(float(a['age'].median()),1),'age_range':[int(a['age'].min()),int(a['age'].max())],
 'male_pct':round(100*a['male'].mean(),1),'grade':{k:int(v) for k,v in a['grade'].value_counts().items()},
 'subtype':{k:int(v) for k,v in a['subtype'].value_counts().items()},
 'idh_mut':int((a['idhwt']==0).sum()),'idh_wt':int((a['idhwt']==1).sum()),
 'median_fu_months':round(fu(a),1),'median_os_low':round(float(KaplanMeierFitter().fit(a[a['组']=='低风险']['t'],a[a['组']=='低风险']['e']).median_survival_time_),1),
 'median_os_high':round(float(KaplanMeierFitter().fit(a[a['组']=='高风险']['t'],a[a['组']=='高风险']['e']).median_survival_time_),1)}
print('TCGA baseline:',json.dumps(out['tcga'],ensure_ascii=False,default=str))

print('\n--- TCGA stratified log-rank ---')
for nm,sub in [('IDHmut',a[a['idhwt']==0]),('IDHwt',a[a['idhwt']==1]),('G2',a[a['g3']==0]),('G3',a[a['g3']==1])]:
    r=multivariate_logrank_test(sub['t'],sub['组'],sub['e'])
    print(f'  {nm}: n={len(sub)} deaths={int(sub["e"].sum())} p={r.p_value:.3e}')

# CGGA
c=pd.read_csv('cgga_dataset.csv',encoding='utf-8-sig')
cl=c[c['Grade'].isin(['WHO II','WHO III'])].copy()
print('\n--- CGGA LGG stratified log-rank ---')
for nm,sub in [('IDHmut',cl[cl['idhwt']==0]),('IDHwt',cl[cl['idhwt']==1]),('WHO II',cl[cl['grade_hi']==0]),('WHO III',cl[cl['grade_hi']==1])]:
    if sub['ev'].sum()<3 or len(sub)<10: print(f'  {nm}: n={len(sub)} too small'); continue
    s=sub.copy(); s['grp']=np.where(s['risk']>cl['risk'].median(),'high','low')
    if s['grp'].nunique()<2: print(f'  {nm}: no split'); continue
    r=multivariate_logrank_test(s['t'],s['grp'],s['ev'])
    print(f'  {nm}: n={len(sub)} deaths={int(sub["ev"].sum())} p={r.p_value:.3e}')
print('\nCGGA: median risk IDHmut vs IDHwt (LGG): %.3f vs %.3f'%(cl[cl['idhwt']==0]['risk'].median(),cl[cl['idhwt']==1]['risk'].median()))
print('CGGA LGG median follow-up: %.1f mo'%fu(cl,'ev'))
print('CGGA all: n=%d deaths=%d median fu=%.1f'%(len(c),int(c['ev'].sum()),fu(c,'ev')))

# TCGA: multivariable model B recompute for exact values
d=a[['t','e','risk_z','age','male','idhwt','g3']].dropna()
cB=CoxPHFitter().fit(d,'t','e')
print('\n--- TCGA Model B exact ---')
print(cB.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4).to_string())
print('C=',round(concordance_index(d['t'],-cB.predict_partial_hazard(d),d['e']),3))
out['tcga_modelB']={k:{'HR':round(float(v['exp(coef)']),3),'lo':round(float(v['exp(coef) lower 95%']),3),'hi':round(float(v['exp(coef) upper 95%']),3),'p':float(v['p'])} for k,v in cB.summary.iterrows()}
json.dump(out,open('results_digest.json','w'),ensure_ascii=False,indent=1)
print('\nsaved results_digest.json')
