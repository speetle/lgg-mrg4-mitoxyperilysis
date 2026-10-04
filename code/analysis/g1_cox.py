# -*- coding: utf-8 -*-
import pandas as pd, numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
pd.set_option('display.width',200)

m=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
m=m.dropna(subset=['risk','t','e','age','sex','idh','grade']).copy()
m['high']=(m['组']=='高风险').astype(int)
m['idhwt']=(m['idh']=='IDHwt').astype(int)
m['g3']=(m['grade']=='G3').astype(int)
m['male']=(m['sex'].str[0]=='M').astype(int)
m['codel']=(m['codel']=='codel').astype(int); m['noncodel']=(m['codel']=='non-codel').astype(int)
m['risk_z']=(m['risk']-m['risk'].mean())/m['risk'].std()
print('analysis n =',len(m),'| deaths =',int(m['e'].sum()),'| EPV(3 var)=',round(125/3,1))
print('median follow-up (reverse KM) approx:',round(m.loc[m['e']==0,'t'].median(),1),'months')

def cox(df,cols,label):
    d=df[['t','e']+cols].dropna()
    c=CoxPHFitter().fit(d,'t','e')
    s=c.summary[['coef','exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']]
    print(f'\n--- {label} (n={len(d)}, events={int(d["e"].sum())}) ---')
    print(s.round(4).to_string())
    ci=concordance_index(d['t'],-c.predict_partial_hazard(d),d['e'])
    print('  Harrell C =',round(ci,3),'| AIC=',round(c.AIC_partial_,1),'| log-likelihood ratio p=',f'{c.log_likelihood_ratio_test().p_value:.3e}')
    return c

print('\n================ UNIVARIABLE ================')
for v in ['risk_z','idhwt','g3']:
    cox(m,[v],f'univariable: {v}')

print('\n================ MODEL A: risk + age + sex  (= the paper\'s claim) ================')
cA=cox(m,['risk_z','age','male'],'risk_z + age + sex(male)')

print('\n================ MODEL B: + IDH + GRADE  (G1 GATE) ================')
cB=cox(m,['risk_z','age','male','idhwt','g3'],'risk_z + age + sex + IDHwt + G3')

print('\n================ MODEL C: + 3-level molecular subtype ================')
cC=cox(m,['risk_z','age','male','idhwt','codel'],'risk_z + age + sex + IDHwt + codel(1p/19q)')

print('\n================ MODEL D: risk GROUP (binary) instead of score ================')
cD=cox(m,['high','age','male','idhwt','g3'],'high-risk group + age + sex + IDHwt + G3')

# stratified KM
print('\n================ RISK-STRATIFIED LOG-RANK ================')
for name,sub in [('All',m),('IDHmut',m[m['idh']=='IDHmut']),('IDHwt',m[m['idh']=='IDHwt']),
                 ('Grade2',m[m['grade']=='G2']),('Grade3',m[m['grade']=='G3'])]:
    if sub['组'].nunique()<2 or len(sub)<10: 
        print(f'  {name}: n={len(sub)} — insufficient'); continue
    r=multivariate_logrank_test(sub['t'],sub['组'],sub['e'])
    km_lo=sub[sub['组']=='低风险']['e'].sum()/len(sub[sub['组']=='低风险'])
    km_hi=sub[sub['组']=='高风险']['e'].sum()/len(sub[sub['组']=='高风险'])
    print(f'  {name:7s} n={len(sub):3d}  death-rate low={km_lo:.3f} high={km_hi:.3f}  log-rank p={r.p_value:.3e}')

# risk score distribution by IDH
from scipy.stats import mannwhitneyu
u,p=mannwhitneyu(m[m['idh']=='IDHwt']['risk'],m[m['idh']=='IDHmut']['risk'])
print('\nrisk score IDHwt vs IDHmut: median %.3f vs %.3f, Mann-Whitney p=%.3e'%(m[m['idh']=='IDHwt']['risk'].median(),m[m['idh']=='IDHmut']['risk'].median(),p))
print('\nmedian OS by group (KM) :')
kmf=KaplanMeierFitter()
for g in ['低风险','高风险']:
    s=m[m['组']==g]; kmf.fit(s['t'],s['e']); print(' ',g,'median OS =',kmf.median_survival_time_,'months')
