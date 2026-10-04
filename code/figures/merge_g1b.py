# -*- coding: utf-8 -*-
import json
import pandas as pd, numpy as np
base='/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/'
risk=pd.read_csv(base+'风险评分.csv',encoding='utf-8-sig')

raw=json.load(open('lgg_clinical_raw.json'))
byk={}
for r in raw:
    byk.setdefault(r['uniquePatientKey'],{'pid':r['patientId']})[r['clinicalAttributeId']]=r['value']
C=pd.DataFrame([dict(upk=k,pid=v.get('pid'),subtype=v.get('SUBTYPE'),
                     os_months=float(v['OS_MONTHS']) if v.get('OS_MONTHS') else np.nan,
                     os_status=v.get('OS_STATUS'), age_cb=int(v['AGE']) if v.get('AGE') else np.nan,
                     sex_cb=v.get('SEX'), rad=v.get('RADIATION_THERAPY')) for k,v in byk.items()])
sg=pd.read_json('lgg_sample_clin.json')
g=sg[sg['clinicalAttributeId']=='GRADE'][['uniquePatientKey','value']].rename(columns={'uniquePatientKey':'upk','value':'grade'}).drop_duplicates('upk')
C=C.merge(g,on='upk',how='left')

m=risk.merge(C,left_on='pat',right_on='upk',how='left')
print('merged rows:',len(m),'| matched:',int(m['pid'].notna().sum()),'| unmatched:',int(m['pid'].isna().sum()))
chk=m.dropna(subset=['os_months']).copy(); chk['d']=chk['t']-chk['os_months']
print('\n=== UNIT CHECK local t vs OS_MONTHS ===')
print(' n:',len(chk),'| mean diff:',round(chk['d'].mean(),4),'| max abs:',round(chk['d'].abs().max(),4))
print(' within 0.1 mo:',int((chk['d'].abs()<0.1).sum()),'/',len(chk))
print(' event agreement:',int((chk['e']==chk['os_status'].str[0].astype(int)).sum()),'/',len(chk))
print('\nage local vs cbio equal:',int((m['age']==m['age_cb']).sum()),'/',int(m['age_cb'].notna().sum()))
print('sex local vs cbio equal:',int((m['sex'].str[0]==m['sex_cb'].str[0]).sum()),'/',int(m['sex_cb'].notna().sum()))
print('\nsubtype x risk group:')
print(pd.crosstab(m['subtype'].fillna('NA'),m['组']))
m['idh']=m['subtype'].map(lambda s:'IDHmut' if isinstance(s,str) and 'IDHmut' in s else ('IDHwt' if s=='LGG_IDHwt' else None))
m['codel']=m['subtype'].map(lambda s:'codel' if s=='LGG_IDHmut-codel' else ('non-codel' if s=='LGG_IDHmut-non-codel' else ('wt' if s=='LGG_IDHwt' else None)))
print('\nIDH x risk group:'); print(pd.crosstab(m['idh'].fillna('NA'),m['组']))
print('grade x risk group:'); print(pd.crosstab(m['grade'].fillna('NA'),m['组']))
from scipy.stats import chi2_contingency
for v in ['idh','grade']:
    t=pd.crosstab(m[v],m['组'])
    chi2,p,dof,_=chi2_contingency(t); print(f'chi2 {v} vs risk group: chi2={chi2:.1f} p={p:.3e}')
m.to_csv('lgg_g1_dataset.csv',index=False,encoding='utf-8-sig')
print('\nsaved')
