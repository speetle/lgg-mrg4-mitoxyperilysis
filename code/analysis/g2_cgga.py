# -*- coding: utf-8 -*-
import io,csv
import pandas as pd, numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from lifelines.utils import concordance_index
from scipy.stats import spearmanr, mannwhitneyu

D='cgga/'
txt=open(D+'CGGA.mRNAseq_325_clinical.20200506.txt',encoding='utf-8-sig',errors='ignore').read().replace('\r','\n')
clin=pd.read_csv(io.StringIO(txt),sep='\t')
clin.columns=[c.strip() for c in clin.columns]
print('clinical:',clin.shape,'| cols:',list(clin.columns))
print('Grade:',dict(clin['Grade'].value_counts(dropna=False)))
print('IDH:',dict(clin['IDH_mutation_status'].value_counts(dropna=False)))

expr=pd.read_csv(D+'CGGA.mRNAseq_325.RSEM-genes.20200506.txt',sep='\t',index_col=0)
print('expr:',expr.shape)
coef={r['基因']:float(r['系数']) for r in csv.DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))}
present=[g for g in coef if g in expr.index]
print('model genes present in CGGA:',len(present),'/',len(coef),'| missing:',[g for g in coef if g not in expr.index])

E=expr.loc[present].T
E=E[~E.index.duplicated()]
Zs=E.apply(lambda c:(c-c.mean())/c.std())
risks=(Zs*pd.Series(coef)[present]).sum(axis=1)

C=clin.set_index('CGGA_ID')
d=pd.DataFrame({'risk':risks})
d=d.join(C[['Grade','IDH_mutation_status','Gender','Age','OS','Censor (alive=0; dead=1)','Histology']],how='inner')
d=d.rename(columns={'Censor (alive=0; dead=1)':'ev','IDH_mutation_status':'idh','Gender':'sex','Age':'age'})
d['ev']=pd.to_numeric(d['ev'],errors='coerce')
d['t']=pd.to_numeric(d['OS'],errors='coerce')/30.4375
d['age']=pd.to_numeric(d['age'],errors='coerce')
d=d.dropna(subset=['risk','t','ev'])
print('\nmerged analysable:',len(d),'| events:',int(d['ev'].sum()))
print('Grade dist:',dict(d['Grade'].value_counts()))
print('IDH dist:',dict(d['idh'].value_counts(dropna=False)))
d['idhwt']=(d['idh'].astype(str).str.lower().str.startswith('wild')).astype(float)
d['male']=(d['sex'].astype(str).str[0]=='M').astype(float)
d['grade_hi']=d['Grade'].isin(['WHO III','WHO IV']).astype(float)

def cox(df,cols,label):
    dd=df[['t','ev']+cols].dropna()
    c=CoxPHFitter().fit(dd,'t','ev')
    s=c.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4)
    ci=concordance_index(dd['t'],-c.predict_partial_hazard(dd),dd['ev'])
    print(f'\n[{label}] n={len(dd)} events={int(dd["ev"].sum())} C={ci:.3f}')
    print(s.to_string())
    return c

print('\n########## CGGA-325 LGG ONLY (WHO II + III) ##########')
lgg=d[d['Grade'].isin(['WHO II','WHO III'])].copy()
print('LGG n=',len(lgg),'events=',int(lgg['ev'].sum()))
print('risk vs IDH (LGG):', 'median %.3f vs %.3f'%(lgg[lgg['idhwt']==1]['risk'].median(),lgg[lgg['idhwt']==0]['risk'].median()),
      'p=%.2e'%mannwhitneyu(lgg[lgg['idhwt']==1]['risk'],lgg[lgg['idhwt']==0]['risk'])[1])
cox(lgg,['risk'],'LGG: risk score alone')
cox(lgg,['idhwt'],'LGG: IDHwt alone')
cox(lgg,['risk','age','male','idhwt','grade_hi'],'LGG: risk + age + sex + IDH + grade')

print('\n########## CGGA-325 ALL GRADES (sensitivity) ##########')
cox(d,['risk'],'ALL: risk score alone')
cox(d,['risk','age','male','idhwt','grade_hi'],'ALL: risk + age + sex + IDH + grade')

# median split KM
for name,sub in [('LGG',lgg),('ALL',d)]:
    s=sub.copy(); s['grp']=np.where(s['risk']>s['risk'].median(),'high','low')
    r=multivariate_logrank_test(s['t'],s['grp'],s['ev'])
    kmf=KaplanMeierFitter(); med={}
    for g in ['low','high']:
        x=s[s['grp']==g]; kmf.fit(x['t'],x['ev']); med[g]=kmf.median_survival_time_
    print(f'\n{name} median-split log-rank p={r.p_value:.3e} | median OS low={med["low"]:.1f}mo high={med["high"]:.1f}mo')
d.to_csv('cgga_dataset.csv',index=False,encoding='utf-8-sig')
