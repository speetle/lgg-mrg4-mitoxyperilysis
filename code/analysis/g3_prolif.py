# -*- coding: utf-8 -*-
import json,csv,urllib.request
import pandas as pd, numpy as np
from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from scipy.stats import spearmanr,pearsonr
B='https://www.cbioportal.org/api'
def post(url,payload):
    req=urllib.request.Request(url,data=json.dumps(payload).encode(),headers={'Content-Type':'application/json','Accept':'application/json'})
    return json.load(urllib.request.urlopen(req,timeout=120))

G=json.load(open('gene_ids.json'))
entrez=[v for v in G['model_entrez'].values()]+[v for v in G['prolif_entrez'].values()]
prof='lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores'
sl='lgg_tcga_pan_can_atlas_2018_rna_seq_v2_mrna'
data=post(f'{B}/molecular-profiles/{prof}/molecular-data/fetch?projection=SUMMARY',{'entrezGeneIds':entrez,'sampleListId':sl})
print('records:',len(data),'| samples:',len(set(d['sampleId'] for d in data)))
inv={v:k for k,v in {**G['model_entrez'],**G['prolif_entrez']}.items()}
df=pd.DataFrame([{'sampleId':d['sampleId'],'patientId':d['patientId'],'gene':inv.get(d['entrezGeneId']),'z':d['value']} for d in data])
piv=df.pivot_table(index=['patientId'],columns='gene',values='z')
print('matrix:',piv.shape)

coef={r['基因']:float(r['系数']) for r in csv.DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))}
present=[g for g in coef if g in piv.columns]
print('model genes present in expr matrix:',len(present),'/',len(coef))
Zs=piv[present].apply(lambda c:(c-c.mean())/c.std())
risk_hat=(Zs*[coef[g] for g in present]).sum(axis=1)

pl=[g for g in G['prolif'] if g in piv.columns]
Zp=piv[pl].apply(lambda c:(c-c.mean())/c.std())
prolif=Zp.mean(axis=1)
print('proliferation markers used:',len(pl))

R=pd.DataFrame({'patientId':piv.index,'risk_hat':risk_hat.values,'prolif':prolif.values})
# join survival via sample->patient using clinical
raw=json.load(open('lgg_clinical_raw.json'))
byk={}
for r in raw:
    byk.setdefault(r['uniquePatientKey'],{'pid':r['patientId']})[r['clinicalAttributeId']]=r['value']
C=pd.DataFrame([dict(pid=v.get('pid'),t=float(v['OS_MONTHS']) if v.get('OS_MONTHS') else np.nan,
                     e=1 if v.get('OS_STATUS','').startswith('1') else (0 if v.get('OS_STATUS','').startswith('0') else np.nan),
                     age=float(v['AGE']) if v.get('AGE') else np.nan, sex=v.get('SEX'),
                     subtype=v.get('SUBTYPE')) for k,v in byk.items()])
R=R.merge(C,left_on='patientId',right_on='pid',how='left')
R['idhwt']=(R['subtype']=='LGG_IDHwt').astype(float); R.loc[R['subtype'].isna(),'idhwt']=np.nan
sg=pd.read_json('lgg_sample_clin.json'); g=sg[sg['clinicalAttributeId']=='GRADE'][['patientId','value']].rename(columns={'value':'grade'}).drop_duplicates('patientId')
R=R.merge(g,on='patientId',how='left'); R['g3']=(R['grade']=='G3').astype(float); R.loc[R['grade'].isna(),'g3']=np.nan
R['male']=(R['sex'].str[0]=='M').astype(float)

loc=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')[['pid','risk','组']]
R=R.merge(loc,left_on='patientId',right_on='pid',how='left')
ok=R.dropna(subset=['risk','risk_hat'])
print('\n=== REPRODUCTION of the published risk score from frozen coefficients ===')
print(' n:',len(ok),'| Pearson r =',round(pearsonr(ok['risk'],ok['risk_hat'])[0],4),'| Spearman rho =',round(spearmanr(ok['risk'],ok['risk_hat'])[0],4))

d=R.dropna(subset=['risk_hat','prolif','t','e','age','male','idhwt','g3'])
print('\n=== IS IT JUST PROLIFERATION? (n=%d, events=%d) ==='%(len(d),int(d['e'].sum())))
print(' Spearman risk_hat vs proliferation score: rho=%.3f p=%.2e'%spearmanr(d['risk_hat'],d['prolif']))
def cox(df,cols):
    c=CoxPHFitter().fit(df[['t','e']+cols],'t','e')
    s=c.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4)
    ci=concordance_index(df['t'],-c.predict_partial_hazard(df),df['e'])
    print(f'\n -- {cols} | C={ci:.3f}\n{s.to_string()}')
for cols in [['prolif'],['risk_hat'],['risk_hat','prolif'],['risk_hat','prolif','idhwt','g3','age','male']]:
    cox(d,cols)
d.to_csv('g3_dataset.csv',index=False,encoding='utf-8-sig')
