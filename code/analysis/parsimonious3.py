# -*- coding: utf-8 -*-
"""简约模型 v3：LASSO CV曲线 + 稳定性选择核心模型 + CGGA冻结系数外部验证"""
import json, io
import numpy as np, pandas as pd
from sksurv.util import Surv
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.metrics import concordance_index_censored
from sklearn.model_selection import KFold
from lifelines import CoxPHFitter
np.random.seed(42); pd.set_option('display.width',210)
C=lambda t,e,r: concordance_index_censored(np.asarray(e).astype(bool),np.asarray(t),np.asarray(r))[0]

E=pd.read_csv('panel_expr_tcga.csv',index_col=0,encoding='utf-8-sig')
G=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
d=G.merge(E,left_on='pat',right_index=True,how='inner').rename(columns={'e':'ev','t':'time'})
genes=list(E.columns); T=d['time'].values; Ev=d['ev'].values
X=d[genes].values.astype(float); y=Surv.from_arrays(Ev.astype(bool),T)
mu,sd=X.mean(0),X.std(0)+1e-9; Xs=(X-mu)/sd
def fit_at(Xa,ya,a,ratio=1.0):
    return CoxnetSurvivalAnalysis(l1_ratio=ratio,alphas=[a],max_iter=100000).fit(Xa,ya).coef_[:,0]
alphas=np.array(CoxnetSurvivalAnalysis(l1_ratio=1.0,alpha_min_ratio=0.01,max_iter=100000).fit(Xs,y).alphas_)
print('[TCGA] n=%d genes=%d events=%d | 25-gene apparent C=%.3f'%(len(d),len(genes),int(Ev.sum()),C(T,Ev,d['risk'].values)))

# ---- CV 曲线 (10-fold) ----
kf=KFold(n_splits=10,shuffle=True,random_state=42); S=np.zeros((10,len(alphas)))
for i,(tr,te) in enumerate(kf.split(Xs)):
    m=CoxnetSurvivalAnalysis(l1_ratio=1.0,alphas=alphas,max_iter=100000).fit(Xs[tr],y[tr])
    for j,a in enumerate(alphas):
        try: r=m.predict(Xs[te],alpha=a)
        except Exception: r=np.zeros(len(te))
        S[i,j]=C(T[te],Ev[te],r)
cvmean=S.mean(0); cvse=S.std(0,ddof=1)/np.sqrt(10)
ks=np.array([int((np.abs(fit_at(Xs,y,a))>1e-8).sum()) for a in alphas])
b=int(cvmean.argmax()); thr=cvmean[b]-cvse[b]; i1=int(np.where(cvmean>=thr)[0].min())
print('\n=== CV 曲线选点 ===')
print('lambda.min: alpha=%.5f k=%2d CVC=%.3f'%(alphas[b],ks[b],cvmean[b]))
print('lambda.1se: alpha=%.5f k=%2d CVC=%.3f  (1se阈值=%.3f)'%(alphas[i1],ks[i1],cvmean[i1],thr))
# 稀疏度对照表
print('\n稀疏度对照（alpha -> k, CV C）:')
seen=set()
for j in range(len(alphas)):
    if ks[j] not in seen and ks[j]<=30:
        seen.add(ks[j]); print('   k=%2d alpha=%.5f CVC=%.3f+-%.3f'%(ks[j],alphas[j],cvmean[j],cvse[j]))

# ---- bootstrap 稳定性选择 ----
B=200; cnt={g:0 for g in genes}; ks_b=[]
for _ in range(B):
    idx=np.random.choice(len(Xs),len(Xs),replace=True)
    if Ev[idx].sum()<10: continue
    try:
        c=fit_at(Xs[idx],y[idx],alphas[i1]); nz=np.where(np.abs(c)>1e-8)[0]
        ks_b.append(len(nz))
        for j in nz: cnt[genes[j]]+=1
    except Exception: pass
stab=sorted([(v/B,g) for g,v in cnt.items() if v>0],reverse=True)
print('\n[Bootstrap %d] 中位入选数=%d'%(len(ks_b),np.median(ks_b)))
print('选择频率 top15:',[('%s %.0f%%'%(g,100*v)) for v,g in stab[:15]])

# ---- 稳定性选择核心模型 (pi>=0.75) ----
core=[g for v,g in stab if v>=0.75]
print('\n[核心模型] pi>=0.75 -> %d 基因: %s'%(len(core),core))
Zc_tc=(d[core].values-mu[[genes.index(g) for g in core]])/sd[[genes.index(g) for g in core]]
tc=d.copy()
for i,g in enumerate(core): tc['z_'+g]=Zc_tc[:,i]
mod=CoxPHFitter().fit(tc[['time','ev']+['z_'+g for g in core]],'time','ev')
core_coef={g:float(mod.params_['z_'+g]) for g in core}
print(mod.summary[['coef','exp(coef)','p']].round(4).to_string())
tc['lp_core']=sum(tc['z_'+g]*core_coef[g] for g in core)
core_C_app=C(tc['time'],tc['ev'],tc['lp_core'])
# 核心模型诚实 CV：在每折内重拟合核心模型
hon=[]
for tr,te in KFold(10,shuffle=True,random_state=7).split(tc):
    f=CoxPHFitter().fit(tc.iloc[tr][['time','ev']+['z_'+g for g in core]],'time','ev')
    hon.append(C(tc.iloc[te]['time'],tc.iloc[te]['ev'],f.predict_partial_hazard(tc.iloc[te])))
hon=np.array(hon)
print('核心模型 apparent C=%.3f | 10折CV C=%.3f+-%.3f | EPV=%.1f'%(core_C_app,hon.mean(),hon.std(ddof=1)/np.sqrt(10),Ev.sum()/len(core)))

# ---- CGGA 外部验证 ----
cexpr=pd.read_csv('cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt',sep='\t',index_col=0)
_raw=open('cgga/CGGA.mRNAseq_325_clinical.20200506.txt',encoding='utf-8-sig',errors='ignore').read().replace('\r','\n')
clin=pd.read_csv(io.StringIO(_raw),sep='\t'); clin.columns=[c.strip() for c in clin.columns]
clin=clin.set_index('CGGA_ID')
clin['ev']=pd.to_numeric(clin['Censor (alive=0; dead=1)'],errors='coerce')
clin['mos']=pd.to_numeric(clin['OS'],errors='coerce')/30.4375
clin['age']=pd.to_numeric(clin['Age'],errors='coerce')
clin['male']=(clin['Gender'].astype(str).str[0]=='M').astype(float)
clin['idhwt']=(clin['IDH_mutation_status'].astype(str).str.lower().str.startswith('wild')).astype(float)
clin['grade_hi']=clin['Grade'].isin(['WHO III','WHO IV']).astype(float)
print('\nCGGA clinical n=%d (Grade: %s)'%(len(clin),dict(clin['Grade'].value_counts())))
# 两种标准化：A=TCGA冻结参数(严格冻结)  B=队列内z化(尺度自适应)
out={}
for tag,standardize in [('frozen_TCGA',True),('cohort_z',False)]:
    def score(genelist,coefmap):
        gg=[g for g in genelist if g in cexpr.index]
        M=cexpr.loc[gg].T; M=M[~M.index.duplicated()]
        if standardize:
            Z=(M.values-mu[[genes.index(g) for g in gg]])/sd[[genes.index(g) for g in gg]]
        else:
            Z=(M.values-M.values.mean(0))/(M.values.std(0)+1e-9)
        return pd.Series(Z@np.array([coefmap[g] for g in gg]),index=M.index), gg
    lp25,_=score(list(dict.fromkeys([r for r in json.load(open('gene_ids.json'))['model']])),
                 {g:float(v) for g,v in zip([r['基因'] for r in __import__('csv').DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))],
                   [float(r['系数']) for r in __import__('csv').DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))])})
    lpc,gc=score(core,{**core_coef, **{'z_'+k:v for k,v in core_coef.items()}})
    lpc,_=score(core,core_coef)
    lpl1,_=score([g for g in genes if g in cexpr.index], {g:fit_at(Xs,y,alphas[i1])[genes.index(g)] for g in genes if g in cexpr.index})
    T_=clin['mos']; Evv=clin['ev']
    sub_lgg=clin['Grade'].isin(['WHO II','WHO III'])
    row={}
    for nm,s in [('25gene',lp25),('LASSO_1se',lpl1),('core',lpc)]:
        v=pd.concat([clin,s.rename('lp')],axis=1).dropna(subset=['lp','mos','ev'])
        vL=v[v['Grade'].isin(['WHO II','WHO III'])]
        row[nm]={'C_all':round(C(v['mos'],v['ev'],v['lp']),3),'C_lgg':round(C(vL['mos'],vL['ev'],vL['lp']),3),
                 'n_lgg':int(len(vL)),'deaths_lgg':int(vL['ev'].sum())}
        if nm=='core':
            vv=vL.copy(); vv['lp_z']=(vv['lp']-vL['lp'].mean())/vL['lp'].std()
            f=CoxPHFitter().fit(vv[['mos','ev','lp_z','age','male','idhwt','grade_hi']].dropna(),'mos','ev')
            row['core_adj']={k:{'HR':round(float(r['exp(coef)']),3),'lo':round(float(r['exp(coef) lower 95%']),3),
                                'hi':round(float(r['exp(coef) upper 95%']),3),'p':float(r['p'])} for k,r in f.summary.iterrows()}
    out[tag]=row
    print('\n=== CGGA 外部验证 [标准化=%s] ==='%tag)
    for nm in ['25gene','LASSO_1se','core']:
        print('  %-9s LGG n=%d deaths=%d  C(LGG)=%.3f  C(全分级)=%.3f'%(nm,row[nm]['n_lgg'],row[nm]['deaths_lgg'],row[nm]['C_lgg'],row[nm]['C_all']))
    print('  core 校正后:'); [print('     %-10s HR=%.3f (%.3f-%.3f) p=%.4g'%(k,v['HR'],v['lo'],v['hi'],v['p'])) for k,v in row['core_adj'].items()]

res={'cv':{'alphas':alphas.tolist(),'mean':cvmean.tolist(),'se':cvse.tolist(),'k':ks.tolist(),
           'i_min':b,'i_1se':i1,'alpha_min':float(alphas[b]),'alpha_1se':float(alphas[i1]),
           'k_min':int(ks[b]),'k_1se':int(ks[i1]),'cvC_min':float(cvmean[b]),'cvC_1se':float(cvmean[i1])},
     'lasso_1se':{'alpha':float(alphas[i1]),'genes':genes,'coef':{g:float(fit_at(Xs,y,alphas[i1])[genes.index(g)]) for g in genes}},
     'stability':{g:round(v,3) for v,g in stab},'boot_k_median':float(np.median(ks_b)),
     'core':{'genes':core,'coef':core_coef,'pi_threshold':0.75,'apparent_C':float(core_C_app),
             'cv_C':float(hon.mean()),'cv_C_se':float(hon.std(ddof=1)/np.sqrt(10)),'EPV':round(Ev.sum()/len(core),1)},
     'external':out,'epv':{'events':int(Ev.sum()),'k_25':25,'epv_25':round(Ev.sum()/25,1),
             'k_1se':int(ks[i1]),'epv_1se':round(Ev.sum()/ks[i1],1)},
     'corr_with_25gene':float(np.corrcoef(tc['lp_core'],tc['risk'])[0,1])}
json.dump(res,open('parsimonious_results.json','w'),ensure_ascii=False,indent=1)
print('\n核心评分 vs 25基因评分 r=%.3f'%res['corr_with_25gene'])
print('EPV: 25基因=%.1f | LASSO1se(%d)=%.1f | 核心(%d)=%.1f'%(res['epv']['epv_25'],ks[i1],res['epv']['epv_1se'],len(core),res['core']['EPV']))
print('\nsaved parsimonious_results.json')
