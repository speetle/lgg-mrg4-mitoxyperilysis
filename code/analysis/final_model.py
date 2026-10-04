# -*- coding: utf-8 -*-
"""最终 MRG-4 模型：锁定系数 + 内外验证 + 出图 12/13/14"""
import json, io, csv
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from sksurv.metrics import concordance_index_censored
from sklearn.model_selection import KFold
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.9,'figure.dpi':300})
CB='#1f5fa8'; CR='#c0392b'; CG='#666666'
C=lambda t,e,r: concordance_index_censored(np.asarray(e).astype(bool),np.asarray(t),np.asarray(r))[0]

# ---------- 数据 ----------
E=pd.read_csv('panel_expr_tcga.csv',index_col=0,encoding='utf-8-sig')
G=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
d=G.merge(E,left_on='pat',right_index=True,how='inner').rename(columns={'e':'ev','t':'time'})
genes=list(E.columns); mu,sd=E.values.mean(0),E.values.std(0)+1e-9
R=json.load(open('parsimonious_results.json'))
CORE=R['core']['genes']                      # 稳定性选择 pi>=0.75 (5基因)
MRG4=[g for g in CORE if g!='OSGIN2']        # 去非显著成员 -> 主模型
print('MRG-4 主模型:',MRG4,'| 5基因敏感性:',CORE)

# ---------- 拟合 MRG-4 ----------
def zmat(df,glist):
    return (df[glist].values-mu[[genes.index(g) for g in glist]])/sd[[genes.index(g) for g in glist]]
T=d.copy()
for i,g in enumerate(MRG4): T['z_'+g]=zmat(d,MRG4)[:,i]
f4=CoxPHFitter().fit(T[['time','ev']+['z_'+g for g in MRG4]],'time','ev')
coef4={g:float(f4.params_['z_'+g]) for g in MRG4}
T['mrg4']=sum(T['z_'+g]*coef4[g] for g in MRG4)
print('\nMRG-4 系数（TCGA 拟合，冻结）:'); print(f4.summary[['coef','exp(coef)','p']].round(4).to_string())
app4=C(T['time'],T['ev'],T['mrg4'])
print('TCGA apparent C=%.3f | EPV=%.1f'%(app4,125/len(MRG4)))
# 内部 CV（每折重拟合系数）
cv4=[]
for tr,te in KFold(10,shuffle=True,random_state=7).split(T):
    ff=CoxPHFitter().fit(T.iloc[tr][['time','ev']+['z_'+g for g in MRG4]],'time','ev')
    cv4.append(C(T.iloc[te]['time'],T.iloc[te]['ev'],ff.predict_partial_hazard(T.iloc[te])))
cv4=np.array(cv4); print('MRG-4 10折CV C=%.3f±%.3f'%(cv4.mean(),cv4.std(ddof=1)/np.sqrt(10)))
# 5基因版对照
T5=T.copy(); T5['mrg5']=sum(T5['z_'+g]*R['core']['coef'][g] for g in CORE if g!='OSGIN2')
T5['z_OSGIN2']=zmat(d,['OSGIN2'])[:,0]; T5['mrg5']=T5['mrg4']+T5['z_OSGIN2']*R['core']['coef']['OSGIN2']
app5=C(T5['time'],T5['ev'],T5['mrg5']); print('5基因版 apparent C=%.3f'%app5)

# ---------- TCGA 校正（IDH+grade） ----------
T['idhwt']=(T['idh']=='IDHwt').astype(int); T['g3']=(T['grade']=='G3').astype(int)
T['male']=(T['sex'].str[0]=='M').astype(int); T['mrg4z']=(T['mrg4']-T['mrg4'].mean())/T['mrg4'].std()
t2=T.dropna(subset=['subtype','grade'])
fT=CoxPHFitter().fit(t2[['time','ev','mrg4z','age','male','idhwt','g3']],'time','ev')
print('\n--- TCGA MRG-4 校正 IDH+grade (n=%d, deaths=%d) ---'%(len(t2),int(t2['ev'].sum())))
print(fT.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4).to_string())

# ---------- CGGA 外部验证 ----------
cexpr=pd.read_csv('cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt',sep='\t',index_col=0)
_raw=open('cgga/CGGA.mRNAseq_325_clinical.20200506.txt',encoding='utf-8-sig',errors='ignore').read().replace('\r','\n')
cl=pd.read_csv(io.StringIO(_raw),sep='\t'); cl.columns=[c.strip() for c in cl.columns]; cl=cl.set_index('CGGA_ID')
cl['ev']=pd.to_numeric(cl['Censor (alive=0; dead=1)'],errors='coerce'); cl['mos']=pd.to_numeric(cl['OS'],errors='coerce')/30.4375
cl['age']=pd.to_numeric(cl['Age'],errors='coerce'); cl['male']=(cl['Gender'].astype(str).str[0]=='M').astype(float)
cl['idhwt']=(cl['IDH_mutation_status'].astype(str).str.lower().str.startswith('wild')).astype(float)
cl['grade_hi']=cl['Grade'].isin(['WHO III','WHO IV']).astype(float)
M=cexpr.loc[MRG4].T; M=M[~M.index.duplicated()]
res={}
for tag,std in [('frozen_TCGA',True),('cohort_z',False)]:
    if std: Z=(M.values-mu[[genes.index(g) for g in MRG4]])/sd[[genes.index(g) for g in MRG4]]
    else: Z=(M.values-M.values.mean(0))/(M.values.std(0)+1e-9)
    lp=pd.Series(Z@np.array([coef4[g] for g in MRG4]),index=M.index)
    v=pd.concat([cl,lp.rename('lp')],axis=1).dropna(subset=['lp','mos','ev'])
    vL=v[v['Grade'].isin(['WHO II','WHO III'])].copy()
    vL['lp_z']=(vL['lp']-vL['lp'].mean())/vL['lp'].std()
    ff=CoxPHFitter().fit(vL[['mos','ev','lp_z','age','male','idhwt','grade_hi']].dropna(),'mos','ev')
    res[tag]={'C_lgg':float(C(vL['mos'],vL['ev'],vL['lp'])),'C_all':float(C(v['mos'],v['ev'],v['lp'])),
              'n_lgg':int(len(vL)),'deaths_lgg':int(vL['ev'].sum()),
              'adj':{k:{'HR':round(float(r['exp(coef)']),3),'lo':round(float(r['exp(coef) lower 95%']),3),
                        'hi':round(float(r['exp(coef) upper 95%']),3),'p':float(r['p'])} for k,r in ff.summary.iterrows()}}
    print('\n--- CGGA[%s] MRG-4 LGG n=%d deaths=%d C=%.3f (全分级 %.3f) ---'%(
          tag,len(vL),int(vL['ev'].sum()),res[tag]['C_lgg'],res[tag]['C_all']))
    print(ff.summary[['exp(coef)','exp(coef) lower 95%','exp(coef) upper 95%','p']].round(4).to_string())
    if tag=='cohort_z':
        vL.to_csv('cgga_mrg4_perpatient.csv',encoding='utf-8-sig')
        # 中位分割 KM
        vL['grp']=np.where(vL['lp']>vL['lp'].median(),'high','low')
        r=multivariate_logrank_test(vL['mos'],vL['grp'],vL['ev'])
        kmf=KaplanMeierFitter(); med={}
        for g in ['low','high']:
            s=vL[vL['grp']==g]; kmf.fit(s['mos'],s['ev']); med[g]=float(kmf.median_survival_time_)
        res['km']={'p':float(r.p_value),'median_low':med['low'],'median_high':med['high'],
                   'n_low':int((vL['grp']=='low').sum()),'n_high':int((vL['grp']=='high').sum())}
        print(' 中位分割 log-rank P=%.3e 中位OS 低=%.1f 高=%.1f 月'%(r.p_value,med['low'],med['high']))
        # 似然比比较：MRG-4 vs 25基因
        import csv as _csv
        co25={r_['基因']:float(r_['系数']) for r_ in _csv.DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))}
        g25=[g for g in co25 if g in cexpr.index]
        M25=cexpr.loc[g25].T; M25=M25[~M25.index.duplicated()]
        Z25=(M25.values-M25.values.mean(0))/(M25.values.std(0)+1e-9)
        lp25=pd.Series(Z25@np.array([co25[g] for g in g25]),index=M25.index)
        vv=pd.concat([vL.reset_index(drop=True),lp25.rename('lp25').reset_index(drop=True)],axis=1)
        vv=vv.dropna(subset=['mos','ev','lp_z','lp25'])
        ll4=CoxPHFitter().fit(vv[['mos','ev','lp_z']].astype(float),'mos','ev'); ll25=CoxPHFitter().fit(vv[['mos','ev','lp25']].astype(float),'mos','ev')
        res['lr_test']={'MRG4_loglik':float(ll4.log_likelihood_),'MRG4_AIC':float(ll4.AIC_partial_),
                        'g25_loglik':float(ll25.log_likelihood_),'g25_AIC':float(ll25.AIC_partial_),'n':int(len(vv))}
        print(' 外部队列 LR: MRG-4 loglik=%.2f AIC=%.1f | 25基因 loglik=%.2f AIC=%.1f (n=%d)'%(
              ll4.log_likelihood_,ll4.AIC_partial_,ll25.log_likelihood_,ll25.AIC_partial_,len(vv)))
res['tcga']={'apparent_C':app4,'cv_C':float(cv4.mean()),'cv_C_se':float(cv4.std(ddof=1)/np.sqrt(10)),
             'EPV':round(125/len(MRG4),1),'apparent_C_5gene':app5,
             'adj':{k:{'HR':round(float(r['exp(coef)']),3),'lo':round(float(r['exp(coef) lower 95%']),3),
                       'hi':round(float(r['exp(coef) upper 95%']),3),'p':float(r['p'])} for k,r in fT.summary.iterrows()},
             'n_adj':int(len(t2)),'deaths_adj':int(t2['ev'].sum())}
res['model']={'MRG4':MRG4,'coef':coef4,'CORE5':CORE,'coef5':R['core']['coef']}
T['grp4']=np.where(T['mrg4']>T['mrg4'].median(),'high','low')
rk=multivariate_logrank_test(T['time'],T['grp4'],T['ev']); kmf=KaplanMeierFitter(); md={}
for g in ['low','high']:
    ss=T[T['grp4']==g]; kmf.fit(ss['time'],ss['ev']); md[g]=float(kmf.median_survival_time_)
res['tcga_km']={'p':float(rk.p_value),'median_low':md['low'],'median_high':md['high']}
print('\nTCGA MRG-4 中位分割 log-rank P=%.3e 中位OS 低=%.1f 高=%.1f 月'%(rk.p_value,md['low'],md['high']))
json.dump(res,open('mrg4_results.json','w'),ensure_ascii=False,indent=1)
T[['pat','mrg4','risk','time','ev','age','sex','subtype','grade']].to_csv('mrg4_tcga_perpatient.csv',index=False,encoding='utf-8-sig')
print('\nsaved mrg4_results.json')
