# -*- coding: utf-8 -*-
"""MRG-4 vs 25基因模型：NRI/IDI、时间依赖AUC、决策曲线(DCA)"""
import json, io, csv
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from sksurv.util import Surv
from sksurv.metrics import cumulative_dynamic_auc, concordance_index_censored
np.random.seed(42); pd.set_option('display.width',210)
C=lambda t,e,r: concordance_index_censored(np.asarray(e).astype(bool),np.asarray(t),np.asarray(r))[0]
M4=json.load(open('mrg4_results.json')); MRG4=M4['model']['MRG4']; coef4=M4['model']['coef']

# ---- 表达与临床 ----
E=pd.read_csv('panel_expr_tcga.csv',index_col=0,encoding='utf-8-sig')
genes=list(E.columns); mu,sd=E.values.mean(0),E.values.std(0)+1e-9
co25={r['基因']:float(r['系数']) for r in csv.DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))}
g25=[g for g in co25 if g in genes]
cexpr=pd.read_csv('cgga/CGGA.mRNAseq_325.RSEM-genes.20200506.txt',sep='\t',index_col=0)
_raw=open('cgga/CGGA.mRNAseq_325_clinical.20200506.txt',encoding='utf-8-sig',errors='ignore').read().replace('\r','\n')
cl=pd.read_csv(io.StringIO(_raw),sep='\t'); cl.columns=[c.strip() for c in cl.columns]; cl=cl.set_index('CGGA_ID')
cl['ev']=pd.to_numeric(cl['Censor (alive=0; dead=1)'],errors='coerce'); cl['mos']=pd.to_numeric(cl['OS'],errors='coerce')/30.4375
G=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
tc=G.merge(E,left_on='pat',right_index=True,how='inner').rename(columns={'e':'ev','t':'mos'})

def zc2(df,glist,cohort=False,mu_=None,sd_=None):
    V=df[glist].values
    if not cohort: V=(V-mu_[ [genes.index(g) for g in glist]])/sd_[[genes.index(g) for g in glist]]
    else: V=(V-V.mean(0))/(V.std(0)+1e-9)
    return V

# 标准化（队列内，与前面口径一致）
Z_tc=lambda gl: (tc[gl].values-tc[gl].values.mean(0))/(tc[gl].values.std(0)+1e-9)
Z_cg=lambda gl: (cexpr.loc[gl].T.values-cexpr.loc[gl].T.values.mean(0))/(cexpr.loc[gl].T.values.std(0)+1e-9)
lgg_ids=cl[cl['Grade'].isin(['WHO II','WHO III'])].index
pos=[list(cexpr.columns).index(i) for i in lgg_ids if i in set(cexpr.columns)]
ok=[i for i in lgg_ids if i in set(cexpr.columns)]

X_tc4=Z_tc(MRG4); X_tc25=Z_tc(g25)
allg=MRG4+[g for g in g25 if g not in MRG4]
Mc=cexpr.loc[allg].T.loc[ok]
X_cg_all=(Mc.values-Mc.values.mean(0))/(Mc.values.std(0)+1e-9)
X_cg4=X_cg_all[:,[allg.index(g) for g in MRG4]]
X_cg25=X_cg_all[:,[allg.index(g) for g in g25]]

# ---- 在 TCGA 拟合两个 Cox（用于预测外部风险）----
f4=CoxPHFitter().fit(pd.DataFrame(X_tc4,columns=MRG4).assign(mos=tc['mos'].values,ev=tc['ev'].values),'mos','ev')
f25=CoxPHFitter().fit(pd.DataFrame(X_tc25,columns=g25).assign(mos=tc['mos'].values,ev=tc['ev'].values),'mos','ev')

cg=pd.DataFrame(X_cg_all,columns=allg,index=ok); cg['mos']=cl.loc[ok,'mos']; cg['ev']=cl.loc[ok,'ev']
cg['lp4']=f4.predict_partial_hazard(cg[MRG4]); cg['lp25']=f25.predict_partial_hazard(cg[g25])
_n0=len(cg); cg=cg.dropna(subset=['lp4','lp25','mos','ev']); print('CGGA LGG 可评估: %d/%d (缺失基因值剔除 %d)'%(len(cg),_n0,_n0-len(cg)))
res={}
def risk_at(f,df,gl,t):
    sf=f.predict_survival_function(df[gl],times=[t])
    return 1-sf.iloc[0].values
for T_ in [36,60]:
    cg['r4']=risk_at(f4,cg,MRG4,T_); cg['r25']=risk_at(f25,cg,g25,T_)
    d=cg.dropna(subset=['r4','r25','mos','ev']).copy()
    d['status']=np.where(d['mos']<=T_,d['ev'],0)
    d=d[(d['mos']>T_)|(d['ev']==1)].copy()   # 排除 T 前删失且未死
    y=d['status'].astype(int)
    up=(d['r4']>d['r25']); dn=(d['r4']<d['r25'])
    n_e=int(y.sum()); n_n=int((1-y).sum())
    nri_e=(up[y==1].sum()-dn[y==1].sum())/n_e
    nri_n=(dn[y==0].sum()-up[y==0].sum())/n_n
    nri=nri_e+nri_n
    isl4=d.loc[y==1,'r4'].mean()-d.loc[y==0,'r4'].mean(); isl25=d.loc[y==1,'r25'].mean()-d.loc[y==0,'r25'].mean()
    idi=isl4-isl25
    auc4=C(d['mos'],d['ev'],d['r4']); auc25=C(d['mos'],d['ev'],d['r25'])
    res['NRI_IDI_%dmo'%T_]={'n':int(len(d)),'events':n_e,'NRI':float(nri),'NRI_events':float(nri_e),'NRI_nonevents':float(nri_n),
                            'IDI':float(idi),'AUC4':float(auc4),'AUC25':float(auc25)}
    print('[%d 月] n=%d 事件=%d | 连续NRI=%+.3f（事件%+.3f / 非事件%+.3f）| IDI=%+.4f | 时点AUC: MRG-4 %.3f vs 25基因 %.3f'%(
          T_,len(d),n_e,nri,nri_e,nri_n,idi,auc4,auc25))
# ---- 自举 NRI/IDI CI ----
B=500; bs={36:[],60:[]}
for T_ in [36,60]:
    cg['r4']=risk_at(f4,cg,MRG4,T_); cg['r25']=risk_at(f25,cg,g25,T_)
    d=cg.dropna(subset=['r4','r25','mos','ev']).copy()
    d['status']=np.where(d['mos']<=T_,d['ev'],0); d=d[(d['mos']>T_)|(d['ev']==1)].copy()
    for _ in range(B):
        s=d.sample(len(d),replace=True); y=s['status'].astype(int)
        up=(s['r4']>s['r25']); dn=(s['r4']<s['r25'])
        if y.sum()<3 or (1-y).sum()<3: continue
        bs[T_].append(((up[y==1].sum()-dn[y==1].sum())/y.sum())+((dn[y==0].sum()-up[y==0].sum())/(1-y).sum()))
for T_ in [36,60]:
    a=np.array(bs[T_]); lo,hi=np.percentile(a,[2.5,97.5])
    res['NRI_IDI_%dmo'%T_]['NRI_CI']=[float(lo),float(hi)]
    print('  %d月 NRI 95%%CI (bootstrap %d): %.3f 到 %.3f'%(T_,len(a),lo,hi))

# ---- 时间依赖 AUC（CGGA LGG） ----
tr=Surv.from_arrays(tc['ev'].astype(bool).values,tc['mos'].values)
te=Surv.from_arrays(cg['ev'].astype(bool).values,cg['mos'].values)
times=np.arange(12,121,6)
estimates=np.vstack([cg['lp4'].values,cg['lp25'].values]).T*-1  # AUC 需要风险越高值越大 -> sksurv 用 estimate 为风险
auc,mean_auc=cumulative_dynamic_auc(tr,te,cg['lp4'].values,times)
auc25,mean25=cumulative_dynamic_auc(tr,te,cg['lp25'].values,times)
res['td_auc']={'times':times.tolist(),'MRG4':auc.tolist(),'g25':auc25.tolist(),'mean_MRG4':float(mean_auc),'mean_g25':float(mean25)}
print('\n时间依赖AUC(CGGA LGG): MRG-4 均值=%.3f | 25基因 均值=%.3f'%(mean_auc,mean25))

# ---- DCA（3年死亡） ----
def dca(d,r,col):
    out=[]
    for pt in np.arange(0.2,0.81,0.02):
        pr=(d[r].values>=pt).astype(int); y=d['status'].values
        tp=((pr==1)&(y==1)).sum(); fp=((pr==1)&(y==0)).sum(); n=len(d)
        out.append((tp/n)-(fp/n)*(pt/(1-pt)))
    return np.array(out)
T_=36
cg['r4']=risk_at(f4,cg,MRG4,T_); cg['r25']=risk_at(f25,cg,g25,T_)
d=cg.dropna(subset=['r4','r25']).copy(); d['status']=np.where(d['mos']<=T_,d['ev'],0)
d=d[(d['mos']>T_)|(d['ev']==1)].copy()
pts=np.arange(0.2,0.81,0.02)
nb4=dca(d,'r4',None); nb25=dca(d,'r25',None)
prev=d['status'].mean()
nb_all=prev-(1-prev)*(pts/(1-pts))
res['dca']={'pt':pts.tolist(),'MRG4':nb4.tolist(),'g25':nb25.tolist(),'treat_all':nb_all.tolist(),'prevalence':float(prev)}
# ---- 类别型 NRI（各自模型三分位；对校准不敏感）----
d36=d.copy()
for c_ in ['r4','r25']:
    d36[c_+'_cat']=pd.qcut(d36[c_],3,labels=[0,1,2]).astype(int)
y=d36['status'].astype(int)
lo_=lambda c: ((d36[c]==0))
hi_=lambda c: ((d36[c]==2))
cat_nri=((hi_('r4_cat')&(d36['r25_cat']<d36['r4_cat']))[y==1].sum()-(lo_('r4_cat')&(d36['r25_cat']>d36['r4_cat']))[y==1].sum())/y.sum() \
      + ((lo_('r4_cat')&(d36['r25_cat']>d36['r4_cat']))[y==0].sum()-(hi_('r4_cat')&(d36['r25_cat']<d36['r4_cat']))[y==0].sum())/(1-y).sum()
res['cat_NRI_36mo']=float(cat_nri); print('\n类别型NRI(三分位, 36月) = %+.3f'%cat_nri)
print('\nDCA(36月) 阈值内 MRG-4 净获益均值=%.3f | 25基因=%.3f | 全治疗=%.3f'%(nb4.mean(),nb25.mean(),nb_all.mean()))
json.dump(res,open('nri_idi_auc.json','w'),ensure_ascii=False,indent=1)
print('\nsaved nri_idi_auc.json')
