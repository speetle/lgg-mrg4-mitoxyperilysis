# -*- coding: utf-8 -*-
import json, io, csv
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from sksurv.metrics import concordance_index_censored
from sklearn.model_selection import KFold
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.9,'figure.dpi':300})
CB='#1f5fa8'; CR='#c0392b'; CG='#888888'
C=lambda t,e,r: concordance_index_censored(np.asarray(e).astype(bool),np.asarray(t),np.asarray(r))[0]
import sys
ONLY14='--only14' in sys.argv
P=json.load(open('parsimonious_results.json')); M4=json.load(open('mrg4_results.json'))
if ONLY14: sys.argv=sys.argv
cv=P['cv']; alphas=np.array(cv['alphas']); ks=np.array(cv['k']); mean=np.array(cv['mean']); se=np.array(cv['se'])

# ================= FIG 12: LASSO CV + 稳定性选择 =================
fig=plt.figure(figsize=(9.2,3.4)); gs=gridspec.GridSpec(1,2,width_ratios=[1.15,1],wspace=0.28)
ax=fig.add_subplot(gs[0])
x=np.log10(alphas); m=ks<=60
ax.errorbar(x[m],mean[m],yerr=se[m],fmt='o',ms=2.6,lw=0.8,color=CB,capsize=1.6,alpha=0.85)
ax.axvline(np.log10(cv['alpha_1se']),color=CR,ls='--',lw=1.0)
ax.axvline(np.log10(cv['alpha_min']),color='#2e7d32',ls=':',lw=1.0)
ax.axhline(0.5,color=CG,ls='-',lw=0.7,alpha=0.6)
ax.text(np.log10(cv['alpha_1se']),0.53,'  $\\lambda_{1se}$ (%d genes)'%cv['k_1se'],color=CR,fontsize=7.5,rotation=90,va='bottom',ha='right')
ax.text(np.log10(cv['alpha_min']),0.53,'$\\lambda_{min}$',color='#2e7d32',fontsize=7.5,va='bottom',ha='left')
ax.set_xlabel('log$_{10}(\\lambda)$'); ax.set_ylabel('10-fold CV C-index')
ax.set_title('A  LASSO-Cox cross-validation',fontsize=9.5)
ax.set_ylim(0.48,0.85); ax.grid(alpha=0.15,lw=0.5); ax.spines[['top','right']].set_visible(False)
axt=ax.twiny(); axt.set_xlim(ax.get_xlim())
ticks=[np.log10(alphas[np.argmin(np.abs(ks-k))]) for k in [5,10,20,30]]
axt.set_xticks(ticks); axt.set_xticklabels(['5','10','20','30'],fontsize=8)
axt.set_xlabel('Number of genes in model',fontsize=8)

ax2=fig.add_subplot(gs[1])
st=sorted(P['stability'].items(),key=lambda kv:-kv[1])[:14][::-1]
cols=[CR if v>=0.75 else CG for _,v in st]
ax2.barh([g for g,_ in st],[v*100 for _,v in st],color=cols,alpha=0.85,height=0.68)
ax2.axvline(75,color='black',ls='--',lw=0.9)
ax2.text(76,0.1,'$\\pi$ = 0.75',fontsize=7.5,va='bottom')
ax2.set_xlabel('Selection frequency (%) in 200 bootstraps'); ax2.set_xlim(0,100)
ax2.set_title('B  Stability selection',fontsize=9.5); ax2.tick_params(labelsize=8)
ax2.grid(axis='x',alpha=0.15,lw=0.5); ax2.spines[['top','right']].set_visible(False)
plt.savefig('Fig12_LASSO与稳定性选择.png',bbox_inches='tight'); plt.savefig('Fig12_LASSO与稳定性选择.pdf',bbox_inches='tight'); plt.close()
print('Fig12 written')

# ================= 计算 25 基因内部 CV（供 Fig13B） =================
E=pd.read_csv('panel_expr_tcga.csv',index_col=0,encoding='utf-8-sig')
G=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
d=G.merge(E,left_on='pat',right_index=True,how='inner').rename(columns={'e':'ev','t':'time'})
genes=list(E.columns); mu,sd=E.values.mean(0),E.values.std(0)+1e-9
co25={r['基因']:float(r['系数']) for r in csv.DictReader(open('/Users/lianbin/Documents/workbuddy/科研工作台/02-数据/LGG/多因素Cox.csv',encoding='utf-8-sig'))}
g25=[g for g in co25 if g in genes]
Z=d[g25].apply(lambda c:(c-c.mean())/c.std()); Z['time']=d['time']; Z['ev']=d['ev']
cv25=[]
for tr,te in KFold(10,shuffle=True,random_state=7).split(Z):
    f=CoxPHFitter().fit(Z.iloc[tr][g25+['time','ev']],'time','ev')
    cv25.append(C(Z.iloc[te]['time'],Z.iloc[te]['ev'],f.predict_partial_hazard(Z.iloc[te])))
cv25=np.array(cv25)
l26=[g for g in genes if abs(P['lasso_1se']['coef'].get(g,0))>1e-8]
Zl=d[l26].apply(lambda c:(c-c.mean())/c.std()); Zl['time']=d['time']; Zl['ev']=d['ev']
f26=CoxPHFitter().fit(Zl[l26+['time','ev']],'time','ev')
app26=C(Zl['time'],Zl['ev'],f26.predict_partial_hazard(Zl))
ext=P['external']['cohort_z']; extT=P['external']['frozen_TCGA']
bars={'25-gene\noriginal':[0.821,cv25.mean(),None,ext['25gene']['C_lgg']],
      'LASSO\n(%d genes)'%len(l26):[app26,0.777,None,ext['LASSO_1se']['C_lgg']],
      'MRG-4\n(this study)':[M4['tcga']['apparent_C'],M4['tcga']['cv_C'],None,ext['core']['C_lgg']]}
print('25gene CV=%.3f | LASSO26 apparent=%.3f CV=0.777 | ext C: 25=%.3f l26=%.3f mrg4=%.3f'%(
      cv25.mean(),app26,ext['25gene']['C_lgg'],ext['LASSO_1se']['C_lgg'],ext['core']['C_lgg']))

# ================= FIG 13: 外部验证 KM + C-index 对比 =================
vL=pd.read_csv('cgga_mrg4_perpatient.csv')
vL['grp']=np.where(vL['lp']>vL['lp'].median(),'high','low')
fig,axes=plt.subplots(1,2,figsize=(9.0,3.4))
ax=axes[0]
for g,col,lab in [('low',CB,'Low MRG-4'),('high',CR,'High MRG-4')]:
    s=vL[vL['grp']==g]; k=KaplanMeierFitter(); k.fit(s['mos'],s['ev'],label=lab)
    ax.step(k.survival_function_.index,k.survival_function_.iloc[:,0],where='post',color=col,lw=1.6)
r=multivariate_logrank_test(vL['mos'],vL['grp'],vL['ev'])
ax.set_xlabel('Time (months)'); ax.set_ylabel('Overall survival'); ax.set_ylim(0,1.02); ax.set_xlim(0,None)
ax.set_title('A  CGGA-325 WHO II-III, external validation\\n(n = %d, %d deaths)'%(len(vL),int(vL['ev'].sum())),fontsize=9)
ax.text(0.97,0.05,'log-rank $P$ = %.1e'%r.p_value,transform=ax.transAxes,ha='right',fontsize=8.5,
        bbox=dict(facecolor='white',edgecolor='none',alpha=0.9,pad=1.5))
ax.legend(frameon=False,fontsize=8,loc='upper right'); ax.grid(alpha=0.15,lw=0.5)
ax.spines[['top','right']].set_visible(False)

ax=axes[1]
labels=list(bars.keys()); xpos=np.arange(len(labels)); w=0.27
series=[('TCGA apparent',[bars[k][0] for k in labels],CB),
        ('TCGA 10-fold CV',[bars[k][1] for k in labels],'#7fb3d5'),
        ('CGGA external',[bars[k][3] for k in labels],CR)]
for i,(nm,vals,col) in enumerate(series):
    b=ax.bar(xpos+(i-1)*w,vals,w,label=nm,color=col,alpha=0.88)
    for r_ in b: ax.text(r_.get_x()+r_.get_width()/2,r_.get_height()+0.008,'%.3f'%r_.get_height(),ha='center',fontsize=6.6)
ax.set_xticks(xpos); ax.set_xticklabels(labels,fontsize=8)
ax.set_ylabel('C-index'); ax.set_ylim(0,0.93); ax.axhline(0.5,color=CG,ls='--',lw=0.8)
ax.set_title('B  Discrimination: parsimony improves\\ngeneralisation',fontsize=9)
ax.legend(frameon=False,fontsize=7.5,ncol=1,loc='upper left',bbox_to_anchor=(0.0,0.99))
ax.grid(axis='y',alpha=0.15,lw=0.5); ax.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('Fig13_外部验证与区分度.png',bbox_inches='tight'); plt.savefig('Fig13_外部验证与区分度.pdf',bbox_inches='tight'); plt.close()
print('Fig13 written')

# ================= FIG 14: MRG-4 森林图（TCGA + CGGA 校正） =================
rows=[('H','TCGA-LGG (n = %d, %d deaths)'%(M4['tcga']['n_adj'],M4['tcga']['deaths_adj']),None,None,None,None)]
for k,lab in [('mrg4z','MRG-4 score (per SD)'),('age','Age (per year)'),('male','Male sex'),('idhwt','IDH-wildtype'),('g3','WHO grade 3')]:
    r=M4['tcga']['adj'][k]; rows.append(('R',lab,r['HR'],r['lo'],r['hi'],r['p']))
cl=M4['cohort_z']['adj']
rows.append(('H','CGGA-325 WHO II-III (n = %d)'%M4['cohort_z']['n_lgg'],None,None,None,None))
for k,lab in [('lp_z','MRG-4 score (per SD)'),('age','Age (per year)'),('male','Male sex'),('idhwt','IDH-wildtype'),('grade_hi','WHO grade III')]:
    r=cl[k]; rows.append(('R',lab,r['HR'],r['lo'],r['hi'],r['p']))
ys=[];y=0
for r_ in rows:
    y-= 1.35 if r_[0]=='H' else 1.0; ys.append(y)
ymin,ymax=min(ys)-1.0,max(ys)+1.5
fig=plt.figure(figsize=(8.4,4.4)); gs=gridspec.GridSpec(1,3,width_ratios=[1.15,1.0,0.62],wspace=0.05)
axl=fig.add_subplot(gs[0]); axf=fig.add_subplot(gs[1],sharey=axl); axr=fig.add_subplot(gs[2],sharey=axl)
for a in (axl,axf,axr): a.set_ylim(ymin,ymax); a.set_yticks([])
for a in (axl,axr):
    for s_ in ['top','right','left','bottom']: a.spines[s_].set_visible(False)
    a.set_xticks([]); a.tick_params(bottom=False,left=False)
axf.spines[['top','right']].set_visible(False)
axf.axvline(1,color='#999',ls='--',lw=0.9); axf.set_xscale('log'); axf.set_xlim(0.5,22)
tr=axl.get_yaxis_transform()
for (kind,lab,hr,lo,hi,p),yy in zip(rows,ys):
    if kind=='H':
        axl.text(1.0,yy,lab,ha='right',va='center',fontsize=8.2,fontweight='bold',transform=tr); continue
    axl.text(0.985,yy,lab,ha='right',va='center',fontsize=8.2,transform=tr)
    col=CR if (lo>1 or hi<1) else CG
    axf.plot([lo,hi],[yy,yy],color=col,lw=1.4,solid_capstyle='round'); axf.plot([hr],[yy],'s',color=col,ms=4.6)
    axr.text(0.02,yy,'%.2f (%.2f-%.2f)'%(hr,lo,hi),ha='left',va='center',fontsize=7.7,transform=axr.get_yaxis_transform())
    axr.text(0.99,yy,('%.1e'%p) if p<0.001 else ('%.3f'%p),ha='right',va='center',fontsize=7.7,transform=axr.get_yaxis_transform())
axl.text(1.0,ymax-0.05,'Covariate',ha='right',va='center',fontsize=8,style='italic',color='#444',transform=tr)
axr.text(0.02,ymax-0.05,'HR (95% CI)',ha='left',va='center',fontsize=8,style='italic',color='#444',transform=axr.get_yaxis_transform())
axr.text(0.99,ymax-0.05,'$P$',ha='right',va='center',fontsize=8,style='italic',color='#444',transform=axr.get_yaxis_transform())
axf.set_xticks([0.5,1,2,4,8,16]); axf.set_xticklabels(['0.5','1.0','2','4','8','16'])
axf.set_xlabel('Hazard ratio (log scale)'); axf.grid(axis='x',alpha=0.12,lw=0.5)
axf.set_title('MRG-4 adjusted for IDH status and grade',fontsize=9.5,pad=14)
plt.savefig('Fig14_MRG4森林图.png',bbox_inches='tight'); plt.savefig('Fig14_MRG4森林图.pdf',bbox_inches='tight'); plt.close()
print('Fig14 written')
json.dump({'bars':bars,'cv25':float(cv25.mean()),'app26':float(app26),'cv26':0.777,
           'ext_frozen':{k:{'lgg':extT[k]['C_lgg'],'all':extT[k]['C_all']} for k in extT},
           'ext_cohort':{k:{'lgg':ext[k]['C_lgg'],'all':ext[k]['C_all']} for k in ext}},
          open('fig13_data.json','w'),indent=1)
