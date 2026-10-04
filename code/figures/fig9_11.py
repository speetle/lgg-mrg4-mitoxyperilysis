# -*- coding: utf-8 -*-
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test
from scipy.stats import spearmanr
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.9,'figure.dpi':300})
CB='#1f5fa8'; CR='#c0392b'
m=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
m['idhwt']=(m['idh']=='IDHwt').astype(int); m['g3']=(m['grade']=='G3').astype(int)
a=m.dropna(subset=['subtype','grade']).copy()
c=pd.read_csv('cgga_dataset.csv',encoding='utf-8-sig'); cl=c[c['Grade'].isin(['WHO II','WHO III'])].copy()
cl['grp']=np.where(cl['risk']>cl['risk'].median(),'高风险','低风险')

def km(ax,df,t,e,grp,pval,title):
    for g,col,lab in [('低风险',CB,'Low risk'),('高风险',CR,'High risk')]:
        s=df[df[grp]==g]; k=KaplanMeierFitter(); k.fit(s[t],s[e],label=lab)
        ax.step(k.survival_function_.index,k.survival_function_.iloc[:,0],where='post',color=col,lw=1.6)
    ax.set_ylim(0,1.02); ax.set_xlabel('Time (months)'); ax.set_ylabel('Overall survival')
    ax.set_title(title,fontsize=9); ax.text(0.97,0.05,f'log-rank $P$ = {pval}',transform=ax.transAxes,ha='right',fontsize=8.5,bbox=dict(facecolor='white',edgecolor='none',alpha=0.88,pad=1.6))
    ax.legend(frameon=False,fontsize=8,loc='upper right'); ax.grid(alpha=0.15,lw=0.5)

fig,axes=plt.subplots(1,3,figsize=(9.2,3.0))
sub=a[a['idhwt']==0]; r=multivariate_logrank_test(sub['t'],sub['组'],sub['e'])
km(axes[0],sub,'t','e','组',f'{r.p_value:.1e}','TCGA-LGG, IDH-mutant (n=%d)'%len(sub))
sub=a[a['g3']==1]; r=multivariate_logrank_test(sub['t'],sub['组'],sub['e'])
km(axes[1],sub,'t','e','组',f'{r.p_value:.1e}','TCGA-LGG, WHO grade 3 (n=%d)'%len(sub))
r=multivariate_logrank_test(cl['t'],cl['grp'],cl['ev'])
km(axes[2],cl,'t','ev','grp',f'{r.p_value:.1e}','CGGA-325, WHO II-III (n=%d)'%len(cl))
plt.tight_layout()
plt.savefig('Fig9_分层与外部验证KM.png',bbox_inches='tight'); plt.savefig('Fig9_分层与外部验证KM.pdf',bbox_inches='tight'); plt.close()
print('Fig9 written')

# Fig11
gd=pd.read_csv('g3_dataset.csv',encoding='utf-8-sig').dropna(subset=['risk_hat','prolif'])
rho,pv=spearmanr(gd['risk_hat'],gd['prolif'])
fig,axes=plt.subplots(1,2,figsize=(8.0,3.2))
axes[0].scatter(gd['prolif'],gd['risk_hat'],s=5,alpha=0.45,color=CB,lw=0)
axes[0].set_xlabel('Proliferation score (mean $z$)'); axes[0].set_ylabel('MRG risk score')
axes[0].set_title(f'A  Spearman $\\rho$ = {rho:.2f}, $P$ = {pv:.1e}',fontsize=9)
bp=axes[1].boxplot([a[a['idhwt']==0]['risk'],a[a['idhwt']==1]['risk']],patch_artist=True,widths=0.55,
                   medianprops=dict(color='black'),tick_labels=['IDH-mutant','IDH-wildtype'])
for p,col in zip(bp['boxes'],[CB,CR]): p.set_facecolor(col); p.set_alpha(0.55)
axes[1].set_ylabel('MRG risk score'); axes[1].set_title('B  Risk score by IDH status',fontsize=9)
for ax in axes: ax.grid(alpha=0.15,lw=0.5); ax.spines[['top','right']].set_visible(False)
plt.tight_layout()
plt.savefig('Fig11_非增殖指数验证.png',bbox_inches='tight'); plt.savefig('Fig11_非增殖指数验证.pdf',bbox_inches='tight'); plt.close()
print('Fig11 written (B: risk by IDH)')
