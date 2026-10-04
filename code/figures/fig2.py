# -*- coding: utf-8 -*-
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, matplotlib.gridspec as gridspec
from lifelines import CoxPHFitter
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'figure.dpi':300})
CR='#c0392b'; CG='#555555'

m=pd.read_csv('lgg_g1_dataset.csv',encoding='utf-8-sig')
m['idhwt']=(m['idh']=='IDHwt').astype(int); m['g3']=(m['grade']=='G3').astype(int); m['male']=(m['sex'].str[0]=='M').astype(int)
m['risk_z']=(m['risk']-m['risk'].mean())/m['risk'].std()
a=m.dropna(subset=['subtype','grade']).copy()
c=pd.read_csv('cgga_dataset.csv',encoding='utf-8-sig'); cl=c[c['Grade'].isin(['WHO II','WHO III'])].copy()
gd=pd.read_csv('g3_dataset.csv',encoding='utf-8-sig').dropna(subset=['risk_hat','prolif','t','e','age','male','idhwt','g3'])

cB=CoxPHFitter().fit(a[['t','e','risk_z','age','male','idhwt','g3']].dropna(),'t','e')
d=cl.dropna(subset=['risk','t','ev','age','idhwt','grade_hi'])
cC=CoxPHFitter().fit(d[['t','ev','risk','age','male','idhwt','grade_hi']],'t','ev')
cD=CoxPHFitter().fit(gd[['t','e','risk_hat','prolif','idhwt','g3','age','male']],'t','e')

rows=[]
def add(gname,items,fit):
    rows.append(('H',gname,None,None,None,None))
    for lab,k in items:
        s=fit.summary.loc[k]
        rows.append(('R',lab,float(s['exp(coef)']),float(s['exp(coef) lower 95%']),float(s['exp(coef) upper 95%']),float(s['p'])))
add('TCGA-LGG, multivariable (n = 502)',[('Risk score (per SD)','risk_z'),('Age (per year)','age'),('Male sex','male'),('IDH-wildtype','idhwt'),('WHO grade 3','g3')],cB)
add('CGGA-325, WHO II-III (n = 172)',[('Risk score (per SD)','risk'),('Age (per year)','age'),('Male sex','male'),('IDH-wildtype','idhwt'),('WHO grade III','grade_hi')],cC)
add('TCGA-LGG + proliferation score (n = 505)',[('Risk score (per SD)','risk_hat'),('Proliferation score (per SD)','prolif')],cD)

ys=[];y=0
for r in rows:
    y-= 1.35 if r[0]=='H' else 1.0
    ys.append(y)
ymin=min(ys)-1.0; ymax=max(ys)+1.5

fig=plt.figure(figsize=(8.4,4.6))
gs=gridspec.GridSpec(1,3,width_ratios=[1.15,1.0,0.62],wspace=0.05)
axl=fig.add_subplot(gs[0]); axf=fig.add_subplot(gs[1],sharey=axl); axr=fig.add_subplot(gs[2],sharey=axl)
for ax in (axl,axf,axr): ax.set_ylim(ymin,ymax); ax.set_yticks([])
for s in ['top','right','left','bottom']: axl.spines[s].set_visible(False); axr.spines[s].set_visible(False)
axf.spines[['top','right']].set_visible(False)
for ax in (axl,axr):
    ax.set_xticks([]); ax.tick_params(bottom=False,left=False)
tr=axl.get_yaxis_transform()
axf.axvline(1,color='#999',ls='--',lw=0.9); axf.set_xscale('log'); axf.set_xlim(0.6,40)

for (kind,lab,hr,lo,hi,p),yy in zip(rows,ys):
    if kind=='H':
        axl.text(1.0,yy,lab,ha='right',va='center',fontsize=8.2,fontweight='bold',color='#222',transform=tr)
        continue
    axl.text(0.985,yy,lab,ha='right',va='center',fontsize=8.2,transform=tr)
    col=CR if (lo>1 or hi<1) else CG
    axf.plot([lo,hi],[yy,yy],color=col,lw=1.4,solid_capstyle='round'); axf.plot([hr],[yy],'s',color=col,ms=4.6)
    axr.text(0.02,yy,f'{hr:.2f} ({lo:.2f}-{hi:.2f})',ha='left',va='center',fontsize=7.7,transform=axr.get_yaxis_transform())
    axr.text(0.99,yy,('%.1e'%p) if p<0.001 else ('%.3f'%p),ha='right',va='center',fontsize=7.7,transform=axr.get_yaxis_transform())
axl.text(1.0,ymax-0.05,'Covariate',ha='right',va='center',fontsize=8,style='italic',color='#444',transform=tr)
axr.text(0.02,ymax-0.05,'HR (95% CI)',ha='left',va='center',fontsize=8,style='italic',color='#444',transform=axr.get_yaxis_transform())
axr.text(0.99,ymax-0.05,'$P$',ha='right',va='center',fontsize=8,style='italic',color='#444',transform=axr.get_yaxis_transform())
axf.set_xticks([1,1.5,2,3,5,8,15]); axf.set_xticklabels(['1.0','1.5','2','3','5','8','15'])
axf.set_xlabel('Hazard ratio (log scale)'); axf.grid(axis='x',alpha=0.12,lw=0.5)
axf.set_title('Multivariable Cox models',fontsize=10,pad=14)
plt.savefig('Fig10_多因素森林图.png',bbox_inches='tight'); plt.savefig('Fig10_多因素森林图.pdf',bbox_inches='tight'); print('Fig10 written')
