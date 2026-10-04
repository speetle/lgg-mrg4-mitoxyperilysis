# -*- coding: utf-8 -*-
"""Figure 8 + Supplementary Table S8 — corrected somatic mutation analysis (n=255/255)."""
import json, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.linewidth':0.9,'figure.dpi':300})
CB='#1f5fa8'; CR='#c0392b'; CG='#888888'

S = json.load(open('mut_final.json'))
rows = S['rows']
nh, nl = S['n_hi'], S['n_lo']
genes = [r['gene'] for r in rows]
order = sorted(range(len(genes)), key=lambda i: -(rows[i]['ph'] - rows[i]['pl']))
genes = [genes[i] for i in order]; rows = [rows[i] for i in order]

hi = np.array([r['ph'] for r in rows]); lo = np.array([r['pl'] for r in rows])
pv = [r['p'] for r in rows]

def star(p):
    return '***' if p < 1e-3 else '**' if p < 1e-2 else '*' if p < 0.05 else 'n.s.'

# ---------------- Figure 8 ----------------
fig, ax = plt.subplots(figsize=(7.0, 4.2))
y = np.arange(len(genes)); h = 0.38
ax.barh(y + h/2, hi, height=h, color=CR, label='High-risk (n = %d)' % nh, edgecolor='none')
ax.barh(y - h/2, lo, height=h, color=CB, label='Low-risk (n = %d)' % nl, edgecolor='none')
for i, (a, b) in enumerate(zip(hi, lo)):
    ax.text(a + 1.2, i + h/2, '%.1f' % a, va='center', fontsize=7.2, color=CR)
    ax.text(b + 1.2, i - h/2, '%.1f' % b, va='center', fontsize=7.2, color=CB)
    ax.text(104, i, star(pv[i]), va='center', ha='left', fontsize=7.8,
            color='#222222' if pv[i] < 0.05 else CG)
ax.text(104, -0.95, 'Fisher $P$', va='center', ha='left', fontsize=7.8, color='#222222')
ax.set_yticks(y); ax.set_yticklabels(genes, fontsize=8.6)
ax.set_ylim(len(genes) - 0.5, -1.35)
ax.set_xlim(0, 115); ax.set_xticks([0, 20, 40, 60, 80, 100])
ax.set_xlabel('Patients with a non-synonymous mutation (%)')
ax.set_title('Somatic mutation frequencies by risk group (TCGA-LGG, n = 510)', fontsize=9.4)
ax.grid(axis='x', alpha=0.15, lw=0.5); ax.set_axisbelow(True)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(frameon=False, fontsize=8.2, loc='lower left', bbox_to_anchor=(0.30, -0.235), ncol=2)
fig.text(0.99, 0.005, 'All 510 patients had evaluable calls for the ten genes; *** $P$ < 0.001, ** $P$ < 0.01, * $P$ < 0.05, n.s. not significant.',
         ha='right', va='bottom', fontsize=6.4, color=CG)
fig.tight_layout(rect=(0, 0.045, 1, 1))
for ext in ('png', 'pdf'):
    fig.savefig('Fig8_突变谱_校正版.%s' % ext, dpi=300, bbox_inches='tight')
print('OK  Fig8_突变谱_校正版.png/.pdf')

# ---------------- Supplementary Table S8 ----------------
mut = json.load(open('lgg_mut_target.json'))
NS = {'Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins',
      'In_Frame_Del','In_Frame_Ins','Splice_Site','Nonstop_Mutation','Translation_Start_Site'}
from collections import defaultdict
pat = defaultdict(set)
for g, recs in mut.items():
    for m in recs:
        if m['mutationType'] in NS:
            pat[g].add(m['patientId'])

COH = {}
with open('lgg_g1_dataset.csv', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f):
        COH[r['pid']] = r

G = ["IDH1", "TP53", "ATRX", "EGFR", "PTEN", "NF1", "CIC", "FUBP1", "NOTCH1", "PIK3CA"]
with open('Table_S8_mutation_per_patient.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Patient_ID'] + G + ['N_genes_mutated', 'Risk_group', 'ModelA_score_25gene'])
    for p, v in COH.items():
        flags = [1 if p in pat[g] else 0 for g in G]
        w.writerow([p] + flags + [sum(flags), v['组'], round(float(v['risk']), 6)])
print('OK  Table_S8_mutation_per_patient.csv  rows=%d' % len(COH))

with open('Table_S9_mutation_summary.csv', 'w', newline='', encoding='utf-8-sig') as f:
    w = csv.writer(f)
    w.writerow(['Gene', 'High_n', 'High_total', 'High_pct', 'Low_n', 'Low_total', 'Low_pct',
                'Fisher_P', 'Direction'])
    for r in rows:
        d = 'enriched in high-risk' if r['ph'] > r['pl'] else 'enriched in low-risk'
        if r['p'] >= 0.05:
            d = 'no difference'
        w.writerow([r['gene'], r['h'], r['nh'], '%.1f' % r['ph'], r['l'], r['nl'],
                    '%.1f' % r['pl'], '%.3g' % r['p'], d])
print('OK  Table_S9_mutation_summary.csv')
