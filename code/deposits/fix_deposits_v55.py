# -*- coding: utf-8 -*-
"""v5.5 沉积件修正：③ 旧表号归一 / ④ Table 3→4 / ⑤ RMST 显示精度 /
⑥ 0.870 出处披露 / ⑨ PROBAST 理由 / ⑪ S25 说明 / ⑫ 两个 SD。

与 v5.4 的 purge_versions_v54.py 同一纪律：
  · 逐条按 (主键列) 定位行，不按行号（行号会随他处编辑漂移）
  · 每条替换都断言命中数，任何一条对不上就整体不写出
  · 写出到 _stage55/，由 fix_supplement_v55.py 再装入交付目录
"""
import csv
import io
import os
import shutil

# ⚠️ 改过的六张读**归档的 v5.4 原件**，不是交付目录里的现件。
#    本脚本是对 v5.4 那一次状态做定点修正；若读现件，第二次运行就会因
#    「该改的已经改过了」而断言失败（实测：S13 题名在行内命中 0 次）。
#    没被本轮改过的（S2/S4/S6 等）归档里没有，回落到交付目录现件。
LIVE = ('/Users/lianbin/Desktop/Mitoxyperilysis课题总包/04-学生论文框架/'
        '学生1_LGG_交付_20261004/补充材料')
OLD = ('/Users/lianbin/Desktop/Mitoxyperilysis课题总包/04-学生论文框架/'
       '学生1_LGG_交付_20261004/历史版本_v5.4/补充材料_被修订')
STAGE = '/tmp/mitoxy/_stage55'
os.makedirs(STAGE, exist_ok=True)

CHANGED = {}


def load(name):
    p = os.path.join(OLD, name)
    if not os.path.exists(p):
        p = os.path.join(LIVE, name)
    with io.open(p, encoding='utf-8-sig', newline='') as f:
        rr = list(csv.reader(f))
    return [r for r in rr if r]


def save(name, rows):
    with io.open(os.path.join(STAGE, name), 'w', encoding='utf-8-sig',
                 newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerows(rows)
    CHANGED[name] = len(rows)


def cell(row, old, new, label):
    """在一行的若干列里做唯一替换，未命中或多次命中都报错。"""
    hits = [i for i, c in enumerate(row) if old in c]
    if len(hits) != 1:
        raise SystemExit('!! %s：在行内命中 %d 次（要求 1 次）\n   %r'
                         % (label, len(hits), row))
    i = hits[0]
    row[i] = row[i].replace(old, new)


# ═══════════════════════════════════════════════════════ S13 TRIPOD
S13 = 'Table_S13_TRIPOD清单.csv'
rows = load(S13)
assert rows[0] == ['Section', 'Item', 'Checklist item', 'Location in this manuscript']
for r in rows[1:]:
    sec, item = r[0], r[1]
    if item == '1':
        cell(r, 'score is an externally reproduced prognostic factor',
             'score is an independent, externally reproduced prognostic factor',
             'S13 题名')
    if item == '3b':
        cell(r, 'relative to [12]', 'relative to [17]', 'S13 [12]→[17]')
    if item == '7':
        cell(r, 'Table 2 and Table 7 (coefficients)',
             'Table 2 (MRG-4 coefficients) and Supplementary Table S12 '
             '(Model A coefficients)', 'S13 row7')
    if item == '11a':
        cell(r, 'Table 9, Table 10, Table 6.',
             'Table 3, Table 4 and Supplementary Table S30.', 'S13 row11a')
    if item == '11b':
        cell(r, 'Figure 6C,D', 'Figure 5C,D', 'S13 row11b')
    if item == '11c':
        cell(r, 'Table 6 (stratified log-rank',
             'Supplementary Table S29 (stratified log-rank', 'S13 row11c')
    if item == '23':
        cell(r, 'Supplementary Tables S1-S15; Figures 1-9 (composite, 28 panels)',
             'Supplementary Tables S1\u2013S30; Figures 1\u20138 (composite, '
             '27 panels)', 'S13 row23')
save(S13, rows)

# ═══════════════════════════════════════════════════════ S25 Model A 一致性
S25 = 'Table_S25_ModelA定义一致性.csv'
rows = load(S25)
for r in rows[1:]:
    if r[0] == 'TCGA apparent C-index' and r[1] == 'archived LP':
        cell(r, 'Table 9 现用值 0.821', 'Table 4 现用值 0.821', 'S25 Table 9→4')
rows += [
    ['TCGA apparent C-index', 'MRG-4, recomputed from the deposited per-patient scores',
     '0.7938', 'Table 4 现用值 0.794；可由 Table_S4 的 MRG4_score 直接重算'],
    ['MRG-4 apparent C-index, second internal object',
     'parsimonious_results (not deposited)',
     '0.7926',
     '同一四基因拟合的第二个结果对象报 0.7926，与 mrg4_results 的 0.7943 差 0.0017，'
     '与 S4 重算值 0.7938 差 0.0012；正文统一引用 0.794（Table 4 与 3.3 节），'
     '任何结论都不依赖这 0.0017'],
    ['SD of the Model A score', 'TCGA-LGG complete annotation, n = 502',
     '1.061719', 'Table 3 的 per-SD 分母；由 Table_S2 的 MRG_risk_score 重算，'
     '且 exp(coef × SD) = 2.9266 与 Table 3 未校正行的 2.927 吻合'],
    ['SD of the Model A score', 'TCGA-LGG all, n = 510', '1.057266',
     '同一列在全 510 例上的 SD'],
    ['SD of the MRG-4 score', 'TCGA-LGG, n = 510', '0.903031',
     'Table 3 的 per-SD 分母；与 Table_S17 的 0.903029 一致'],
]
save(S25, rows)

# ═══════════════════════════════════════════════════════ S27 PROBAST
S27 = 'Table_S27_PROBAST自评.csv'
rows = load(S27)
Q = {}
for r in rows[1:]:
    Q[(r[0], r[1])] = r

# ⑨-a  Q4.1：原答案「No for Model A / Yes for MRG-4」把对照模型的 EPV 混进主模型的域评级
r = Q[('4 Analysis', '4.1')]
r[3] = 'Yes'
r[4] = ('MRG-4: 125 deaths for 4 predictors, 31.2 events per variable. This signalling '
        'question is answered for the model the paper reports as its main result. The '
        'Model A comparator has 125 deaths for 25 predictors (5.0 events per variable), '
        'below the conventional minimum, but a comparator that is never used for '
        'prediction does not carry this domain; its EPV is reported with the comparison '
        'in Section 3.2 and Table 4 rather than folded into the rating here.')

# ⑨-b  Q4.5：sign-flip 数在 Table 4，不在 Table 2
r = Q[('4 Analysis', '4.5')]
cell(r, '(Table 2)', '(Table 4)', 'S27 Q4.5 表号')

# ⑨-c  Q4.6：终点是总生存，本无竞争风险；原理由站不住
r = Q[('4 Analysis', '4.6')]
r[3] = 'Yes'
r[4] = ('Censoring is handled by the Cox partial likelihood throughout and by the '
        'concordance, time-dependent AUC (inverse-probability-of-censoring weighted, '
        'censoring estimated in the cohort being evaluated) and RMST estimators in '
        'validation. Competing risks are not a complexity of this endpoint: overall '
        'survival is a single all-cause event, so by construction there is no competing '
        'event. No absolute-risk claim is made from the frozen model.')

# ⑨-d  Q4.7：原文写「calibration 与 clinical utility 未评估」，与正文和 S27 图注相矛盾
r = Q[('4 Analysis', '4.7')]
r[3] = 'Probably yes'
r[4] = ('Discrimination was evaluated thoroughly (apparent, non-nested cross-validated, '
        'fully nested and external C-indices; time-dependent AUC with a patient-level '
        'paired bootstrap; RMST differences). Calibration and clinical utility were also '
        'evaluated: the calibration slope and the tertile calibration plot are in Figure '
        '5C, the integrated Brier score in Figure 5D and Table 4, and net benefit and '
        'decision curves in Figure 5B, against a four-variable clinical model refitted in '
        'the same patients. The one comparison not made is against an externally published '
        'LGG model of these MRGs, for which no comparable published model existed; that is '
        'declared in the Limitations.')

# ⑨-e  Q4.9：四个 MRG-4 系数在 S17，Model A 的 25 个系数在 S12
r = Q[('4 Analysis', '4.9')]
cell(r, 'deposited (Table S12)',
     'deposited (the four MRG-4 coefficients in Supplementary Table S17, the 25 Model A '
     'coefficients in Supplementary Table S12)', 'S27 Q4.9 表号')

# ⑨-f  Domain 4 判语：Q4.1 与 Q4.7 已改判，只剩 Q4.5
r = Q[('4 Analysis', 'JUDGEMENT')]
cell(r, 'Q4.1 No for Model A, Q4.5 No, Q4.7 No.', 'Q4.5 No.', 'S27 域 4 判语')
save(S27, rows)

# ═══════════════════════════════════════════════════════ S28 随机种子
S28 = 'Table_S28_随机种子与软件环境.csv'
rows = load(S28)
for r in rows[1:]:
    body = ' | '.join(r)
    if 'k-means subtypes of the 25 Model A genes' in body:
        cell(r, 'Figure 8A', 'Figure 7A', 'S28 k-means 图号')
    if '200-bootstrap stability selection' in body:
        cell(r, '(Table 7, Table S6)', '(Table 2, Table S6)', 'S28 稳定性表号')
    if 'Randomised-panel benchmark' in body:
        cell(r, '(Table 10)', '(Table 4)', 'S28 随机面板表号')
    if 'Time-dependent AUC and RMST, paired bootstrap' in body:
        cell(r, '(Figures 6C, 6D; Table S23, Table S24)',
             '(Figure 5A; Table S23, Table S24)', 'S28 AUC/RMST 图号')
    if 'Single-cell aggregation' in body:
        cell(r, '(Figures 9A, 9B; Tables S10, S11, S16)',
             '(Figures 8A, 8B; Tables S10, S11, S16)', 'S28 单细胞图号')
save(S28, rows)

# ═══════════════════════════════════════════════════════ S30 不确定性
S30 = 'Table_S30_不确定性.csv'
rows = load(S30)
for r in rows[1:]:
    if r[0].startswith('Non-nested 10-fold CV'):
        cell(r, 'Table 3', 'Table 4', 'S30 Table 3→4')
for r in rows[1:]:
    if r[0].startswith('26-gene LASSO apparent C-index'):
        raise SystemExit('!! S30 已有 26-gene 行，重复')
idx = next(i for i, r in enumerate(rows) if r[0].startswith('Non-nested 10-fold CV'))
rows.insert(idx + 1, [
    '26-gene LASSO apparent C-index, TCGA', '0.870',
    'quoted from the derivation run and printed in Figure 3B; the only value in Table 4 '
    'that cannot be rebuilt from the deposited objects, because the 26-gene coefficient '
    'vector at the 1-SE alpha was not retained. Its cross-validated (0.777) and external '
    '(0.726) counterparts are rebuildable and are the values used in the text'])
idx = next(i for i, r in enumerate(rows) if r[0].startswith('26-gene LASSO apparent'))
rows.insert(idx + 1, [
    'Standard deviation of each score, TCGA-LGG',
    'MRG-4 0.903 (n = 510); Model A 1.062 (n = 502)',
    'the per-SD denominators of every Table 3 hazard ratio; both recomputed from the '
    'deposited per-patient scores (Supplementary Tables S2 and S4); the Model A value '
    'is 1.061719'])
save(S30, rows)

# ═══════════════════════════════════════════════════════ S24 RMST
#  用与 round4_stats.py 相同的估计量（KM 阶梯函数的梯形积分）重算四位小数，
#  断言一位小数与交付件逐格一致，再把四位小数写成新列。
import numpy as np                                       # noqa: E402
from lifelines import KaplanMeierFitter                  # noqa: E402


def rmst(t, e, tau):
    k = KaplanMeierFitter().fit(t, e)
    sf = k.survival_function_
    x = np.concatenate([[0], sf.index.values])
    y = np.concatenate([[1], sf.iloc[:, 0].values])
    keep = x <= tau
    x, y = x[keep], y[keep]
    if x[-1] < tau:
        x = np.append(x, tau)
        y = np.append(y, y[-1])
    return float(np.trapezoid(y, x))


tcga = load('Table_S4_TCGA_LGG_MRG4逐患者评分.csv')
h4 = tcga[0]
grp = {}
for r in tcga[1:]:
    if not r or not r[0].strip():
        continue
    grp.setdefault(r[h4.index('MRG4_risk_group')], []).append(
        (float(r[h4.index('OS_months')]), int(float(r[h4.index('Death')]))))
cgga = load('Table_S5_CGGA325_MRG4逐患者评分.csv')
h5 = cgga[0]
grp2 = {}
for r in cgga[1:]:
    if not r or not r[0].strip():
        continue
    if r[h5.index('Included_in_primary_validation')].strip() != 'True':
        continue
    grp2.setdefault(r[h5.index('MRG4_risk_group_cohort_z')], []).append(
        (float(r[h5.index('OS_months')]), int(float(r[h5.index('Death')]))))

old = load('Table_S24_RMST.csv')
hdr_old = old[0]
assert hdr_old == ['Cohort', 'tau_months', 'RMST_low', 'RMST_high', 'RMST_difference',
                   'diff_lo95', 'diff_hi95'], hdr_old
new_rows = [hdr_old + ['RMST_low_4dp', 'RMST_high_4dp', 'RMST_difference_4dp']]
for r in old[1:]:
    coh, tau = r[0], int(r[1])
    g = grp if coh.startswith('TCGA') else grp2
    lo = rmst([t for t, e in g['low']], [e for t, e in g['low']], tau)
    hi = rmst([t for t, e in g['high']], [e for t, e in g['high']], tau)
    # 断言：重算值的一位小数必须与交付件完全一致，否则说明估计量不同，停下
    for got, shown in ((lo, r[2]), (hi, r[3]), (lo - hi, r[4])):
        if '%.1f' % round(got, 2) != shown:
            raise SystemExit('!! S24 重算 %.4f 与交付值 %s 不符（%s τ=%d）'
                             % (got, shown, coh, tau))
    new_rows.append(r + ['%.4f' % lo, '%.4f' % hi, '%.4f' % (lo - hi)])
save('Table_S24_RMST.csv', new_rows)

print('=== 已 staged ===')
for k, v in sorted(CHANGED.items()):
    print('  %-46s %d 行' % (k, v))
