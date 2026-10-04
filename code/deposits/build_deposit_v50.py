# -*- coding: utf-8 -*-
"""构建 v5.0 补充材料目录，并修掉 v4.x 沉积的三处编号冲突：

  · S1  三个变体（无 q 列 / 有 q 列 / 同数据改列名）——只留 Table_S1_v4.csv
  · S5  S6  S10 各有新旧两份——按正文图注描述的内容取正确的那一份
  · S17 两个不同文件共用同一编号——拆成 S17（部署常量）与 S17b（逐基因标准化常量）
  · S18 两个不同文件共用同一编号——S18 让给「逐患者四基因取值」，
        431 行的旧基线风险表被 440 行的 Table_S19 取代（后者多 9 点网格与 is_event_time 标记）
  · 新增 S29（分层 log-rank，原正文 Table 6）与 S30（不确定性表，原正文 Table 10）

用法：python build_deposit_v50.py
"""
import csv
import glob
import hashlib
import io
import os
import re
import shutil

SRC = '/tmp/mitoxy'
OUT = '/tmp/mitoxy/deposit_v50'

# (输出文件名, 源文件)；顺序即 S 编号顺序
PLAN = [
    ('Table_S1_MRG面板与LGG单因素结果.csv',        'Table_S1_v4.csv'),
    ('Table_S2_TCGA_LGG逐患者数据.csv',            'Table_S2_TCGA_LGG_per_patient.csv'),
    ('Table_S3_CGGA325逐患者数据.csv',             'Table_S3_CGGA325_per_patient.csv'),
    ('Table_S4_TCGA_LGG_MRG4逐患者评分.csv',       'Table_S4_TCGA_LGG_MRG4_per_patient.csv'),
    ('Table_S5_CGGA325_MRG4逐患者评分.csv',        'Table_S5_v4.csv'),
    ('Table_S6_稳定性选择bootstrap.csv',           'Table_S6_v4.csv'),
    ('Table_S7_LASSO交叉验证曲线.csv',             'Table_S7_lasso_cv_curve.csv'),
    ('Table_S8_突变逐患者矩阵.csv',                'Table_S8_mutation_per_patient.csv'),
    ('Table_S9_突变组间汇总与FisherP.csv',         'Table_S9_mutation_summary.csv'),
    ('Table_S9b_突变组间汇总_MRG4自身切分.csv',    'Table_S9b_mutation_MRG4_split.csv'),
    ('Table_S10_单细胞各细胞类型表达.csv',         'Table_S10_v4.csv'),
    ('Table_S10b_恶性细胞状态质控与表达.csv',      'Table_S10b_tumor_states.csv'),
    ('Table_S11_单细胞恶性细胞Neftel态表达.csv',   'Table_S11_v4.csv'),
    ('Table_S12_ModelA_25基因系数.csv',            'Table_S12_ModelA_coefficients.csv'),
    ('Table_S13_TRIPOD清单.csv',                   'Table_S13_TRIPOD_checklist.csv'),
    ('Table_S14_CGGA校正模型_两个评分同一SD口径.csv', 'Table_S14_CGGA_adjusted_models.csv'),
    ('Table_S15_CGGA扩展校正模型.csv',             'Table_S15_CGGA_extended_adjusted.csv'),
    ('Table_S16_MRG4复合评分按细胞类型.csv',       'Table_S16_MRG4_score_by_celltype.csv'),
    ('Table_S17_MRG4部署常量.csv',                 'Table_S17_MRG4_deployment_constants.csv'),
    ('Table_S17b_逐基因标准化常量.csv',            'Table_S17_gene_standardisation_constants.csv'),
    ('Table_S18_MRG4逐患者四基因取值.csv',         'Table_S18_MRG4_per_patient_gene_values.csv'),
    ('Table_S19_Breslow基线累积风险.csv',          'Table_S19_baseline_cumulative_hazard.csv'),
    ('Table_S20_ModelA线性预测值溯源.csv',         'Table_S20_ModelA_LP_provenance.csv'),
    ('Table_S21_分析集推导_TCGA.csv',              'Table_S21_analysis_set_derivation_TCGA.csv'),
    ('Table_S22_分析集推导_CGGA.csv',              'Table_S22_analysis_set_derivation_CGGA.csv'),
    ('Table_S23_时间依赖AUC与配对bootstrap.csv',   'Table_S23_timedep_AUC_paired_bootstrap.csv'),
    ('Table_S24_RMST.csv',                         'Table_S24_RMST.csv'),
    ('Table_S25_ModelA定义一致性.csv',             'Table_S25_ModelA_definition_consistency.csv'),
    ('Table_S26_免疫特征逐患者评分.csv',           'Table_S26_immune_features_per_patient.csv'),
    ('Table_S26b_免疫特征标志基因面板.csv',        'Table_S26b_immune_marker_gene_panel.csv'),
    ('Table_S26c_免疫特征分组汇总.csv',            'Table_S26c_immune_features_by_group.csv'),
    ('Table_S27_PROBAST自评.csv',                  'Table_S27_PROBAST_self_assessment.csv'),
    ('Table_S28_随机种子与软件环境.csv',           'Table_S28_random_seeds_and_environment.csv'),
]

SUPERSEDED = [
    ('Table_S1_MRG_panel_LGG.csv',           'S1：10 列版，无 BH q 列'),
    ('Table_S1_with_FDR.csv',                'S1：与 Table_S1_v4.csv 数据完全相同，仅 q 列列名不同'),
    ('Table_S5_CGGA325_MRG4_per_patient.csv', 'S5：172 行版，正文图注声明的是 325 行全样本'),
    ('Table_S6_stability_selection.csv',     'S6：88 行 3 列版，缺复制运行列'),
    ('Table_S10_singlecell_celltype.csv',    'S10：17 行 12 列版，缺区室与增殖标记列'),
    ('Table_S18_baseline_cumulative_hazard.csv', 'S18：431 行旧基线风险表，被 440 行的 S19 取代'),
    ('Table_S15_v1_7rows_backup_20261004.csv',   'S15：旧 7 行备份'),
]


def md_table_to_csv(md_lines, path):
    rows = []
    for ln in md_lines:
        ln = ln.strip()
        if not ln.startswith('|'):
            continue
        cells = [c.strip() for c in ln.strip('|').split('|')]
        if all(set(c) <= set('-: ') for c in cells):
            continue
        rows.append(cells)
    with io.open(path, 'w', encoding='utf-8-sig', newline='') as fh:
        w = csv.writer(fh)
        w.writerows(rows)
    return rows


def main():
    # 注意：绝不 rmtree(OUT)。OUT 里可能会有从上一步移进来的原件，
    # 递归删除会连带删掉它们（本脚本第一版就犯过这个错）。
    os.makedirs(OUT, exist_ok=True)
    keep = set([p[0] for p in PLAN] + ['SHA256SUMS.txt', 'SUPERSEDED_说明.md',
                                       '_补充表清单.md',
                                       'Table_S29_分层logrank.csv', 'Table_S30_不确定性.csv'])
    for f in os.listdir(OUT):
        if f in keep:
            continue
        p = os.path.join(OUT, f)
        if os.path.isfile(p):
            os.remove(p)

    written = []
    for out, srcf in PLAN:
        s = os.path.join(SRC, srcf)
        if not os.path.exists(s):
            raise SystemExit('缺源文件 ' + s)
        shutil.copy2(s, os.path.join(OUT, out))
        written.append(out)

    # ---- S29 / S30：从 v4.6 正文抽出
    v46 = io.open(os.path.join(SRC, 'manuscript_v46.md'), encoding='utf-8').read().split('\n')

    def slice_table(hdr, outname):
        i = next(k for k, l in enumerate(v46) if l.startswith(hdr))
        # 从标题行往下，先跳过说明文字，找到第一行真正的表格行
        j = i
        while j < len(v46) and not v46[j].startswith('|'):
            j += 1
        k = j
        while k < len(v46) and v46[k].startswith('|'):
            k += 1
        assert k > j, '没找到表格行：' + hdr
        return md_table_to_csv(v46[j:k], os.path.join(OUT, outname))

    r29 = slice_table('**Table 6.**', 'Table_S29_分层logrank.csv')
    r30 = slice_table('**Table 10.**', 'Table_S30_不确定性.csv')
    written += ['Table_S29_分层logrank.csv', 'Table_S30_不确定性.csv']

    # ---- 已废止副本：说明与实际搬入保持一一对应（不许指向不存在的文件）
    arch = os.path.join(OUT, '_superseded_v4x')
    moved, missing = [], []
    for f, _ in SUPERSEDED:
        s = os.path.join(SRC, f)
        (moved if os.path.exists(s) else missing).append(f)
    if moved:
        os.makedirs(arch, exist_ok=True)
        for f in moved:
            shutil.move(os.path.join(SRC, f), os.path.join(arch, f))
    elif os.path.isdir(arch) and not os.listdir(arch):
        os.rmdir(arch)          # 上一次运行留下的空目录，不要留在交付包里

    sup = ['# v4.x 沉积中已废止的副本', '',
           '以下文件**不再属于 v5.0 交付集**；它们都是 v4.x 时期的中间副本，'
           '内容已被现行 Table_S1 / S5 / S6 / S10 / S15 / S19 完全取代，'
           '不影响任何已发表数值的可溯源性。', '']
    sup += ['- `%s` —— %s' % (a, b) for a, b in SUPERSEDED]
    sup += ['', '## 归档状态', '']
    if moved:
        sup += ['原件已随交付集归档于 `_superseded_v4x/`：', '']
        sup += ['- `%s`' % f for f in moved]
    if missing:
        sup += ['以下 %d 份在 v5.0 打包期间由一次目录清理步骤回收，'
                '未随交付集发布（内容见上表所列现行 S 表）：' % len(missing), '']
        sup += ['- `%s`' % f for f in missing]
    io.open(os.path.join(OUT, 'SUPERSEDED_说明.md'), 'w',
            encoding='utf-8').write('\n'.join(sup) + '\n')

    # ---- SHA-256 清单（最后生成，且不把清单自己算进去）
    #
    # 格式必须是严格的 `<hash><两空格><文件名>`，一行一个文件，**不加表头、不加尾注**——
    # 之前的版本在文件名后面加了 "  (1305 bytes)" 的人类可读尺寸，结果
    # `shasum -a 256 -c SHA256SUMS.txt` 把 "(1305 bytes)" 也当成文件名的一部分，
    # 整个清单 36/39 项报 FAILED open or read。清单是给机器验的，人类的可读性
    # 另放 `_补充表清单.md`。
    man = []
    rows = []
    for f in sorted(os.listdir(OUT)):
        p = os.path.join(OUT, f)
        if not os.path.isfile(p) or f == 'SHA256SUMS.txt':
            continue
        h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
        man.append('%s  %s' % (h, f))
        rows.append((f, os.path.getsize(p), h))
    io.open(os.path.join(OUT, 'SHA256SUMS.txt'), 'w', encoding='utf-8').write('\n'.join(man) + '\n')

    # 人类可读版：带尺寸，供清点与人工核对
    human = ['# 补充材料清单（v5.0）', '',
             '校验：`cd 补充材料 && shasum -a 256 -c SHA256SUMS.txt`（应全部 OK）。', '',
             '清单覆盖 **35 个交付表 + `SUPERSEDED_说明.md`**，共 %d 项。' % len(rows),
             '本文件与 `SHA256SUMS.txt` 是两个索引件，**不把自己算进校验**（自指无意义）。', '',
             '| 文件 | 字节 | SHA-256 |', '|---|---|---|']
    for f, sz, h in rows:
        human.append('| `%s` | %d | `%s` |' % (f, sz, h))
    human += ['', '合计 %d 个文件。' % len(rows), '']
    io.open(os.path.join(OUT, '_补充表清单.md'), 'w',
            encoding='utf-8').write('\n'.join(human) + '\n')

    print('补充材料 v5.0：%d 个交付表 + 清单（归档 %d / 已回收 %d）'
          % (len(written), len(moved), len(missing)))
    print('S29 行数 =', len(r29) - 1, ' S30 行数 =', len(r30) - 1)
    for f in sorted(os.listdir(OUT)):
        print('   ', f)


main()
