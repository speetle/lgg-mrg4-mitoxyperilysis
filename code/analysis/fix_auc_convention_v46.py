# -*- coding: utf-8 -*-
"""送审前第三次修补：把时间依赖 AUC 的"平均"口径写清。

第四轮正文写 "Time-dependent AUC averaged 0.819 ... versus 0.711 ..."，但 0.819/0.711
是 sksurv cumulative_dynamic_auc 返回的**事件加权**均值；同一 19 点数组的算术均值是
0.791 / 0.682。评审人只要对 Table S23 的列取一次平均就会得到 0.791，必须自行披露口径，
否则看起来像又一处对不上。
"""
import re, hashlib

P = 'manuscript_v46.md'
t = open(P, encoding='utf-8').read()
h0 = hashlib.md5(t.encode()).hexdigest()

E = [
    # Results 正文
    ('Time-dependent AUC in CGGA-325 averaged 0.819 for MRG-4 versus 0.711 for Model A, and '
     'MRG-4 was higher at every one of the 19 time points (24 months: 0.875 vs 0.759; 120 months: '
     '0.695 vs 0.618) (Figure 6A).',
     'Time-dependent AUC in CGGA-325, averaged over the 19 horizons from 12 to 120 months with the '
     'event-weighted scheme of the estimator, was 0.819 for MRG-4 versus 0.711 for Model A, and MRG-4 '
     'was higher at every one of the 19 time points (24 months: 0.875 vs 0.759; 120 months: 0.695 vs '
     '0.618) (Figure 6A). The convention matters for the reading of a single number and is therefore '
     'stated: the unweighted mean of the same 19 values is 0.791 against 0.682, and no comparison in '
     'this paper depends on which of the two is quoted.'),

    # Table 9 行
    ('| Time-dependent AUC, mean (12–120 mo), cohort z | 0.711 | — | **0.819** |',
     '| Time-dependent AUC, event-weighted mean over the 19 horizons of 12–120 mo, cohort z '
     '(unweighted mean 0.682 vs 0.791) | 0.711 | — | **0.819** |'),

    # Figure 6A 图注
    ('mean 0.819 for MRG-4 versus 0.711 for the 25-gene model; MRG-4 is higher at every one of the 19 '
     'time points',
     'event-weighted mean 0.819 for MRG-4 versus 0.711 for the 25-gene model (unweighted means 0.791 '
     'and 0.682); MRG-4 is higher at every one of the 19 time points'),
]

for i, (old, new) in enumerate(E, 1):
    c = t.count(old)
    assert c == 1, 'E%d 锚点命中 %d 次' % (i, c)
    t = t.replace(old, new)
    print('E%d 已替换' % i)

open(P, 'w', encoding='utf-8').write(t)
print('\nmd5 %s -> %s' % (h0, hashlib.md5(t.encode()).hexdigest()))
print('行数 %d | 字节 %d | 文献 %d | 占位符 %d'
      % (t.count('\n') + 1, len(t.encode()),
         len(re.findall(r'^\s*\d+\.\s.+$', t[t.find('## References'):], re.M)),
         t.count('{{CITE:')))
