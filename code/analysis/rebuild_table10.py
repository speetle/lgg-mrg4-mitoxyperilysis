# -*- coding: utf-8 -*-
"""从 manuscript_v43.md 的 Table 10 段落重建 Table_10_v4.csv，保证沉积表与正文逐行一致。"""
import csv, re, io

MS = '/tmp/mitoxy/manuscript_v43.md'
OUT = '/tmp/mitoxy/Table_10_v4.csv'

txt = open(MS, encoding='utf-8').read()
m10 = txt.split('**Table 10.**')[1].split('**Supplementary Table S1.**')[0]

rows = []
for line in m10.split('\n'):
    s = line.strip()
    if not s.startswith('|'):
        continue
    cells = [c.strip() for c in s.strip('|').split('|')]
    if set(''.join(cells)) <= set('-: '):
        continue
    rows.append(cells)

assert rows, 'no table rows parsed'
hdr = rows[0]
body = rows[1:]
assert all(len(r) == len(hdr) for r in body), [r for r in body if len(r) != len(hdr)]

with open(OUT, 'w', encoding='utf-8-sig', newline='') as f:
    w = csv.writer(f)
    w.writerow(hdr)
    w.writerows(body)

print('header:', hdr)
print('rows written:', len(body))
for r in body[:3]:
    print('  sample:', r[0][:60], '|', r[1][:40])
print()
# 关键量是否落到沉积表里
flat = '\n'.join(','.join(r) for r in body)
must = ['0.722', '0.812', '0.68%', '2,287', '97.8%', '348', '49.3%', '15.2', '0.5713',
        'n = 160', '410', '+0.440', 'n = 156']
missing = [m for m in must if m not in flat]
print('关键量缺失:', missing if missing else '无')
# 不该再出现的旧值
stale = ['v4.0', 'n = 411', '0.639', '75% of', '1.36', 'LR chi2 = 39.3', 'matrix score']
present = [s for s in stale if s in flat]
print('残留旧值:', present if present else '无')
