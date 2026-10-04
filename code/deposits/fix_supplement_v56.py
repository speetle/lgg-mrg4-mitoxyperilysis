# -*- coding: utf-8 -*-
"""补充材料（沉积件）v5.6 归位 —— 本轮**无沉积修订**，故只点货、索引、自检。

为什么这一版不像 v5.5 那样搬移归档：sir 2026-10-05 的四项答复**只涉及元数据与声明**
（作者块、代码可用性措辞、目标刊、另两篇只读扫描），**沉积件一个字未改**。
既无修订，就无「被修订件」需要搬进 `历史版本_v5.5/`；若照搬 v5.5 的搬移逻辑，
会在同名同内容时凭空产生 `_旧2` 一类中间件，把归档弄脏。

三步：
  1. 点货 —— 断言 35 张沉积表齐备；
  2. 重生成两份索引件（`_补充表清单.md`、`SHA256SUMS.txt`），幂等；
  3. `shasum -a 256 -c` 自检。

用法：python fix_supplement_v56.py
"""
import hashlib
import io
import os
import subprocess
import sys

DST = ('/Users/lianbin/Desktop/Mitoxyperilysis课题总包/04-学生论文框架/'
       '学生1_LGG_交付_20261004')
SUPP = os.path.join(DST, '补充材料')
IDX = '_补充表清单.md'


def sha256(p):
    """⚠️ 索引件里必须是真 SHA-256（64 位十六进制），不是 MD5：
    写 32 位进去会让 `shasum -a 256 -c` 直接报
    "no properly formatted SHA checksum lines found"。"""
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


print('=' * 74)
print('补充材料归位（v5.6）—— 本轮无沉积修订，只点货 + 索引 + 自检')
print('=' * 74)

tables = sorted(f for f in os.listdir(SUPP) if f.endswith('.csv'))
assert tables, '补充材料里一个 .csv 都没有？'
assert len(tables) == 35, (
    '沉积表应是 35 张（S1–S30 + S9b/S10b/S17b/S26b/S26c 分表），实际 %d 张 —— '
    '少一张说明装配没装全，多一张说明混进了非沉积件' % len(tables))
for f in tables:
    if not f.startswith('Table_S'):
        raise SystemExit('出现非 Table_S 开头的沉积件：%s' % f)
print('\n[1] 点货：沉积表 %d 张，全部以 Table_S 开头' % len(tables))

# ------------------------------------------------------------------ 2. 索引件
_buf = io.StringIO()
_buf.write('# 补充材料清单\n\n')
_buf.write('校验：`cd 补充材料 && shasum -a 256 -c SHA256SUMS.txt`（应全部 OK）。\n\n')
_buf.write('清单覆盖 **%d 个沉积表**，共 %d 项。\n' % (len(tables), len(tables) + 1))
_buf.write('本文件与 `SHA256SUMS.txt` 是两个索引件，**不把自己算进校验**（自指无意义）。\n\n')
_buf.write('| 文件 | 字节 | SHA-256 |\n|---|---|---|\n')
for fn in tables:
    p = os.path.join(SUPP, fn)
    _buf.write('| `%s` | %d | `%s` |\n' % (fn, os.path.getsize(p), sha256(p)))
IDX_TXT = _buf.getvalue()

_p = os.path.join(SUPP, IDX)
_old = io.open(_p, encoding='utf-8').read() if os.path.exists(_p) else None
if _old == IDX_TXT:
    print('[2] %s 已是最新，跳过（幂等）' % IDX)
else:
    io.open(_p, 'w', encoding='utf-8').write(IDX_TXT)
    print('[2] 已重写 %s' % IDX)

lst = tables + [IDX]
_buf2 = io.StringIO()
for fn in lst:
    _buf2.write('%s  %s\n' % (sha256(os.path.join(SUPP, fn)), fn))
SUM_TXT = _buf2.getvalue()
_ps = os.path.join(SUPP, 'SHA256SUMS.txt')
if os.path.exists(_ps) and io.open(_ps, encoding='utf-8').read() == SUM_TXT:
    print('    SHA256SUMS.txt 已是最新，跳过（幂等）')
else:
    io.open(_ps, 'w', encoding='utf-8').write(SUM_TXT)
    print('    已重写 SHA256SUMS.txt（%d 项）' % len(lst))

# ------------------------------------------------------------------ 3. 自检
r = subprocess.run(['shasum', '-a', '256', '-c', 'SHA256SUMS.txt'],
                   cwd=SUPP, capture_output=True, text=True)
bad = [l for l in (r.stdout or '').splitlines() if not l.endswith(': OK')]
print('[3] shasum -c：%d 行，失败 %d 行 %s'
      % (len((r.stdout or '').splitlines()), len(bad), bad[:3]))
if bad:
    raise SystemExit('!! 校验清单与实际不符')
left = [f for f in os.listdir(SUPP) if f.endswith('.md') and f != IDX]
print('    补充材料里的 .md（应只剩索引件）：%s' % left)
print('    文件总数：%d（%d 表 + 2 索引 + %d 其他）'
      % (len(os.listdir(SUPP)), len(tables), len(os.listdir(SUPP)) - len(tables) - 2))
print('\n完成（零删除、零搬移：本轮无沉积修订）')
sys.exit(0)
