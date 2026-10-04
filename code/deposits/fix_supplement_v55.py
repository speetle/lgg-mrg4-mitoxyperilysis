# -*- coding: utf-8 -*-
"""补充材料（沉积件）归位：把 v5.5 改过的六张沉积件装入。

sir 2026-10-04「根据你的建议来做，选项是 abcaa」（八项小修）
--------------------------------------------------------------------------
三件事：

1. **装入改过的沉积件**。本轮改的是**六张**：
   · S13（TRIPOD 清单）：表号 Table 9/10/6 → Table 3/4/S30、图号 6C,D → 5C,D、
     `Figure 1-9 (28 panels)` → `Figures 1–8 (27 panels)`、作者序引用 [12] → [17]；
   · S24（RMST）：补三位小数落四位（`RMST_low_4dp` 等三列）+ 取整约定说明；
   · S25（Model A 定义一致性）：补 0.7938 重算值、`parsimonious_results` 0.7926 与
     `mrg4_results` 0.7943 两个结果对象、三个 SD（含 Model A 的 1.061719）；
   · S27（PROBAST）：Q4.1 → Yes、Q4.6 → Yes、Q4.7 → Probably yes、
     Q4.5 表号、Q4.9 出处点名 S17/S12，Domain 4 判断收敛为「Q4.5 No.」；
   · S28（随机种子与软件环境）：五处图号/表号引用改现行号（8A→7A 等）；
   · S30（不确定性）：Non-nested CV 行的表号 → Table 4，补 0.870 出处行与两倍 SD 行。

   **正文与沉积表是同一句话**，只改正文不改沉积件就会分叉 —— 而审稿人手里拿的是
   沉积件。旧件先搬进 `历史版本_v5.4/补充材料_被修订/`，逐文件 MD5 前后比对，
   不删除。

2. **确认 `SUPERSEDED_说明.md` 不在投稿包里**（上一轮已搬出，本步幂等跳过）。
   它写的是「v4.x 沉积中已废止的副本」，是发布管理笔记，不是科学内容；投给期刊的
   补充材料里不该出现中文的版本回收说明。搬进 `历史版本_v5.4/` —— 它记录的正是
   「被废止的副本」，归档目录才是它的位置。**搬移不是删除。**

3. **重生成两份索引件**：`_补充表清单.md` 与 `SHA256SUMS.txt`。
   两份都按磁盘现扫重算，保证与实物一致。沿用原约定：两个索引件不把自己算进校验
   （自指无意义）。⚠️ SHA-256 必须是**真 64 位**，写 md5 进去会让
   `shasum -a 256 -c` 直接报 "no properly formatted SHA checksum lines found"。

用法：python fix_supplement_v55.py
"""
import hashlib
import io
import os
import shutil
import sys

WD = '/tmp/mitoxy'
STAGE = os.path.join(WD, '_stage55')
DST = ('/Users/lianbin/Desktop/Mitoxyperilysis课题总包/04-学生论文框架/'
       '学生1_LGG_交付_20261004')
SUPP = os.path.join(DST, '补充材料')
ARCH = os.path.join(DST, '历史版本_v5.4')
ARCH_SUPP = os.path.join(ARCH, '补充材料_被修订')

LOG = []


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def sha256(p):
    """⚠️ 索引件里必须是真 SHA-256（64 位十六进制）。
    第一次写这个脚本时误把 md5（32 位）填进了 SHA256SUMS.txt，
    `shasum -a 256 -c` 直接报 "no properly formatted SHA checksum lines found"。"""
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def move(src, dstdir):
    os.makedirs(dstdir, exist_ok=True)
    dst = os.path.join(dstdir, os.path.basename(src))
    if os.path.exists(dst):
        stem, ext = os.path.splitext(os.path.basename(src))
        n = 2
        while os.path.exists(os.path.join(dstdir, '%s_旧%d%s' % (stem, n, ext))):
            n += 1
        dst = os.path.join(dstdir, '%s_旧%d%s' % (stem, n, ext))
    b = md5(src)
    shutil.move(src, dst)
    a = md5(dst)
    if a != b:
        raise SystemExit('MD5 不一致（搬移失败）：%s' % src)
    LOG.append(('move', src, dst, os.path.getsize(dst), b))
    return dst


def copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    b = md5(src)
    if os.path.exists(dst) and md5(dst) == b:
        LOG.append(('same', src, dst, os.path.getsize(dst), b))
        return dst
    shutil.copy2(src, dst)
    a = md5(dst)
    if a != b:
        raise SystemExit('MD5 不一致（复制失败）：%s' % src)
    LOG.append(('copy', src, dst, os.path.getsize(dst), b))
    return dst


print('=' * 74)
print('补充材料归位（v5.5）')
print('=' * 74)

# ---------------------------------------------------------------- 1. 装入改过的沉积件
print('\n[1] 装入改过的沉积件（旧件搬进 %s）'
      % os.path.relpath(ARCH_SUPP, DST))
for f in sorted(os.listdir(STAGE)):
    if not f.endswith('.csv'):
        continue
    src = os.path.join(STAGE, f)
    dst = os.path.join(SUPP, f)
    if not os.path.exists(dst):
        raise SystemExit('沉积件不在补充材料里：%s' % dst)
    if md5(dst) == md5(src):
        print('    置同跳过 %s' % f)
        continue
    move(dst, ARCH_SUPP)
    copy(src, dst)
    print('    已更新 %s  （旧件已归档、MD5 前后一致）' % f)

# ---------------------------------------------------------------- 2. SUPERSEDED 归位
print('\n[2] SUPERSEDED_说明.md 移出投稿包')
p = os.path.join(SUPP, 'SUPERSEDED_说明.md')
if os.path.isfile(p):
    tgt = move(p, ARCH)
    print('    → %s' % os.path.relpath(tgt, DST))
else:
    print('    不在补充材料里（已归位过），跳过')

# ---------------------------------------------------------------- 3. 重生成两份索引件
print('\n[3] 重生成索引件（按磁盘现扫）')
tables = sorted(f for f in os.listdir(SUPP) if f.endswith('.csv'))
assert tables, '补充材料里一个 .csv 都没有？'
assert len(tables) == 35, (
    '沉积表应是 35 张（S1–S30 + S9b/S10b/S17b/S26b/S26c 分表），实际 %d 张 —— '
    '少一张说明装配没装全，多一张说明混进了非沉积件' % len(tables))
for f in tables:
    if not f.startswith('Table_S'):
        raise SystemExit('出现非 Table_S 开头的沉积件：%s' % f)
print('    沉积表 %d 个' % len(tables))

IDX = '_补充表清单.md'
# 索引件是**可再生的派生件**，不该每跑一次就归档一个副本。
# 故先算出新内容，与现盘比对：内容相同就跳过（幂等），不同才搬旧件再写新的。
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
    print('     %s 已是最新，跳过（幂等）' % IDX)
else:
    if _old is not None:
        move(_p, ARCH_SUPP)
    io.open(_p, 'w', encoding='utf-8').write(IDX_TXT)

lst = tables + [IDX]
_buf2 = io.StringIO()
for fn in lst:
    _buf2.write('%s  %s\n' % (sha256(os.path.join(SUPP, fn)), fn))
SUM_TXT = _buf2.getvalue()
_ps = os.path.join(SUPP, 'SHA256SUMS.txt')
if os.path.exists(_ps) and io.open(_ps, encoding='utf-8').read() == SUM_TXT:
    print('     SHA256SUMS.txt 已是最新，跳过（幂等）')
else:
    io.open(_ps, 'w', encoding='utf-8').write(SUM_TXT)
print('    已重写 %s（%d 表）与 SHA256SUMS.txt（%d 项）'
      % (IDX, len(tables), len(lst)))

# ---------------------------------------------------------------- 4. 自检
print('\n[4] 自检')
import subprocess                                                   # noqa: E402
r = subprocess.run(['shasum', '-a', '256', '-c', 'SHA256SUMS.txt'],
                   cwd=SUPP, capture_output=True, text=True)
bad = [l for l in (r.stdout or '').splitlines() if not l.endswith(': OK')]
print('    shasum -c：%d 行，失败 %d 行 %s'
      % (len((r.stdout or '').splitlines()), len(bad), bad[:3]))
if bad:
    raise SystemExit('!! 校验清单与实际不符')
left = [f for f in os.listdir(SUPP)
        if f.endswith('.md') and f not in (IDX,)]
print('    补充材料里的 .md（应只剩索引件）：%s' % left)
print('    文件总数：%d（%d 表 + 2 索引 + %d 其他）'
      % (len(os.listdir(SUPP)), len(tables), len(os.listdir(SUPP)) - len(tables) - 2))

print('\n操作明细：')
for op, s, d, sz, h in LOG:
    print('  %-5s %-46s → %-46s %8d  %s'
          % (op, os.path.basename(s), os.path.relpath(d, DST), sz, h[:12]))
print('\n完成（零删除）')
sys.exit(0)
