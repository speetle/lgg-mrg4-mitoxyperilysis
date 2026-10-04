# -*- coding: utf-8 -*-
"""由 layout.json 生成 Adobe Illustrator 的 ExtendScript（.jsx）。

设计要点
--------
* make_figs_v5.py 导出的每一张面板 PDF 都是**整幅画布尺寸**（例如 7.4 × 7.7 in 的一页），
  因此 JSX 里只要把每张都放在原点、100% 缩放，拼出来的版面就与 matplotlib 整图逐像素一致，
  不需要任何坐标换算。
* 面板字母用 **活字文本**（不是位图），位置由 layout.json 里的图幅比例换算成
  Illustrator 的左上原点坐标。字体 Helvetica-Bold，字号 = matplotlib 的 pts。
* 每张图输出 .ai（可继续手改）+ PNG(600 dpi) + PDF。
"""
import json
import os
import re

LAYOUT = json.load(open('/tmp/mitoxy/figs_v5/layout.json', encoding='utf-8'))
SRC = '/tmp/mitoxy/figs_v5/panels'
OUT = '/tmp/mitoxy/figs_v5/ai'

# v5.0（sir 2026-10-04）：**取消补充图**。原 FigS1 / FigS2 升为正文图。
# v5.3（sir 2026-10-04）：**图号按正文首次引用顺序重排**（旧→新 1→1 2→2 4→3 3→4
#   5→5 7→6 8→7 6→8），与 make_figs_v5.py 的 PLAN 保持一致。
# 交付集 = 8 张正文图，无补充图。
#
# ⚠️ ORDER / LETTER_ORDER 一律**从 layout.json 推导**，不再硬编码：
#   图号一旦重排（v5.3 就重排过一次），任何手抄的名字表都会立刻变成悬挂引用，
#   而报错只会在 Illustrator 跑的时候才出现。推导则永远不会错位。
ORDER = sorted(LAYOUT, key=lambda n: int(re.match(r'Figure(\d+)', n).group(1)))
LETTER_ORDER = {n: sorted(LAYOUT[n]['panels']) for n in ORDER}
_expected = ['Figure%d' % i for i in range(1, 9)]
if [re.match(r'Figure\d+', n).group(0) for n in ORDER] != _expected:
    raise SystemExit('!! layout.json 的图名不是 Figure1..Figure8：%s' % ORDER)


def esc(s):
    return s.replace('\\', '\\\\').replace('"', '\\"')


figs = []
for name in ORDER:
    rec = LAYOUT[name]
    w, h = rec['figsize_in']
    panels = []
    for L in LETTER_ORDER[name]:
        p = rec['panels'][L]
        panels.append({'L': L, 'file': '%s_%s.pdf' % (name, L),
                       'lx': p.get('letter_x', p['x0'] - 0.04),
                       'ty': p.get('letter_ytop', p['y1'] + 0.03),
                       'fs': p.get('letter_fs', 10.4)})
    figs.append({'name': name, 'w': w, 'h': h, 'panels': panels})

lines = []
A = lines.append
A('#target illustrator')
A('// 由 /tmp/mitoxy/make_ai_jsx.py 自动生成 —— 请勿手改。')
A('// 每张面板 PDF 都是整幅画布尺寸，全部放在原点 100% 缩放即可精确复原版面。')
A('')
A('var SRC = "%s/";' % esc(SRC))
A('var OUT = "%s/";' % esc(OUT))
A('var FIGURES = [')
for f in figs:
    A('  {name:"%s", w:%s, h:%s, panels:[' % (esc(f['name']), f['w'], f['h']))
    for p in f['panels']:
        A('    {L:"%s", file:"%s", lx:%s, ty:%s, fs:%s},'
          % (p['L'], esc(p['file']), p['lx'], p['ty'], p['fs']))
    A('  ]},')
A('];')
A('')
A('function makeFolder(p) {')
A('  var f = new Folder(p); if (!f.exists) f.create(); return f;')
A('}')
A('')
A('function compose(f) {')
A('  // v5.3 关键修正：**不要对点值取整**。')
A('  // 原来写 Math.round(f.w * 72)，7.4 in → 533 pt（本该 532.8），600 dpi 导出后')
A('  // 画布变成 4442 px 而 matplotlib 整幅是 4440 px，宽多 2、高少 3 —— 于是')
A('  // 「Illustrator 拼版件 vs matplotlib 参照件」的逐像素对照永远差几个像素，')
A('  // 分不清是拼版错了还是导出取整。用实数点值后两边尺寸严格相同（4440×4620），')
A('  // 对照才有意义。Illustrator 的 artboardRect 接受小数。')
A('  var W = f.w * 72, H = f.h * 72;')
A('  var doc = app.documents.add(DocumentColorSpace.RGB, W, H);')
A('  doc.artboards[0].artboardRect = [0, 0, W, -H];')
A('  // ---- 放置面板（全部在原点，页面尺寸一致） ----')
A('  for (var i = 0; i < f.panels.length; i++) {')
A('    var p = f.panels[i];')
A('    var pl = doc.placedItems.add();')
A('    pl.file = new File(SRC + p.file);')
A('    pl.width = W; pl.height = H;')
A('    pl.position = [0, 0];')
A('    pl.name = p.L + " — " + p.file;')
A('  }')
A('  // ---- 面板字母：活字文本 ----')
A('  var FONTS = ["Helvetica-Bold", "Arial-BoldMT", "HelveticaNeue-Bold", "ArialMT"];')
A('  var theFont = null;')
A('  for (var m = 0; m < FONTS.length; m++) {')
A('    try { theFont = app.textFonts.getByName(FONTS[m]); break; } catch (e) {}')
A('  }')
A('  var black = new RGBColor(); black.red = 0; black.green = 0; black.blue = 0;')
A('  for (var j = 0; j < f.panels.length; j++) {')
A('    var q = f.panels[j];')
A('    var tf = doc.textFrames.add();')
A('    tf.contents = q.L;')
A('    var ca = tf.textRange.characterAttributes;')
A('    if (theFont) { try { ca.textFont = theFont; } catch (e) {} }')
A('    ca.size = q.fs;')
A('    ca.fillColor = black;   // 中文版 AI 的色板名不是 "Black"，直接给 RGBColor 最稳')
A('    // 坐标换算：matplotlib 图幅比例（原点左下）→ Illustrator（原点左上，y 向下为负）')
A('    //   q.lx  = 文字左沿的图幅比例 x')
A('    //   q.ty  = 文字**墨迹上沿**的图幅比例 y')
A('    //   tf.position 是文本框左上角。Illustrator 把活字排在文本框顶端下面，')
A('    //   框沿到字形墨迹上沿有 0.1461 em 的间距（Helvetica-Bold @10.4 pt 实测），')
A('    //   而 matplotlib 的 bbox 上沿到墨迹上沿只有约 0.0346 em —— 两者之差')
A('    //   0.1115 em 就是这里要补的量。')
A('    //   ⚠️ 原来写的是 0.25 em（凭"约 0.25 em"的估），实测字母整体偏高 12–13 px')
A('    //      （600 dpi 下 1.5 pt），8 张图 27 个字母无一例外，是系统偏差不是抖动。')
A('    //      改 0.1115 em 后与 matplotlib 的落点一致（x 偏移 ≤2 px、y 偏移 ≤2 px）。')
A('    var GAP_EM = 0.1115;')
A('    var x = q.lx * W;')
A('    var yTop = -(1 - q.ty) * H + GAP_EM * q.fs;')
A('    tf.position = [x, yTop];')
A('    tf.name = "letter " + q.L;')
A('  }')
A('  makeFolder(OUT);')
A('  var so = new IllustratorSaveOptions();')
A('  so.pdfCompatible = true;')
A('  so.embedLinkedFiles = true;')
A('  doc.saveAs(new File(OUT + f.name + ".ai"), so);')
A('  // PNG 600 dpi：文档是 72 dpi，故缩放 = 600/72*100')
A('  var po = new ExportOptionsPNG24();')
A('  po.artBoardClipping = true;')
A('  po.antiAliasing = true;')
A('  po.transparency = true;')
A('  po.horizontalScale = 600 / 72 * 100;')
A('  po.verticalScale = 600 / 72 * 100;')
A('  doc.exportFile(new File(OUT + f.name + ".png"), ExportType.PNG24, po);')
A('  var pf = new PDFSaveOptions();')
A('  pf.preserveEditability = true;')
A('  doc.saveAs(new File(OUT + f.name + ".pdf"), pf);')
A('  doc.close(SaveOptions.DONOTSAVECHANGES);')
A('  return "ok";')
A('}')
A('')
A('function main() {')
A('  app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;')
A('  var log = [];')
A('  for (var i = 0; i < FIGURES.length; i++) {')
A('    try { log.push(FIGURES[i].name + " -> " + compose(FIGURES[i])); }')
A('    catch (e) { log.push(FIGURES[i].name + " -> ERROR: " + e.message + " @line " + e.line); }')
A('  }')
A('  var fo = new File(OUT + "compose_log.txt");')
A('  fo.encoding = "UTF-8"; fo.open("w"); fo.write(log.join("\\n")); fo.close();')
A('  return log.join("\\n");')
A('}')
A('main();')

os.makedirs(OUT, exist_ok=True)
with open('/tmp/mitoxy/compose_ai.jsx', 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(lines) + '\n')
print('wrote /tmp/mitoxy/compose_ai.jsx  (%d figures, %d lines)'
      % (len(figs), len(lines)))
for f in figs:
    print('  %-34s %.1f x %.1f in  panels=%s' % (f['name'], f['w'], f['h'],
                                                  ','.join(p['L'] for p in f['panels'])))
