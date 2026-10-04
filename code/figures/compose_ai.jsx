#target illustrator
// 由 /tmp/mitoxy/make_ai_jsx.py 自动生成 —— 请勿手改。
// 每张面板 PDF 都是整幅画布尺寸，全部放在原点 100% 缩放即可精确复原版面。

var SRC = "/tmp/mitoxy/figs_v5/panels/";
var OUT = "/tmp/mitoxy/figs_v5/ai/";
var FIGURES = [
  {name:"Figure1_发现与ModelA", w:7.4, h:7.7, panels:[
    {L:"A", file:"Figure1_发现与ModelA_A.pdf", lx:0.08839061519146266, ty:0.9047590187590188, fs:10.4},
    {L:"B", file:"Figure1_发现与ModelA_B.pdf", lx:0.44937617702448207, ty:0.9047590187590188, fs:10.4},
    {L:"C", file:"Figure1_发现与ModelA_C.pdf", lx:0.45170044726930314, ty:0.6444192362198439, fs:10.4},
    {L:"D", file:"Figure1_发现与ModelA_D.pdf", lx:0.44937617702448207, ty:0.38407945368066904, fs:10.4},
  ]},
  {name:"Figure2_ModelA与随机面板基准", w:7.4, h:5.5, panels:[
    {L:"A", file:"Figure2_ModelA与随机面板基准_A.pdf", lx:0.07902542372881355, ty:0.9122626262626262, fs:10.4},
    {L:"B", file:"Figure2_ModelA与随机面板基准_B.pdf", lx:0.525635593220339, ty:0.9122626262626262, fs:10.4},
    {L:"C", file:"Figure2_ModelA与随机面板基准_C.pdf", lx:0.07902542372881355, ty:0.4478181818181818, fs:10.4},
    {L:"D", file:"Figure2_ModelA与随机面板基准_D.pdf", lx:0.525635593220339, ty:0.4478181818181818, fs:10.4},
  ]},
  {name:"Figure3_MRG4外部验证", w:7.4, h:7.4, panels:[
    {L:"A", file:"Figure3_MRG4外部验证_A.pdf", lx:0.07979166666666665, ty:0.9055195195195195, fs:10.4},
    {L:"B", file:"Figure3_MRG4外部验证_B.pdf", lx:0.5318750000000001, ty:0.9055195195195195, fs:10.4},
    {L:"C", file:"Figure3_MRG4外部验证_C.pdf", lx:0.08029032258064514, ty:0.5385122955773111, fs:10.4},
  ]},
  {name:"Figure4_MRG4推导与稳定性", w:7.4, h:5.4, panels:[
    {L:"A", file:"Figure4_MRG4推导与稳定性_A.pdf", lx:0.08294573643410852, ty:0.9127489711934156, fs:10.4},
    {L:"B", file:"Figure4_MRG4推导与稳定性_B.pdf", lx:0.5575581395348836, ty:0.9127489711934156, fs:10.4},
    {L:"C", file:"Figure4_MRG4推导与稳定性_C.pdf", lx:0.08294573643410852, ty:0.42377086900363464, fs:10.4},
    {L:"D", file:"Figure4_MRG4推导与稳定性_D.pdf", lx:0.5575581395348836, ty:0.42377086900363464, fs:10.4},
  ]},
  {name:"Figure5_头对头与决策分析", w:7.4, h:5.5, panels:[
    {L:"A", file:"Figure5_头对头与决策分析_A.pdf", lx:0.07902542372881355, ty:0.9122626262626262, fs:10.4},
    {L:"B", file:"Figure5_头对头与决策分析_B.pdf", lx:0.525635593220339, ty:0.9122626262626262, fs:10.4},
    {L:"C", file:"Figure5_头对头与决策分析_C.pdf", lx:0.07902542372881355, ty:0.4407122386657271, fs:10.4},
    {L:"D", file:"Figure5_头对头与决策分析_D.pdf", lx:0.525635593220339, ty:0.4407122386657271, fs:10.4},
  ]},
  {name:"Figure6_分层生存", w:7.4, h:2.6, panels:[
    {L:"A", file:"Figure6_分层生存_A.pdf", lx:0.095, ty:0.9415555555555557, fs:10.4},
    {L:"B", file:"Figure6_分层生存_B.pdf", lx:0.385625, ty:0.9415555555555557, fs:10.4},
    {L:"C", file:"Figure6_分层生存_C.pdf", lx:0.67625, ty:0.9415555555555557, fs:10.4},
  ]},
  {name:"Figure7_分子相关性", w:7.4, h:7.6, panels:[
    {L:"A", file:"Figure7_分子相关性_A.pdf", lx:0.113, ty:0.9412228147272478, fs:10.4},
    {L:"B", file:"Figure7_分子相关性_B.pdf", lx:0.6092230215827338, ty:0.9217192982456143, fs:10.4},
    {L:"C", file:"Figure7_分子相关性_C.pdf", lx:0.039999999999999994, ty:0.5042915622389307, fs:10.4},
  ]},
  {name:"Figure8_单细胞定位", w:7.8, h:4.9, panels:[
    {L:"A", file:"Figure8_单细胞定位_A.pdf", lx:0.13626465107892047, ty:0.9154784580498866, fs:10.4},
    {L:"B", file:"Figure8_单细胞定位_B.pdf", lx:0.6551175091663162, ty:0.9154784580498866, fs:10.4},
  ]},
];

function makeFolder(p) {
  var f = new Folder(p); if (!f.exists) f.create(); return f;
}

function compose(f) {
  // v5.3 关键修正：**不要对点值取整**。
  // 原来写 Math.round(f.w * 72)，7.4 in → 533 pt（本该 532.8），600 dpi 导出后
  // 画布变成 4442 px 而 matplotlib 整幅是 4440 px，宽多 2、高少 3 —— 于是
  // 「Illustrator 拼版件 vs matplotlib 参照件」的逐像素对照永远差几个像素，
  // 分不清是拼版错了还是导出取整。用实数点值后两边尺寸严格相同（4440×4620），
  // 对照才有意义。Illustrator 的 artboardRect 接受小数。
  var W = f.w * 72, H = f.h * 72;
  var doc = app.documents.add(DocumentColorSpace.RGB, W, H);
  doc.artboards[0].artboardRect = [0, 0, W, -H];
  // ---- 放置面板（全部在原点，页面尺寸一致） ----
  for (var i = 0; i < f.panels.length; i++) {
    var p = f.panels[i];
    var pl = doc.placedItems.add();
    pl.file = new File(SRC + p.file);
    pl.width = W; pl.height = H;
    pl.position = [0, 0];
    pl.name = p.L + " — " + p.file;
  }
  // ---- 面板字母：活字文本 ----
  var FONTS = ["Helvetica-Bold", "Arial-BoldMT", "HelveticaNeue-Bold", "ArialMT"];
  var theFont = null;
  for (var m = 0; m < FONTS.length; m++) {
    try { theFont = app.textFonts.getByName(FONTS[m]); break; } catch (e) {}
  }
  var black = new RGBColor(); black.red = 0; black.green = 0; black.blue = 0;
  for (var j = 0; j < f.panels.length; j++) {
    var q = f.panels[j];
    var tf = doc.textFrames.add();
    tf.contents = q.L;
    var ca = tf.textRange.characterAttributes;
    if (theFont) { try { ca.textFont = theFont; } catch (e) {} }
    ca.size = q.fs;
    ca.fillColor = black;   // 中文版 AI 的色板名不是 "Black"，直接给 RGBColor 最稳
    // 坐标换算：matplotlib 图幅比例（原点左下）→ Illustrator（原点左上，y 向下为负）
    //   q.lx  = 文字左沿的图幅比例 x
    //   q.ty  = 文字**墨迹上沿**的图幅比例 y
    //   tf.position 是文本框左上角。Illustrator 把活字排在文本框顶端下面，
    //   框沿到字形墨迹上沿有 0.1461 em 的间距（Helvetica-Bold @10.4 pt 实测），
    //   而 matplotlib 的 bbox 上沿到墨迹上沿只有约 0.0346 em —— 两者之差
    //   0.1115 em 就是这里要补的量。
    //   ⚠️ 原来写的是 0.25 em（凭"约 0.25 em"的估），实测字母整体偏高 12–13 px
    //      （600 dpi 下 1.5 pt），8 张图 27 个字母无一例外，是系统偏差不是抖动。
    //      改 0.1115 em 后与 matplotlib 的落点一致（x 偏移 ≤2 px、y 偏移 ≤2 px）。
    var GAP_EM = 0.1115;
    var x = q.lx * W;
    var yTop = -(1 - q.ty) * H + GAP_EM * q.fs;
    tf.position = [x, yTop];
    tf.name = "letter " + q.L;
  }
  makeFolder(OUT);
  var so = new IllustratorSaveOptions();
  so.pdfCompatible = true;
  so.embedLinkedFiles = true;
  doc.saveAs(new File(OUT + f.name + ".ai"), so);
  // PNG 600 dpi：文档是 72 dpi，故缩放 = 600/72*100
  var po = new ExportOptionsPNG24();
  po.artBoardClipping = true;
  po.antiAliasing = true;
  po.transparency = true;
  po.horizontalScale = 600 / 72 * 100;
  po.verticalScale = 600 / 72 * 100;
  doc.exportFile(new File(OUT + f.name + ".png"), ExportType.PNG24, po);
  var pf = new PDFSaveOptions();
  pf.preserveEditability = true;
  doc.saveAs(new File(OUT + f.name + ".pdf"), pf);
  doc.close(SaveOptions.DONOTSAVECHANGES);
  return "ok";
}

function main() {
  app.userInteractionLevel = UserInteractionLevel.DONTDISPLAYALERTS;
  var log = [];
  for (var i = 0; i < FIGURES.length; i++) {
    try { log.push(FIGURES[i].name + " -> " + compose(FIGURES[i])); }
    catch (e) { log.push(FIGURES[i].name + " -> ERROR: " + e.message + " @line " + e.line); }
  }
  var fo = new File(OUT + "compose_log.txt");
  fo.encoding = "UTF-8"; fo.open("w"); fo.write(log.join("\n")); fo.close();
  return log.join("\n");
}
main();
