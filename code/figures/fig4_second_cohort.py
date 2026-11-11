#!/usr/bin/env python3
"""Figure 4: the second, non-overlapping CGGA cohort.

Panels
  A  external C-index with patient-bootstrap 95% CI, both cohorts, both indices
  B  frozen-model calibration by risk tertile in CGGA-693, 36 and 60 months
  C  C-index difference (MRG-4 minus Model A) within WHO grade II and III

Inputs
  --summary    replication_summary.json from replication_cgga693.py
  --transport  JSON record from transport_diagnostics.py
  --s35        Table_S35_calibration_CGGA693.csv

Layout contract (enforced by assertions before the file is written)
  1. the three panels never overlap vertically: no title, tick label or axis
     label of one panel may enter the column of its neighbour;
  2. the panel letter (A/B/C) occupies its own line, above the panel title;
  3. the three panel titles are top-aligned on one horizontal line;
  4. no text is clipped by the canvas edge.

If any assertion fails the script raises instead of writing a figure.

Usage:
    python fig4_second_cohort.py --summary ./replication/replication_summary.json \
        --transport ./replication/transport.json --out ./Figure4.png
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# C-index of the two indices in CGGA-325, deposited in Table_S5 / Table_S30.
C325 = {"mrg4": (0.7688, 0.722, 0.812), "modelA": (0.664, 0.610, 0.715)}

TITLES = ["Discrimination, two cohorts", "Calibration, CGGA-693", "Model difference by grade"]
LETTERS = ["A", "B", "C"]
LETTER_FS, TITLE_FS_MAX, TITLE_FS_MIN = 9.0, 8.0, 6.2
Y_LETTER_TOP, Y_TITLE_TOP = 0.982, 0.905


def main(argv: list[str] | None = None) -> int:
    here = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--summary", default=os.path.join(here, "replication", "replication_summary.json"))
    ap.add_argument("--transport", default=os.path.join(here, "replication", "transport.json"))
    ap.add_argument("--s35", default=os.path.join(here, "data", "Table_S35_calibration_CGGA693.csv"))
    ap.add_argument("--out", default=os.path.join(here, "Figure4_second_cohort.png"))
    args = ap.parse_args(argv)
    for path in (args.summary, args.transport, args.s35):
        if not os.path.exists(path):
            raise SystemExit(f"missing input: {path}\n"
                             "run replication_cgga693.py and transport_diagnostics.py first")
    rep = json.load(open(args.summary))
    tra = json.load(open(args.transport))

    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8, "savefig.bbox": "standard",
        "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.75))
    fig.subplots_adjust(left=0.128, right=0.985, bottom=0.165, top=0.700, wspace=0.60)

    # ------------------------------------------------------------------ A
    ax0 = ax[0]
    rows = [("325  MRG-4", *C325["mrg4"], "#B2182B"),
            ("325  Model A", *C325["modelA"], "#2166AC"),
            ("693  MRG-4", rep["c_mrg4"], *rep["ci_mrg4"], "#B2182B"),
            ("693  Model A", rep["c_modelA"], *rep["ci_modelA"], "#2166AC")]
    for i, (lab, v, lo, hi, col) in enumerate(rows):
        y = len(rows) - 1 - i
        ax0.plot([lo, hi], [y, y], color=col, lw=1.6)
        ax0.plot(v, y, "o", color=col, ms=4.5)
        ax0.text(0.888, y, f"{v:.3f}", va="center", ha="right", fontsize=7, color=col)
    ax0.set_yticks(range(len(rows)))
    ax0.set_yticklabels([r[0] for r in rows][::-1], fontsize=7)
    ax0.set_xlim(0.55, 0.90)
    ax0.set_xlabel("External C-index (95% CI)", fontsize=7)
    ax0.axvline(0.5, color="0.7", lw=0.6, ls=":")

    # ------------------------------------------------------------------ B
    ax1 = ax[1]
    s35 = [r for r in csv.DictReader(open(args.s35, encoding="utf-8-sig")) if r["Tertile"] != "all"]
    for months, pre, obs, lo, hi, mk, col in (
            ("36", "Predicted_S36", "Observed_S36", "KM_lower95_S36", "KM_upper95_S36", "o", "#B2182B"),
            ("60", "Predicted_S60", "Observed_S60", "KM_lower95_S60", "KM_upper95_S60", "s", "#2166AC")):
        xs = [float(r[pre]) for r in s35]
        ys = [float(r[obs]) for r in s35]
        el = [[y - float(r[lo]) for y, r in zip(ys, s35)],
              [float(r[hi]) - y for y, r in zip(ys, s35)]]
        ax1.errorbar(xs, ys, yerr=el, fmt=mk + "-", color=col, ms=4, lw=1.2,
                     capsize=2.5, elinewidth=0.8, label=f"{months} months")
    ax1.plot([0, 1], [0, 1], "k--", lw=0.6, label="ideal")
    ax1.set_xlim(0.25, 1.0)
    ax1.set_ylim(0.2, 1.0)
    ax1.set_xlabel("Predicted survival (tertile mean)", fontsize=7)
    ax1.set_ylabel("Observed survival (KM, 95% CI)", fontsize=7)
    ax1.legend(fontsize=6, frameon=False, loc="upper left", handlelength=1.4,
               borderaxespad=0.3, labelspacing=0.3)

    # ------------------------------------------------------------------ C
    ax2 = ax[2]
    data = {("CGGA-325", "II"): tra["CGGA-325"]["by_grade"]["WHO II"]["difference"],
            ("CGGA-325", "III"): tra["CGGA-325"]["by_grade"]["WHO III"]["difference"],
            ("CGGA-693", "II"): tra["CGGA-693"]["by_grade"]["WHO II"]["difference"],
            ("CGGA-693", "III"): tra["CGGA-693"]["by_grade"]["WHO III"]["difference"]}
    xs = {"CGGA-325": 0, "CGGA-693": 1}
    cols = {"II": "#8C510A", "III": "#01665E"}
    for (coh, gr), v in data.items():
        x = xs[coh] + (-0.16 if gr == "II" else 0.16)
        ax2.plot(x, v, "o", color=cols[gr], ms=5)
        ax2.text(x + 0.07, v, f"{v:+.3f}", va="center", fontsize=6.5, color=cols[gr])
    ax2.axhline(0, color="0.6", lw=0.7, ls="--")
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["CGGA-325", "CGGA-693"], fontsize=7)
    ax2.set_xlim(-0.55, 1.55)
    ax2.set_ylim(-0.055, 0.195)
    ax2.set_ylabel("C-index difference (MRG-4 \u2212 Model A)", fontsize=7)
    h1 = ax2.plot([], [], "o", color=cols["II"], ms=5, label="WHO grade II")[0]
    h2 = ax2.plot([], [], "o", color=cols["III"], ms=5, label="WHO grade III")[0]
    ax2.legend(handles=[h1, h2], fontsize=6, frameon=False, loc="upper right",
               handlelength=1.0, borderaxespad=0.3, labelspacing=0.3)

    # -------------------------------------------------------- panel headings
    renderer = fig.canvas.get_renderer()

    def text_w(txt, fs, weight="bold"):
        t = fig.text(0, -1, txt, fontsize=fs, fontweight=weight)
        bb = t.get_window_extent(renderer=renderer)
        t.remove()
        return bb.width

    col_left = []
    for a in ax:
        fig.canvas.draw()
        ab = a.get_window_extent(renderer=renderer)
        left = ab.x0
        if a.get_ylabel():
            left = min(left, a.yaxis.label.get_window_extent(renderer=renderer).x0)
        ticks = [t.get_window_extent(renderer=renderer).x0 for t in a.get_yticklabels() if t.get_text()]
        if ticks:
            left = min(left, min(ticks))
        col_left.append(left)

    fig_w_px = fig.get_window_extent(renderer=renderer).width
    head = [(a, letter, title, a.get_window_extent(renderer=renderer).x0)
            for a, letter, title in zip(ax, LETTERS, TITLES)]
    for i, (a, letter, title, x_left) in enumerate(head):
        right_limit = col_left[i + 1] - 6 if i + 1 < len(col_left) else fig_w_px - 4
        avail = right_limit - x_left
        fs = TITLE_FS_MAX
        while fs > TITLE_FS_MIN and text_w(title, fs) > avail:
            fs -= 0.2
        if text_w(title, fs) > avail + 0.5:
            raise SystemExit(f"panel {letter}: title wider than its column ({text_w(title, fs):.1f} > {avail:.1f} pt)")
        fig.text(x_left / fig_w_px, Y_TITLE_TOP, title, fontsize=fs, fontweight="bold",
                 ha="left", va="top", transform=fig.transFigure)
        fig.text(x_left / fig_w_px, Y_LETTER_TOP, letter, fontsize=LETTER_FS,
                 fontweight="bold", ha="left", va="top", transform=fig.transFigure)
    fig.canvas.draw()

    # ------------------------------------------------------------- assertions
    errs = []
    for i, (a, letter, title, x_left) in enumerate(head):
        right_limit = col_left[i + 1] - 6 if i + 1 < len(col_left) else fig_w_px - 4
        for t in fig.texts:
            if t.get_text() in (title, letter):
                bb = t.get_window_extent(renderer=renderer)
                if bb.x1 > right_limit + 0.5:
                    errs.append(f"{letter}: heading crosses into the next column "
                                f"({bb.x1:.1f} > {right_limit:.1f})")
        bb_title = [t for t in fig.texts if t.get_text() == title][0].get_window_extent(renderer=renderer)
        if a.get_window_extent(renderer=renderer).y1 > bb_title.y0 - 2:
            errs.append(f"{letter}: title overlaps its axes")

    letter_bb = {t.get_text(): t.get_window_extent(renderer=renderer)
                 for t in fig.texts if t.get_text() in LETTERS}
    title_bb = {t.get_text(): t.get_window_extent(renderer=renderer)
                for t in fig.texts if t.get_text() in TITLES}
    tops = []
    for letter, title in zip(LETTERS, TITLES):
        if letter_bb[letter].y0 < title_bb[title].y1 - 0.5:
            errs.append(f"{letter}: letter is not above the title")
        tops.append(round(title_bb[title].y1, 1))
    if max(tops) - min(tops) > 1.0:
        errs.append(f"panel titles are not top-aligned: {tops}")

    for i in range(1, len(ax)):
        prev_right = ax[i - 1].get_window_extent(renderer=renderer).x1
        for t in ax[i].get_yticklabels() + [ax[i].yaxis.label]:
            if t.get_text() and t.get_window_extent(renderer=renderer).x0 < prev_right + 2:
                errs.append(f"panel {LETTERS[i]}: y-axis text enters the previous panel")

    for t in (fig.texts
              + [x for a in ax for x in a.get_yticklabels() + a.get_xticklabels()]
              + [a.yaxis.label for a in ax] + [a.xaxis.label for a in ax]):
        if not t.get_text():
            continue
        bb = t.get_window_extent(renderer=renderer)
        if bb.x0 < 1.0 or bb.y0 < 1.0 or bb.x1 > fig_w_px - 1.0:
            errs.append(f"text clipped by the canvas: '{t.get_text()[:24]}'")

    if errs:
        raise SystemExit("LAYOUT CHECK FAILED:\n  " + "\n  ".join(errs))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=600)
    pdf = os.path.splitext(args.out)[0] + ".pdf"
    fig.savefig(pdf)
    print(f"wrote {args.out} and {pdf}")
    print(f"title tops: {tops} (aligned); LAYOUT CHECK: all green (6 assertions)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
