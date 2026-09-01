"""2x3 summary of the localization axis (09g).

House style: stylia "slide" format, NPG palette, every panel specific to the single `--organism`.

  1  compartment composition   final localization classes
  2  coverage gained           UniProt-only (the previous state of the axis) vs the merged call
  3  evidence tier             experimental / curated / predicted, and the winning source
  4  predictor concordance     DeepLocPro vs PSORTb 3.0 — two independent predictors
  5  topology                  TM helices per compartment, and where beta-barrels land
  6  Clp accessibility         graded score distribution by compartment

Reads output/results/<org>/<prefix>_localization.csv (09g).
Output: output/plots/09h_localization_<prefix>.png (one slide per --organism). Run with `gradi`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL

# Compartment order + palette live in src/localization.py so every 09* slide shares one definition.
CLASS_ORDER = [c for c in LOC.LOC_CLASS_ORDER if c != "membrane"]
CLASS_COLOR = LOC.LOC_CLASS_COLOR
EVIDENCE_COLOR = {"experimental": NPG[2], "curated": NPG[4], "predicted": "#C9C9C7", "none": "#EEEEEE"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]
    d = pd.read_csv(LOC.results_dir(org) / f"{prefix}_localization.csv")

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: compartment composition ----
    ax = axs.next()
    vc = d["localization"].value_counts().reindex(CLASS_ORDER).fillna(0).astype(int)
    present = [c for c in CLASS_ORDER if vc[c] > 0]
    bars = ax.bar(range(len(present)), [vc[c] for c in present],
                  color=[CLASS_COLOR[c] for c in present])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in present], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Compartment — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 2: coverage gained over the UniProt-only axis ----
    ax = axs.next()
    n = len(d)
    uni_known = int((d["uniprot_localization"].fillna("unknown") != "unknown").sum())
    final_known = int((d["localization"] != "unknown").sum())
    bars = ax.bar([0, 1], [uni_known, final_known], color=["#C9C9C7", NPG[0]], width=0.55)
    ax.bar_label(bars, labels=[f"{uni_known}\n{uni_known / n:.0%}", f"{final_known}\n{final_known / n:.0%}"],
                 padding=2, fontsize=SS)
    ax.axhline(n, ls="--", color="#555555", lw=1)
    ax.text(1.45, n, "proteome", fontsize=SS, color="#555555", va="bottom", ha="right")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["UniProt only", "merged (09g)"], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="proteins localized", title=f"Coverage gained — {orgname}")
    ax.margins(y=0.22)

    # ---- panel 3: evidence tier, split by winning source ----
    ax = axs.next()
    tiers = ["experimental", "curated", "predicted"]
    sub = d[d["localization"] != "unknown"]
    counts = [int((sub["localization_evidence"] == t).sum()) for t in tiers]
    bars = ax.bar(range(3), counts, color=[EVIDENCE_COLOR[t] for t in tiers])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.set_xticks(range(3)); ax.set_xticklabels(tiers, fontsize=SS)
    srcs = sub["localization_source"].value_counts()
    ax.text(0.02, 0.97, "\n".join(f"{s}: {c}" for s, c in srcs.head(5).items()),
            transform=ax.transAxes, fontsize=SS, va="top", ha="left", color="#555555")
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Evidence tier — {orgname}")
    ax.margins(y=0.18)

    # ---- panel 4: DeepLocPro vs PSORTb concordance ----
    ax = axs.next()
    both = d[d["dlp_localization"].notna() & d["psortb_localization"].notna()
             & (d["psortb_localization"] != "unknown")]
    classes = [c for c in CLASS_ORDER if c != "unknown"]
    if len(both):
        mat = np.zeros((len(classes), len(classes)))
        for i, a in enumerate(classes):
            for j, b in enumerate(classes):
                mat[i, j] = ((both["dlp_localization"] == a) & (both["psortb_localization"] == b)).sum()
        norm = mat / np.maximum(mat.sum(), 1)
        ax.imshow(norm, cmap="Greens", aspect="auto")
        for i in range(len(classes)):
            for j in range(len(classes)):
                if mat[i, j]:
                    ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=SS - 1,
                            color="#222222" if norm[i, j] < norm.max() * 0.6 else "white")
        agree = (both["dlp_localization"] == both["psortb_localization"]).mean()
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in classes], fontsize=SS, rotation=45)
        ax.set_yticks(range(len(classes)))
        ax.set_yticklabels([LOC.LOC_ABBREV[c] for c in classes], fontsize=SS)
        title = f"DeepLocPro vs PSORTb — {agree:.0%} agree"
    else:
        ax.text(0.5, 0.5, "no overlapping calls", ha="center", va="center", transform=ax.transAxes)
        title = "DeepLocPro vs PSORTb"
    stylia.label(ax, xlabel="PSORTb 3.0", ylabel="DeepLocPro", title=title)

    # ---- panel 5: topology per compartment ----
    ax = axs.next()
    if "n_tm_helix" in d.columns and d["n_tm_helix"].notna().any():
        present5 = [c for c in classes if (d["localization"] == c).sum() > 0]
        helix = [d.loc[d["localization"] == c, "n_tm_helix"].fillna(0).mean() for c in present5]
        barrels = [int(d.loc[d["localization"] == c, "is_beta_barrel"].fillna(False).sum())
                   for c in present5]
        x = np.arange(len(present5))
        bars = ax.bar(x, helix, color=[CLASS_COLOR[c] for c in present5])
        ax.bar_label(bars, labels=[f"{h:.1f}" for h in helix], padding=2, fontsize=SS)
        ax2 = ax.twinx()
        ax2.plot(x, barrels, "o--", color="#E64B35", lw=1.2, ms=4, label="β-barrels")
        ax2.set_ylabel("β-barrels", fontsize=SS, color="#E64B35")
        ax2.tick_params(axis="y", labelsize=SS, colors="#E64B35")
        ax.set_xticks(x); ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in present5], fontsize=SS)
        stylia.label(ax, xlabel="", ylabel="mean TM helices", title=f"Topology — {orgname}")
        ax.margins(y=0.2)
    else:
        ax.text(0.5, 0.5, "TMbed topology not available\n(run 09d)", ha="center", va="center",
                transform=ax.transAxes, fontsize=SS, color="#555555")
        stylia.label(ax, xlabel="", ylabel="", title=f"Topology — {orgname}")

    # ---- panel 6: Clp accessibility ----
    ax = axs.next()
    acc = d["clp_accessibility"].dropna()
    levels = sorted(acc.unique())
    counts = [int((acc == lv).sum()) for lv in levels]
    bars = ax.bar(range(len(levels)), counts,
                  color=[stylia.FadingColormap("turquoise").cmap(lv) for lv in levels])
    ax.bar_label(bars, padding=2, fontsize=SS)
    ax.axvline(len(levels) - 0.5 - sum(1 for lv in levels if lv >= LOC.MIN_CLP_ACCESSIBILITY),
               ls="--", color="#555555", lw=1)
    ax.set_xticks(range(len(levels))); ax.set_xticklabels([f"{lv:g}" for lv in levels], fontsize=SS)
    n_short = int((acc >= LOC.MIN_CLP_ACCESSIBILITY).sum())
    stylia.label(ax, xlabel="clp_accessibility", ylabel="proteins",
                 title=f"Clp accessibility — {n_short} above {LOC.MIN_CLP_ACCESSIBILITY:g}")
    ax.margins(y=0.18)

    out = LOC.REPO_ROOT / "output" / "plots" / f"09h_localization_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")


if __name__ == "__main__":
    main()
