"""Evidence and predictor behaviour behind the localization call (docs §5.1).

Where each call actually came from, how many independent opinions back it, and how confident the
merge is entitled to be. This is the slide the axis defends itself with: 09h says *what* the
compartments are, this says *why you should believe them*.

  1  source x compartment      which source supplies which compartment (row-normalised)
  2  independent opinions      # sources per protein, stacked by evidence tier
  3  predictor agreement       fraction of UniProt/DeepLocPro/PSORTb agreeing with the final call
  4  confidence by tier        localization_confidence distribution per evidence tier
  5  DeepLocPro margin         top1-top2 probability gap per predicted class — where the model hedges
  6  transferred evidence      the 09f track: source x evidence, and ortholog donor agreement

Reads output/results/<org>/<prefix>_localization.csv.
Output: output/plots/09j_evidence_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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

EVIDENCE_ORDER = ["experimental", "curated", "predicted"]
EVIDENCE_COLOR = {"experimental": "#00A087", "curated": "#FCBF49", "predicted": "#C9C9C7"}

# Source labels short enough for an axis tick.
SOURCE_LABEL = {
    "uniprot_experimental": "UniProt exp.",
    "uniprot_curated": "UniProt cur.",
    "stepdb": "STEPdb",
    "stepdb_curated": "STEPdb cur.",
    "ortholog_transfer": "Ortholog",
    "deeplocpro": "DeepLocPro",
    "psortb": "PSORTb",
    "none": "none",
}
CLASSES = [c for c in LOC.LOC_CLASS_ORDER if c not in ("membrane", "unknown")]


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

    # ---- panel 1: source x compartment, row-normalised ----
    ax = axs.next()
    sources = [s for s in LOC.SOURCE_PRECEDENCE if (d["localization_source"] == s).any()]
    mat = np.zeros((len(sources), len(CLASSES)))
    for i, s in enumerate(sources):
        sub = d[d["localization_source"] == s]
        for j, c in enumerate(CLASSES):
            mat[i, j] = (sub["localization"] == c).sum()
    row_tot = mat.sum(axis=1, keepdims=True)
    norm = np.divide(mat, np.maximum(row_tot, 1))
    ax.imshow(norm, cmap="Greens", aspect="auto", vmin=0, vmax=1)
    for i in range(len(sources)):
        for j in range(len(CLASSES)):
            if mat[i, j]:
                ax.text(j, i, int(mat[i, j]), ha="center", va="center", fontsize=SS - 1.5,
                        color="white" if norm[i, j] > 0.55 else "#222222")
    ax.set_xticks(range(len(CLASSES)))
    ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in CLASSES], fontsize=SS)
    ax.set_yticks(range(len(sources)))
    ax.set_yticklabels([f"{SOURCE_LABEL.get(s, s)} ({int(row_tot[i][0]):,})"
                        for i, s in enumerate(sources)], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="", title=f"Which source calls what — {orgname}")

    # ---- panel 2: number of independent sources per protein ----
    ax = axs.next()
    levels = sorted(d["n_localization_sources"].dropna().unique())
    bottom = np.zeros(len(levels))
    for tier in EVIDENCE_ORDER:
        vals = np.array([int(((d["n_localization_sources"] == lv)
                              & (d["localization_evidence"] == tier)).sum()) for lv in levels],
                        dtype=float)
        ax.bar(range(len(levels)), vals, bottom=bottom, color=EVIDENCE_COLOR[tier],
               label=tier, width=0.7)
        bottom += vals
    for i, tot in enumerate(bottom):
        ax.text(i, tot, f"{int(tot):,}", ha="center", va="bottom", fontsize=SS, color="#555555")
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"{int(lv)}" for lv in levels], fontsize=SS)
    ax.legend(fontsize=SS, frameon=False, loc="upper left")
    ax.margins(y=0.18)
    stylia.label(ax, xlabel="independent sources with an opinion", ylabel="proteins",
                 title=f"Corroboration depth — {orgname}")

    # ---- panel 3: predictor agreement ----
    ax = axs.next()
    ag = d["predictor_agreement"].dropna()
    levels = sorted(ag.unique())
    bottom = np.zeros(len(levels))
    for tier in EVIDENCE_ORDER:
        vals = np.array([int(((d["predictor_agreement"] == lv)
                              & (d["localization_evidence"] == tier)).sum()) for lv in levels],
                        dtype=float)
        ax.bar(range(len(levels)), vals, bottom=bottom, color=EVIDENCE_COLOR[tier], width=0.7)
        bottom += vals
    for i, tot in enumerate(bottom):
        ax.text(i, tot, f"{int(tot):,}", ha="center", va="bottom", fontsize=SS, color="#555555")
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"{lv:.2f}".rstrip("0").rstrip(".") for lv in levels], fontsize=SS)
    ax.margins(y=0.18)
    stylia.label(ax, xlabel="fraction of opinions agreeing with the final call", ylabel="proteins",
                 title=f"Predictor agreement — mean {ag.mean():.2f}")

    # ---- panel 4: confidence by evidence tier ----
    ax = axs.next()
    data, labels, colors = [], [], []
    for tier in EVIDENCE_ORDER:
        v = d.loc[d["localization_evidence"] == tier, "localization_confidence"].dropna()
        if len(v):
            data.append(v.to_numpy()); labels.append(f"{tier}\n({len(v):,})")
            colors.append(EVIDENCE_COLOR[tier])
    if data:
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False,
                        medianprops=dict(color="#2B2333", lw=1.4))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_edgecolor("#8A8A96"); patch.set_alpha(0.9)
        ax.set_xticklabels(labels, fontsize=SS)
    ax.set_ylim(0, 1.05)
    stylia.label(ax, xlabel="", ylabel="localization_confidence",
                 title=f"Confidence by evidence tier — {orgname}")

    # ---- panel 5: DeepLocPro margin per predicted class ----
    ax = axs.next()
    present = [c for c in CLASSES if (d["dlp_localization"] == c).sum() >= 5]
    data = [d.loc[d["dlp_localization"] == c, "dlp_margin"].dropna().to_numpy() for c in present]
    if data:
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False,
                        medianprops=dict(color="#2B2333", lw=1.4))
        for patch, c in zip(bp["boxes"], present):
            patch.set_facecolor(LOC.LOC_CLASS_COLOR[c]); patch.set_edgecolor("#8A8A96")
            patch.set_alpha(0.9)
        ax.set_xticklabels([f"{LOC.LOC_ABBREV[c]}\n({int((d['dlp_localization'] == c).sum()):,})"
                            for c in present], fontsize=SS)
    ax.set_ylim(0, 1.05)
    stylia.label(ax, xlabel="", ylabel="top1 − top2 probability",
                 title=f"DeepLocPro decisiveness — {orgname}")

    # ---- panel 6: what DeepLocPro conflates with what ----
    # Mean probability vector per predicted class. The diagonal is decisiveness; off-diagonal mass
    # is the compartment the model would have picked instead — a within-model view, distinct from
    # 09h's between-predictor confusion matrix.
    ax = axs.next()
    prob_cols = [f"dlp_p_{c}" for c in CLASSES if f"dlp_p_{c}" in d.columns]
    pred = [c for c in CLASSES if (d["dlp_localization"] == c).sum() >= 5]
    if prob_cols and pred:
        mat = np.array([[d.loc[d["dlp_localization"] == r, pc].mean() for pc in prob_cols]
                        for r in pred])
        ax.imshow(mat, cmap="Greens", aspect="auto", vmin=0, vmax=1)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                if mat[i, j] >= 0.01:
                    ax.text(j, i, f"{mat[i, j]:.2f}".lstrip("0"), ha="center", va="center",
                            fontsize=SS - 1.5, color="white" if mat[i, j] > 0.55 else "#222222")
        ax.set_xticks(range(len(prob_cols)))
        ax.set_xticklabels([LOC.LOC_ABBREV[pc.replace("dlp_p_", "")] for pc in prob_cols],
                           fontsize=SS)
        ax.set_yticks(range(len(pred)))
        ax.set_yticklabels([LOC.LOC_ABBREV[c] for c in pred], fontsize=SS)
    else:
        ax.text(0.5, 0.5, "no DeepLocPro probabilities", ha="center", va="center",
                transform=ax.transAxes, fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="mean probability assigned", ylabel="predicted class",
                 title=f"What DeepLocPro conflates — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"09j_evidence_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")


if __name__ == "__main__":
    main()
