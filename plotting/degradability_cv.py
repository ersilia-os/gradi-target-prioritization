"""How well does the degradability model actually work? Cross-validated, on the measured labels.

    degradability_cv.png   1  ROC, both activators
                           2  precision-recall, both activators
                           3  against the baseline that matters and the yardstick that is not 1.0
                           4  the out-of-fold score distribution, split by the measured label

**AUROC IS NEVER REPORTED ALONE HERE, and that is a project rule with a measurement behind it.**
On the grid that chose this estimator, TabPFN's gain over the hand-set forest was ~4x larger on
PR-AUC than on AUROC, and AUROC alone called 4 of 6 arms "no difference". The gain is concentrated
at the top of the ranking, which is the part a shortlist consumes. So the PR panel is not decoration
beside the ROC panel -- it is the panel that discriminates.

**The numbers quoted are the 5-SEED MEAN +/- SD from the stage's CV table, NOT the area under the
curve drawn beside them.** The curves are one out-of-fold realisation and come out slightly
different -- ADEP4 0.8761 against a 5-seed 0.8738, ONC212 0.7706 against 0.7671, both inside one SD.
A single-seed estimate carries about +/-0.004 of pure arbitrariness, which is larger than several
effects this stage tested and rejected, so the seed-averaged value is the one that may be quoted.

**The baseline is LENGTH, not 0.5.** Short proteins are degraded more, so a model that learned only
size would already score 0.775 / 0.680. That is the number the model has to beat, and the ROC
diagonal is not a meaningful reference for this problem.

**`CROSS_ASSAY_AUROC` (0.877 / 0.815) is a YARDSTICK, NOT A CEILING.** It is how well one
activator's measured labels predict the other's -- i.e. how reproducible the assay itself is. The
model sits just under it here and exceeds it at three of five label cutoffs. Do not draw it as a
maximum.

**Cluster-grouped CV**: folds hold out whole sequence clusters, so a protein cannot be scored by its
own near-duplicate. Measured `leakage_gap` is +0.001 (ADEP4) and -0.005 (ONC212) -- unusually small,
which says this label set has little near-duplicate structure, not that grouping is unnecessary.

**The two activators are not two experiments.** Their probabilities correlate at rho 0.89 while the
measured labels agree at only 0.52. Two panels, one opinion.

Everything is read from the stage's own tables through `src/degradability.py`; the curves are
rendered from the stored out-of-fold predictions and no metric is recomputed for display.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/degradability_cv.py
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import stylia  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from sklearn.metrics import precision_recall_curve, roc_curve  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import degradability as D  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
ACT_COLOR = {"adep4": PAL.PRIMARY, "onc212": PAL.SECONDARY}
ACT_LABEL = {"adep4": "ADEP4", "onc212": "ONC212"}
OOF_COLUMN = "oof_esmc"


def _curves(activator: str) -> tuple[pd.DataFrame, pd.Series]:
    oof = D.load_oof(activator)
    cv = D.load_cv(activator).iloc[0]
    return oof, cv


def plot_roc(ax, data: dict) -> None:
    """ROC, with the LENGTH baseline marked rather than the diagonal.

    Deliberately NOT `set_aspect("equal")`: a square ROC is conventional on its own, but in a
    four-panel row it shrinks this axes box and leaves every title at a different height.

    The diagonal is drawn because readers expect it, but it is the wrong reference: a model that
    learned nothing but protein size already reaches 0.775 / 0.680 on these labels."""
    for act, (oof, cv) in data.items():
        fpr, tpr, _ = roc_curve(oof["hit"], oof[OOF_COLUMN])
        ax.plot(fpr, tpr, color=ACT_COLOR[act], lw=2.2,
                label=f"{ACT_LABEL[act]}   {cv.roc_auc_clustered:.3f} ± {cv.roc_auc_clustered_sd:.3f}")
    ax.plot([0, 1], [0, 1], color=PAL.MUTED, lw=1.0, ls="--")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=SS * 0.85, frameon=False, loc="lower right", handletextpad=0.5)
    stylia.label(ax, xlabel="False positive rate", ylabel="True positive rate",
                 title="ROC (5-seed mean ± SD)")


def plot_pr(ax, data: dict) -> None:
    """Precision-recall, with each activator's BASE RATE as its own floor.

    The floors differ -- 0.137 for ADEP4, 0.246 for ONC212 -- so the two PR curves are NOT on a
    common scale and the higher curve is not automatically the better model. That is exactly why
    the base rate is drawn."""
    for act, (oof, cv) in data.items():
        prec, rec, _ = precision_recall_curve(oof["hit"], oof[OOF_COLUMN])
        ax.plot(rec, prec, color=ACT_COLOR[act], lw=2.2,
                label=f"{ACT_LABEL[act]}   {cv.pr_auc_clustered:.3f} ± {cv.pr_auc_clustered_sd:.3f}")
        ax.axhline(cv.base_rate, color=ACT_COLOR[act], lw=1.0, ls=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(fontsize=SS * 0.85, frameon=False, loc="upper right", handletextpad=0.5)
    stylia.label(ax, xlabel="Recall", ylabel="Precision",
                 title="Precision-recall (dotted = base rate)")


def plot_benchmark(ax, data: dict) -> None:
    """The model against the baseline it must beat and the yardstick it is measured against.

    Three marks per activator: length-only (what size alone buys), the model (+/- SD over 5 seeds),
    and the cross-assay AUROC -- how well one activator's labels predict the other's, i.e. the
    assay's own reproducibility. The last is drawn hollow because it is a REFERENCE, not a ceiling:
    the model exceeds it at three of five label cutoffs."""
    acts = list(data)
    y = np.arange(len(acts))[::-1]
    for yi, act in zip(y, acts):
        cv = data[act][1]
        ax.plot([D.V1_LENGTH_ONLY_AUROC[act], cv.roc_auc_clustered], [yi, yi],
                color="#CFCFCB", lw=1.6, zorder=1)
        ax.scatter([D.V1_LENGTH_ONLY_AUROC[act]], [yi], s=60, color=PAL.MUTED,
                   linewidths=0, zorder=3)
        ax.errorbar([cv.roc_auc_clustered], [yi], xerr=[cv.roc_auc_clustered_sd], fmt="o",
                    color=ACT_COLOR[act], markersize=9, capsize=3, lw=1.4, zorder=4,
                    markeredgecolor="white", markeredgewidth=0.6)
        ax.scatter([D.CROSS_ASSAY_AUROC[act]], [yi], s=70, facecolors="white",
                   edgecolors=PAL.INK, linewidths=1.3, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([ACT_LABEL[a] for a in acts], fontsize=SS)
    ax.set_ylim(-0.7, len(acts) - 0.3)
    ax.set_xlim(0.6, 0.95)
    ax.grid(axis="x", color="#EFEFEC", lw=0.8)
    ax.set_axisbelow(True)
    handles = [
        Line2D([], [], marker="o", ls="", markersize=6, color=PAL.MUTED, label="length only"),
        Line2D([], [], marker="o", ls="", markersize=7, color=PAL.INK, label="TabPFN (± SD)"),
        Line2D([], [], marker="o", ls="", markersize=7, markerfacecolor="white",
               markeredgecolor=PAL.INK, color="none", label="cross-assay yardstick"),
    ]
    ax.legend(handles=handles, fontsize=SS * 0.8, frameon=False, loc="lower right",
              handletextpad=0.3, labelspacing=0.3)
    stylia.label(ax, xlabel="AUROC", ylabel="", title="Beating size, not beating 0.5")


def plot_distribution(ax, data: dict) -> None:
    """Out-of-fold score, split by the MEASURED label. This is what the two areas summarise.

    Violins rather than histograms because the comparison is between four distributions on one
    scale, and four overlaid histograms is unreadable. Filled = measured substrate, hollow = not,
    the same convention the enrichment figure uses."""
    pos = 0
    ticks, tick_labels = [], []
    for act, (oof, cv) in data.items():
        for hit in (0, 1):
            v = oof.loc[oof["hit"] == hit, OOF_COLUMN].dropna().to_numpy()
            parts = ax.violinplot([v], positions=[pos], widths=0.75, showextrema=False)
            for body in parts["bodies"]:
                body.set_alpha(1.0)
                if hit:
                    body.set_facecolor(ACT_COLOR[act])
                    body.set_edgecolor("white")
                else:
                    body.set_facecolor("white")
                    body.set_edgecolor(ACT_COLOR[act])
                body.set_linewidth(1.1)
            ax.scatter([pos], [np.median(v)], s=14, color=PAL.INK, zorder=4)
            ticks.append(pos)
            tick_labels.append(f"{'substrate' if hit else 'not'}\nn={len(v):,}")
            pos += 1
        # the cut this axis actually uses, over that activator's pair only
        ax.plot([pos - 2.45, pos - 0.55], [D.BASE_RATE_THRESHOLD[act]] * 2,
                color=ACT_COLOR[act], lw=1.2, ls="--")
        pos += 0.6

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels, fontsize=SS * 0.78)
    ax.set_ylim(0, 1)
    # Activator names go BELOW the tick labels, in axes-fraction y, which is the one place they
    # cannot collide with either the violins or the title.
    trans = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for act, xc in zip(data, (0.5, 3.1)):
        ax.text(xc, -0.115, ACT_LABEL[act], ha="center", va="top", fontsize=SS,
                color=ACT_COLOR[act], transform=trans)
    stylia.label(ax, xlabel="", ylabel="Out-of-fold probability",
                 title="Dashed = the cut this axis uses")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  degradability cross-validation -- S. aureus measured labels")
    data = {act: _curves(act) for act in D.ACTIVATORS}

    fig, axs = stylia.create_figure(1, 4, width_ratios=[1, 1, 1.05, 1.25], width=1.0, height=0.42)
    plot_roc(axs.next(), data)
    plot_pr(axs.next(), data)
    plot_benchmark(axs.next(), data)
    plot_distribution(axs.next(), data)
    out = OUT_DIR / "degradability_cv.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say(f"\n  {D.DEFAULT_ESTIMATOR} {D.TABPFN_MODEL} on ESM-C, "
        f"{D.N_SPLITS}-fold cluster-grouped CV x {D.N_CV_SEEDS} seeds")
    say(f"    {'':<10} {'n':>6} {'pos':>5} {'base':>7} {'AUROC':>17} {'PR-AUC':>17} {'length':>8} {'yardstick':>10}")
    for act, (oof, cv) in data.items():
        say(f"    {ACT_LABEL[act]:<10} {int(cv.n):>6,} {int(cv.n_pos):>5} {cv.base_rate:>7.3f}"
            f" {cv.roc_auc_clustered:>9.4f} ± {cv.roc_auc_clustered_sd:<5.4f}"
            f" {cv.pr_auc_clustered:>9.4f} ± {cv.pr_auc_clustered_sd:<5.4f}"
            f" {D.V1_LENGTH_ONLY_AUROC[act]:>8.3f} {D.CROSS_ASSAY_AUROC[act]:>10.3f}")

    say("\n  THE CURVE IS NOT THE NUMBER")
    for act, (oof, cv) in data.items():
        from sklearn.metrics import average_precision_score, roc_auc_score

        say(f"    {ACT_LABEL[act]:<10} curve drawn: AUROC {roc_auc_score(oof['hit'], oof[OOF_COLUMN]):.4f}"
            f"  AP {average_precision_score(oof['hit'], oof[OOF_COLUMN]):.4f}"
            f"   vs 5-seed {cv.roc_auc_clustered:.4f} / {cv.pr_auc_clustered:.4f}")
    say("    One out-of-fold realisation against the seed-averaged value. Both inside one SD.")
    say("    Quote the 5-seed number; a single seed carries ±0.004 of arbitrariness.")

    say("\n  LEAKAGE CHECK (clustered minus plain AUROC)")
    for act, (oof, cv) in data.items():
        say(f"    {ACT_LABEL[act]:<10} clustered {cv.roc_auc_clustered:.4f}   "
            f"plain {cv.roc_auc_plain:.4f}   gap {cv.leakage_gap:+.4f}")
    say("    Unusually small: this label set has little near-duplicate structure. It does NOT")
    say("    mean grouping is unnecessary -- it means grouping cost nothing here.")

    say("\n  CAVEATS")
    say("    - PR-AUC is reported beside AUROC, never AUROC alone: the estimator's gain was ~4x")
    say("      larger on PR-AUC, and AUROC alone called 4 of 6 arms 'no difference'.")
    say("    - The base rates differ (0.137 vs 0.246), so the two PR curves are NOT on a common")
    say("      scale. A higher PR curve is not automatically a better model.")
    say("    - The cross-assay AUROC is a yardstick, NOT a ceiling: the model exceeds it at")
    say("      three of five label cutoffs.")
    say("    - These are S. aureus measurements. Nothing here validates the Kp or Ec columns,")
    say("      which are extrapolations priced by nn_similarity in degradability.png.")


if __name__ == "__main__":
    main()
