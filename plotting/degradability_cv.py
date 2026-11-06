"""How well does the degradability model actually work? Cross-validated, on the measured labels.

    degradability_cv.png   1  ROC, both activators, with the spread across CV folds
                           2  precision-recall, same
                           3  where the scores land in each proteome -- and seen vs unseen

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

**THE BAND IS THE SPREAD ACROSS THE 5 CV FOLDS, not a confidence interval.** Each fold's held-out
curve is interpolated onto a common grid and the shaded envelope is their min-to-max. It answers
"how much does this curve move if you resample the data", which is the honest thing to show beside
a single line. It is NOT a standard error and must not be read as one -- with 46 positives per fold
(ADEP4) and ~51 (ONC212), a single fold's PR curve is genuinely noisy, and the band is wide at high
recall for that reason alone.

**The baseline is LENGTH, not 0.5.** Short proteins are degraded more, so a model that learned only
size would already score 0.775 / 0.680. That and the cross-assay level (0.877 / 0.815) are drawn on
the ROC panel as **iso-AUROC reference curves** -- binormal curves whose area is exactly that value.
An AUROC is an area, so a scalar cannot be drawn on an ROC axis any other way. They are LEVELS, not
measurements: the real length-only classifier traces some other path with the same area, and these
curves must never be read as "what the baseline did".

**The third panel is the deployment picture, not a validation.** S. aureus is the only organism with
measurements, so it is split into the proteins the model was fitted on (`seen` -- scored
out-of-fold) and the rest of its proteome (`unseen`). K. pneumoniae and E. coli are unseen in
their entirety. **`seen` is a different set for each activator** -- 1,677 proteins for ADEP4 against
1,045 for ONC212 -- because the two screens measured different numbers of proteins.

**`CROSS_ASSAY_AUROC` is a YARDSTICK, NOT A CEILING.** It is how well one activator's measured
labels predict the other's -- how reproducible the assay itself is. The model sits just under it
here and exceeds it at three of five label cutoffs. Never draw it as a maximum.

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
from scipy.stats import norm  # noqa: E402
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


#: How much wider the deployment panel is than each square curve panel. Used BOTH in
#: `width_ratios` and in that panel's `set_box_aspect`, so the three axes boxes end up exactly the
#: same height -- which is what keeps the three titles on one line.
VIOLIN_RATIO = 1.3

GRID = np.linspace(0.0, 1.0, 201)


def _iso_auroc(auc: float) -> np.ndarray:
    """A smooth ROC curve whose area is exactly `auc`, for drawing a reference LEVEL as a line.

    An AUROC is an area, so a scalar baseline cannot be drawn on an ROC axis directly. This is the
    standard binormal form -- equal-variance, `TPR = Phi(Phi^-1(FPR) + d)` with `d = sqrt(2)
    Phi^-1(auc)` -- which integrates to `auc` by construction.

    **It is a LEVEL, not a measurement.** The real length-only classifier would trace some other
    path with the same area, so these curves must never be read as "what the baseline did", only as
    "where a model of this quality sits". Labelled that way in the legend."""
    f = np.clip(GRID, 1e-6, 1 - 1e-6)
    return norm.cdf(norm.ppf(f) + np.sqrt(2.0) * norm.ppf(auc))


def _fold_band(oof: pd.DataFrame, which: str) -> tuple[np.ndarray, np.ndarray]:
    """Min-to-max envelope of the per-fold curves on a common grid.

    Each CV fold is a held-out set, so its curve is a complete, independent estimate; interpolating
    them onto one grid is the only way to stack them. A fold whose curve does not reach the end of
    the grid is extended with its last value rather than dropped -- dropping would silently narrow
    the band exactly where it should be widest."""
    curves = []
    for _, g in oof.groupby("fold"):
        y, sc = g["hit"].to_numpy(), g[OOF_COLUMN].to_numpy()
        if y.sum() == 0 or y.sum() == len(y):
            continue
        if which == "roc":
            x, v, _ = roc_curve(y, sc)
        else:
            v, x, _ = precision_recall_curve(y, sc)
            x, v = x[::-1], v[::-1]          # recall ascending
        curves.append(np.interp(GRID, x, v, left=v[0], right=v[-1]))
    a = np.vstack(curves)
    return a.min(axis=0), a.max(axis=0)


def plot_roc(ax, data: dict) -> None:
    """ROC, pooled line with the across-fold envelope behind it.

    **Square, via `set_box_aspect(1)`.** A ROC read on a stretched axis is misleading -- the
    distance of the curve from the diagonal is the whole visual argument, and it is only
    proportional to the effect when x and y are on the same scale. `set_box_aspect` squares the
    AXES BOX directly, where `set_aspect("equal")` squares it by shrinking it inside its slot and
    so drops the title relative to its neighbours. The deployment panel's box aspect is set from
    the same constant, so all three boxes come out the same height."""
    for act, (oof, cv) in data.items():
        lo, hi = _fold_band(oof, "roc")
        ax.fill_between(GRID, lo, hi, color=ACT_COLOR[act], alpha=0.18, linewidth=0, zorder=1)
        fpr, tpr, _ = roc_curve(oof["hit"], oof[OOF_COLUMN])
        ax.plot(fpr, tpr, color=ACT_COLOR[act], lw=2.2, zorder=3,
                label=f"{ACT_LABEL[act]}   {cv.roc_auc_clustered:.3f} ± {cv.roc_auc_clustered_sd:.3f}")
    ax.plot([0, 1], [0, 1], color=PAL.MUTED, lw=1.0, ls="--", zorder=0)

    # The two reference LEVELS, as iso-AUROC curves in each activator's colour: what the model has
    # to beat (length alone) and what the assay reproduces (cross-assay). Thin and translucent so
    # they read as grid, not as data.
    for act in data:
        ax.plot(GRID, _iso_auroc(D.V1_LENGTH_ONLY_AUROC[act]), color=ACT_COLOR[act],
                lw=1.1, ls=":", alpha=0.75, zorder=2)
        ax.plot(GRID, _iso_auroc(D.CROSS_ASSAY_AUROC[act]), color=ACT_COLOR[act],
                lw=1.1, ls=(0, (5, 2)), alpha=0.75, zorder=2)

    # Named ON the curve, not in the legend: a reference level is a position on this axis, and a
    # legend entry makes the reader hunt for which line it means. Annotated on the ADEP4 pair only
    # -- the ONC212 pair repeats the line styles, so two labels carry four curves. The white box
    # keeps the text legible where a model curve passes behind it.
    x_at = 0.10
    for auc_map, name in ((D.V1_LENGTH_ONLY_AUROC, "length only"),
                          (D.CROSS_ASSAY_AUROC, "cross-assay")):
        y_at = float(np.interp(x_at, GRID, _iso_auroc(auc_map["adep4"])))
        ax.text(x_at + 0.03, y_at, name, fontsize=SS * 0.8, color=PAL.MUTED,
                va="center", ha="left",
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.2, "alpha": 0.85})

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)
    ax.legend(fontsize=SS * 0.82, frameon=False, loc="lower right", handletextpad=0.5)
    stylia.label(ax, xlabel="False positive rate", ylabel="True positive rate",
                 title="ROC (band = 5 CV folds)")


def plot_pr(ax, data: dict) -> None:
    """Precision-recall, with each activator's BASE RATE as its own floor.

    The floors differ -- 0.137 for ADEP4, 0.246 for ONC212 -- so the two curves are NOT on a common
    scale and the higher one is not automatically the better model. That is why the floor is drawn."""
    for act, (oof, cv) in data.items():
        lo, hi = _fold_band(oof, "pr")
        ax.fill_between(GRID, lo, hi, color=ACT_COLOR[act], alpha=0.18, linewidth=0, zorder=1)
        prec, rec, _ = precision_recall_curve(oof["hit"], oof[OOF_COLUMN])
        ax.plot(rec, prec, color=ACT_COLOR[act], lw=2.2, zorder=3,
                label=f"{ACT_LABEL[act]}   {cv.pr_auc_clustered:.3f} ± {cv.pr_auc_clustered_sd:.3f}")
        ax.axhline(cv.base_rate, color=ACT_COLOR[act], lw=1.0, ls=":", zorder=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1)
    ax.legend(fontsize=SS * 0.82, frameon=False, loc="upper right", handletextpad=0.5)
    stylia.label(ax, xlabel="Recall", ylabel="Precision",
                 title="Precision-recall (dotted = base rate)")


#: The four groups of the deployment panel, as (tick label, species, seen-or-None).
GROUPS = [
    ("K. pneumoniae\nunseen", "kpneumoniae", None),
    ("E. coli\nunseen", "ecoli", None),
    ("S. aureus\nseen", "saureus", True),
    ("S. aureus\nunseen", "saureus", False),
]


#: Centre-to-centre spacing of the four groups. Wider than 1.0 because each group is a SPLIT
#: violin: at spacing 1.0 the right half of one group overlaps the left half of the next, which
#: reads as a single malformed shape rather than as two activators.
GROUP_STEP = 1.7
VIOLIN_WIDTH = 1.35


def _half_violin(ax, values: np.ndarray, pos: float, side: str, color, filled: bool) -> None:
    """One side of a split violin.

    matplotlib has no split violin, so the body is drawn whole and then its path is clipped to one
    side of `pos`. Clipping the vertices rather than halving the kernel keeps both halves on the
    same density estimate, which is the only way the two are comparable."""
    parts = ax.violinplot([values], positions=[pos], widths=VIOLIN_WIDTH, showextrema=False)
    for b in parts["bodies"]:
        v = b.get_paths()[0].vertices
        v[:, 0] = np.clip(v[:, 0], -np.inf, pos) if side == "left" else np.clip(v[:, 0], pos, np.inf)
        b.set_alpha(1.0)
        b.set_facecolor(color if filled else "white")
        b.set_edgecolor(color)
        b.set_linewidth(1.0)


def plot_species(ax, species_scores: dict) -> None:
    """Where the predicted scores actually land, per proteome, split by activator.

    Left half ADEP4, right half ONC212. The S. aureus `seen` group is the labelled set scored
    out-of-fold; everything else is a model prediction on an organism the model never saw. The two
    `seen` sets are NOT the same proteins -- 1,677 for ADEP4, 1,045 for ONC212 -- because the two
    screens measured different numbers of proteins."""
    for i, (label, _sp, _seen) in enumerate(GROUPS):
        x = i * GROUP_STEP
        for act, side in (("adep4", "left"), ("onc212", "right")):
            v = species_scores[(label, act)]
            _half_violin(ax, v, x, side, ACT_COLOR[act], filled=(side == "left"))
        for act, dx in (("adep4", -0.22), ("onc212", 0.22)):
            ax.scatter([x + dx], [np.median(species_scores[(label, act)])], s=12,
                       color=PAL.INK, zorder=5)

    for act, ls in (("adep4", "--"), ("onc212", ":")):
        ax.axhline(D.BASE_RATE_THRESHOLD[act], color=ACT_COLOR[act], lw=1.1, ls=ls, zorder=0)

    ax.set_xticks([i * GROUP_STEP for i in range(len(GROUPS))])
    ax.set_xticklabels([g[0] for g in GROUPS], fontsize=SS * 0.8)
    ax.set_xlim(-VIOLIN_WIDTH / 2 - 0.2, (len(GROUPS) - 1) * GROUP_STEP + VIOLIN_WIDTH / 2 + 0.2)
    ax.set_ylim(0, 1)
    ax.set_box_aspect(1 / VIOLIN_RATIO)
    handles = [
        Line2D([], [], marker="s", ls="", markersize=7, color=ACT_COLOR["adep4"], label="ADEP4"),
        Line2D([], [], marker="s", ls="", markersize=7, markerfacecolor="white",
               markeredgecolor=ACT_COLOR["onc212"], color="none", label="ONC212"),
    ]
    ax.legend(handles=handles, fontsize=SS * 0.82, frameon=False, loc="upper right",
              handletextpad=0.4, labelspacing=0.3)
    stylia.label(ax, xlabel="", ylabel="Predicted probability",
                 title="Where the scores land (lines = the cuts)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  degradability cross-validation -- S. aureus measured labels")
    data = {act: _curves(act) for act in D.ACTIVATORS}

    # The deployment panel reads the DELIVERABLE, not the OOF table: `<act>_prob` is already the
    # out-of-fold value wherever a protein was measured, so the two are on one scale by design.
    species_scores: dict[tuple[str, str], np.ndarray] = {}
    counts: dict[tuple[str, str], int] = {}
    for label, sp, seen in GROUPS:
        df = D.load(sp)
        for act in D.ACTIVATORS:
            sub = df
            if seen is not None:
                has = df["uniprot_ac"].map(D.measured(act)).notna()
                sub = df[has] if seen else df[~has]
            v = pd.to_numeric(sub[f"{act}_prob"], errors="coerce").dropna().to_numpy()
            species_scores[(label, act)] = v
            counts[(label, act)] = len(v)

    fig, axs = stylia.create_figure(1, 3, width_ratios=[1, 1, VIOLIN_RATIO], width=1.0, height=0.46)
    plot_roc(axs.next(), data)
    plot_pr(axs.next(), data)
    plot_species(axs.next(), species_scores)
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

    say("\n  PER-FOLD SPREAD (the shaded band; min-to-max over 5 folds, NOT a standard error)")
    from sklearn.metrics import average_precision_score, roc_auc_score

    for act, (oof, cv) in data.items():
        per = [(roc_auc_score(g["hit"], g[OOF_COLUMN]), average_precision_score(g["hit"], g[OOF_COLUMN]))
               for _, g in oof.groupby("fold") if 0 < g["hit"].sum() < len(g)]
        r = [x[0] for x in per]
        a = [x[1] for x in per]
        say(f"    {ACT_LABEL[act]:<10} AUROC {min(r):.3f}-{max(r):.3f}   "
            f"PR {min(a):.3f}-{max(a):.3f}   over {len(per)} folds")

    say("\n  THE CURVE IS NOT THE NUMBER")
    for act, (oof, cv) in data.items():
        say(f"    {ACT_LABEL[act]:<10} pooled curve: AUROC {roc_auc_score(oof['hit'], oof[OOF_COLUMN]):.4f}"
            f"  AP {average_precision_score(oof['hit'], oof[OOF_COLUMN]):.4f}"
            f"   vs 5-seed {cv.roc_auc_clustered:.4f} / {cv.pr_auc_clustered:.4f}")
    say("    One out-of-fold realisation against the seed-averaged value. Both inside one SD.")
    say("    Quote the 5-seed number; a single seed carries ±0.004 of arbitrariness.")

    say("\n  WHERE THE SCORES LAND   (median, and % above that activator's cut)")
    for label, _sp, _seen in GROUPS:
        cells = []
        for act in D.ACTIVATORS:
            v = species_scores[(label, act)]
            above = (v >= D.BASE_RATE_THRESHOLD[act]).mean() * 100
            cells.append(f"{ACT_LABEL[act]} n={counts[(label, act)]:>5,} med {np.median(v):.3f} "
                         f"{above:>5.1f}%")
        say(f"    {label.replace(chr(10), ' '):<24} {'   '.join(cells)}")

    say("\n  CAVEATS")
    say("    - The band is the MIN-TO-MAX across 5 CV folds, not a confidence interval. With ~46")
    say("      positives per fold it is wide at high recall for sample-size reasons alone.")
    say("    - PR-AUC is reported beside AUROC, never AUROC alone: the estimator's gain was ~4x")
    say("      larger on PR-AUC, and AUROC alone called 4 of 6 arms 'no difference'.")
    say("    - The base rates differ (0.137 vs 0.246), so the two PR curves are NOT on a common")
    say("      scale. A higher PR curve is not automatically a better model.")
    say("    - `seen` is a DIFFERENT protein set per activator (1,677 vs 1,045).")
    say("    - Only the S. aureus `seen` group is validated. Kp, Ec and the Sa `unseen` group are")
    say("      extrapolations; nothing in panels 1-2 licenses them.")


if __name__ == "__main__":
    main()
