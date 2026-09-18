"""Figures for the stage-04 degradability predictions, across all three species.

`plots/degradability.py` answers *how good is the model* on the one species that has labels.
This answers the different question: **what did it predict everywhere else, and how much should you
trust it there.**

Three figures:

  predictions.png    the score distributions -- per species, per activator, against the base-rate
                     threshold; and, for *S. aureus*, split by measured hit status, which is the
                     only place a predicted score can be read next to a real one
  agreement.png      adep4 vs onc212 predicted probability, per species. They correlate at
                     rho 0.83-0.86 while the underlying LABELS agree at only rho 0.52 / Jaccard
                     0.32 -- so "both activators agree" is close to one opinion, not two
  extrapolation.png  what the Kp/Ec numbers are worth: distance to the nearest training protein,
                     how each proteome distributes across the distance bands, and the measured
                     out-of-fold AUROC in the only two bands that carry one

*E. coli* and *K. pneumoniae* have **no measured labels at all** -- `_hit` is empty for every row --
so every point plotted for them is a ranking hypothesis. Nothing here is refitted; the tables and
`evidence/domain_bands.tsv` are read as written.

Run with the `gradi` env:
    python scripts/plots/degradability_predictions.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "degradability"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
SPECIES_COLOR = {"kpneumoniae": NC.plum, "ecoli": NC.orange, "saureus": NC.mint}

# The band edges stage 04 uses to price the extrapolation, from `evidence/domain_bands.tsv`.
BAND_EDGES = (0.80, 0.90, 0.95)


def violins(ax, groups: list[np.ndarray], colors: list, positions=None) -> None:
    """Violins with a mean marker, coloured per group. Shared by several panels."""
    pos = list(range(1, len(groups) + 1)) if positions is None else positions
    parts = ax.violinplot(groups, positions=pos, showextrema=False, widths=0.8)
    for body, color in zip(parts["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(1.0)
    ax.scatter(pos, [g.mean() for g in groups],
               color=NC.white, edgecolor=NC.black, zorder=3)


# ---------------------------------------------------------------- figure 1


def plot_score_by_species(ax, data: dict, activator: str, abc: str) -> None:
    """Predicted probability per species, against the base-rate threshold."""
    species = list(data)
    groups = [data[s][f"{activator}_prob"].to_numpy(float) for s in species]
    violins(ax, groups, [SPECIES_COLOR[s] for s in species])
    thr = D.BASE_RATE_THRESHOLD[activator]
    ax.axhline(thr, color=NC.gray, linestyle="--")
    ax.set_xticks(range(1, len(species) + 1))
    ax.set_xticklabels([LABELS[s] for s in species], style="italic")
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel=f"{activator} probability",
                 title=f"{activator}  ·  base-rate cut {thr}", abc=abc)


def plot_score_by_measured(ax, sa: pd.DataFrame, activator: str, abc: str) -> None:
    """*S. aureus* only: the score split by what the screens actually measured.

    The two measured groups are out-of-fold, so this is an honest read of the separation the model
    achieves -- and it puts the unmeasured proteins, which is all Kp and Ec are, on the same axis.
    """
    hit = sa[f"{activator}_hit"]
    prob = sa[f"{activator}_prob"]
    groups, names, colors = [], [], []
    for name, mask, color in (
        ("measured\nhit", hit == 1, NC.plum),
        ("measured\nnon-hit", hit == 0, NC.blue),
        ("unmeasured", hit.isna(), NC.gray),
    ):
        vals = prob[mask].to_numpy(float)
        if len(vals) < 2:
            continue
        groups.append(vals)
        names.append(f"{name}\nn={len(vals):,}")
        colors.append(color)
    violins(ax, groups, colors)
    ax.axhline(D.BASE_RATE_THRESHOLD[activator], color=NC.gray, linestyle="--")
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels(names)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel=f"{activator} probability",
                 title=f"{activator}  ·  S. aureus, by measured status", abc=abc)


# ---------------------------------------------------------------- figure 2


def plot_agreement(ax, df: pd.DataFrame, species: str, abc: str) -> float:
    """adep4 vs onc212 predicted probability. Hexbin, because n is thousands.

    `n` is in the title on purpose: this is EVERY row, so for *S. aureus* it mixes out-of-fold
    scores (the measured proteins) with pure predictions. Stage 04's own printed rho is over the
    predicted rows only, which for Sa is a different population (n=1,212, rho +0.866) and therefore
    a different number. Kp and Ec have no measured rows, so the two agree exactly there.
    """
    x = df["adep4_prob"].to_numpy(float)
    y = df["onc212_prob"].to_numpy(float)
    rho = float(spearmanr(x, y).statistic)
    ax.hexbin(x, y, gridsize=45, cmap=stylia.FadingColormap("plum").cmap, mincnt=1,
              linewidths=0)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="adep4 probability", ylabel="onc212 probability",
                 title=f"{LABELS[species]}  ·  rho {rho:+.3f}  (n={len(df):,})", abc=abc)
    return rho


# ---------------------------------------------------------------- figure 3


def plot_nn_similarity(ax, data: dict, abc: str) -> None:
    """Distance to the nearest training protein, per species -- the extrapolation's raw material."""
    species = list(data)
    groups = [data[s]["nn_similarity"].to_numpy(float) for s in species]
    violins(ax, groups, [SPECIES_COLOR[s] for s in species])
    for edge in BAND_EDGES:
        ax.axhline(edge, color=NC.gray, linestyle="--")
    ax.set_xticks(range(1, len(species) + 1))
    ax.set_xticklabels([LABELS[s] for s in species], style="italic")
    stylia.label(ax, xlabel="", ylabel="ESM-C cosine to nearest training protein",
                 title="Distance to the training set", abc=abc)


def plot_band_share(ax, bands: pd.DataFrame, activator: str, abc: str) -> dict:
    """How each proteome distributes across the bands, and how much of it carries no AUROC at all.

    A band with a `roc_auc` was measured on the labeled *S. aureus* out-of-fold set; the two nearest
    bands were too small to estimate one, so proteins there are **unpriced** -- the reweighted
    figure simply does not cover them.
    """
    sub = bands[(bands.activator == activator) & (bands.band != "REWEIGHTED")]
    priced = sub[sub.scope == "saureus_labeled_oof"].set_index("band")["roc_auc"]
    band_order = [b for b in priced.index]
    species = [s for s in LABELS if s in set(sub.scope)]

    bottom = np.zeros(len(species))
    uncovered = {}
    # One legend entry per *category*, not per band: there are four bands but only two kinds, and
    # labelling two of the four would imply the other two do not exist.
    labelled = set()
    for band in band_order:
        vals = np.array([
            float(sub[(sub.scope == s) & (sub.band == band)]["share"].sum()) for s in species
        ])
        has_auc = bool(np.isfinite(priced[band]))
        kind = "AUROC estimated" if has_auc else "no estimate"
        ax.bar([LABELS[s] for s in species], 100 * vals, bottom=100 * bottom,
               color=NC.plum if has_auc else NC.gray,
               label=None if kind in labelled else kind)
        labelled.add(kind)
        if not has_auc:
            for s, v in zip(species, vals):
                uncovered[s] = uncovered.get(s, 0.0) + v
        bottom += vals
    ax.set_ylim(0, 100)
    # The bars fill the axes, so any in-axes legend covers data; put it outside.
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2)
    stylia.label(ax, xlabel="", ylabel="% of proteome",
                 title=f"{activator}  ·  band composition", abc=abc)
    return uncovered


def plot_band_auroc(ax, bands: pd.DataFrame, abc: str) -> None:
    """The measured out-of-fold AUROC per band, with CI. Only two bands ever carry one."""
    oof = bands[(bands.scope == "saureus_labeled_oof") & bands.roc_auc.notna()]
    activators = list(D.ACTIVATORS)
    colors = {"adep4": NC.plum, "onc212": NC.orange}
    # Band labels are set once, outside the loop: the two activators have DIFFERENT n per band
    # (adep4 202/1468, onc212 118/923), so writing n into a shared tick label would show only
    # whichever activator happened to be drawn last. The n's go to stdout instead.
    bands_shown = list(oof[oof.activator == activators[0]]["band"])
    for i, act in enumerate(activators):
        sub = oof[oof.activator == act].set_index("band").reindex(bands_shown)
        x = np.arange(len(bands_shown)) + (i - 0.5) * 0.18
        lo = sub.roc_auc - sub.roc_auc_lo
        hi = sub.roc_auc_hi - sub.roc_auc
        ax.errorbar(x, sub.roc_auc, yerr=[lo, hi], fmt="o", color=colors[act],
                    capsize=4, label=act)
    ax.set_xticks(np.arange(len(bands_shown)))
    ax.set_xticklabels(bands_shown)
    ax.set_xlim(-0.5, len(bands_shown) - 0.5)
    ax.axhline(0.5, color=NC.gray, linestyle="--")
    ax.legend(loc="lower left")
    stylia.label(ax, xlabel="", ylabel="Out-of-fold AUROC",
                 title="Measured AUROC by band", abc=abc)


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(D.SPECIES), choices=list(D.SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    order = [s for s in LABELS if s in args.species]
    data = {s: D.load(s) for s in order}
    bands = D.load_domain_bands()
    activators = list(D.ACTIVATORS)

    # Figure 1 -- the scores themselves
    fig, axs = stylia.create_figure(len(activators), 2)
    for i, act in enumerate(activators):
        plot_score_by_species(axs.next(), data, act, abc="AC"[i])
        if "saureus" in data:
            plot_score_by_measured(axs.next(), data["saureus"], act, abc="BD"[i])
        else:
            axs.next().axis("off")
    out = OUT_DIR / "predictions.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")
    for act in activators:
        thr = D.BASE_RATE_THRESHOLD[act]
        for s in order:
            n_hit = len(D.hits(data[s], act))
            say(f"     {act:<7} {s:<14} median prob {data[s][f'{act}_prob'].median():.3f}"
                f"   at base-rate cut {thr}: {n_hit:>5} / {len(data[s]):>5} "
                f"({100 * n_hit / len(data[s]):>4.1f}%)")

    # Figure 2 -- the two activators are not independent evidence
    fig, axs = stylia.create_figure(1, len(order))
    for abc, s in zip("ABC", order):
        rho = plot_agreement(axs.next(), data[s], s, abc)
        say(f"     rho(adep4_prob, onc212_prob)  {s:<14} {rho:+.3f}")
    out = OUT_DIR / "agreement.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    # Figure 3 -- what the extrapolation is worth
    fig, axs = stylia.create_figure(1, 3)
    plot_nn_similarity(axs.next(), data, abc="A")
    uncovered = plot_band_share(axs.next(), bands, activators[0], abc="B")
    plot_band_auroc(axs.next(), bands, abc="C")
    out = OUT_DIR / "extrapolation.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")
    for s in order:
        say(f"     {s:<14} nn_similarity median {data[s]['nn_similarity'].median():.4f}"
            f"   unpriced share {100 * uncovered.get(s, 0.0):>5.1f}%")
    oof = bands[(bands.scope == "saureus_labeled_oof")]
    for _, r in oof.iterrows():
        auc = "no estimate" if not np.isfinite(r.roc_auc) else \
            f"AUROC {r.roc_auc:.4f} [{r.roc_auc_lo:.4f}, {r.roc_auc_hi:.4f}]"
        say(f"     band {r.activator:<7} {r.band:<12} n={int(r.n):>5} "
            f"pos={int(r.n_pos):>4}  {auc}")
    rew = bands[bands.band == "REWEIGHTED"]
    for _, r in rew.iterrows():
        say(f"     {r.activator:<7} {r.scope:<14} expected AUROC {r.expected_roc_auc:.4f}"
            f"   over {100 * r.share:.1f}% of the proteome")


if __name__ == "__main__":
    main()
