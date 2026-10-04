"""Unexplored AND essential: the quadrant the collaboration is actually looking for.

    novelty_vs_essentiality.png   A  every Kp protein: literature against essentiality
                                  B  the confound, priced -- essentiality by studiedness tier

**The confound is real and it is printed, not hidden.** Essentiality and studiedness correlate:
`geptop_ess` runs rho 0.38-0.50 against the literature counts, and alone among the axes it SURVIVES
stratification (0.35-0.39). People have studied essential genes more. So the top-left quadrant is
partly an artifact of that correlation, and panel B shows how much by splitting essentiality across
the five studiedness evidence tiers: if the quadrant were pure artifact, the tiers would separate
completely.

**A 0 on the x-axis is ambiguous in the shipped table** -- it is `no_hit` (nothing among 575,748
curated entries resembles this protein, the strongest novelty claim this project makes) or
`below_floor` (an in-scope hit under the 40% identity floor). Only `load_transfer()` distinguishes
them, so this figure joins it back rather than reading the zero at face value.

The x-axis is a PAPER COUNT -- curated SwissProt references, nothing scaled or blended. A 0-1
composite shipped once and was removed; do not reintroduce one.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
`stylia.set_style("ersilia")` MUST precede `NamedColors()`, or NC.plum raises AttributeError.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/novelty.py
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

from plotting import filters as F  # noqa: E402
from src import interest as I  # noqa: E402
from src import studiedness as ST  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

TIER_ORDER = ["swissprot_direct", "swissprot_close", "swissprot_homolog", "below_floor", "no_hit"]
TIER_COLOR = {
    "swissprot_direct": "#00A087", "swissprot_close": "#4DBBA5",
    "swissprot_homolog": "#9AD5C8", "below_floor": "#F39B7F", "no_hit": "#E64B35",
}
PAPERS = "n_papers_uniprot_prokaryotic"


def plot_scatter(ax, df, surv, panel, cut, abc: str) -> None:
    """The map. x is log1p papers because the distribution spans four orders of magnitude and a
    linear axis would put 95% of the proteome in one pixel column."""
    x = np.log1p(pd.to_numeric(df[PAPERS], errors="coerce").fillna(0))
    ax.scatter(x, df["screens_ess_mean"], s=1.6, color=NC.gray, linewidths=0, rasterized=True)

    xp = np.log1p(pd.to_numeric(panel[PAPERS], errors="coerce").fillna(0))
    ax.scatter(xp, panel["screens_ess_mean"], s=11, color=NC.blue, linewidths=0,
               label=f"consortium panel (n={len(panel)})")
    xs = np.log1p(pd.to_numeric(surv[PAPERS], errors="coerce").fillna(0))
    ax.scatter(xs, surv["screens_ess_mean"], s=16, color=NC.plum, linewidths=0,
               edgecolors="none", label=f"shortlist (n={len(surv)})")

    ax.axhline(cut, color=NC.black, lw=1.0, ls="--")
    ax.axvline(np.log1p(5), color=NC.black, lw=1.0, ls="--")
    ax.text(0.06, 0.965, "unexplored + essential", transform=ax.transAxes, fontsize=SS,
            color=NC.black, va="top")
    ticks = [0, 1, 5, 20, 100, 500]
    ax.set_xticks(np.log1p(ticks))
    ax.set_xticklabels([str(t) for t in ticks], fontsize=SS)
    ax.legend(fontsize=SS, frameon=False, loc="lower right")
    stylia.label(ax, xlabel="curated papers (prokaryotic donor)", ylabel="essentiality (screens)",
                 title="Where the opportunity is", abc=abc)


def plot_confound(ax, df, abc: str) -> None:
    """Essentiality by studiedness tier. If novelty were pure artifact these would separate
    cleanly; they overlap heavily, which is why the quadrant survives as a lead."""
    present = [t for t in TIER_ORDER if (df["evidence"] == t).sum() >= 30]
    data = [df.loc[df["evidence"] == t, "screens_ess_mean"].dropna() for t in present]
    bp = ax.boxplot(data, vert=True, patch_artist=True, widths=0.6, showfliers=False)
    for patch, t in zip(bp["boxes"], present):
        patch.set_facecolor(TIER_COLOR[t])
        patch.set_edgecolor(NC.black)
        patch.set_linewidth(0.8)
    for el in ("medians", "whiskers", "caps"):
        for ln in bp[el]:
            ln.set_color(NC.black)
            ln.set_linewidth(0.9)
    ax.set_xticklabels([f"{t}\nn={(df['evidence'] == t).sum():,}" for t in present],
                       fontsize=SS * 0.8, rotation=25, ha="right")
    stylia.label(ax, xlabel="", ylabel="essentiality (screens)",
                 title="The confound, priced", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  novelty vs essentiality -- {LABELS[sp]}")
    df = F.load_joined(sp)
    tr = ST.load_transfer(sp)[["uniprot_ac", "evidence"]]
    df = df.merge(tr, on="uniprot_ac", how="left")
    surv, _steps = F.cascade(df)
    panel = I.annotate(df)
    panel = panel[panel["is_interest"]]
    cut = float(df["screens_ess_mean"].quantile(F.ESSENTIAL_PERCENTILE / 100.0))

    fig, axs = stylia.create_figure(1, 2, width_ratios=[3, 2], width=1.0, height=0.40)
    plot_scatter(axs.next(), df, surv, panel, cut, abc="A")
    plot_confound(axs.next(), df, abc="B")
    out = OUT_DIR / "novelty_vs_essentiality.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    quad = df[(df["screens_ess_mean"] >= cut) & (pd.to_numeric(df[PAPERS], errors="coerce").fillna(0) <= 5)]
    say(f"\n  UNEXPLORED + ESSENTIAL   {len(quad):,} proteins (<=5 papers, essentiality top decile)")
    say(f"    of which on the shortlist  {len(set(quad['uniprot_ac']) & set(surv['uniprot_ac'])):,}")
    rho = (np.log1p(pd.to_numeric(df[PAPERS], errors="coerce").fillna(0))
           .corr(df["screens_ess_mean"], method="spearman"))
    say(f"\n  CONFOUND  rho(papers, essentiality) = {rho:.3f}")
    say("    Essentiality is the UNPRICED confound on this axis. Note the column matters:")
    say("    geptop_ess runs rho 0.38-0.50 against studiedness and SURVIVES stratification")
    say("    (0.35-0.39), while screens_ess_mean -- the column this deck ranks on -- is far")
    say("    weaker, measured above. Ranking on screens_ess_mean therefore buys a less")
    say("    circular quadrant, but stacking novelty and essentiality still double-counts.")
    say("\n  CAVEAT")
    say("    A 0 on the x-axis is `no_hit` (real novelty) OR `below_floor` (unknown). The shipped")
    say("    table cannot tell them apart; panel B joins load_transfer() back to do so.")


if __name__ == "__main__":
    main()
