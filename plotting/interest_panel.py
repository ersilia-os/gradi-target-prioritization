"""The consortium's own 43 proteins of interest, scored on every axis this project built.

    interest_panel.png   A  which of the 43 panel symbols are findable in K. pneumoniae
                         B  where the panel sits, as a percentile of the Kp proteome, per axis
                         C  the panel's degradability against the proteome's

**This is the slide for this audience, and it is not flattering.** The panel is excellent biology
and the wrong chemistry for a degrader: essential (far right tail), already well studied (NOT
novel), and trending DEPLETED on predicted degradability -- because it is almost entirely envelope,
and envelope proteins are either membrane-inserted co-translationally or exported, so activated
ClpP reaches them only through a pre-export window, if at all.

**`src/interest.py` is an EXPANSION OF PROSE, flagged as a draft.** There is no machine-readable
consortium list in this repository; the targets are stated in prose in two legacy documents, and
the kick-off's "targets in the periplasm" conflicts with the v5 proposal's cytosolic GyrA/GyrB --
never settled. Panel A exists so the naming gap is visible rather than silent: matching is by gene
symbol and Kp `gene_name` covers only 63.4% of the proteome, so a missing symbol is usually a
NAMING gap, not a biological absence.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`, which returns the ersilia plum/orange/mint set.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME:
    python plotting/interest_panel.py
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

from plotting import palette as PAL  # noqa: E402

from plotting import filters as F  # noqa: E402
from src import interest as I  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: ersilia. The STYLE sets typography and grid; the COLOURS come from
# `plotting/palette.py` -- stylia's npg palette, not stylia's ersilia palette.
stylia.set_format("slide")
stylia.set_style("ersilia")

SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

# (column, display, higher-is-more-of-the-thing). Percentiles are against the whole proteome.
AXES = [
    ("screens_ess_mean", "essentiality", True),
    ("n_papers_uniprot_prokaryotic", "studiedness", True),
    ("adep4_prob", "degradability (ADEP4)", True),
    ("onc212_prob", "degradability (ONC212)", True),
    ("bacterial_panel_orthologs", "conservation", True),
    ("p2rank_score", "pocket score", True),
]


def plot_found(ax, cov: pd.DataFrame, abc: str) -> None:
    """Found vs missing per family. Missing is a naming gap far more often than an absence."""
    y = np.arange(len(cov))[::-1]
    ax.barh(y, cov["n_found"], color=PAL.PRIMARY, height=0.6, label="found by symbol")
    ax.barh(y, cov["n_panel"] - cov["n_found"], left=cov["n_found"], color=PAL.MUTED,
            height=0.6, label="not matched")
    for yi, f_, n_ in zip(y, cov["n_found"], cov["n_panel"]):
        ax.text(n_ + 0.25, yi, f"{f_}/{n_}", va="center", fontsize=SS, color=PAL.INK)
    ax.set_yticks(y)
    ax.set_yticklabels(cov["family"], fontsize=SS * 0.95)
    ax.set_xlim(0, cov["n_panel"].max() * 1.25)
    ax.legend(fontsize=SS, frameon=False, loc="lower right")
    stylia.label(ax, xlabel="Panel members", ylabel="", title="Is the panel findable?", abc=abc)


def plot_percentiles(ax, panel: pd.DataFrame, full: pd.DataFrame, abc: str) -> None:
    """Where the panel sits on each axis, as a percentile of the proteome. The median is the dot;
    the bar is the interquartile range, so a wide family is visibly wide."""
    rows = []
    for col, name, _h in AXES:
        ref = pd.to_numeric(full[col], errors="coerce")
        vals = pd.to_numeric(panel[col], errors="coerce").dropna()
        if vals.empty:
            continue
        pct = vals.apply(lambda v: (ref < v).mean() * 100.0)
        rows.append((name, pct.quantile(0.25), pct.median(), pct.quantile(0.75)))

    y = np.arange(len(rows))[::-1]
    for yi, (_n, q1, med, q3) in zip(y, rows):
        color = PAL.SECONDARY if med < 50 else PAL.PRIMARY
        ax.plot([q1, q3], [yi, yi], color=color, lw=4, solid_capstyle="round", alpha=0.55)
        ax.plot([med], [yi], "o", color=color, markersize=7)
        ax.text(med, yi + 0.33, f"{med:.0f}", fontsize=SS, color=PAL.INK, ha="center")
    ax.axvline(50, color=PAL.INK, lw=1.0, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=SS)
    ax.set_xlim(0, 100)
    stylia.label(ax, xlabel="Percentile within proteome", ylabel="",
                 title="Right biology, wrong chemistry", abc=abc)


def plot_degradability(ax, panel: pd.DataFrame, full: pd.DataFrame, abc: str) -> None:
    """The panel's ADEP4 distribution against the proteome's, as survival curves.

    Survival rather than a histogram because the question is "how far up the ranking does this set
    reach", which a cumulative curve answers directly and a binned count does not."""
    for data, name, color in (
        (full["adep4_prob"], f"proteome (n={len(full):,})", PAL.MUTED),
        (panel["adep4_prob"], f"consortium panel (n={len(panel)})", PAL.PRIMARY),
    ):
        v = np.sort(pd.to_numeric(data, errors="coerce").dropna().to_numpy())[::-1]
        ax.step(v, np.arange(1, len(v) + 1) / len(v) * 100, where="post", color=color, lw=2,
                label=name)
    from src import degradability as D

    ax.axvline(D.BASE_RATE_THRESHOLD["adep4"], color=PAL.SECONDARY, lw=1.3, ls="--")
    ax.text(D.BASE_RATE_THRESHOLD["adep4"], 86, " base-rate cut", fontsize=SS, color=PAL.SECONDARY)
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    stylia.label(ax, xlabel="ADEP4 probability", ylabel="% of set at or above",
                 title="The panel is depleted, not enriched", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  consortium panel -- {LABELS[sp]}")
    full = F.load_joined(sp)
    cov = I.coverage(full)
    ann = I.annotate(full)
    panel = ann[ann["is_interest"]]
    say(f"  {cov['n_found'].sum()} of {cov['n_panel'].sum()} panel symbols matched; "
        f"{len(panel):,} proteome rows carry one")

    fig, axs = stylia.create_figure(1, 3, width_ratios=[2.6, 2.6, 2.8], width=1.0, height=0.40)
    plot_found(axs.next(), cov, abc="A")
    plot_percentiles(axs.next(), panel, full, abc="B")
    plot_degradability(axs.next(), panel, full, abc="C")
    out = OUT_DIR / "interest_panel.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  COVERAGE BY FAMILY")
    for _, r in cov.iterrows():
        say(f"    {r['family']:<26} {r['n_found']}/{r['n_panel']}"
            + (f"   missing: {r['missing']}" if r["missing"] else ""))

    say("\n  PANEL PERCENTILE WITHIN PROTEOME (median)")
    for col, name, _h in AXES:
        ref = pd.to_numeric(full[col], errors="coerce")
        vals = pd.to_numeric(panel[col], errors="coerce").dropna()
        if vals.empty:
            continue
        med = vals.apply(lambda v: (ref < v).mean() * 100.0).median()
        say(f"    {name:<24} {med:>5.1f}")

    say("\n  CAVEATS")
    say("    - PANEL is an EXPANSION OF PROSE and a draft, not the consortium's own file. Every")
    say("      judgement is in src/interest.py:EXPANSION_NOTES.")
    say("    - Matching is by gene symbol; Kp gene_name is 63.4%, so 'missing' is usually a")
    say("      NAMING gap, not a biological absence.")
    say("    - The panel is almost entirely ENVELOPE, which is why it is degradability-depleted:")
    say("      membrane proteins never present a soluble cytoplasmic chain.")


if __name__ == "__main__":
    main()
