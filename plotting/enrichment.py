"""What kind of protein does activated ClpP reach? Function and compartment, enriched and depleted.

    enrichment.png   left   COG functional categories, ordered by effect
                     right  compartment and membrane topology

Fisher enrichment of the TOP DECILE of each proteome's predictions, against the rest of that
proteome. Both activators are shown because agreement between them is the only internal check this
figure has -- but they are NOT two opinions: the probabilities correlate at rho 0.89 while the
measured labels agree at only 0.52.

**The compartment panel carries the finding.** `wholly cytoplasmic` is the strongest enrichment in
the whole table, and every membrane and export feature is depleted: `TM helix`, `signal peptide`,
`beta-barrel`, inner membrane, periplasm. But `extracellular` is ENRICHED (Kp 3.80 under ONC212,
Ec 2.23, Sa 2.08) -- and that is not a contradiction, it is the mechanism. A secreted protein
transits the cytoplasm unfolded and IS reachable; a membrane protein is inserted co-translationally
and never presents a soluble chain. So the right filter is **"not membrane"**, never
"cytoplasm only", which would discard the compartment with the highest measured ONC212 hit rate.

**Read `n=` before reading any ratio.** The largest effect in the table is COG `B`, chromatin
structure, at 27x under ADEP4 and 64x under ONC212 -- on **8 members**, 6 and 7 of which are hits.
Every row therefore carries its group size in the tick label, and dot area is proportional to it,
so a category with 8 members cannot look like one with 800.

**Categories with NO significant result under either activator are not drawn** -- they carry no
evidence in either direction, and on a slide they cost the height that makes the rest readable
(Kp: 6 of 25, `A D O R V W`). `--all` keeps them, and the run log names them with their numbers
either way, so nothing is hidden, only undrawn.

**Categories below `MIN_ENRICHMENT_GROUP` members are dropped for a different and harder reason,
and that one is a correctness fix, not tidying.** `log2_or_adj` applies a Haldane +0.5 continuity correction so that groups with a zero
cell still have a finite logarithm; on a group of ONE with zero hits the correction outvotes the
data and the sign flips -- COG `Z` (cytoskeleton, 1 member, 0 hits, raw odds ratio 0) comes out at
**+1.58, reading as enriched**. No row with 5 or more members flips. Dropped rows are named in the
run log with their counts.

**THE TWO LARGEST ENRICHMENTS ARE BOTH SHORT-PROTEIN CATEGORIES, AND THAT IS NOT A COINCIDENCE.**
Measured on Kp, 2026-10-04: `spearman(length, adep4_prob) = -0.677`, and the median length is
**88 aa for `unclassified` and 97 aa for COG `B`, against 297 aa for a COG-classified protein**.
Stratifying by length deciles, `unclassified`'s top-decile rate ratio falls from **3.85 to 2.64** --
attenuated by a third, so a real residual survives, but a third of that bar is protein size wearing
a category's name. COG assignment itself fails more often on short proteins, so `unclassified` is
partly a length bin. Treat the top two rows of the function panel as the weakest in it, despite
being the largest.

**These are enrichments of PREDICTIONS.** For K. pneumoniae and E. coli they inherit the whole
S. aureus -> Gram-negative extrapolation. The yardstick is the MEASURED S. aureus hit rate by
compartment -- cytoplasm 0.180/0.275, membrane 0.031/0.105, extracellular 0.049/0.346 -- and the
right question is whether the predictions track those, not whether they track the previous model.

Everything is read from the stage's own table through `src/degradability.py`; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`, which returns the ersilia plum/orange/mint set.
`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/enrichment.py
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

from plotting import palette as PAL  # noqa: E402

from src import degradability as D  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
ACT_COLOR = {"adep4": PAL.PRIMARY, "onc212": PAL.SECONDARY}

# Dot area in points^2 for the smallest and largest group. Area, not radius, scales with group
# size -- radius would make a 10x bigger category look 100x more important.
AREA_MIN, AREA_MAX = 14.0, 230.0

# Short, readable names for the compartment/topology rows. The table's own labels are column
# expressions (`n_tm_helix > 0`) which are precise and unreadable on a slide.
LOC_PRETTY = {
    "Cyt": "cytoplasm",
    "CM": "inner membrane",
    "Peri": "periplasm",
    "OM": "outer membrane",
    "CW": "cell wall",
    "Ext": "extracellular",
    "signal peptide": "signal peptide",
    "TM helix": "any TM helix",
    "beta-barrel": "beta-barrel (8+ strands)",
    "wholly cytoplasmic": "wholly cytoplasmic",
}


def _areas(n: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Dot area on a sqrt scale between the smallest and largest group in the figure.

    sqrt rather than linear because group sizes span 8 to 3,828 here: linear would collapse
    everything below ~500 members into the same invisible dot."""
    if hi <= lo:
        return np.full_like(n, (AREA_MIN + AREA_MAX) / 2, dtype=float)
    f = (np.sqrt(n) - np.sqrt(lo)) / (np.sqrt(hi) - np.sqrt(lo))
    return AREA_MIN + f * (AREA_MAX - AREA_MIN)


def plot_block(ax, df: pd.DataFrame, order: list[str], pretty: dict[str, str],
               title: str, size_lo: float, size_hi: float, show_legend: bool) -> None:
    """One dot plot: a row per category, a dot per activator.

    Filled = significant after BH correction within its own test family; hollow = not, which is a
    statement about evidence and not about direction. A connecting line joins the two activators so
    that a category where they DISAGREE is visible as a long segment rather than as two dots the
    eye has to pair up itself."""
    y = {k: i for i, k in enumerate(order)}

    # the segment first, so dots sit on top of it
    for key in order:
        s = df[df["key"] == key]
        if len(s) == 2:
            ax.plot(s["log2_or_adj"], [y[key]] * 2, color="#CFCFCB", lw=1.4, zorder=1)

    for act, color in ACT_COLOR.items():
        s = df[df["activator"] == act]
        if s.empty:
            continue
        yy = np.array([y[k] for k in s["key"]])
        a = _areas(s["group_n"].to_numpy(float), size_lo, size_hi)
        sig = s["significant"].to_numpy(bool)
        ax.scatter(s["log2_or_adj"][sig], yy[sig], s=a[sig], color=color,
                   linewidths=0.6, edgecolors="white", zorder=3)
        ax.scatter(s["log2_or_adj"][~sig], yy[~sig], s=a[~sig], facecolors="white",
                   edgecolors=color, linewidths=1.2, zorder=3)

    ax.axvline(0, color=PAL.INK, lw=1.0, ls="--", zorder=0)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(
        [f"{pretty.get(k, k)}   n={int(df.loc[df['key'] == k, 'group_n'].iloc[0]):,}"
         for k in order],
        fontsize=SS * 0.82,
    )
    ax.set_ylim(-0.8, len(order) - 0.2)
    ax.grid(axis="x", color="#EFEFEC", lw=0.8, zorder=0)
    ax.set_axisbelow(True)

    if show_legend:
        handles = [
            Line2D([], [], marker="o", ls="", markersize=6, color=ACT_COLOR["adep4"], label="ADEP4"),
            Line2D([], [], marker="o", ls="", markersize=6, color=ACT_COLOR["onc212"], label="ONC212"),
            Line2D([], [], marker="o", ls="", markersize=6, markerfacecolor="white",
                   markeredgecolor=PAL.INK, color="none", label="not significant"),
        ]
        ax.legend(handles=handles, fontsize=SS * 0.85, frameon=False, loc="lower right",
                  handletextpad=0.3, labelspacing=0.3, borderpad=0.1)

    stylia.label(ax, xlabel="Depleted   <-   log2 odds ratio   ->   Enriched", ylabel="",
                 title=title)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("--min-group", type=int, default=D.MIN_ENRICHMENT_GROUP,
                    help="drop categories smaller than this (see the Haldane sign flip)")
    ap.add_argument("--all", action="store_true",
                    help="also draw categories with no significant result under either activator")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  degradability enrichment -- {LABELS[sp]}")

    cog_all = D.load_enrichment(kind="cog", species=sp)
    loc_all = D.load_enrichment(kind="loc", species=sp)
    cog = cog_all[cog_all["group_n"] >= args.min_group]
    loc = loc_all[loc_all["group_n"] >= args.min_group]

    dropped = cog_all[cog_all["group_n"] < args.min_group]
    if len(dropped):
        say(f"\n  DROPPED below {args.min_group} members (Haldane correction outvotes the data)")
        for _, r in dropped.drop_duplicates("key").iterrows():
            say(f"    COG {r['key']:<14} n={int(r['group_n'])}  hits={int(r['a_group_top'])}"
                f"  raw OR={r['odds_ratio']:.2f}  but log2_or_adj={r['log2_or_adj']:+.2f}")

    if not args.all:
        keep = cog.groupby("key")["significant"].any()
        mute = sorted(keep[~keep].index)
        if mute:
            say(f"\n  NOT DRAWN -- no significant result under either activator ({len(mute)} of "
                f"{keep.size} COG categories; --all keeps them)")
            for k in mute:
                r = cog[(cog["key"] == k) & (cog["activator"] == "adep4")].iloc[0]
                o = cog[(cog["key"] == k) & (cog["activator"] == "onc212")].iloc[0]
                say(f"    COG {k:<3} n={int(r['group_n']):<5} hits={int(r['a_group_top']):<4}"
                    f" log2 OR {r['log2_or_adj']:+.2f} / {o['log2_or_adj']:+.2f}   {r['label'][:44]}")
        cog = cog[cog["key"].isin(keep[keep].index)]
        lkeep = loc.groupby("key")["significant"].any()
        loc = loc[loc["key"].isin(lkeep[lkeep].index)] if lkeep.any() else loc

    top_n = int(cog["top_n"].iloc[0])
    say(f"\n  hits = top {top_n:,} of the proteome (top decile), against the rest")

    # Order each block by mean effect so the eye finds the extremes; enriched at the TOP, so the
    # reading direction matches the x-axis arrow.
    cog_order = cog.groupby("key")["log2_or_adj"].mean().sort_values().index.tolist()
    loc_order = loc.groupby("key")["log2_or_adj"].mean().sort_values().index.tolist()

    lo = min(cog["group_n"].min(), loc["group_n"].min())
    hi = max(cog["group_n"].max(), loc["group_n"].max())

    cog_pretty = {k: f"{k}  {v}" if len(k) == 1 else v
                  for k, v in zip(cog["key"], cog["label"])}

    fig, axs = stylia.create_figure(1, 2, width_ratios=[3.1, 2.3], width=1.0, height=0.46)
    plot_block(axs.next(), cog, cog_order, cog_pretty,
               "Function", lo, hi, show_legend=False)
    plot_block(axs.next(), loc, loc_order, LOC_PRETTY,
               "Compartment and topology", lo, hi, show_legend=True)
    out = OUT_DIR / "enrichment.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    for name, block, pretty in (("COMPARTMENT", loc, LOC_PRETTY), ("FUNCTION", cog, cog_pretty)):
        say(f"\n  {name}  (log2 OR, * = significant)")
        piv = block.pivot_table(index="key", columns="activator",
                                values=["log2_or_adj", "significant", "group_n", "a_group_top"],
                                aggfunc="first")
        order = loc_order if name == "COMPARTMENT" else cog_order
        for k in reversed(order):
            n = int(piv[("group_n", "adep4")][k])
            hits = f"{int(piv[('a_group_top', 'adep4')][k]):>4}/{n:<5}"
            cells = "  ".join(
                f"{piv[('log2_or_adj', a)][k]:+6.2f}{'*' if piv[('significant', a)][k] else ' '}"
                for a in ("adep4", "onc212")
            )
            say(f"    {pretty.get(k, k)[:44]:<46} {hits} {cells}")

    say("\n  CAVEATS")
    say("    - These are enrichments of PREDICTIONS, not of measurements. For Kp and Ec they")
    say("      inherit the S. aureus -> Gram-negative extrapolation.")
    say("    - The two activators are NOT two opinions: rho 0.89 on the probabilities against")
    say("      0.52 on the measured labels. Agreement here is weaker evidence than it looks.")
    say("    - Judge against the MEASURED hit rates (cytoplasm 0.180/0.275, membrane")
    say("      0.031/0.105, extracellular 0.346 under ONC212), never against the previous model.")
    say("    - Both estimators track protein LENGTH far more strongly than the labels do")
    say("      (-0.61 against -0.33), so some of this is length wearing a category's name.")
    say("    - CONCRETELY, for the two largest enrichments (measured on Kp, 2026-10-04):")
    say("      median length is 88 aa for `unclassified` and 97 aa for COG B, against 297 aa")
    say("      for a classified protein, and spearman(length, adep4_prob) = -0.677. Within")
    say("      length deciles `unclassified`'s rate ratio falls 3.85 -> 2.64. Real, but a third")
    say("      of it is size. They are the LARGEST bars and the WEAKEST evidence in the panel.")


if __name__ == "__main__":
    main()
