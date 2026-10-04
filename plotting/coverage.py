"""What each axis can actually say about K. pneumoniae -- and where a zero means "we do not know".

    axis_coverage.png   A  per axis: informative / structurally-zero / unknown, as % of proteome
                        B  the ambiguity that matters most -- studiedness's five evidence tiers

The honesty slide, and the one most likely to earn the room's trust. Every axis in this project
ships ONE COMPLETE ROW PER PROTEIN, so no axis has missing rows -- but a 0 is not the same claim in
each. In function a 0 means "nothing is annotated"; in studiedness it is either `no_hit` (the
strongest novelty claim the axis makes) or `below_floor` (an in-scope hit under 40% identity), and
the deliverable cannot tell them apart. In pockets an NA means "could not look", never "looked and
found nothing" -- so never `fillna(0)` there.

Panel A therefore splits each bar three ways rather than drawing a single coverage number, because
a single number would be the exact overstatement this slide exists to prevent.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.
`stylia.set_style("ersilia")` MUST precede `NamedColors()`, or NC.plum raises AttributeError.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/coverage.py
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

from src import degradability as D  # noqa: E402
from src import essentiality as E  # noqa: E402
from src import function as FN  # noqa: E402
from src import ligandability as LG  # noqa: E402
from src import localization as LOC  # noqa: E402
from src import orthology as O  # noqa: E402
from src import pockets as PK  # noqa: E402
from src import studiedness as ST  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

# The five studiedness tiers, ordered strongest-evidence first. Two of them score 0 meaning
# different things, which is the whole point of panel B.
TIER_ORDER = ["swissprot_direct", "swissprot_close", "swissprot_homolog", "below_floor", "no_hit"]
TIER_COLOR = {
    "swissprot_direct": "#00A087",
    "swissprot_close": "#4DBBA5",
    "swissprot_homolog": "#9AD5C8",
    "below_floor": "#F39B7F",
    "no_hit": "#E64B35",
}


def axis_coverage(sp: str) -> pd.DataFrame:
    """One row per axis: how many proteins it is informative about, and how many it cannot speak to.

    'informative' is deliberately strict -- it is the count the axis would stand behind, not the
    count of non-null cells. Every axis has 100% non-null by construction."""
    rows = []

    loc = LOC.load(sp)
    rows.append(("localization", len(loc), 0, 0,
                 "DeepLocPro always returns a call: 100% is a property of the method"))

    fn = FN.load(sp)
    has_go = (fn["goslim_terms"].astype(str).str.len() > 0).sum()
    has_cog = (fn["cog_categories"].astype(str).str.len() > 0).sum()
    rows.append(("function (GO-slim)", has_go, 0, len(fn) - has_go,
                 "an empty list is NOT ANNOTATED, never 'ruled out'"))
    rows.append(("function (COG)", has_cog, 0, len(fn) - has_cog,
                 "NCBI's own curators reach 81.6% of E. coli; there is no headroom"))

    ort = O.load(sp)
    rows.append(("orthology", len(ort), 0, 0,
                 "OrthoFinder runs de novo on our FASTAs, so a 0 is a MEASURED 0"))

    ess = E.load(sp)
    gep = E.load_geptop(sp)
    no_orth = (gep["geptop_evidence"] == "no_orthologs").sum() if "geptop_evidence" in gep else 0
    rows.append(("essentiality", len(ess) - no_orth, 0, no_orth,
                 "geptop 0 is 'orthologs, none essential' (evidence) vs 'no orthologs' (absence)"))

    deg = D.load(sp)
    rows.append(("degradability", len(deg), 0, 0,
                 "every protein scored, but Kp rows are HYPOTHESES extrapolated from S. aureus"))

    st = ST.load(sp)
    tr = ST.load_transfer(sp)
    unknown = tr["evidence"].isin(["no_hit", "below_floor"]).sum()
    rows.append(("studiedness", len(st) - unknown, 0, unknown,
                 "a 0 is `no_hit` (novelty) or `below_floor` (unknown) -- the table cannot tell"))

    lg = LG.load(sp)
    assayed = (pd.to_numeric(lg["n_assayed_bacterial"], errors="coerce") > 0).sum()
    rows.append(("ligands", assayed, 0, len(lg) - assayed,
                 "0 potent against 0 assayed is an OPEN QUESTION, not a negative"))

    pk = PK.load(sp)
    modelled = pd.to_numeric(pk["af_plddt"], errors="coerce").notna().sum()
    rows.append(("pockets", modelled, 0, len(pk) - modelled,
                 "NA means 'could not look', never 'looked and found nothing' -- never fillna(0)"))

    df = pd.DataFrame(rows, columns=["axis", "informative", "structural", "unknown", "note"])
    df["n"] = df[["informative", "structural", "unknown"]].sum(axis=1)
    return df


def plot_coverage(ax, cov: pd.DataFrame, abc: str) -> None:
    """Stacked bars, informative vs unknown. Stacked rather than a single coverage number because a
    single number is exactly the overstatement this panel exists to prevent."""
    y = np.arange(len(cov))[::-1]
    pct_i = cov["informative"] / cov["n"] * 100
    pct_u = cov["unknown"] / cov["n"] * 100

    ax.barh(y, pct_i, color=NC.plum, height=0.62, label="axis can speak to this protein")
    ax.barh(y, pct_u, left=pct_i, color=NC.gray, height=0.62, label="a 0 here means UNKNOWN")
    for yi, pi in zip(y, pct_i):
        ax.text(101, yi, f"{pi:.0f}%", va="center", fontsize=SS, color=NC.black)

    ax.set_yticks(y)
    ax.set_yticklabels(cov["axis"], fontsize=SS)
    ax.set_xlim(0, 112)
    ax.legend(fontsize=SS, frameon=False, loc="lower left", bbox_to_anchor=(0, -0.22), ncol=2)
    stylia.label(ax, xlabel="% of proteome", ylabel="", title="What each axis can say", abc=abc)


def plot_tiers(ax, tr: pd.DataFrame, abc: str) -> None:
    """Studiedness's five tiers. The last two both read 0 in the shipped table, and only one of
    them is a novelty claim -- this is the axis where the ambiguity costs the most."""
    counts = tr["evidence"].value_counts()
    order = [t for t in TIER_ORDER if t in counts.index]
    vals = [counts[t] for t in order]
    y = np.arange(len(order))[::-1]

    ax.barh(y, vals, color=[TIER_COLOR[t] for t in order], height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v + max(vals) * 0.02, yi, f"{v:,}", va="center", fontsize=SS, color=NC.black)
    ax.set_yticks(y)
    ax.set_yticklabels(order, fontsize=SS)
    ax.set_xlim(0, max(vals) * 1.2)
    stylia.label(ax, xlabel="proteins", ylabel="",
                 title="Studiedness: both red tiers read 0", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  axis coverage -- {LABELS[sp]}")
    cov = axis_coverage(sp)
    tr = ST.load_transfer(sp)

    fig, axs = stylia.create_figure(1, 2, width_ratios=[3, 2], width=1.0, height=0.40)
    plot_coverage(axs.next(), cov, abc="A")
    plot_tiers(axs.next(), tr, abc="B")
    out = OUT_DIR / "axis_coverage.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  COVERAGE (informative / unknown, of proteome)")
    for _, r in cov.iterrows():
        say(f"    {r['axis']:<20} {r['informative']:>6,} / {r['unknown']:>6,}"
            f"   {r['informative'] / r['n'] * 100:>5.1f}%")
        say(f"    {'':<20} {r['note']}")
    say("\n  CAVEAT")
    say("    Every axis ships one complete row per protein -- no axis has MISSING rows. What")
    say("    differs is what a 0 claims. Panel A is about meaning, not about completeness.")


if __name__ == "__main__":
    main()
