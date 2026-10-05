"""Essentiality agreement at three levels: between screens, across species, between predictors.

    essentiality_agreement.png   1  two measured CRISPRi screens against each other
                                 2  the cross-species essentialome, 447 proteins x 12 genomes
                                 3  the two independent predictors against each other

**Panel 1 -- do two experimental screens agree?** Rousset 2018 and Wang 2018 CRISPRi, both
genome-wide, both continuous, both E. coli, coloured by Keio's independent arrayed-knockout call.

**ONE CORRELATION NUMBER WOULD MISLEAD HERE, which is the point of the panel.** Over all 3,719
shared proteins Pearson is **+0.908** and Spearman only **+0.393**. That gap is the finding, not an
artifact to explain away: split by Keio, Spearman is **+0.726 among the 274 essential genes** and
**+0.252 among the 3,445 dispensable ones**. The screens agree strongly about what is ESSENTIAL and
barely at all about how dispensable the dispensable genes are -- and the essential tail (373
proteins below log2FC -2) is what anchors the Pearson line. Quote the stratified Spearmans; never
the Pearson alone.

**Panel 2 -- the cross-species essentialome.** Of 3,170 K. pneumoniae proteins with a measured call
in the 12-genome Enterobacteriaceae compendium, **447 are essential in at least one genome**: 159 in
all twelve, 65 in eleven, and 89 in exactly one. Sorted by conservation it separates into a solid
core block and a sparse species-specific tail, which is the shape that makes "essential" an
organism-specific word rather than a universal one -- and the reason this project never merges
screens into a single label.

**Panel 3 -- do the two independent predictors agree?** `geptop_ess` (orthology to a curated
reference set of essential genes) against `proteomelm_ess` (a protein language model head). They
share no inputs, no training data and no method.

**`geptop_ess` IS 66.3% TIES AT EXACTLY 0 on Kp**, which is the dense wall at x=0 and is what holds
its Spearman down. A 0 there means `orthologs_none_essential` -- a confident NON-essential call,
58.5% of the proteome -- or `no_orthologs`, 7.8%; the two are separable only in `geptop_<sp>.tsv`.

**PROVENANCE: PANELS 1 AND 2 READ LEGACY v1 TABLES, deliberately and exceptionally.**
`output/results/ecoli/ec_ess_experimental.csv` and
`output/results/kpneumoniae/kp_ess_publications.csv` (both 2026-07, from `legacy/scripts/07n_*` and
`07l_*`) are the only places in this repository where two experimental screens share a continuous
readout on one key, and where the 12-genome compendium is joined to the anchor proteome. v2's
training sets carry exactly ONE continuous column between them (Goodall's `insertion_index`) and no
cross-species matrix. Neither table is on `legacy/HISTORY.md`'s do-not-trust list. Panel 3 is v2.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/essentiality_agreement.py
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

from matplotlib.patches import Patch  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import essentiality as E  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# The two LEGACY (v1) inputs -- see the module docstring for why they are the only route.
CRISPRI_TABLE = REPO_ROOT / "output" / "results" / "ecoli" / "ec_ess_experimental.csv"
CROSS_TABLE = REPO_ROOT / "output" / "results" / "kpneumoniae" / "kp_ess_publications.csv"
CR_X, CR_Y = "ecoli_crispri_rousset18_log2fc", "ecoli_crispri_wang18_fitness"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL

#: Shortened genome labels for the matrix columns, in the compendium's own order.
GENOME_SHORT = {
    "pub_ess__K. pneumoniae ECL8": "Kp ECL8",
    "pub_ess__K. pneumoniae RH201207": "Kp RH201207",
    "pub_ess__E. coli BW25113": "Ec BW25113",
    "pub_ess__E. coli EC958": "Ec EC958",
    "pub_ess__E. coli NCTC13441": "Ec NCTC13441",
    "pub_ess__C. rodentium ICC168": "C. rodentium",
    "pub_ess__S. Typhi Ty2": "S. Typhi",
    "pub_ess__S. Tm A130": "S. Tm A130",
    "pub_ess__S. Tm D23580": "S. Tm D23580",
    "pub_ess__S. Tm SL3261": "S. Tm SL3261",
    "pub_ess__S. Tm SL1344": "S. Tm SL1344",
    "pub_ess__S. Enteritidis P125109": "S. Enteritidis",
}


def plot_screens(ax) -> dict:
    """Two continuous CRISPRi screens, coloured by Keio's independent call.

    Reported as TWO stratified Spearmans rather than one correlation: a single number either
    overstates the agreement (Pearson, anchored by the essential tail) or understates it where it
    matters (Spearman over a bulk that is mostly noise around zero)."""
    d = pd.read_csv(CRISPRI_TABLE).dropna(subset=[CR_X, CR_Y])
    keio = d["ecoli_keio_essential"].fillna(0).astype(bool)
    x, y = d[CR_X], d[CR_Y]

    ax.scatter(x[~keio], y[~keio], s=3.0, color=PAL.MUTED, alpha=0.45, linewidths=0,
               rasterized=True, label=f"Keio non-essential  {int((~keio).sum()):,}")
    ax.scatter(x[keio], y[keio], s=8.0, color=PAL.SECONDARY, alpha=0.85, linewidths=0,
               rasterized=True, label=f"Keio essential  {int(keio.sum()):,}")

    out = {"n": len(d), "spearman": x.corr(y, method="spearman"), "pearson": x.corr(y),
           "rho_ess": x[keio].corr(y[keio], method="spearman"),
           "rho_non": x[~keio].corr(y[~keio], method="spearman"),
           "n_ess": int(keio.sum()), "n_non": int((~keio).sum())}
    ax.text(0.03, 0.97, f"ρ {out['rho_ess']:+.2f}  essential\nρ {out['rho_non']:+.2f}  rest",
            transform=ax.transAxes, fontsize=SS * 0.82, color=PAL.INK, va="top")
    ax.legend(fontsize=SS * 0.74, frameon=False, loc="lower right", handletextpad=0.4,
              labelspacing=0.3, markerscale=2.2)
    stylia.label(ax, xlabel="Rousset 2018 CRISPRi (log2FC)",
                 ylabel="Wang 2018 CRISPRi (fitness)", title="Two screens, where they agree")
    return out


def plot_cross_species(ax) -> dict:
    """The 12-genome essentialome: every protein essential somewhere, sorted by conservation.

    All 447 rows rather than a hand-picked subset. The gradient from a solid core block to a sparse
    species-specific tail IS the result, and choosing 15 representative genes would be choosing the
    shape. No row labels at this height; the landmark genes go to the run log instead."""
    d = pd.read_csv(CROSS_TABLE)
    cols = [c for c in GENOME_SHORT if c in d.columns]
    cov = d[d["pub_covered"].fillna(0).astype(bool)]
    m = cov[cols].fillna(0).astype(int)
    n_ess = m.sum(axis=1)
    keep = n_ess > 0
    order = n_ess[keep].sort_values(ascending=False).index
    mat = m.loc[order].to_numpy()

    ax.imshow(mat, aspect="auto", cmap=PAL.SEQUENTIAL, interpolation="nearest", vmin=0, vmax=1)
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([GENOME_SHORT[c] for c in cols], rotation=55, ha="right",
                       fontsize=SS * 0.72)
    ax.set_yticks([])
    # A binary image with no key is a guessing game: name both states explicitly.
    ax.legend(handles=[Patch(facecolor=PAL.PRIMARY, label="essential"),
                       Patch(facecolor="white", edgecolor=PAL.MUTED, label="not essential")],
              fontsize=SS * 0.72, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.225), ncol=2, handletextpad=0.4, handlelength=1.2)
    counts = n_ess[keep].value_counts()
    stylia.label(ax, xlabel="", ylabel=f"{len(mat):,} proteins, sorted by conservation",
                 title="Essential where?")
    return {"covered": len(cov), "any": int(len(mat)), "all12": int(counts.get(len(cols), 0)),
            "one": int(counts.get(1, 0)), "dist": counts.sort_index().to_dict()}


def plot_predictors(ax, species: str) -> dict:
    """The two independent predictors against each other.

    They share no inputs, no training data and no method, which is why their agreement is worth a
    panel at all. Geptop's tie block at exactly 0 is the wall at x=0, annotated rather than hidden:
    it is two thirds of the anchor proteome and it is what holds the first rho down."""
    d = E.load(species)
    x = pd.to_numeric(d["geptop_ess"], errors="coerce")
    y = pd.to_numeric(d["proteomelm_ess"], errors="coerce")
    ax.scatter(x, y, s=3.0, color=PAL.PRIMARY, alpha=0.30, linewidths=0, rasterized=True)

    rho = x.corr(y, method="spearman")
    ties = (x == 0).mean() * 100
    nz = x > 0
    rho_nz = x[nz].corr(y[nz], method="spearman")
    ax.text(0.03, 0.97, f"ρ {rho:+.2f}  all\nρ {rho_nz:+.2f}  Geptop > 0",
            transform=ax.transAxes, fontsize=SS * 0.82, color=PAL.INK, va="top")
    ax.text(0.03, 0.82, f"{ties:.1f}% of Geptop is exactly 0", transform=ax.transAxes,
            fontsize=SS * 0.74, color=PAL.INK, va="top",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5, "alpha": 0.85})
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="Geptop (orthology)", ylabel="ProteomeLM-Ess (language model)",
                 title="Two predictors, no shared inputs")
    return {"rho": rho, "rho_nonzero": rho_nz, "ties": ties, "n_nonzero": int(nz.sum())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=["kpneumoniae", "ecoli"],
                    help="species for the predictor panel (the other two are fixed by their data)")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  essentiality agreement -- screens, species, predictors")
    fig, axs = stylia.create_figure(1, 3, width_ratios=[1.15, 1.0, 1.1], width=1.0, height=0.42)
    cr = plot_screens(axs.next())
    cs = plot_cross_species(axs.next())
    pr = plot_predictors(axs.next(), args.species)
    out = OUT_DIR / "essentiality_agreement.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  TWO CRISPRi SCREENS  (E. coli, both continuous, legacy v1 table)")
    say(f"    n={cr['n']:,} shared proteins")
    say(f"    Pearson, everything       {cr['pearson']:+.3f}   <- never quote this alone")
    say(f"    Spearman, everything      {cr['spearman']:+.3f}")
    say(f"    Spearman, Keio essential  {cr['rho_ess']:+.3f}  (n={cr['n_ess']:,})")
    say(f"    Spearman, the rest        {cr['rho_non']:+.3f}  (n={cr['n_non']:,})")
    say("    They agree about what is ESSENTIAL, not about how dispensable the rest is.")

    say("\n  CROSS-SPECIES ESSENTIALOME  (12 Enterobacteriaceae genomes, legacy v1 table)")
    say(f"    {cs['covered']:,} Kp proteins with a call   {cs['any']:,} essential in >=1 genome")
    say(f"    essential in all 12: {cs['all12']:,}      in exactly one: {cs['one']:,}")
    say("    n genomes -> proteins: " + "  ".join(f"{k}:{v}" for k, v in cs["dist"].items()))
    say("    A solid core block and a sparse species-specific tail: 'essential' is an")
    say("    organism-specific word, which is why this project never merges screens.")

    say(f"\n  TWO PREDICTORS  ({args.species})")
    say(f"    Geptop ~ ProteomeLM-Ess   rho {pr['rho']:+.3f} over everything")
    say(f"    {'':<26} rho {pr['rho_nonzero']:+.3f} where Geptop > 0 (n={pr['n_nonzero']:,})")
    say(f"    Geptop is {pr['ties']:.1f}% ties at exactly 0, which is what holds the first rho down.")

    say("\n  CAVEATS")
    say("    - Panels 1 and 2 read LEGACY v1 tables (output/results/, 2026-07). They are the only")
    say("      route to two continuous screens on one key, and to the 12-genome matrix. Panel 3")
    say("      is v2 throughout. Neither legacy table is on HISTORY.md's do-not-trust list.")
    say("    - A `geptop_ess` of 0 is `orthologs_none_essential` (a confident NON-essential call,")
    say("      58.5% of Kp) or `no_orthologs` (7.8%); only geptop_<sp>.tsv separates them.")
    say("    - The cross-species matrix is binary per genome and covers the 3,170 Kp proteins")
    say("      with a compendium call, not the whole proteome.")


if __name__ == "__main__":
    main()
