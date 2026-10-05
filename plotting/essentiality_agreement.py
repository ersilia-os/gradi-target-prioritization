"""Essentiality agreement at three levels: between screens, across species, between predictors.

    essentiality_agreement.png   1  two measured CRISPRi screens against each other
                                 2  the cross-species essentialome, 447 proteins x 12 genomes
                                 3  ProteomeLM-Ess ranked, with the top genes named

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

**Panel 3 -- what the language model actually ranks first.** Every K. pneumoniae protein ordered by
`proteomelm_ess`, with the top of the list named. The head is exactly what a reader should want it
to be: aminoacyl-tRNA synthetases (`leuS` 0.9997, `metG`, `aspS`, `argS`, `valS`, `glyS`), RNA
polymerase (`rpoB`, `rpoC`, `rpoD`), peptidoglycan (`murG`) and gyrase (`gyrB`). **Three of the top
twelve are consortium panel targets** -- `lpxL`, `lptG` and `gyrB` -- and they are marked, because
that is the one thing this panel says that a list of ribosomal genes would not.

**The score means something different in each species** and Kp is the hardest case:
`proteomelm_ess` on K. pneumoniae is `unseen_species` -- no *Klebsiella* is in the authors' 89
genomes. It is comparable WITHIN a species and never across one.

**`geptop_ess` is not drawn here but is still measured** -- the two predictors share no inputs, no
training data and no method, and they agree at rho +0.439 over the whole proteome and +0.663 where
Geptop is above zero. Geptop is 66.3% ties at exactly 0 on Kp, which is what holds the first number
down; a 0 there is `orthologs_none_essential` (a confident NON-essential call, 58.5%) or
`no_orthologs` (7.8%), separable only in `geptop_<sp>.tsv`. Both numbers are in the run log.

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
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli"}

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


#: How many of the top-ranked proteins to name on panel 3.
N_LABEL = 12


def plot_proteomelm_rank(ax, species: str, n_label: int) -> dict:
    """Every protein ranked by `proteomelm_ess`, with the head of the list named.

    Log x because the proteins worth naming are the first 0.2% of the rank order; on a linear axis
    the head is one pixel column. Labels are a leader-line ladder rather than text pinned to each
    marker -- the top dozen differ by thousandths (0.9997 down to 0.9946) and sit on top of one
    another, so anchored text smears into an unreadable block.

    Consortium panel members in the head are marked, via `src.interest.annotate()`. Matching is by
    gene symbol and Kp `gene_name` covers 63.4% of the proteome, so this flags what it can find and
    is not a claim that nothing else in the head is of interest."""
    from src import interest as I
    from src import proteomes as P

    d = E.load(species).merge(P.load(species)[["uniprot_ac", "gene_name"]], on="uniprot_ac")
    d = I.annotate(d)
    d["p"] = pd.to_numeric(d["proteomelm_ess"], errors="coerce")
    d = d.dropna(subset=["p"]).sort_values("p", ascending=False).reset_index(drop=True)
    d["rank"] = np.arange(1, len(d) + 1)

    ax.plot(d["rank"], d["p"], color=PAL.PRIMARY, lw=2.0, zorder=2)
    ax.fill_between(d["rank"], 0, d["p"], color=PAL.PRIMARY, alpha=0.18, linewidth=0, zorder=1)

    top = d.head(n_label)
    panel = top["is_interest"].to_numpy(bool)
    ax.scatter(top.loc[~panel, "rank"], top.loc[~panel, "p"], s=20, color=PAL.PRIMARY,
               zorder=4, linewidths=0.5, edgecolors="white")
    ax.scatter(top.loc[panel, "rank"], top.loc[panel, "p"], s=34, color=PAL.ACCENT,
               zorder=5, linewidths=0.6, edgecolors="white")

    y_top, y_bot = float(top["p"].iloc[0]), float(top["p"].iloc[-1]) - 0.30
    ladder = np.linspace(y_top, y_bot, len(top))
    for (_, r), y_lab in zip(top.iterrows(), ladder):
        name = r["gene_name"] if isinstance(r["gene_name"], str) and r["gene_name"] else r["uniprot_ac"]
        hit = bool(r["is_interest"])
        ax.annotate(
            f"{name}{' *' if hit else ''}  {r['p']:.4f}",
            xy=(r["rank"], r["p"]), xytext=(len(d) ** 0.40, y_lab), textcoords="data",
            ha="left", va="center", fontsize=SS * 0.70,
            color=PAL.ACCENT if hit else PAL.INK,
            arrowprops={"arrowstyle": "-", "lw": 0.6, "color": PAL.MUTED,
                        "shrinkA": 0, "shrinkB": 2},
        )

    ax.set_xscale("log")
    ax.set_xlim(1, len(d) * 1.05)
    ax.set_ylim(0, 1.04)
    ax.text(0.97, 0.06, "* consortium panel target", transform=ax.transAxes, ha="right",
            fontsize=SS * 0.72, color=PAL.ACCENT)
    stylia.label(ax, xlabel=f"Rank within {LABELS[species]} (log)",
                 ylabel="ProteomeLM-Ess p(essential)",
                 title="What the language model ranks first")
    return {"n": len(d), "top": top[["gene_name", "uniprot_ac", "p", "is_interest"]],
            "above_0_9": int((d["p"] > 0.9).sum()), "median": float(d["p"].median())}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=["kpneumoniae", "ecoli"],
                    help="species for the predictor panel (the other two are fixed by their data)")
    ap.add_argument("--label", type=int, default=N_LABEL, help="how many top proteins to name")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  essentiality agreement -- screens, species, predictors")
    fig, axs = stylia.create_figure(1, 3, width_ratios=[1.15, 1.0, 1.1], width=1.0, height=0.42)
    cr = plot_screens(axs.next())
    cs = plot_cross_species(axs.next())
    pr = plot_proteomelm_rank(axs.next(), args.species, args.label)
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

    say(f"\n  PROTEOMELM-ESS, TOP {args.label}  ({args.species})")
    for i, (_, r) in enumerate(pr["top"].iterrows(), 1):
        name = r["gene_name"] if isinstance(r["gene_name"], str) and r["gene_name"] else r["uniprot_ac"]
        say(f"    {i:>3}. {name:<10} {r['p']:.4f}"
            + ("   <- consortium panel target" if r["is_interest"] else ""))
    say(f"    median {pr['median']:.4f}   above 0.9: {pr['above_0_9']:,} of {pr['n']:,}")

    say(f"\n  THE OTHER PREDICTOR, measured but not drawn  ({args.species})")
    d = E.load(args.species)
    gx = pd.to_numeric(d["geptop_ess"], errors="coerce")
    gy = pd.to_numeric(d["proteomelm_ess"], errors="coerce")
    nz = gx > 0
    say(f"    Geptop ~ ProteomeLM-Ess   rho {gx.corr(gy, method='spearman'):+.3f} over everything")
    say(f"    {'':<26} rho {gx[nz].corr(gy[nz], method='spearman'):+.3f} where Geptop > 0 "
        f"(n={int(nz.sum()):,})")
    say(f"    Geptop is {100 * (gx == 0).mean():.1f}% ties at exactly 0 -- what holds the first down.")

    say("\n  CAVEATS")
    say("    - Panels 1 and 2 read LEGACY v1 tables (output/results/, 2026-07). They are the only")
    say("      route to two continuous screens on one key, and to the 12-genome matrix. Panel 3")
    say("      is v2 throughout. Neither legacy table is on HISTORY.md's do-not-trust list.")
    say("    - A `geptop_ess` of 0 is `orthologs_none_essential` (a confident NON-essential call,")
    say("      58.5% of Kp) or `no_orthologs` (7.8%); only geptop_<sp>.tsv separates them.")
    say("    - proteomelm_ess on Kp is `unseen_species` -- no Klebsiella is in the authors' 89")
    say("      genomes. Comparable WITHIN a species, never across one.")
    say("    - Panel-target matching is by gene symbol and Kp gene_name is 63.4%, so the marks")
    say("      show what can be found, not that nothing else in the head is of interest.")
    say("    - The cross-species matrix is binary per genome and covers the 3,170 Kp proteins")
    say("      with a compendium call, not the whole proteome.")


if __name__ == "__main__":
    main()
