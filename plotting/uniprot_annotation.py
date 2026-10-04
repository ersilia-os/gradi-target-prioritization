"""E. coli is annotated in UniProt; K. pneumoniae largely is not. The biology is the same.

    uniprot_annotation.png   left   `proteomes_evidence`, UniProt's identity ladder, 100% stacked
                             right  the fraction with a human ortholog

**The column is UniProt's own identity ladder**, built by `src/proteomes.py:identity_evidence()`:

    3  the entry carries its OWN identity -- SwissProt-reviewed, or a gene symbol on the anchor
       entry itself -- AND a specific protein name
    2  one of those two, not both
    1  neither: no own identity and no specific name

Measured: **Kp 2,205 / 2,469 / 1,054** against **Ec 0 / 648 / 3,755**. E. coli is **85.3%** at the
top tier where K. pneumoniae is **18.4%**, and K. pneumoniae has **2,205 proteins at tier 1 where
E. coli has none at all**.

**Why this column and not `is_reviewed`.** The blunt version of this fact is that E. coli is 100%
SwissProt-reviewed and K. pneumoniae is **7 proteins of 5,728**. That is a more dramatic number and
a worse column: it is degenerate per species -- all-or-nothing -- which is exactly why
`src/proteomes.py` replaced it with this three-level ladder. The ladder shows the same contrast and
still separates proteins *within* K. pneumoniae, which a boolean that is false 99.9% of the time
cannot. The raw boolean survives in `evidence/proteome_full_<species>.tsv` via `load_full()`.

**THE SECOND PANEL IS THE CONTROL, and it is what makes the first one a statement about curation
rather than about proteomes.** If K. pneumoniae were simply a stranger organism, its biology should
look different too. It does not: **16.6% of Kp has a human ortholog against 19.0% of Ec**, at a
median identity of **33.5% against 33.0%** -- the same fraction of the same kind of protein. So the
4.6x gap in annotation tier 3 (18.4% against 85.3%) is a fact about how much work has been done,
not about what is there to find.

Human orthology is the union of OrthoFinder and RBH, and the error direction is stated in
`selectivity.png`: under-detecting human homology makes a target look more selective than it is.

**This is about CURATION, not about biology.** A tier 1 protein is not a worse protein; it is a
protein nobody has written an identity for. K. pneumoniae HS11286 is a dark TrEMBL proteome and
E. coli K-12 MG1655 is the most annotated bacterial genome there is, so this panel measures the
history of who studied what, which is the same confound `studiedness_essentiality.png` prices.

100% stacked rather than raw counts because the two proteomes differ in size (5,728 vs 4,403) and
the claim is compositional; the counts are annotated inside the bars so the absolute numbers are
still readable.

Everything is read from the stage's own table through `src/proteomes.py`; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/uniprot_annotation.py
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

from src import orthology as O  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
SPECIES = ["kpneumoniae", "ecoli"]
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli"}

#: Tier -> (legend label, colour). Same ramp as the studiedness tiers in
#: `studiedness_essentiality.png`, so a reader who has seen one bar reads this one the same way.
TIERS = {
    1: ("1  no own identity, no specific name", PAL.NPG[0]),
    2: ("2  one of the two", PAL.NPG[1]),
    3: ("3  own identity and a specific name", PAL.NPG[4]),
}


def plot_annotation(ax, data: dict) -> None:
    """The identity ladder, 100% stacked, counts inside the bars.

    Drawn bottom-up weakest-first so the eye reads the bar in the same direction as the ladder,
    and so the two proteomes' tier-3 blocks share a top edge and can be compared at a glance."""
    x = np.arange(len(SPECIES))
    bottom = np.zeros(len(SPECIES))
    for tier, (label, color) in TIERS.items():
        pct, counts = [], []
        for sp in SPECIES:
            n = int((data[sp]["proteomes_evidence"] == tier).sum())
            counts.append(n)
            pct.append(n / len(data[sp]) * 100)
        ax.bar(x, pct, bottom=bottom, width=0.6, color=color, label=label,
               edgecolor="white", linewidth=0.9)
        for xi, v, b, n in zip(x, pct, bottom, counts):
            if v >= 5:
                ax.text(xi, b + v / 2, f"{n:,}", ha="center", va="center",
                        fontsize=SS * 0.9, color="white")
        bottom += np.array(pct)

    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[sp] for sp in SPECIES], fontsize=SS, style="italic")
    ax.set_xlim(-0.6, len(SPECIES) - 0.4)
    ax.set_ylim(0, 100)
    ax.set_box_aspect(1)
    ax.legend(fontsize=SS * 0.78, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.12), handletextpad=0.4, labelspacing=0.3)
    stylia.label(ax, xlabel="", ylabel="% of proteome",
                 title="UniProt identity, by evidence tier")


#: The comparability metrics, in reading order: how connected the protein is to other bacteria,
#: then to human, then to itself. Each is (label, how to compute it from the two frames).
CONSERVATION = [
    ("In an orthogroup", lambda d, n: n["in_orthogroup"].astype(bool)),
    ("Conserved in >=50% of the panel",
     lambda d, n: pd.to_numeric(d["bacterial_panel_orthologs"], errors="coerce") >= 0.5),
    ("Core: >=90% of the panel",
     lambda d, n: pd.to_numeric(d["bacterial_panel_orthologs"], errors="coerce") >= 0.9),
    ("Has a human ortholog", lambda d, n: d["has_human_ortholog"].astype(bool)),
    ("Has a paralog", lambda d, n: pd.to_numeric(n["n_paralogs"], errors="coerce") > 0),
]


def plot_conservation(ax, ortho: dict, dense: dict) -> dict:
    """Five ways of asking "is this a normal bacterial proteome", for each species.

    The control for the panel beside it. If K. pneumoniae were simply a stranger organism its
    biology should look different too; across all five measures it does not, while its annotation
    differs 4.6-fold. Horizontal bars because the labels are sentences.

    **`bacterial_panel_orthologs` IS A FRACTION, NOT A COUNT**, despite the name -- 0 to 1 over the
    28 bacterial proteomes (26 tier-C comparators plus the three anchors, minus this protein's own
    species). It counts SPECIES, never proteins, so a paralog pair cannot inflate it."""
    y = np.arange(len(CONSERVATION))[::-1]
    height = 0.36
    vals: dict[str, list[float]] = {}
    for i, sp in enumerate(SPECIES):
        v = [100 * f(ortho[sp], dense[sp]).mean() for _label, f in CONSERVATION]
        vals[sp] = v
        ax.barh(y + (0.5 - i) * height, v, height=height, color=PAL.SPECIES_COLOR[sp],
                label=LABELS[sp])
        for yi, x in zip(y + (0.5 - i) * height, v):
            ax.text(x + 1.4, yi, f"{x:.1f}", va="center", fontsize=SS * 0.78, color=PAL.INK)

    ax.set_yticks(y)
    ax.set_yticklabels([label for label, _f in CONSERVATION], fontsize=SS * 0.82)
    ax.set_xlim(0, 112)
    ax.set_box_aspect(1)
    ax.legend(fontsize=SS * 0.82, frameon=False, loc="lower right", handletextpad=0.5)
    stylia.label(ax, xlabel="% of proteome", ylabel="", title="Orthology analysis")
    return vals


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  UniProt annotation -- K. pneumoniae against E. coli")
    data = {sp: P.load(sp) for sp in SPECIES}
    ortho = {sp: O.load(sp) for sp in SPECIES}

    # Narrow: two square panels, not a slide's width. `set_box_aspect(1)` on each squares the axes
    # box directly rather than shrinking it inside its slot.
    dense = {sp: O.load_dense(sp) for sp in SPECIES}
    fig, axs = stylia.create_figure(1, 2, width=0.72, height=0.40)
    plot_annotation(axs.next(), data)
    cons = plot_conservation(axs.next(), ortho, dense)
    out = OUT_DIR / "uniprot_annotation.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  IDENTITY TIERS  (proteomes_evidence)")
    for sp in SPECIES:
        d = data[sp]
        c = d["proteomes_evidence"].value_counts().sort_index()
        say(f"    {LABELS[sp]:<16} n={len(d):>5,}   " + "  ".join(
            f"{t}:{int(c.get(t, 0)):>5,} ({100 * c.get(t, 0) / len(d):>5.1f}%)" for t in (1, 2, 3)))

    say("\n  THE SAME FACT IN THREE OTHER COLUMNS")
    full = {sp: P.load_full(sp) for sp in SPECIES}
    for sp in SPECIES:
        d, f = data[sp], full[sp]
        rev = f["is_reviewed"].astype(str).str.lower().isin(["true", "1"]).mean() * 100
        gn = (d["gene_name"].astype(str).str.len() > 0).mean() * 100
        say(f"    {LABELS[sp]:<16} SwissProt-reviewed {rev:>5.1f}%   "
            f"gene symbol {gn:>5.1f}%   tier 3 "
            f"{100 * (d['proteomes_evidence'] == 3).mean():>5.1f}%")
    say("    `is_reviewed` is the blunter number (Ec 100% against Kp 7 of 5,728) and the worse")
    say("    column: degenerate per species, which is why proteomes_evidence replaced it.")

    say("\n  CONSERVATION  (the control: same biology, different curation)")
    say(f"    {'':<34} {'Kp':>8} {'Ec':>8}   diff")
    for i, (label, _f) in enumerate(CONSERVATION):
        k, e = cons["kpneumoniae"][i], cons["ecoli"][i]
        say(f"    {label:<34} {k:>7.1f}% {e:>7.1f}%   {e - k:+5.1f} pp")
    for sp in SPECIES:
        bp = pd.to_numeric(ortho[sp]["bacterial_panel_orthologs"], errors="coerce")
        idp = pd.to_numeric(
            dense[sp].loc[dense[sp]["has_human_ortholog"].astype(bool), "best_identity_human"],
            errors="coerce")
        say(f"    {LABELS[sp]:<16} median panel fraction {bp.median():.3f}   "
            f"median identity to the human ortholog {idp.median():.1f}%")
    say("    Every measure within ~9 pp, against a 4.6x gap in annotation tier 3. The annotation")
    say("    gap is about how much work has been done, not about what is there to find.")
    say("    `bacterial_panel_orthologs` is a FRACTION over 28 bacterial proteomes, not a count.")

    say("\n  CAVEATS")
    say("    - This is CURATION, not biology. A tier 1 protein is not a worse protein; it is one")
    say("      nobody has written an identity for.")
    say("    - Kp HS11286 is a dark TrEMBL proteome and Ec K-12 MG1655 is the best-annotated")
    say("      bacterial genome there is, so this measures the history of who studied what.")
    say("    - It is the same confound studiedness_essentiality.png prices: annotation depth")
    say("      tracks essentiality, so 'unannotated' is not a clean proxy for 'novel target'.")


if __name__ == "__main__":
    main()
