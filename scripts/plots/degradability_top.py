"""The top of the degradability ranking, named: which proteins, and what they are.

`plots/degradability_predictions.py` shows the score *distributions*. This shows the **proteins**,
because a ranking nobody can read is not a deliverable.

Three figures:

  top_adep4.png       the 20 highest-ranked proteins per species, labelled and coloured by COG group
  top_onc212.png      the same for the second activator
  top_composition.png what the top of the list is made of, against the proteome it came from, and
                      how much the two activators' shortlists actually share

Labels resolve `gene_name` (stage 00) -> eggNOG `preferred_name` (stage 02) -> accession. That order
matters on the anchor: Kp `gene_name` is 63.4% and eggNOG's `preferred_name` 67.0%, so neither alone
names the list. In practice the top of the ranking is well annotated -- high-scoring proteins are
small, conserved and abundant -- so accessions are rare there.

**Ranking, not thresholding.** The model was validated on AUROC, so the top-N is the defensible way
to read it; `src.degradability.hits()` and its base-rate cut are for when a cut is unavoidable.

For *S. aureus* a `*` marks a protein the screens actually measured as a hit: there, a high rank is
confirmation rather than prediction. Kp and Ec have no measured labels at all, so every bar is a
hypothesis.

Run with the `gradi` env:
    python scripts/plots/degradability_top.py
    python scripts/plots/degradability_top.py --top 30
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402
from src import function as F  # noqa: E402
from src import localization as LOC  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "degradability"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
SPECIES_COLOR = {"kpneumoniae": NC.plum, "ecoli": NC.orange, "saureus": NC.mint}

# COG's four top-level groups, plus the unclassified. Short keys because they are legend entries.
COG_GROUPS = (
    "INFORMATION STORAGE AND PROCESSING",
    "CELLULAR PROCESSES AND SIGNALING",
    "METABOLISM",
    "POORLY CHARACTERIZED",
)
COG_GROUP_SHORT = {
    "INFORMATION STORAGE AND PROCESSING": "Information",
    "CELLULAR PROCESSES AND SIGNALING": "Cellular processes",
    "METABOLISM": "Metabolism",
    "POORLY CHARACTERIZED": "Poorly characterised",
    "": "unclassified",
}
COG_GROUP_COLOR = {
    "INFORMATION STORAGE AND PROCESSING": NC.plum,
    "CELLULAR PROCESSES AND SIGNALING": NC.orange,
    "METABOLISM": NC.mint,
    "POORLY CHARACTERIZED": NC.blue,
    "": NC.gray,
}


def annotated(species: str) -> pd.DataFrame:
    """The prediction table joined to identity, function and compartment.

    Every join is on `uniprot_ac`, which every stage keys on, and every one is a left join from the
    predictions so a row can never be dropped by a missing annotation.
    """
    pred = D.load(species)
    ident = P.load(species)[["uniprot_ac", "gene_name", "protein_name"]]
    egg = F.load_eggnog(species)[["uniprot_ac", "preferred_name"]]
    cog = F.load_cog(species)[["uniprot_ac", "cog_category", "cog_group", "cog_name"]]
    loc = LOC.load_deeplocpro(species)[["uniprot_ac", "localization"]]
    df = (pred.merge(ident, on="uniprot_ac", how="left")
              .merge(egg, on="uniprot_ac", how="left")
              .merge(cog, on="uniprot_ac", how="left")
              .merge(loc, on="uniprot_ac", how="left"))
    # Empty string, not NaN, is how stages 00 and 02 record "no value" -- see docs/function.md.
    for col in ("gene_name", "preferred_name", "cog_group", "cog_name"):
        df[col] = df[col].fillna("").astype(str).str.strip()
    df["label"] = np.where(df.gene_name != "", df.gene_name,
                           np.where(df.preferred_name != "", df.preferred_name, df.uniprot_ac))
    return df


def top(df: pd.DataFrame, activator: str, n: int) -> pd.DataFrame:
    """The n highest-scoring proteins, best first."""
    return df.nlargest(n, f"{activator}_prob")


# ---------------------------------------------------------------- figure 1 / 2


def plot_top(ax, df: pd.DataFrame, species: str, activator: str, n: int, abc: str) -> pd.DataFrame:
    """Horizontal bars for the top n, labelled and coloured by COG group."""
    t = top(df, activator, n).iloc[::-1]  # best at the top of the axis
    y = np.arange(len(t))
    colors = [COG_GROUP_COLOR.get(g, NC.gray) for g in t.cog_group]
    ax.barh(y, t[f"{activator}_prob"], color=colors)

    measured_hit = t[f"{activator}_hit"] == 1
    labels = [f"{lab} *" if hit else lab for lab, hit in zip(t.label, measured_hit)]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.7, len(t) - 0.3)
    ax.set_xlim(0, float(t[f"{activator}_prob"].max()) * 1.08)
    stylia.label(ax, xlabel=f"{activator} probability", ylabel="",
                 title=LABELS[species], abc=abc)
    return t.iloc[::-1]


def _legend(ax, pairs, ncol: int) -> None:
    """Legend for a top-vs-proteome panel, including the tick that marks the proteome share.

    Without that last handle the black dashes are unexplained, which is the difference between a
    reader seeing an enrichment and a reader seeing decoration.
    """
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=c, label=lab) for lab, c in pairs]
    handles.append(Line2D([], [], color=NC.black, marker="_", linestyle="none",
                          label="proteome"))
    # Below the panel, not inside it: the tallest bars here reach 94%, so an in-axes legend covers
    # exactly the data the panel exists to show. Two columns keeps each legend inside its own
    # panel's width, so neighbouring panels' legends cannot collide.
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=ncol)


def cog_legend(ax) -> None:
    """One shared legend, drawn under a figure's panels."""
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=COG_GROUP_COLOR[g], label=COG_GROUP_SHORT[g])
               for g in (*COG_GROUPS, "")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=5)


# ---------------------------------------------------------------- figure 3


def plot_group_composition(ax, data: dict, activator: str, n: int, abc: str) -> None:
    """COG group share in the top n, against the same share across the whole proteome."""
    species = list(data)
    groups = (*COG_GROUPS, "")
    x = np.arange(len(species))
    for i, g in enumerate(groups):
        top_share, all_share = [], []
        for s in species:
            df = data[s]
            t = top(df, activator, n)
            top_share.append((t.cog_group == g).mean())
            all_share.append((df.cog_group == g).mean())
        off = (i - (len(groups) - 1) / 2) * 0.16
        ax.bar(x + off, 100 * np.array(top_share), 0.15,
               color=COG_GROUP_COLOR[g], label=COG_GROUP_SHORT[g])
        # The proteome share as a tick, so enrichment is read as bar-vs-tick, not bar-vs-bar.
        ax.scatter(x + off, 100 * np.array(all_share), marker="_",
                   color=NC.black, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in species], style="italic")
    _legend(ax, [(COG_GROUP_SHORT[g], COG_GROUP_COLOR[g]) for g in groups], ncol=2)
    stylia.label(ax, xlabel="", ylabel=f"% of top {n}",
                 title=f"{activator} · function", abc=abc)


def plot_localization_composition(ax, data: dict, activator: str, n: int, abc: str) -> None:
    """Compartment share in the top n, against the proteome. Clp has to reach these."""
    species = list(data)
    classes = [c for c in LOC.LOC_CLASSES]
    x = np.arange(len(species))
    for i, c in enumerate(classes):
        top_share, all_share = [], []
        for s in species:
            df = data[s]
            t = top(df, activator, n)
            top_share.append((t.localization == c).mean())
            all_share.append((df.localization == c).mean())
        off = (i - (len(classes) - 1) / 2) * 0.14
        ax.bar(x + off, 100 * np.array(top_share), 0.13,
               color=LOC.LOC_CLASS_COLOR[c], label=LOC.LOC_CLASS_ABBREV[c])
        ax.scatter(x + off, 100 * np.array(all_share), marker="_",
                   color=NC.black, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([LABELS[s] for s in species], style="italic")
    _legend(ax, [(LOC.LOC_CLASS_ABBREV[c], LOC.LOC_CLASS_COLOR[c]) for c in classes], ncol=2)
    stylia.label(ax, xlabel="", ylabel=f"% of top {n}",
                 title=f"{activator} · compartment", abc=abc)


def plot_overlap(ax, data: dict, abc: str) -> dict:
    """How much the two activators' shortlists share, as a function of list length.

    The dashed line is chance: drawing two independent top-N lists from M proteins gives an expected
    overlap of N^2/M. Rank correlation is global; this is what happens at the sharp end, which is
    where a shortlist is actually taken.
    """
    ns = np.unique(np.round(np.logspace(np.log10(10), np.log10(1000), 25)).astype(int))
    out = {}
    for s, df in data.items():
        a_rank = df.nlargest(len(df), "adep4_prob").uniprot_ac.to_numpy()
        o_rank = df.nlargest(len(df), "onc212_prob").uniprot_ac.to_numpy()
        frac = []
        for n in ns:
            frac.append(len(set(a_rank[:n]) & set(o_rank[:n])) / n)
        ax.plot(ns, 100 * np.array(frac), color=SPECIES_COLOR[s], label=LABELS[s])
        out[s] = len(set(a_rank[:100]) & set(o_rank[:100]))
    m = float(np.mean([len(df) for df in data.values()]))
    ax.plot(ns, 100 * ns / m, color=NC.gray, linestyle="--", label="chance")
    ax.set_xscale("log")
    ax.set_ylim(0, 100)
    ax.legend(loc="upper left")
    stylia.label(ax, xlabel="Top N", ylabel="% shared by both activators",
                 title="Shortlist overlap", abc=abc)
    return out


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(D.SPECIES), choices=list(D.SPECIES))
    ap.add_argument("--top", type=int, default=20, help="proteins per species in the named figures")
    ap.add_argument("--composition-top", type=int, default=100,
                    help="list length for the composition and enrichment panels")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    order = [s for s in LABELS if s in args.species]
    data = {s: annotated(s) for s in order}
    activators = list(D.ACTIVATORS)

    # Figures 1 and 2 -- the named top of each ranking
    for act in activators:
        fig, axs = stylia.create_figure(1, len(order))
        last = None
        for abc, s in zip("ABC", order):
            last = axs.next()
            t = plot_top(last, data[s], s, act, args.top, abc)
            named = int((t.label != t.uniprot_ac).sum())
            # How many of the top N are MEASURED hits rather than predictions. Only Sa can score
            # here, and it is the closest thing to a validation of the top of the list.
            confirmed = int((t[f"{act}_hit"] == 1).sum())
            measured = int(t[f"{act}_hit"].notna().sum())
            say(f"  {act}  {LABELS[s]}  top {args.top}   "
                f"{named}/{len(t)} named   "
                f"measured {measured}/{len(t)}, of which {confirmed} are confirmed hits   "
                f"prob {t[f'{act}_prob'].min():.3f}-{t[f'{act}_prob'].max():.3f}")
            say("      " + ", ".join(t.label.head(10)))
        cog_legend(last)
        out = OUT_DIR / f"top_{act}.png"
        stylia.save_figure(str(out))
        say(f"  -> {out.relative_to(REPO_ROOT)}")

    # Figure 3 -- what the top is made of, and how much the two lists share
    n = args.composition_top
    fig, axs = stylia.create_figure(1, 3)
    plot_group_composition(axs.next(), data, activators[0], n, abc="A")
    plot_localization_composition(axs.next(), data, activators[0], n, abc="B")
    overlap = plot_overlap(axs.next(), data, abc="C")
    out = OUT_DIR / "top_composition.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    for s in order:
        df, t = data[s], top(data[s], activators[0], n)
        say(f"     {s:<14} top{n} cytoplasmic {100 * (t.localization == 'cytoplasm').mean():>5.1f}%"
            f" (proteome {100 * (df.localization == 'cytoplasm').mean():>5.1f}%)"
            f"   information-storage {100 * (t.cog_group == COG_GROUPS[0]).mean():>5.1f}%"
            f" (proteome {100 * (df.cog_group == COG_GROUPS[0]).mean():>5.1f}%)")
    for s, k in overlap.items():
        say(f"     {s:<14} top{n} shared by both activators: {k}/{n}")


if __name__ == "__main__":
    main()
