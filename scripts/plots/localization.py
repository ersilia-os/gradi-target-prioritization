"""Figures for the stage-03 localization predictions.

Three figures, each answering one question the stage raises:

  composition.png   what the three proteomes are made of, by compartment
  cross_check.png   do the two predictors agree? DeepLocPro calls a compartment, TMbed counts
                    residues -- they share no machinery, so their agreement is the stage's real
                    validation, and this is where it is visible
  gram_remap.png    how much of *S. aureus*'s `extracellular` is the Gram-positive remap rather
                    than the model

Nothing is recomputed except the raw-argmax labels in `gram_remap`, which come from the un-remapped
probability vector stage 03 persisted for exactly this purpose.

Compartment colours come from `src.localization.LOC_CLASS_COLOR`, not from a stylia palette: the
class-to-colour mapping is pinned in the vocabulary so figures cannot drift between v1 and v2, or
between figures here. Bars and boxes are ordered `LOC_CLASSES`, i.e. inward -> outward through the
cell envelope, which is also the order Clp reachability falls in.

Run with the `gradi` env:
    python scripts/plots/localization.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import localization as LOC  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "localization"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}


def frame(species: str) -> pd.DataFrame:
    """One species' compartment call joined to its TMbed topology.

    `confidence` and `has_signal_peptide` are no longer in the deliverable -- joined back from
    `load_deeplocpro()` and `load_tmbed()`, because the confidence-by-compartment box and the
    signal-peptide panels below are the whole point of the extracellular cross-check."""
    return (LOC.load(species)
            .merge(LOC.load_deeplocpro(species)[["uniprot_ac", "confidence"]], on="uniprot_ac")
            .merge(LOC.load_tmbed(species)[["uniprot_ac", "has_signal_peptide"]], on="uniprot_ac")
            .merge(LOC.load_topology(species), on="uniprot_ac"))


def present(df: pd.DataFrame) -> list[str]:
    """The compartments actually called, in envelope order.

    A compartment missing here is a structural zero, not missing data: `cell_wall_surface` cannot
    occur in a Gram-negative, and `periplasm`/`outer_membrane` are remapped away in a Gram-positive.
    """
    return [c for c in LOC.LOC_CLASSES if (df["localization"] == c).any()]


def colors_for(classes: list[str]) -> list[str]:
    return [LOC.LOC_CLASS_COLOR[c] for c in classes]


def abbrevs(classes: list[str]) -> list[str]:
    return [LOC.LOC_CLASS_ABBREV[c] for c in classes]


def violin_by_compartment(ax, pooled: pd.DataFrame, column: str, classes: list[str],
                          legend_loc: str = "best") -> None:
    """One violin per compartment, plus a mean marker.

    A violin rather than a box because the distribution that matters most here is **bimodal**:
    `extracellular` has a median near 0.96 and a mean of 0.58, so a box reports a "typical" value
    that describes almost none of the class and hides the split that is the whole finding.
    """
    nc = stylia.NamedColors()
    data = [pooled.loc[pooled["localization"] == c, column].to_numpy() for c in classes]
    parts = ax.violinplot(data, showextrema=False, widths=0.85)
    for body, color in zip(parts["bodies"], colors_for(classes)):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(1.0)
    ax.scatter(range(1, len(classes) + 1), [d.mean() for d in data],
               color=nc.white, edgecolor=nc.black, zorder=3, label="mean")
    ax.set_xticks(range(1, len(classes) + 1))
    ax.set_xticklabels(abbrevs(classes))
    ax.legend(loc=legend_loc)


def box_by_compartment(ax, pooled: pd.DataFrame, column: str, classes: list[str]) -> None:
    """One box per compartment, coloured by compartment."""
    data = [pooled.loc[pooled["localization"] == c, column].to_numpy() for c in classes]
    bp = ax.boxplot(data, patch_artist=True, showfliers=False,
                    medianprops={"color": "white"})
    for patch, color in zip(bp["boxes"], colors_for(classes)):
        patch.set_facecolor(color)
        patch.set_edgecolor(color)
    ax.set_xticks(range(1, len(classes) + 1))
    ax.set_xticklabels(abbrevs(classes))


# ---------------------------------------------------------------- figure 1


def plot_composition(ax, counts: pd.DataFrame, classes: list[str]) -> None:
    """Stacked horizontal bars: what each proteome is made of."""
    species = list(LABELS)
    left = np.zeros(len(species))
    for c in classes:
        vals = np.array([
            counts.loc[(counts.species == s) & (counts.localization == c), "pct"].sum()
            for s in species
        ])
        ax.barh(species, vals, left=left, color=LOC.LOC_CLASS_COLOR[c],
                label=LOC.LOC_CLASS_ABBREV[c])
        left += vals
    ax.set_yticks(range(len(species)))
    ax.set_yticklabels([LABELS[s] for s in species], style="italic")
    ax.set_xlim(0, 100)
    ax.legend(handles=[Patch(facecolor=LOC.LOC_CLASS_COLOR[c],
                             label=LOC.LOC_CLASS_ABBREV[c]) for c in classes],
              loc="center left", bbox_to_anchor=(1.0, 0.5))
    stylia.label(ax, xlabel="% of proteome", ylabel="", title="Compartment composition", abc="A")


def plot_confidence(ax, pooled: pd.DataFrame, classes: list[str]) -> None:
    """DeepLocPro's own confidence, per compartment. Extracellular is the weak one."""
    box_by_compartment(ax, pooled, "confidence", classes)
    ax.set_ylim(0, 1)
    stylia.label(ax, xlabel="", ylabel="DeepLocPro confidence",
                 title="Confidence by compartment", abc="B")


# ---------------------------------------------------------------- figure 2


def plot_cyto_fraction(ax, pooled: pd.DataFrame, classes: list[str]) -> None:
    """TMbed's residue count, per DeepLocPro compartment -- the independent cross-check."""
    violin_by_compartment(ax, pooled, "cytoplasmic_fraction", classes,
                          legend_loc="lower left")
    ax.set_ylim(-0.05, 1.05)
    stylia.label(ax, xlabel="", ylabel="Cytoplasmic fraction",
                 title="TMbed residues vs DeepLocPro call", abc="A")


def plot_signal_peptide(ax, pooled: pd.DataFrame, classes: list[str]) -> None:
    """Share carrying a signal peptide: says *why* a fraction is near zero."""
    pct = [100 * pooled.loc[pooled["localization"] == c, "has_signal_peptide"].mean()
           for c in classes]
    ax.bar(abbrevs(classes), pct, color=colors_for(classes))
    ax.set_ylim(0, 100)
    stylia.label(ax, xlabel="", ylabel="% with signal peptide",
                 title="Export signal by compartment", abc="B")


def plot_topology(ax, pooled: pd.DataFrame, classes: list[str]) -> None:
    """Mean TM helices vs mean TM strands: helices in the membrane, strands in the barrel."""
    nc = stylia.NamedColors()
    x = np.arange(len(classes))
    hel = [pooled.loc[pooled["localization"] == c, "n_tm_helix"].mean() for c in classes]
    strd = [pooled.loc[pooled["localization"] == c, "n_tm_strand"].mean() for c in classes]
    ax.bar(x - 0.2, hel, 0.4, color=nc.plum, label="TM helices")
    ax.bar(x + 0.2, strd, 0.4, color=nc.orange, label="TM strands")
    ax.set_xticks(x)
    ax.set_xticklabels(abbrevs(classes))
    ax.legend()
    stylia.label(ax, xlabel="", ylabel="Mean segments per protein",
                 title="Membrane topology by compartment", abc="C")


# ---------------------------------------------------------------- figure 3


def raw_argmax(species: str) -> pd.Series:
    """The un-remapped argmax label, from the persisted six-vector."""
    p = LOC.load_probabilities(species)
    mat = p[[f"p_{c}" for c in LOC.DLP_CANON]].to_numpy(float)
    return pd.Series([LOC.DLP_CANON[i] for i in mat.argmax(1)],
                     index=p["uniprot_ac"].to_numpy())


def plot_remap_counts(ax, species: str) -> dict:
    """Raw argmax vs shipped label, per compartment. The gap is the remap."""
    nc = stylia.NamedColors()
    raw = raw_argmax(species)
    final = LOC.load_deeplocpro(species).set_index("uniprot_ac")["localization"]
    raw = raw.reindex(final.index)

    classes = [c for c in LOC.LOC_CLASSES if (raw == c).any() or (final == c).any()]
    x = np.arange(len(classes))
    r = [int((raw == c).sum()) for c in classes]
    f = [int((final == c).sum()) for c in classes]
    ax.bar(x - 0.2, r, 0.4, color=nc.gray, label="raw argmax")
    ax.bar(x + 0.2, f, 0.4, color=nc.plum, label="shipped label")
    ax.set_xticks(x)
    ax.set_xticklabels(abbrevs(classes))
    ax.set_yscale("log")
    ax.legend()
    stylia.label(ax, xlabel="", ylabel="Proteins (log)",
                 title="Gram-positive remap: before and after", abc="A")
    return {"changed": int((raw != final).sum()),
            "raw_ext": int((raw == "extracellular").sum()),
            "final_ext": int((final == "extracellular").sum())}


def plot_moved_mass(ax, species: str) -> dict:
    """How much probability mass the remap moved into `extracellular`, per protein."""
    nc = stylia.NamedColors()
    p = LOC.load_probabilities(species)
    moved = (p["p_periplasm"] + p["p_outer_membrane"]).to_numpy(float)
    ax.hist(moved, bins=40, color=nc.plum)
    ax.set_yscale("log")
    ax.axvline(0.1, color=nc.orange, linestyle="--")
    stylia.label(ax, xlabel="Probability mass moved into extracellular",
                 ylabel="Proteins (log)", title="Most of it is near zero", abc="B")
    return {"mean": float(moved.mean()), "over": int((moved > 0.1).sum()),
            "max": float(moved.max())}


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(LOC.SPECIES), choices=list(LOC.SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    frames = {sp: frame(sp) for sp in args.species}
    pooled = pd.concat([f.assign(species=sp) for sp, f in frames.items()], ignore_index=True)
    counts = pd.concat([LOC.load_counts(sp) for sp in args.species], ignore_index=True)
    classes = present(pooled)

    # Figure 1 -- composition and confidence
    fig, axs = stylia.create_figure(1, 2, width_ratios=[3, 2])
    plot_composition(axs.next(), counts, classes)
    plot_confidence(axs.next(), pooled, classes)
    out = OUT_DIR / "composition.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")
    for sp, f in frames.items():
        top = f["localization"].value_counts(normalize=True).mul(100).round(1)
        say(f"     {sp:<14} " + "  ".join(
            f"{LOC.LOC_CLASS_ABBREV[c]} {top.get(c, 0.0):>4.1f}%" for c in present(f)))

    # Figure 2 -- the two predictors, cross-checked
    fig, axs = stylia.create_figure(1, 3)
    plot_cyto_fraction(axs.next(), pooled, classes)
    plot_signal_peptide(axs.next(), pooled, classes)
    plot_topology(axs.next(), pooled, classes)
    out = OUT_DIR / "cross_check.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")
    for c in classes:
        g = pooled[pooled["localization"] == c]
        say(f"     {LOC.LOC_CLASS_ABBREV[c]:<5} n={len(g):>5}  cyto_frac {g.cytoplasmic_fraction.mean():.3f}"
            f"  signal {100 * g.has_signal_peptide.mean():>5.1f}%"
            f"  helices {g.n_tm_helix.mean():>4.1f}  strands {g.n_tm_strand.mean():>4.1f}")
    # The doc's suggested flag for "probably actually cytoplasmic", quantified: no export signal and
    # every residue on the cytoplasmic side is not a protein that got secreted.
    ext = pooled[pooled["localization"] == "extracellular"]
    disputed = ext[(~ext.has_signal_peptide) & (ext.cytoplasmic_fraction > 0.5)]
    say(f"     TMbed disputes {len(disputed)}/{len(ext)} "
        f"({100 * len(disputed) / len(ext):.1f}%) extracellular calls "
        "(no signal peptide, cyto_frac > 0.5)")
    for sp, f in frames.items():
        e = f[f["localization"] == "extracellular"]
        d = e[(~e.has_signal_peptide) & (e.cytoplasmic_fraction > 0.5)]
        say(f"       {sp:<14} {len(d):>4}/{len(e):<4} ({100 * len(d) / len(e):>5.1f}%)")

    # Figure 3 -- the Gram-positive remap, S. aureus only
    if "saureus" in args.species:
        fig, axs = stylia.create_figure(1, 2)
        a = plot_remap_counts(axs.next(), "saureus")
        b = plot_moved_mass(axs.next(), "saureus")
        out = OUT_DIR / "gram_remap.png"
        stylia.save_figure(str(out))
        say(f"  -> {out.relative_to(REPO_ROOT)}")
        say(f"     {a['changed']} labels changed; extracellular {a['raw_ext']} raw -> "
            f"{a['final_ext']} shipped "
            f"(+{100 * (a['final_ext'] - a['raw_ext']) / a['raw_ext']:.1f}%)")
        say(f"     moved mass mean {b['mean']:.4f}, above 0.1 for {b['over']} proteins, "
            f"max {b['max']:.3f}")


if __name__ == "__main__":
    main()
