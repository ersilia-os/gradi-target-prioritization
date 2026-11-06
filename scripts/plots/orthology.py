"""Figures for stage 05: how far orthology reaches, how similar the pairs are, and where the two
ortholog callers disagree.

    coverage.png   what fraction of each proteome has an ortholog in each other proteome, and the
                   orthogroup-size distribution -- the panel that says where transfer can reach
    identity.png   identity per species pair. Kp-Ec near 86% against Kp-human near 37% is the whole
                   selectivity argument in one axis
    methods.png    OrthoFinder against reciprocal best hit, pair by pair

Everything is read from the stage's own tables; nothing is recomputed.

Run with the `gradi` env:
    python scripts/plots/orthology.py
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
from src import orthology as O  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "orthology"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli",
          "saureus": "S. aureus", "human": "H. sapiens"}
SHORT = {"kpneumoniae": "Kp", "ecoli": "Ec", "saureus": "Sa", "human": "Hs"}
METHOD_COLOR = {"orthofinder": NC.plum, "rbh": NC.orange, "either": NC.mint}


def species_present(dense: pd.DataFrame) -> list[str]:
    return [s for s in O.SPECIES if s in set(dense.species)]


# ---------------------------------------------------------------- figure 1


def plot_reach(ax, dense: pd.DataFrame, species: list[str], abc: str) -> pd.DataFrame:
    """% of each proteome with >=1 ortholog in each other proteome, by method.

    The diagonal is skipped: orthology is a between-species relation, so a species has no orthologs
    in itself -- that cell would be a category error, not a zero.
    """
    rows = []
    for q in species:
        d = dense[dense.species == q]
        for t in species:
            if q == t:
                continue
            for method, col in (("orthofinder", f"n_orthologs_of_{t}"),
                                ("rbh", f"n_orthologs_rbh_{t}"),
                                ("either", f"n_orthologs_{t}")):
                rows.append({"query": q, "target": t, "method": method,
                             "pct": 100 * (d[col] > 0).mean()})
    tab = pd.DataFrame(rows)
    pairs = [(q, t) for q in species for t in species if q != t]
    x = np.arange(len(pairs))
    for i, method in enumerate(("orthofinder", "rbh", "either")):
        vals = [tab[(tab["query"] == q) & (tab.target == t)
                    & (tab.method == method)].pct.iloc[0] for q, t in pairs]
        ax.bar(x + (i - 1) * 0.27, vals, 0.26, color=METHOD_COLOR[method], label=method)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{SHORT[q]}→{SHORT[t]}" for q, t in pairs])
    ax.set_ylim(0, 100)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    stylia.label(ax, xlabel="", ylabel="% of proteome with an ortholog",
                 title="How far orthology reaches", abc=abc)
    return tab


def plot_orthogroup_sizes(ax, groups: pd.DataFrame, abc: str) -> None:
    """Orthogroup size distribution, and how many species each group spans."""
    counts = groups.n_species.value_counts().sort_index()
    ax.bar(counts.index.astype(int), counts.values, color=NC.plum)
    ax.set_xticks(counts.index.astype(int))
    ax.set_yscale("log")
    stylia.label(ax, xlabel="species in the orthogroup", ylabel="orthogroups (log)",
                 title=f"{len(groups):,} orthogroups", abc=abc)


# ---------------------------------------------------------------- figure 2


def plot_identity(ax, nb: pd.DataFrame, species: list[str], abc: str) -> None:
    """Identity of the reciprocal-best-hit pairs, per species pair.

    RBH rather than all neighbours: a top-5 list includes distant hits by construction, so its
    identity distribution says more about the cut-off than about the biology.
    """
    pairs = [(q, t) for i, q in enumerate(species) for t in species[i + 1:]]
    data, names = [], []
    for q, t in pairs:
        v = nb[(nb.query_species == q) & (nb.target_species == t) & nb.is_rbh].pident.dropna()
        if len(v) < 5:
            continue
        data.append(v.to_numpy())
        names.append(f"{SHORT[q]}–{SHORT[t]}\nn={len(v):,}")
    parts = ax.violinplot(data, showextrema=False, widths=0.85)
    for body in parts["bodies"]:
        body.set_facecolor(NC.plum)
        body.set_edgecolor(NC.plum)
        body.set_alpha(1.0)
    ax.scatter(range(1, len(data) + 1), [np.median(d) for d in data],
               color=NC.white, edgecolor=NC.black, zorder=3, label="median")
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels(names)
    ax.set_ylim(0, 100)
    ax.legend(loc="upper right")
    stylia.label(ax, xlabel="", ylabel="% identity (RBH pairs)",
                 title="Sequence identity by species pair", abc=abc)


def plot_norm_bitscore(ax, nb: pd.DataFrame, species: list[str], abc: str) -> None:
    """Normalised bitscore against identity, coloured by whether the pair is cross-species.

    `bitscore_norm` divides by the query's self-bitscore, so it is comparable across proteins of
    very different length in a way raw bitscore is not.
    """
    cross = nb[nb.query_species != nb.target_species]
    within = nb[nb.query_species == nb.target_species]
    for frame, color, label in ((within, NC.gray, "same species (paralogs)"),
                                (cross, NC.plum, "cross-species")):
        s = frame.sample(min(len(frame), 20000), random_state=0)
        ax.scatter(s.pident, s.bitscore_norm, color=color, alpha=0.25, linewidths=0,
                   rasterized=True, label=label)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="upper left")
    stylia.label(ax, xlabel="% identity", ylabel="bitscore / self-bitscore",
                 title="Identity is not the whole story", abc=abc)


# ---------------------------------------------------------------- figure 3


def plot_methods(ax, ortho: pd.DataFrame, species: list[str], abc: str) -> pd.DataFrame:
    """Where the two ortholog callers agree, and where each is alone."""
    pairs = [(q, t) for i, q in enumerate(species) for t in species[i + 1:]]
    rows, x = [], np.arange(len(pairs))
    bottom = np.zeros(len(pairs))
    cats = [("both", NC.plum), ("OrthoFinder only", NC.orange), ("RBH only", NC.mint)]
    counts = {c: [] for c, _ in cats}
    for q, t in pairs:
        m = ortho[(ortho.query_species == q) & (ortho.target_species == t)]
        both = int((m.is_ortholog_orthofinder & m.is_rbh).sum())
        of = int((m.is_ortholog_orthofinder & ~m.is_rbh).sum())
        rb = int((~m.is_ortholog_orthofinder & m.is_rbh).sum())
        counts["both"].append(both)
        counts["OrthoFinder only"].append(of)
        counts["RBH only"].append(rb)
        rows.append({"pair": f"{SHORT[q]}-{SHORT[t]}", "both": both,
                     "orthofinder_only": of, "rbh_only": rb})
    for cat, color in cats:
        ax.bar(x, counts[cat], 0.6, bottom=bottom, color=color, label=cat)
        bottom += np.array(counts[cat])
    ax.set_xticks(x)
    ax.set_xticklabels([f"{SHORT[q]}–{SHORT[t]}" for q, t in pairs])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)
    stylia.label(ax, xlabel="", ylabel="ordered ortholog pairs",
                 title="OrthoFinder vs reciprocal best hit", abc=abc)
    return pd.DataFrame(rows)


def plot_agreement(ax, ortho: pd.DataFrame, species: list[str], abc: str) -> None:
    """The same thing as a rate, which is what you would quote."""
    pairs = [(q, t) for i, q in enumerate(species) for t in species[i + 1:]]
    vals, names = [], []
    for q, t in pairs:
        m = ortho[(ortho.query_species == q) & (ortho.target_species == t)]
        if not len(m):
            continue
        vals.append(100 * float((m.is_ortholog_orthofinder & m.is_rbh).sum()) / len(m))
        names.append(f"{SHORT[q]}–{SHORT[t]}")
    ax.bar(names, vals, color=NC.plum)
    ax.set_ylim(0, 100)
    stylia.label(ax, xlabel="", ylabel="% of pairs called by both",
                 title="Agreement rate", abc=abc)


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    dense = O.load_dense_all()   # the per-species count columns live in evidence/ now
    species = species_present(dense)
    nb = O.load_neighbors()
    ortho = O.load_orthologs()
    groups = O.load_orthogroups()

    fig, axs = stylia.create_figure(1, 2, width_ratios=[3, 2])
    reach = plot_reach(axs.next(), dense, species, abc="A")
    plot_orthogroup_sizes(axs.next(), groups, abc="B")
    out = OUT_DIR / "coverage.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")
    for _, r in reach[reach.method == "either"].iterrows():
        say(f"     {r['query']:<14} -> {r.target:<14} {r.pct:>5.1f}% have an ortholog")

    fig, axs = stylia.create_figure(1, 2)
    plot_identity(axs.next(), nb, species, abc="A")
    plot_norm_bitscore(axs.next(), nb, species, abc="B")
    out = OUT_DIR / "identity.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    fig, axs = stylia.create_figure(1, 2, width_ratios=[3, 2])
    tab = plot_methods(axs.next(), ortho, species, abc="A")
    plot_agreement(axs.next(), ortho, species, abc="B")
    out = OUT_DIR / "methods.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")
    say(tab.to_string(index=False))


if __name__ == "__main__":
    main()
