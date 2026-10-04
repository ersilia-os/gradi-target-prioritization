"""The shortlist: five stated rules applied in order, and the K. pneumoniae proteins that survive.

    shortlist.png                 A  the funnel -- each rule, what it removed, what remains
                                  B  what each rule COSTS (the cascade with that one rule dropped)
                                  C  the survivors, by axis, ordered by ADEP4 probability
    shortlist_<species>.tsv       the survivors as a table -- the artifact the consortium receives

**This figure does not rank by a score, and that is the point.** Three axes in this project each
REMOVED a composite (see `plotting/filters.py` for the measurements). The cascade is a conjunction
of single stated comparisons instead: disagree with any rule and panel B shows what it cost.

Panel C orders by `adep4_prob` purely so the rows have an order. That is NOT a priority ranking --
the four axes are different units and the project ships no defensible weighting across them.

Everything is read from the stages' own tables through `src/` loaders; nothing is recomputed.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep, not two. Copying the `parents[2]` line resolves REPO_ROOT to the repo's PARENT, and
the resulting ImportError names `src`, not the path, so it reads as a broken conda env.

`stylia.label(..., xlabel=None)` does NOT mean "no label" -- it writes the placeholder string
"X-axis / Units" onto the figure. Pass "" for no label.

`stylia.create_figure(width=, height=)` takes FRACTIONS OF THE FORMAT SIZE, not inches. Passing
inches is silent: `width=13` asks for thirteen slide-widths and yields a 4-gigapixel, 16 MB PNG,
and the script still exits 0. Keep `width` at or below 1.0.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/shortlist.py
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

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
# set_style MUST precede NamedColors(): without it NamedColors returns ArticleColors and every
# NC.plum below raises AttributeError.
stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
SS = stylia.SLIDE_FONTSIZE_SMALL
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}

# Short labels for the funnel. The rules themselves live in plotting/filters.py.
STEP_LABEL = {
    "not_membrane": "not membrane",
    "selective": "no human ortholog",
    "essential": "essential (top decile)",
    "degradable": "ADEP4 above base rate",
}

# The axes shown in panel C, as (column, display name, higher-is-better).
HEATMAP_COLUMNS = [
    ("adep4_prob", "degradability", True),
    ("screens_ess_mean", "essentiality", True),
    ("n_papers_uniprot_prokaryotic", "novelty", False),
    ("p2rank_score", "pocket", True),
]


def plot_funnel(ax, steps: list[F.Step], n_start: int, abc: str) -> None:
    """The cascade as stepped bars. Horizontal because the rules are text and text reads across."""
    names = ["proteome"] + [STEP_LABEL[s.name] for s in steps]
    counts = [n_start] + [s.n_after for s in steps]
    y = np.arange(len(counts))[::-1]

    colors = [NC.gray] + [NC.plum] * len(steps)
    ax.barh(y, counts, color=colors, height=0.62)
    for yi, c in zip(y, counts):
        ax.text(c + n_start * 0.015, yi, f"{c:,}", va="center", fontsize=SS, color=NC.black)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=SS)
    ax.set_xlim(0, n_start * 1.16)
    stylia.label(ax, xlabel="proteins remaining", ylabel="", title="Four rules, in order", abc=abc)


def plot_cost(ax, loo: dict[str, int], n_final: int, abc: str) -> None:
    """What each rule costs: survivors with that ONE rule dropped, against the full cascade.

    A funnel alone overstates the early rules -- they are conjunctive, so the final count does not
    depend on their order, and this is the panel that says which rule actually does the work."""
    names = [STEP_LABEL[k] for k in loo]
    vals = [loo[k] for k in loo]
    y = np.arange(len(vals))[::-1]

    ax.barh(y, vals, color=NC.orange, height=0.62)
    ax.axvline(n_final, color=NC.black, lw=1.2, ls="--")
    ax.text(
        n_final + max(vals) * 0.02, y.min() - 0.42, f"all four: {n_final}",
        fontsize=SS, color=NC.black, ha="left", va="center",
    )
    for yi, v in zip(y, vals):
        ax.text(v + max(vals) * 0.02, yi, f"{v:,}", va="center", fontsize=SS, color=NC.black)

    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=SS)
    ax.set_xlim(0, max(vals) * 1.18)
    stylia.label(ax, xlabel="survivors without this rule", ylabel="", title="What each rule costs", abc=abc)


def plot_survivors(ax, surv: pd.DataFrame, full: pd.DataFrame, top: int, abc: str) -> None:
    """The survivors as within-proteome percentiles, one column per axis.

    Percentile against the WHOLE proteome, not against the survivors -- the question a reader has
    is "how unusual is this protein", and rescaling within 59 rows would answer a different one.
    Novelty is inverted so that bright always means "more interesting" in every column."""
    show = surv.nlargest(top, "adep4_prob")
    labels = [
        g if isinstance(g, str) and g else a
        for g, a in zip(show["gene_name"], show["uniprot_ac"])
    ]

    mat = np.full((len(show), len(HEATMAP_COLUMNS)), np.nan)
    for j, (col, _name, higher) in enumerate(HEATMAP_COLUMNS):
        ref = pd.to_numeric(full[col], errors="coerce")
        vals = pd.to_numeric(show[col], errors="coerce")
        pct = vals.apply(lambda v: np.nan if pd.isna(v) else (ref < v).mean() * 100.0)
        mat[:, j] = pct if higher else 100.0 - pct

    im = ax.imshow(mat, aspect="auto", cmap="BuPu", vmin=0, vmax=100)
    ax.set_xticks(range(len(HEATMAP_COLUMNS)))
    ax.set_xticklabels([n for _c, n, _h in HEATMAP_COLUMNS], rotation=35, ha="right", fontsize=SS)
    ax.set_yticks(range(len(show)))
    ax.set_yticklabels(labels, fontsize=SS * 0.8)
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.set_label("percentile within proteome", fontsize=SS)
    cb.ax.tick_params(labelsize=SS)
    stylia.label(ax, xlabel="", ylabel="", title=f"Top {len(show)} survivors", abc=abc)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", default="kpneumoniae", choices=sorted(LABELS))
    ap.add_argument("--top", type=int, default=30, help="rows in the panel C heatmap")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    sp = args.species
    say(f"\n  shortlist -- {LABELS[sp]}")
    df = F.load_joined(sp)
    say(f"  joined {len(df):,} proteins x {df.shape[1]} columns")

    surv, steps = F.cascade(df)
    loo = F.leave_one_out(df)

    say("\n  THE CASCADE")
    for s in steps:
        say(f"    {s.name:<14} {s.rule}")
        say(f"    {'':<14} {s.n_before:>6,} -> {s.n_after:>6,}   (removed {s.n_removed:,})")
        say(f"    {'':<14} {s.rationale}\n")

    say("  WHAT EACH RULE COSTS (survivors with that one rule dropped)")
    for k, v in loo.items():
        say(f"    without {k:<14} {v:>6,}   (all five: {len(surv):,})")

    fig, axs = stylia.create_figure(1, 3, width_ratios=[3, 3, 2.6], width=1.0, height=0.40)
    plot_funnel(axs.next(), steps, len(df), abc="A")
    plot_cost(axs.next(), loo, len(surv), abc="B")
    plot_survivors(axs.next(), surv, df, args.top, abc="C")
    out = OUT_DIR / "shortlist.png"
    stylia.save_figure(str(out))
    say(f"\n  -> {out.relative_to(REPO_ROOT)}")

    cols = [
        "uniprot_ac", "gene_name", "protein_name", "localization",
        "screens_ess_mean", "geptop_ess", "proteomelm_ess",
        "adep4_prob", "onc212_prob", "nn_similarity",
        "n_papers_uniprot_prokaryotic", "n_ligands_bacterial", "n_assayed_bacterial",
        "p2rank_score", "n_pdb_structures", "af_plddt", "bacterial_panel_orthologs",
    ]
    tsv = OUT_DIR / f"shortlist_{sp}.tsv"
    surv.sort_values("adep4_prob", ascending=False)[cols].to_csv(tsv, sep="\t", index=False)
    say(f"  -> {tsv.relative_to(REPO_ROOT)}  ({len(surv):,} proteins)")

    named = surv["gene_name"].astype(str).str.len().gt(0).sum()
    say("\n  SUMMARY")
    say(f"    survivors               {len(surv):>6,}  ({len(surv) / len(df) * 100:.2f}% of proteome)")
    say(f"    with a gene symbol      {named:>6,}  (the rest are accession-only -- Kp is a dark proteome)")
    say(f"    median papers           {surv['n_papers_uniprot_prokaryotic'].median():>6.0f}")
    say(f"    with a potent ligand    {(surv['n_ligands_bacterial'] > 0).sum():>6,}")
    say("\n  CAVEATS")
    say("    - Kp degradability rows are ranking hypotheses extrapolated from S. aureus labels,")
    say("      never measurements. Panel C is ordered, not ranked: no weighting across axes exists.")
    say("    - 'novelty' is low paper count, which is confounded with essentiality (rho 0.38-0.50).")


if __name__ == "__main__":
    main()
