"""2x2 overview of a focal organism's orthology mapping (03a output).

Lightweight reader of data/processed/other/orthology/<prefix>_orthologs_long.tsv (no network):

  TL  per-species coverage      # anchor proteins with an ortholog in each species, by tier
  TR  per-tier coverage         # anchor proteins with >=1 ortholog in each tier
  BL  conservation breadth      # bacterial species each protein is found in (24 = universal core)
  BR  human-homolog selectivity % identity to closest human ortholog + selectivity cutoff

Note: pident/coverage exist only for the human tier (DIAMOND); bacterial tiers (OrthoFinder) are
presence/absence, hence the bacterial panels are coverage/breadth-based.

Organism selected with --organism (kpneumoniae default, or ecoli). Output:
  output/plots/03b_general_orthology_<prefix>.png
Styling via stylia. Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import stylia
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[1]
ORTHO_DIR = REPO_ROOT / "data" / "processed" / "other" / "orthology"
HUMAN_SELECTIVITY_PIDENT = 30.0  # < this % identity to human -> "selective" (safety)

# (organism key) -> focal-organism metadata. close_tier is anchor-specific.
ORGANISMS = {
    "kpneumoniae": {
        "prefix": "kp",
        "name": "K. pneumoniae",
        "proteome_tsv": REPO_ROOT
        / "data/raw/kpneumoniae/proteome/UP000007841_HS11286.tsv",
        "close_tier": "klebsiella",
    },
    "ecoli": {
        "prefix": "ec",
        "name": "E. coli K-12",
        "proteome_tsv": REPO_ROOT / "data/raw/ecoli/proteome/UP000000625_EcoliK12.tsv",
        "close_tier": "escherichia",
    },
}

# NPG palette (consistent with 01d/02c). Both close-tier keys map to the same "close" colour.
TIER_COLOR = {
    "klebsiella": "#E64B35",
    "escherichia": "#E64B35",
    "gram_negative": "#4DBBD5",
    "bacteria": "#00A087",
    "human": "#3C5488",
}
TIER_LABEL = {
    "klebsiella": "Klebsiella",
    "escherichia": "Escherichia/Shigella",
    "gram_negative": "Gram-negative",
    "bacteria": "Other bacteria",
    "human": "Human",
}
SHARED_TIERS = ["gram_negative", "bacteria", "human"]


def proteome_size(tsv: Path) -> int:
    with open(tsv) as f:
        return sum(1 for _ in f) - 1  # minus header


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    spec = ORGANISMS[args.organism]
    prefix, org_name = spec["prefix"], spec["name"]
    tier_order = [spec["close_tier"]] + SHARED_TIERS
    ortho_tsv = ORTHO_DIR / f"{prefix}_orthologs_long.tsv"
    plot_path = REPO_ROOT / "output" / "plots" / f"03b_general_orthology_{prefix}.png"

    df = pd.read_csv(ortho_tsv, sep="\t", low_memory=False)
    n_anchor = proteome_size(spec["proteome_tsv"])
    print(f"{org_name}: {len(df)} ortholog rows; proteome = {n_anchor} proteins")
    sp_tier = df.drop_duplicates("species").set_index("species")["tier"].to_dict()

    stylia.set_format("slide")
    # 16:9 slide (13 x 7.3 in); top row taller (1.4:1) for the bar plots
    fig, axs = stylia.create_figure(
        2, 2, width=1.0, height=0.5625, height_ratios=[1.4, 1.0]
    )

    # ---- TL: per-species coverage (distinct anchors), colored by tier ----
    ax = axs.next()
    cov = df.groupby("species")["anchor_uniprot"].nunique().sort_values()
    colors = [TIER_COLOR.get(sp_tier.get(s, ""), "#999999") for s in cov.index]
    bars = ax.barh(range(len(cov)), cov.values, color=colors)
    ax.set_yticks(range(len(cov)))
    ax.set_yticklabels([s.replace("_", " ") for s in cov.index])
    ax.set_ylabel("")
    ax.bar_label(
        bars, labels=[f"{100 * v / n_anchor:.0f}%" for v in cov.values], padding=2
    )
    ax.set_xlabel(f"Number of {org_name} proteins with an ortholog")
    ax.set_title("Ortholog coverage per species")
    ax.set_xlim(0, cov.max() * 1.12)
    ax.legend(
        handles=[
            Patch(facecolor=TIER_COLOR[t], label=TIER_LABEL[t]) for t in tier_order
        ],
        loc="lower right",
        title="tier",
    )

    # ---- TR: per-tier coverage summary ----
    ax = axs.next()
    tier_cov = {t: df[df["tier"] == t]["anchor_uniprot"].nunique() for t in tier_order}
    ys = list(range(len(tier_order)))[::-1]
    vals = [tier_cov[t] for t in tier_order]
    bars = ax.barh(ys, vals, color=[TIER_COLOR[t] for t in tier_order])
    ax.set_yticks(ys)
    ax.set_yticklabels([TIER_LABEL[t] for t in tier_order])
    ax.set_ylabel("")
    ax.bar_label(
        bars,
        labels=[f"{v:,} ({100 * v / n_anchor:.0f}%)" for v in vals],
        label_type="center",
        color="white",
    )
    ax.set_xlabel(f"Number of {org_name} proteins")
    ax.set_title("Coverage per tier")
    ax.set_xlim(0, max(vals) * 1.12)
    print(
        "tier coverage:",
        {t: f"{tier_cov[t]} ({100 * tier_cov[t] / n_anchor:.1f}%)" for t in tier_order},
    )

    # ---- BL: conservation breadth across bacterial species ----
    ax = axs.next()
    bac = df[df["tier"] != "human"]
    n_bac_species = bac["species"].nunique()
    breadth = bac.groupby("anchor_uniprot")["species"].nunique()
    bins = np.arange(0.5, n_bac_species + 1.5, 1)
    ax.hist(breadth.values, bins=bins, color=TIER_COLOR["bacteria"])
    ax.set_xlabel(f"Number of bacterial species (of {n_bac_species})")
    ax.set_ylabel("Number of proteins")
    ax.set_title("Conservation breadth")
    ax.annotate(
        f"{int((breadth == n_bac_species).sum())} core\n({n_bac_species}/{n_bac_species})",
        xy=(n_bac_species, (breadth == n_bac_species).sum()),
        xytext=(-8, 8),
        textcoords="offset points",
        ha="right",
        color=TIER_COLOR["bacteria"],
    )

    # ---- BR: human-homolog selectivity ----
    ax = axs.next()
    h = df[df["tier"] == "human"].dropna(subset=["pident"])
    ax.hist(h["pident"], bins=20, color=TIER_COLOR["human"])
    ax.axvline(HUMAN_SELECTIVITY_PIDENT, color="#E64B35", linestyle="--", linewidth=1.5)
    n_sel = int((h["pident"] < HUMAN_SELECTIVITY_PIDENT).sum())
    ax.set_xlabel("% identity to closest human ortholog")
    ax.set_ylabel("Number of proteins")
    ax.set_title(f"Human-homolog selectivity  ({len(h)} with a human homolog)")
    ax.annotate(
        f"{n_sel} selective\n(<{HUMAN_SELECTIVITY_PIDENT:.0f}%)",
        xy=(HUMAN_SELECTIVITY_PIDENT, 0),
        xytext=(6, 6),
        textcoords="offset points",
        color="#E64B35",
        va="bottom",
    )

    plot_path.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(plot_path))
    print(f"\nWrote {plot_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
