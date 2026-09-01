"""2x2 overview of the popularity / studiedness annotation (05a) for both pathogens.

Lightweight reader of data/processed/<org>/bibliometric/<prefix>_popularity.csv (no network):

  TL  tier composition       dark / studied / well_studied, normalised %, K. pneumoniae vs E. coli
  TR  popularity-score dist   score histogram (density) per organism
  BL  top homolog sources     organisms supplying the best-studied ortholog, K. pneumoniae
  BR  top homolog sources     organisms supplying the best-studied ortholog, E. coli

Output: output/plots/05b_popularity.png. Styling via stylia. Run with the `gradi` env.
"""

from __future__ import annotations

import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)  # stylia rmtree's this on import; ensure it exists
import pandas as pd
import stylia
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_PATH = REPO_ROOT / "output" / "plots" / "05b_popularity.png"
TOP_N = 8

ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}  # NPG per organism (matches 01d/02c/03b)
TIERS = ["dark", "studied", "well_studied"]
TIER_COLOR = {"dark": "#B0B0B0", "studied": "#F2C14E", "well_studied": "#00A087"}  # none/some/lots
ORGANISMS = [
    ("kpneumoniae", "kp", "K. pneumoniae"),
    ("ecoli", "ec", "E. coli K-12"),
]


def load(prefix: str, org: str) -> pd.DataFrame:
    return pd.read_csv(REPO_ROOT / "data" / "processed" / org / "bibliometric" / f"{prefix}_popularity.csv")


def clean_org(s: str) -> str:
    return str(s).split(" (")[0]  # drop "(strain ...)" suffix


def main() -> None:
    data = {prefix: load(prefix, org) for org, prefix, _ in ORGANISMS}
    names = {prefix: name for _, prefix, name in ORGANISMS}

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 2, width=1.0, height=0.5625)  # 16:9 slide

    # ---- TL: tier composition (normalised %) ----
    ax = axs.next()
    prefixes = [p for _, p, _ in ORGANISMS]
    x = range(len(prefixes))
    bottoms = [0.0] * len(prefixes)
    for tier in TIERS:
        fracs = [100 * (data[p]["popularity_tier"] == tier).mean() for p in prefixes]
        bars = ax.bar(x, fracs, bottom=bottoms, color=TIER_COLOR[tier], label=tier.replace("_", " "))
        ax.bar_label(bars, labels=[f"{f:.0f}%" for f in fracs], label_type="center", color="white")
        bottoms = [b + f for b, f in zip(bottoms, fracs)]
    ax.set_xticks(list(x)); ax.set_xticklabels([names[p] for p in prefixes])
    ax.set_xlabel(""); ax.set_ylabel("Share of proteome (%)"); ax.set_ylim(0, 100)
    ax.set_title("Studiedness tier composition")
    ax.legend(loc="lower right", title="tier")

    # ---- TR: popularity-score distribution (density) ----
    ax = axs.next()
    for p in prefixes:
        ax.hist(data[p]["popularity_score"], bins=20, range=(0, 1), density=True,
                histtype="stepfilled", alpha=0.5, color=ORG_COLOR[p], label=names[p])
    ax.set_xlabel("Popularity score"); ax.set_ylabel("Density")
    ax.set_title("Popularity-score distribution")
    ax.legend(loc="upper right")

    # ---- BL / BR: top homolog-source organisms per pathogen ----
    for org, prefix, name in ORGANISMS:
        ax = axs.next()
        df = data[prefix]
        src = (df.loc[df["best_homolog_uniprot"].notna() & (df["best_homolog_organism"].astype(str) != ""),
                      "best_homolog_organism"].map(clean_org).value_counts().head(TOP_N))
        ys = range(len(src))[::-1]
        bars = ax.barh(ys, src.values, color=ORG_COLOR[prefix])
        ax.set_yticks(list(ys)); ax.set_yticklabels(src.index)
        ax.set_ylabel(""); ax.set_xlabel("Number of proteins")
        ax.bar_label(bars, label_type="center", color="white")
        ax.set_xlim(0, src.values.max() * 1.05)
        ax.set_title(f"Top homolog sources — {name}")

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(PLOT_PATH))
    print(f"Wrote {PLOT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
