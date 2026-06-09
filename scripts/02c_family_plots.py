"""2x3 family-overview barplots for K. pneumoniae and E. coli.

Lightweight reader: counts top families/domains straight from the per-protein family CSVs
produced by 02a (InterPro) and 02b (PANTHER) -- no network.

  rows = K. pneumoniae HS11286 / E. coli K-12 MG1655
  cols = top InterPro families | top InterPro domains | top PANTHER families

Each bar = number of proteins carrying that family/domain. Bars are coloured per organism with
the NPG palette (matching 01d). Output: output/plots/02c_family_plots.png
Styling via stylia (ersilia-os/stylia). Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_PATH = REPO_ROOT / "output" / "plots" / "02c_family_plots.png"
TOP_N = 10
NAME_MAX = 26  # truncate long family/domain names

NPG = {"kp": "#E64B35", "ec": "#4DBBD5"}  # one NPG colour per organism (matches 01d)
ORGANISMS = [
    ("kpneumoniae", "kp", "K. pneumoniae HS11286"),
    ("ecoli", "ec", "E. coli K-12 MG1655"),
]


def top_entries(
    df: pd.DataFrame, id_col: str, name_col: str, n: int
) -> list[tuple[str, int]]:
    """Count ids across proteins (';'-joined cells) -> top-n [(label 'ID  name', count)]."""
    counts: Counter[str] = Counter()
    names: dict[str, str] = {}
    for ids_str, names_str in zip(df[id_col].fillna(""), df[name_col].fillna("")):
        ids = [x for x in str(ids_str).split(";") if x]
        nms = [x for x in str(names_str).split(";") if x]
        for j, i in enumerate(ids):
            counts[i] += 1
            if j < len(nms):
                names[i] = nms[j]
    out = []
    for i, c in counts.most_common(n):
        nm = names.get(i, "")
        nm = nm if len(nm) <= NAME_MAX else nm[: NAME_MAX - 1] + "…"
        out.append((f"{i}  {nm}", c))
    return out


def main() -> None:
    stylia.set_format("slide")  # bigger fonts for slides
    fig, axs = stylia.create_figure(2, 3, height=0.5)  # default width; slide format

    for i, (org, prefix, name) in enumerate(ORGANISMS):
        fam_dir = REPO_ROOT / "data" / "processed" / org / "families"
        ip = pd.read_csv(fam_dir / f"{prefix}_interpro_annotation.csv")
        pa = pd.read_csv(fam_dir / f"{prefix}_panther_annotation.csv")
        color = NPG[prefix]

        metrics = [
            (
                "InterPro family",
                top_entries(ip, "interpro_family_ids", "interpro_family_names", TOP_N),
            ),
            (
                "InterPro domain",
                top_entries(ip, "interpro_domain_ids", "interpro_domain_names", TOP_N),
            ),
            (
                "PANTHER family",
                top_entries(pa, "panther_family_ids", "panther_family_names", TOP_N),
            ),
        ]
        for j, (col_title, top) in enumerate(metrics):
            ax = axs.next()
            labels = [t[0] for t in top][::-1]  # largest on top
            vals = [t[1] for t in top][::-1]
            bars = ax.barh(range(len(vals)), vals, color=color)
            ax.set_yticks(range(len(vals)))
            ax.set_yticklabels(labels)
            ax.bar_label(
                bars, label_type="center", color="white"
            )  # counts inside the bars
            ax.set_xlabel("Number of proteins")
            ax.set_ylabel(name if j == 0 else "")
            if i == 0:
                ax.set_title(col_title)
            print(f"{name} | {col_title}: top = {top[0][0]} ({top[0][1]})")

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(PLOT_PATH))
    print(f"\nWrote {PLOT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
