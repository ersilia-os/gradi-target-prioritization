"""Stage-02 figures: the COG functional-category profile of the three bacteria.

Kept separate from `02_function.py` so the data stage needs nothing but pandas and cogclassifier --
matplotlib and stylia are a plotting concern, not a classification one. Same split as v1's `02c`.

    python scripts/plots/function.py     ->  output/plots/function/cog_categories.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import stylia  # noqa: E402
from src import function as F  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "function"
NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL

LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(F.SPECIES), choices=list(F.SPECIES))
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = {sp: F.load_cog_counts(sp) for sp in args.species}
    totals = {sp: len(F.load_cog(sp)) for sp in args.species}

    # COG's own category order (fun-24.tab): grouped, and the grouping is meaningful.
    order = counts[args.species[0]]["cog_category"].tolist()
    descs = dict(zip(counts[args.species[0]]["cog_category"],
                     counts[args.species[0]]["cog_category_desc"]))

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(1, 2, width=1.0, height=0.62, width_ratios=[3, 1])

    # -- panel 1: the profile, as % of proteome so the three species are comparable
    ax = axs[0]
    h = 0.8 / len(args.species)
    ys = range(len(order))
    for i, sp in enumerate(args.species):
        c = counts[sp].set_index("cog_category")["n"]
        pct = [100 * c.get(k, 0) / totals[sp] for k in order]
        ax.barh([y + i * h for y in ys], pct, height=h, color=NPG[i], label=LABELS[sp])
    ax.set_yticks([y + 0.4 - h / 2 for y in ys])
    ax.set_yticklabels([f"{k}  {descs[k]}" for k in order], fontsize=SS)
    ax.invert_yaxis()
    ax.legend(fontsize=SS, frameon=False, loc="lower right")
    stylia.label(ax, xlabel="% of proteome", ylabel="", title="COG functional categories")

    # -- panel 2: how much of each proteome got an answer at all
    ax = axs[1]
    bottoms = [0.0] * len(args.species)
    layers = [
        ("informative", lambda d: (~d["cog_category"].isin(["", "R", "S"])).sum(), NPG[3]),
        ("R + S", lambda d: d["cog_category"].isin(["R", "S"]).sum(), NPG[4]),
        ("no COG", lambda d: (d["cog_category"] == "").sum(), "#d9d9d9"),
    ]
    xs = range(len(args.species))
    for name, fn, col in layers:
        vals = []
        for sp in args.species:
            d = F.load_cog(sp)
            vals.append(100 * fn(d) / len(d))
        ax.bar(xs, vals, bottom=bottoms, color=col, label=name, width=0.6)
        bottoms = [b + v for b, v in zip(bottoms, vals)]
    ax.set_xticks(list(xs))
    ax.set_xticklabels([LABELS[sp] for sp in args.species], fontsize=SS, rotation=30, ha="right")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=SS, frameon=False, loc="lower right")
    stylia.label(ax, xlabel="", ylabel="% of proteome", title="Coverage")

    out = OUT_DIR / "cog_categories.png"
    stylia.save_figure(str(out))
    print(f"wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
