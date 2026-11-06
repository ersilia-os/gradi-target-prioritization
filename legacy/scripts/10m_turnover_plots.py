"""2x3 slide: measured protein turnover and which protease does the degrading (docs §3.3).

Built on `10e`, which ingests the three datasets that measure degradation rather than predicting it:
Nagar 2021 (pulsed-SILAC half-lives), Gupta 2024 (half-lives across 13 conditions plus a
protease-knockout panel) and Niwa 2022 (GroE depletion run in WT, Δlon, ΔclpPX, ΔhslVU).

  1  Gupta vs Nagar             do two independent half-life measurements agree?
  2  half-life by protease      Gupta's KO panel, recomputed — who degrades what
  3  ClpP substrates            what the ClpP-attributed proteins look like vs the rest
  4  ClpP vs Lon stabilisation  the KO panel as a 2D map of protease specificity
  5  Niwa head-to-head          Lon vs ClpXP vs HslUV rescue of GroE-depletion loss
  6  the growth-correction cliff  why 1,149 measured half-lives leave ~50 usable ones

Panel 6 is the one to read before trusting any of the others: most E. coli proteins are diluted out
by growth faster than they are proteolysed, so their measured half-life carries almost no
proteolytic signal once corrected.

Reads output/results/<org>/<prefix>_deg_measured.csv (10e).
Output: output/plots/10m_turnover_<prefix>.png. Run with `gradi`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
from scipy import stats  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
GREY = "#C9C9C7"

PROTEASE_COLOR = {"clpP": "#E64B35", "lon": "#3C5488", "hslV": "#00A087",
                  "additive": "#F39B7F", "redundant": "#8491B4", "unattributed": GREY}
PROTEASE_ORDER = ["clpP", "lon", "hslV", "additive", "redundant", "unattributed"]
NIWA_COLOR = {"clpxp": "#E64B35", "lon": "#3C5488", "hslvu": "#00A087",
              "multiple": "#F39B7F", "none": GREY}
NIWA_ORDER = ["clpxp", "lon", "hslvu", "multiple", "none"]
CLASS_COLOR = {"fast": "#E64B35", "intermediate": "#F39B7F", "stable": "#8491B4"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]

    raw = pd.read_csv(LOC.results_dir(org) / f"{prefix}_deg_measured.csv")
    # Count before the inf sweep: Nagar records `inf` for proteins that never measurably decayed,
    # and losing that distinction would flatten the first step of the panel-6 cliff.
    n_halflife_rows = int(raw["nagar_halflife_min"].notna().sum())
    d = raw.replace([np.inf, -np.inf], np.nan)

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: Gupta vs Nagar half-life ----
    ax = axs.next()
    b = d.dropna(subset=["gupta_halflife_hrs", "nagar_halflife_min"])
    if len(b) > 3:
        x = b["gupta_halflife_hrs"] * 60.0        # hours -> minutes, to compare like with like
        y = b["nagar_halflife_min"]
        ax.scatter(x, y, s=10, alpha=0.5, color=NPG[0], linewidths=0, rasterized=True)
        ax.set_xscale("log"); ax.set_yscale("log")
        rho, _ = stats.spearmanr(x, y)
        lo = min(x.min(), y.min()); hi = max(x.max(), y.max())
        ax.plot([lo, hi], [lo, hi], ls="--", color="#999999", lw=0.9)
        ax.text(0.03, 0.97, f"n = {len(b)}\nρ = {rho:.2f}", transform=ax.transAxes, fontsize=SS,
                va="top", ha="left", color="#2B2333", linespacing=1.4)
    else:
        ax.text(0.5, 0.5, "too few shared proteins", ha="center", va="center",
                transform=ax.transAxes, fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="Gupta 2024 half-life (min)", ylabel="Nagar 2021 half-life (min)",
                 title=f"Two independent measurements — {orgname}")

    # ---- panel 2: half-life by protease attribution ----
    ax = axs.next()
    present = [p for p in PROTEASE_ORDER
               if (d["gupta_protease_attribution"] == p).sum() >= 5]
    data = [d.loc[(d["gupta_protease_attribution"] == p) & d["gupta_halflife_hrs"].notna(),
                  "gupta_halflife_hrs"].to_numpy() for p in present]
    data = [v for v in data if len(v)]
    labs = [p for p, v in zip(present, data) if len(v)]
    if data:
        bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False,
                        medianprops=dict(color="#2B2333", lw=1.4))
        for patch, p in zip(bp["boxes"], labs):
            patch.set_facecolor(PROTEASE_COLOR[p]); patch.set_edgecolor("#8A8A96"); patch.set_alpha(0.9)
        ax.set_xticklabels([f"{p}\n({len(v):,})" for p, v in zip(labs, data)], fontsize=SS - 1)
        med_all = d["gupta_halflife_hrs"].median()
        ax.axhline(med_all, ls="--", color="#555555", lw=1)
        ax.text(0.98, 0.02, f"proteome median {med_all:.1f} h", transform=ax.transAxes,
                fontsize=SS, va="bottom", ha="right", color="#555555")
    stylia.label(ax, xlabel="", ylabel="half-life (h), minimal media",
                 title=f"Who degrades what — {orgname}")

    # ---- panel 3: what ClpP substrates look like ----
    ax = axs.next()
    is_clpp = d["gupta_protease_attribution"] == "clpP"
    other = d["gupta_protease_attribution"].notna() & ~is_clpp
    groups = [("ClpP-attributed", d.loc[is_clpp, "nagar_halflife_min"].dropna(), PROTEASE_COLOR["clpP"]),
              ("other measured", d.loc[other, "nagar_halflife_min"].dropna(), GREY)]
    groups = [(lab, v, c) for lab, v, c in groups if len(v) >= 3]
    if groups:
        bins = np.logspace(np.log10(max(min(v.min() for _, v, _ in groups), 1e-2)),
                           np.log10(max(v.max() for _, v, _ in groups)), 26)
        for lab, v, c in groups:
            ax.hist(v, bins=bins, density=True, alpha=0.6, color=c, label=f"{lab} ({len(v):,})")
        ax.set_xscale("log")
        ax.legend(fontsize=SS, frameon=False, loc="upper left")
        if len(groups) == 2:
            u = stats.mannwhitneyu(groups[0][1], groups[1][1], alternative="two-sided")
            ax.text(0.97, 0.97,
                    f"median {groups[0][1].median():.0f} vs {groups[1][1].median():.0f} min\n"
                    f"Mann-Whitney p = {u.pvalue:.1e}",
                    transform=ax.transAxes, fontsize=SS, va="top", ha="right", color="#555555",
                    linespacing=1.4)
    else:
        ax.text(0.5, 0.5, "not enough overlap", ha="center", va="center", transform=ax.transAxes,
                fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="Nagar half-life (min)", ylabel="density",
                 title=f"Are ClpP substrates short-lived? — {orgname}")

    # ---- panel 4: ClpP vs Lon stabilisation ----
    ax = axs.next()
    k = d.dropna(subset=["gupta_log2_clpP", "gupta_log2_lon"])
    if len(k) > 3:
        att = d.loc[k.index, "gupta_protease_attribution"].fillna("unattributed")
        for p in PROTEASE_ORDER:
            m = att == p
            if not m.any():
                continue
            ax.scatter(k.loc[m, "gupta_log2_clpP"], k.loc[m, "gupta_log2_lon"],
                       s=12 if p in ("clpP", "lon", "hslV") else 7,
                       alpha=0.85 if p in ("clpP", "lon", "hslV") else 0.35,
                       color=PROTEASE_COLOR[p], linewidths=0, rasterized=True, label=p)
        ax.axhline(0.5, ls="--", color="#999999", lw=0.9)
        ax.axvline(0.5, ls="--", color="#999999", lw=0.9)
        ax.legend(fontsize=SS - 1, frameon=False, ncol=2, loc="upper left")
        ax.text(0.97, 0.03, f"n = {len(k):,}", transform=ax.transAxes, fontsize=SS,
                va="bottom", ha="right", color="#555555")
    stylia.label(ax, xlabel="stabilisation in ΔclpP (log2)", ylabel="stabilisation in Δlon (log2)",
                 title=f"Protease specificity — {orgname}")

    # ---- panel 5: Niwa head-to-head ----
    ax = axs.next()
    present = [p for p in NIWA_ORDER if (d["niwa_attribution"] == p).sum() > 0]
    counts = [int((d["niwa_attribution"] == p).sum()) for p in present]
    if counts:
        bars = ax.bar(range(len(present)), counts,
                      color=[NIWA_COLOR[p] for p in present], width=0.62)
        ax.bar_label(bars, padding=2, fontsize=SS)
        ax.set_xticks(range(len(present)))
        ax.set_xticklabels([p for p in present], fontsize=SS, rotation=12)
        ax.margins(y=0.2)
        ax.text(0.97, 0.97,
                "rescue of GroE-depletion loss\nin each protease knockout",
                transform=ax.transAxes, fontsize=SS, va="top", ha="right", color="#555555",
                linespacing=1.4)
    stylia.label(ax, xlabel="", ylabel="proteins",
                 title=f"Lon vs ClpXP vs HslUV (Niwa 2022) — {orgname}")

    # ---- panel 6: the growth-correction cliff ----
    ax = axs.next()
    stages = [
        ("Nagar half-life measured", n_halflife_rows, NPG[0]),
        ("finite (not diluted-out)", int(d["nagar_halflife_min"].notna().sum()), NPG[1]),
        ("positive after growth\ncorrection", int(d["nagar_halflife_proteolytic_min"].notna().sum()),
         "#E64B35"),
    ]
    bars = ax.barh(range(len(stages)), [n for _, n, _ in stages],
                   color=[c for _, _, c in stages], height=0.6)
    ax.bar_label(bars, padding=3, fontsize=SS, fmt="%d")
    ax.set_yticks(range(len(stages)))
    ax.set_yticklabels([lab for lab, _, _ in stages], fontsize=SS - 1)
    ax.invert_yaxis()
    ax.margins(x=0.18)
    ax.text(0.97, 0.05,
            "most proteins are diluted by growth\nfaster than they are proteolysed",
            transform=ax.transAxes, fontsize=SS, va="bottom", ha="right", color="#555555",
            linespacing=1.4)
    stylia.label(ax, xlabel="proteins", ylabel="",
                 title=f"How much turnover is proteolysis? — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"10m_turnover_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")


if __name__ == "__main__":
    main()
