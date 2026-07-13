"""2x3 overview of the BindingDB bioactivity annotation (06b), direct vs orthology-transferred.

Same house style as 06i (the ChEMBL slide) and the task-agnostic deck: stylia "slide" format, NPG
palette, white/dark in-bar labels. Every panel is specific to the single `--organism` of the slide
— there is no content shared between the kp and ec slides. Direct = a bacterial BindingDB target at
>=95% identity; transferred = a lower-identity bacterial homolog (>=40%).

  1  potent targets by route             # proteins with <=1 uM / <=100 nM affinity, direct vs transferred
  2  evidence transfer distance          best BindingDB-target %identity distribution (95% = direct)
  3  top targets                         # potent compounds, coloured by route
  4  evidence source organisms           which organisms supply the evidence
  5  test depth vs hits                  # compounds tested vs # potent (log-log)
  6  data coverage by homolog            proteins reachable via direct / bacterial / human homolog

Reads output/results/<org>/<prefix>_bindingdb.csv (06b). Output: output/plots/06j_bindingdb_<prefix>.png
(one slide per --organism). Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)  # stylia rmtree's this on import; ensure it exists
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

REPO_ROOT = L.REPO_ROOT
TOP_N = 8

NPG = stylia.CategoricalPalette("npg").colors                 # vivid Nature palette (deck-wide)
ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}                 # NPG per organism (matches 01d/02c/05b)
ROUTE_COLOR = {"direct": "#00A087", "ortholog": "#F39B7F"}     # NPG green (direct) / salmon (transferred)
ROUTE_LABEL = {"direct": "Direct (≥95% id)", "ortholog": "Orthology-transferred"}
ORGANISMS = [("kpneumoniae", "kp", "K. pneumoniae"), ("ecoli", "ec", "E. coli K-12")]


def load(prefix: str, org: str) -> pd.DataFrame:
    d = pd.read_csv(L.results_dir(org) / f"{prefix}_bindingdb.csv")
    d["has_direct"] = d["bindingdb_direct_has"].fillna(False).astype(bool)
    d["bact_n_potent"] = d["bindingdb_bact_n_potent"].fillna(0)
    d["bact_n_comp"] = d["bindingdb_bact_n_compounds"].fillna(0)
    d["dir_n_potent"] = d["bindingdb_direct_n_potent"].fillna(0)
    d["human_n_comp"] = d["bindingdb_human_n_compounds"].fillna(0)
    d["has_bact"] = d["bact_n_comp"] > 0
    d["route"] = np.where(d["has_direct"], "direct", "ortholog")
    return d


def short_org(name: str) -> str:
    s = str(name).split(" (")[0].strip()        # drop "(strain ...)"
    parts = s.split()
    if len(parts) >= 2 and parts[0][:1].isupper():
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    data = {p: load(p, o) for o, p, _ in ORGANISMS}
    names = {p: n for _, p, n in ORGANISMS}
    org = args.organism
    _, prefix = L.ORGANISMS[org]
    d = data[prefix]
    genes = L.load_genes(org)
    orgname = names[prefix]

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)  # 16:9 slide, 2x3 panels

    # ---- panel 1: potent targets by route, at two potency cut-offs (this organism) ----
    ax = axs.next()
    dir_aff = d["bindingdb_direct_best_paff"].fillna(0)
    bact_aff = d["bindingdb_bact_best_paff"].fillna(0)
    cats = ["≤1 µM\n(pAff≥6)", "≤100 nM\n(pAff≥7)"]
    xc = range(len(cats))
    direct = [int((d["dir_n_potent"] > 0).sum()),
              int((d["has_direct"] & (dir_aff >= 7)).sum())]
    ortho = [int(((d["bact_n_potent"] > 0) & ~d["has_direct"]).sum()),
             int((~d["has_direct"] & (bact_aff >= 7)).sum())]
    b1 = ax.bar(xc, direct, color=ROUTE_COLOR["direct"], label=ROUTE_LABEL["direct"])
    b2 = ax.bar(xc, ortho, bottom=direct, color=ROUTE_COLOR["ortholog"], label=ROUTE_LABEL["ortholog"])
    ax.bar_label(b1, label_type="center", color="white")          # white reads on green
    ax.bar_label(b2, label_type="center", color="#5a2d1a")        # dark reads on salmon
    ax.bar_label(b2, labels=[f"{d_+o_}" for d_, o_ in zip(direct, ortho)], padding=2)
    ax.set_xticks(list(xc)); ax.set_xticklabels(cats)
    ax.set_xlabel(""); ax.set_ylabel("Number of proteins")
    ax.margins(y=0.15)
    ax.set_title(f"Potent targets by route — {orgname}")
    ax.legend(loc="upper right", title="route")

    # ---- panel 2: best BindingDB-target sequence identity (this organism) ----
    ax = axs.next()
    pid = d.loc[d["has_bact"], "bindingdb_bact_best_pident"].dropna()
    ax.hist(pid, bins=range(25, 101, 5), histtype="stepfilled", alpha=0.75, color=ORG_COLOR[prefix])
    ax.axvline(95, color="#555555", linestyle="--")
    ax.text(94, ax.get_ylim()[1] * 0.92, "direct →", ha="right", color="#555555")
    ax.set_xlabel("Best BindingDB-target % identity"); ax.set_ylabel("Number of proteins")
    ax.set_title(f"Evidence transfer distance — {orgname}")

    # ---- panel 3: top targets by # potent compounds, coloured by route (this organism) ----
    ax = axs.next()
    top = d[d["bact_n_potent"] > 0].sort_values("bact_n_potent", ascending=False).head(TOP_N).iloc[::-1]
    ys = range(len(top))
    colors = [ROUTE_COLOR[r] for r in top["route"]]
    vals = top["bact_n_potent"].to_numpy()
    bars = ax.barh(list(ys), vals, color=colors)
    labels = [genes.get(a) or a for a in top["uniprot_accession"]]
    ax.set_yticks(list(ys)); ax.set_yticklabels(labels)
    # log x — # potent compounds spans orders of magnitude
    if len(vals) and vals.max() > 0:
        ax.set_xscale("log"); ax.set_xlim(1, vals.max() * 2.5)
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_ylabel(""); ax.set_xlabel("# potent compounds (≤1 µM, log)")
    ax.set_title(f"Top targets — {orgname}")
    ax.legend(handles=[Patch(color=ROUTE_COLOR["direct"], label=ROUTE_LABEL["direct"]),
                       Patch(color=ROUTE_COLOR["ortholog"], label=ROUTE_LABEL["ortholog"])],
              loc="lower right")

    # ---- panel 4: evidence source organisms, one vivid colour per organism (this organism) ----
    ax = axs.next()
    sub = d[d["has_bact"]].copy()
    sub["src"] = sub["bindingdb_bact_best_organism"].map(short_org)
    counts = sub["src"].value_counts().sort_values(ascending=True).tail(TOP_N)
    ys = range(len(counts))
    colors = [NPG[i % len(NPG)] for i in range(len(counts))][::-1]  # distinct hue per organism
    bars = ax.barh(list(ys), counts.to_numpy(), color=colors)
    ax.set_yticks(list(ys)); ax.set_yticklabels(counts.index.tolist())
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_xlim(0, counts.max() * 1.12 if len(counts) else 1)
    ax.set_ylabel(""); ax.set_xlabel("Number of proteins")
    ax.set_title(f"Evidence source organisms — {orgname}")

    # ---- panel 5: test depth vs hits — # tested vs # potent (log-log), one organism ----
    ax = axs.next()
    YFLOOR = 0.7  # log scale can't show 0 potent; park zero-hit targets on the axis floor
    sub = d[d["bact_n_comp"] > 0]
    for route in ("ortholog", "direct"):  # draw direct last so it sits on top
        s = sub[sub["route"] == route]
        xs = s["bact_n_comp"].to_numpy()
        ys = np.where(s["bact_n_potent"].to_numpy() > 0, s["bact_n_potent"].to_numpy(), YFLOOR)
        ax.scatter(xs, ys, s=22, alpha=0.7, color=ROUTE_COLOR[route],
                   edgecolors="white", linewidths=0.3)
    xmax = sub["bact_n_comp"].max() if len(sub) else 10
    xline = np.array([1, xmax * 1.1])
    for frac, lab in [(1.0, "100%"), (0.1, "10%"), (0.01, "1%")]:  # hit-rate guide diagonals
        ax.plot(xline, np.clip(xline * frac, YFLOOR, None), color="#999999",
                linestyle="--", linewidth=0.7, zorder=0)
        ax.text(xmax * 1.1, max(xmax * 1.1 * frac, YFLOOR), f" {lab}", va="center",
                ha="left", color="#999999", fontsize="small")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(0.8, xmax * 1.6); ax.set_ylim(YFLOOR * 0.85, None)
    ax.set_xlabel("# compounds tested (log)")
    ax.set_ylabel("# potent (≤1 µM, log; 0 at floor)")
    ax.set_title(f"Test depth vs hits — {orgname}")
    ax.legend(handles=[Patch(color=ROUTE_COLOR["direct"], label=ROUTE_LABEL["direct"]),
                       Patch(color=ROUTE_COLOR["ortholog"], label=ROUTE_LABEL["ortholog"])],
              loc="lower right")

    # ---- panel 6: how many proteins have BindingDB data, and via which homolog (this organism) ----
    ax = axs.next()
    tiers = ["Direct\n(≥95% id)", "Bacterial\nhomolog", "Human\nhomolog"]
    xt = np.arange(len(tiers))
    counts = [int(d["has_direct"].sum()),
              int((d["bact_n_comp"] > 0).sum()),
              int((d["human_n_comp"] > 0).sum())]
    tier_colors = [ROUTE_COLOR["direct"], ROUTE_COLOR["ortholog"], NPG[3]]
    bars = ax.bar(xt, counts, color=tier_colors)
    ax.bar_label(bars, padding=2)
    ax.set_xticks(list(xt)); ax.set_xticklabels(tiers)
    ax.set_ylabel("Proteins with BindingDB data"); ax.set_xlabel("")
    ax.set_title(f"Data coverage by homolog — {orgname}")
    ax.margins(y=0.18)  # headroom for the count labels
    n_none = int((~((d["bact_n_comp"] > 0) | (d["human_n_comp"] > 0))).sum())
    ax.text(0.5, -0.24,
            f"Counts are nested (every Direct protein is also a Bacterial homolog).\n"
            f"The other {n_none:,} of {len(d):,} proteins have no BindingDB data.",
            transform=ax.transAxes, ha="center", va="top", color="#777777", fontsize="small")

    out = REPO_ROOT / "output" / "plots" / f"06j_bindingdb_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
