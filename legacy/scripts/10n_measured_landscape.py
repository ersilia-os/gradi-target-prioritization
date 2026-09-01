"""2x3 slide: how much measured degradability evidence exists, and where (docs §3).

The dataset-level view. Every other slide zooms into one experiment; this one asks what the measured
layer covers in total, how much the datasets overlap, where in protein space the evidence sits, and
whether the computed features predict the measured turnover at all.

  1  coverage per dataset      proteins reached by each measured source
  2  dataset overlap           how much the sources share (2D)
  3  any-evidence funnel       proteome -> any measured evidence -> Clp-specific evidence
  4  evidence on the map       measured coverage across the ESM-C protein universe (2D)
  5  do features predict it?   the 10d disorder AUROCs with bootstrap CIs against the 0.5 null
  6  native vs transferred     for K. pneumoniae, all of it is transferred

Panel 5 is the accountability panel: `disorder_auroc` returns a bootstrap interval precisely so the
weakness of the signal cannot be hidden, and every interval here crosses or nearly crosses 0.5.

Reads <prefix>_deg_measured.csv (10e), <prefix>_clpp_activator.csv (10c), the ESM-C projection and
data/processed/<org>/degradability/<prefix>_disorder_validation.csv (10d).
Output: output/plots/10n_landscape_<prefix>.png. Run with `gradi`.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
GREY = "#C9C9C7"
BG = "#D8D8D6"

# Ec->Kp orthogroup transfer ceiling, from the 03a table.
KP_TRANSFER_CEILING = 3179


def load(org: str) -> pd.DataFrame:
    _, prefix = LOC.ORGANISMS[org]
    r = LOC.results_dir(org)
    d = pd.DataFrame({"uniprot_accession": LOC.load_accessions(org)})
    meas = r / f"{prefix}_deg_measured.csv"
    if meas.exists():
        d = d.merge(pd.read_csv(meas), on="uniprot_accession", how="left")
    act = r / f"{prefix}_clpp_activator.csv"
    if act.exists():
        a = pd.read_csv(act)[["uniprot_accession", "activator_evidence"]]
        d = d.merge(a, on="uniprot_accession", how="left")
    proj = r / f"{prefix}_esmc600m_projection.csv"
    if proj.exists():
        d = d.merge(pd.read_csv(proj)[["uniprot_accession", "tsne_x", "tsne_y"]],
                    on="uniprot_accession", how="left")
    for c in ("tsne_x", "tsne_y", "activator_evidence"):
        if c not in d:
            d[c] = np.nan
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]

    d = load(org)
    n_total = len(d)

    # The measured sources, each reduced to a boolean "this protein was measured here".
    sources = {
        "Nagar 2021\nhalf-life": d.get("nagar_halflife_class").notna() if "nagar_halflife_class" in d else pd.Series(False, index=d.index),
        "Gupta 2024\nhalf-life": d.get("gupta_halflife_hrs").notna() if "gupta_halflife_hrs" in d else pd.Series(False, index=d.index),
        "Gupta 2024\nprotease KO": d.get("gupta_protease_attribution").notna() if "gupta_protease_attribution" in d else pd.Series(False, index=d.index),
        "Niwa 2022\nprotease panel": d.get("niwa_attribution").notna() if "niwa_attribution" in d else pd.Series(False, index=d.index),
        "ADEP4/ONC212\nactivator": d["activator_evidence"].notna(),
    }
    any_evidence = pd.concat(sources.values(), axis=1).any(axis=1)

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: coverage per dataset ----
    ax = axs.next()
    names = list(sources)
    counts = [int(v.sum()) for v in sources.values()]
    colors = [NPG[i % len(NPG)] for i in range(len(names))]
    bars = ax.bar(range(len(names)), counts, color=colors, width=0.62)
    ax.bar_label(bars, labels=[f"{c:,}\n{c / n_total:.0%}" for c in counts], padding=2, fontsize=SS - 1)
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, fontsize=SS - 1)
    ax.axhline(n_total, ls="--", color="#555555", lw=1)
    ax.text(0.02, 0.97, f"proteome {n_total:,}", transform=ax.transAxes, fontsize=SS,
            va="top", ha="left", color="#555555")
    if org == "kpneumoniae":
        ax.axhline(KP_TRANSFER_CEILING, ls=":", color="#C98A1E", lw=1.2)
        ax.text(0.98, KP_TRANSFER_CEILING / ax.get_ylim()[1] + 0.02,
                f"ortholog ceiling {KP_TRANSFER_CEILING:,}", transform=ax.transAxes,
                fontsize=SS, va="bottom", ha="right", color="#C98A1E")
    ax.margins(y=0.26)
    stylia.label(ax, xlabel="", ylabel="proteins measured",
                 title=f"Measured coverage by source — {orgname}")

    # ---- panel 2: dataset overlap ----
    ax = axs.next()
    keys = [k.replace("\n", " ") for k in names]
    mat = np.zeros((len(names), len(names)))
    vals = list(sources.values())
    for i in range(len(names)):
        for j in range(len(names)):
            mat[i, j] = int((vals[i] & vals[j]).sum())
    ax.imshow(mat, cmap="Greens", aspect="auto")
    vmax = mat.max() if mat.max() else 1
    for i in range(len(names)):
        for j in range(len(names)):
            if mat[i, j]:
                ax.text(j, i, f"{int(mat[i, j]):,}", ha="center", va="center", fontsize=SS - 2,
                        color="white" if mat[i, j] > vmax * 0.55 else "#222222")
    short = ["Nagar", "Gupta HL", "Gupta KO", "Niwa", "Activator"]
    ax.set_xticks(range(len(names))); ax.set_xticklabels(short, fontsize=SS - 1, rotation=35)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(short, fontsize=SS - 1)
    stylia.label(ax, xlabel="", ylabel="", title=f"Where the sources overlap — {orgname}")

    # ---- panel 3: any-evidence funnel ----
    ax = axs.next()
    clp_specific = (
        (d.get("gupta_protease_attribution", pd.Series(index=d.index, dtype=object)) == "clpP")
        | (d.get("niwa_attribution", pd.Series(index=d.index, dtype=object)) == "clpxp")
        | (d["activator_evidence"].fillna(0) > 0)
    )
    stages = [("proteome", n_total, GREY),
              ("any measured evidence", int(any_evidence.sum()), NPG[0]),
              ("Clp-specific evidence", int(clp_specific.sum()), "#E64B35")]
    widths = np.sqrt(np.array([n for _, n, _ in stages], dtype=float) / stages[0][1])
    for i, (lab, n, c) in enumerate(stages):
        w = widths[i]
        ax.barh(i, w, height=0.7, left=(1 - w) / 2, color=c)
        ax.text(0.5, i, f"{lab}: {n:,}  ({n / n_total:.0%})", ha="center", va="center",
                fontsize=SS, color="#2B2333")
    ax.set_ylim(-0.6, len(stages) - 0.4); ax.set_xlim(0, 1)
    ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
    for s in ("top", "right", "bottom", "left"):
        ax.spines[s].set_visible(False)
    stylia.label(ax, xlabel="", ylabel="", title=f"How much is measured at all — {orgname}")

    # ---- panel 4: evidence on the ESM-C map ----
    ax = axs.next()
    has_xy = d["tsne_x"].notna() & d["tsne_y"].notna()
    xy = d[has_xy]
    ax.scatter(xy["tsne_x"], xy["tsne_y"], s=4, alpha=0.3, color=BG, linewidths=0, rasterized=True,
               label=f"no measured evidence ({int((~any_evidence[has_xy]).sum()):,})")
    m = has_xy & any_evidence
    ax.scatter(d.loc[m, "tsne_x"], d.loc[m, "tsne_y"], s=7, alpha=0.7, color=NPG[0], linewidths=0,
               rasterized=True, label=f"measured ({int(m.sum()):,})")
    c = has_xy & clp_specific
    ax.scatter(d.loc[c, "tsne_x"], d.loc[c, "tsne_y"], s=20, marker="*", color="#E64B35",
               linewidths=0, rasterized=True, label=f"Clp-specific ({int(c.sum()):,})")
    ax.legend(fontsize=SS - 1, frameon=False, loc="lower left")
    ax.set_xticks([]); ax.set_yticks([])
    stylia.label(ax, xlabel="ESM-C tSNE-1", ylabel="ESM-C tSNE-2",
                 title=f"Where the evidence is — {orgname}")

    # ---- panel 5: do computed features predict measured turnover? ----
    ax = axs.next()
    vpath = D.degradability_processed_dir(org) / f"{prefix}_disorder_validation.csv"
    if vpath.exists():
        v = pd.read_csv(vpath).dropna(subset=["auroc"]).copy()
        v = v.reindex(v["auroc"].sub(0.5).abs().sort_values(ascending=False).index).head(12)
        v = v.iloc[::-1]
        y = np.arange(len(v))
        crosses = (v["auroc_lo"] <= 0.5) & (v["auroc_hi"] >= 0.5)
        ax.hlines(y, v["auroc_lo"], v["auroc_hi"], color="#BBBBBB", lw=2)
        ax.scatter(v["auroc"], y, s=22,
                   color=[GREY if c else NPG[0] for c in crosses], zorder=3, linewidths=0)
        ax.axvline(0.5, ls="--", color="#555555", lw=1)
        ax.set_yticks(y)
        ax.set_yticklabels([str(f)[:26] for f in v["feature"]], fontsize=SS - 2)
        ax.margins(y=0.10)
        ax.text(0.98, 0.97,
                f"{int(crosses.sum())}/{len(v)} intervals cross 0.5",
                transform=ax.transAxes, fontsize=SS, va="top", ha="right", color="#555555")
    else:
        ax.text(0.5, 0.5, "no validation table (run 10d)", ha="center", va="center",
                transform=ax.transAxes, fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="AUROC vs measured turnover (95% CI)", ylabel="",
                 title=f"Do computed features predict it? — {orgname}")

    # ---- panel 6: native vs transferred ----
    ax = axs.next()
    src = d.get("measured_source")
    if src is not None and src.notna().any():
        vc = src.value_counts()
        labels = {"native": "measured in this organism",
                  "ortholog_transfer": "transferred from E. coli"}
        names6 = [labels.get(k, k) for k in vc.index]
        cols = [NPG[0] if k == "native" else "#C98A1E" for k in vc.index]
        bars = ax.bar(range(len(vc)), vc.to_numpy(), color=cols, width=0.5)
        ax.bar_label(bars, padding=2, fontsize=SS, fmt="%d")
        ax.set_xticks(range(len(vc))); ax.set_xticklabels(names6, fontsize=SS - 1)
        ax.margins(y=0.25)
        # The activator layer is transferred for both organisms, from S. aureus.
        n_act = int(d["activator_evidence"].notna().sum())
        ax.text(0.98, 0.97,
                f"plus {n_act:,} from S. aureus\n(activator, cross-phylum RBH)",
                transform=ax.transAxes, fontsize=SS, va="top", ha="right", color="#555555",
                linespacing=1.4)
    else:
        ax.text(0.5, 0.5, "no provenance recorded", ha="center", va="center",
                transform=ax.transAxes, fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Whose measurement is it? — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"10n_landscape_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")
    print(f"[{org}] any measured {int(any_evidence.sum())}/{n_total} "
          f"| Clp-specific {int(clp_specific.sum())}")


if __name__ == "__main__":
    main()
