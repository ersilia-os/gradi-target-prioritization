"""Membrane topology and envelope architecture (docs §5.1d/e).

The structural half of the axis: what TMbed and the signal-peptide track actually see, and how those
observations turn into the graded Clp-accessibility score. Panel 1 is the important one — it is the
evidence that the 0.6 / 0.2 split for inner-membrane proteins is a measurement rather than a guess.

  1  cytoplasm-facing fraction   IM proteins either side of the 0.30 cut that grades them 0.6 vs 0.2
  2  beta-barrel evidence        transmembrane strand counts, the >=8 cut, known OM markers annotated
  3  TM helices per compartment  polytopic inner-membrane proteins vs everything else
  4  architecture per compartment which topological features each compartment is built from
  5  export signals              signal-peptide and Sec/SPII lipoprotein calls, with Lol +2 sorting
  6  topology overrides          where topology overruled a predicted compartment call

Reads output/results/<org>/<prefix>_localization.csv.
Output: output/plots/09k_topology_<prefix>.png (one slide per --organism). Run with the `gradi` env.
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
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL
CLASSES = [c for c in LOC.LOC_CLASS_ORDER if c not in ("membrane", "unknown")]

ACCESSIBLE_C = "#00A087"
BURIED_C = "#F39B7F"
BARREL_C = "#8491B4"

# Canonical Gram-negative OM beta-barrels, for annotating panel 2. Strand counts are the textbook
# values; they are what the panel is checking TMbed against.
BARREL_MARKERS = ["ompA", "ompC", "ompF", "lamB", "btuB", "fhuA", "lptD"]

OVERRIDE_LABEL = {
    "membrane_resolved": "UniProt “Membrane”\nside resolved",
    "beta_barrel": "β-barrel →\nouter membrane",
    "lipoprotein_lol_rule": "lipoprotein →\nLol +2 rule",
}
OVERRIDE_COLOR = {"membrane_resolved": NPG[1], "beta_barrel": BARREL_C, "lipoprotein_lol_rule": "#B05CC8"}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    orgname = LOC.ORG_DISPLAY[org]
    d = pd.read_csv(LOC.results_dir(org) / f"{prefix}_localization.csv")
    genes = LOC.load_genes(org)
    d["gene"] = d["uniprot_accession"].map(genes).fillna("")

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: the 0.30 cytoplasm-facing cut ----
    ax = axs.next()
    im = d[d["localization"] == "inner_membrane"]
    frac = im["cyto_residue_fraction"].dropna()
    cut = LOC.CYTO_DOMAIN_FRACTION_MIN
    bins = np.linspace(0, 1, 41)
    below = frac[frac < cut]; above = frac[frac >= cut]
    ax.hist(below, bins=bins, color=BURIED_C, alpha=0.95, label=f"< {cut:g} → 0.2  ({len(below):,})")
    ax.hist(above, bins=bins, color=ACCESSIBLE_C, alpha=0.95, label=f"≥ {cut:g} → 0.6  ({len(above):,})")
    ax.axvline(cut, ls="--", color="#555555", lw=1.2)
    ax.legend(fontsize=SS, frameon=False, loc="upper right")
    # The pile-up at 1.0 is not an artefact: these are proteins their source calls inner-membrane
    # but in which TMbed finds no membrane-spanning segment at all — peripheral membrane proteins,
    # sitting on the cytoplasmic face. Scoring them 0.6 is the intended behaviour.
    n_peripheral = int(((im["cyto_residue_fraction"] >= 0.999)
                        & (im["n_tm_helix"].fillna(0) == 0)).sum())
    if n_peripheral:
        ax.annotate(f"{n_peripheral} with no TM helix\n(peripheral, membrane-associated)",
                    xy=(1.0, ax.get_ylim()[1] * 0.58), xytext=(0.72, ax.get_ylim()[1] * 0.72),
                    fontsize=SS, ha="center", va="center", color="#555555", linespacing=1.35,
                    arrowprops=dict(arrowstyle="->", color="#BBBBBB", lw=0.8))
    stylia.label(ax, xlabel="fraction of residues facing the cytoplasm", ylabel="inner-membrane proteins",
                 title=f"What grades an IM protein — {orgname}")

    # ---- panel 2: beta-barrel strand counts ----
    ax = axs.next()
    strands = d["n_tm_strand"].dropna()
    has = strands[strands > 0]
    top = int(min(has.max(), 30)) if len(has) else 8
    bins = np.arange(0.5, top + 1.5, 1)
    ax.hist(has[has < 8], bins=bins, color="#C9C9C7", label="< 8 strands (noise)")
    ax.hist(has[has >= 8], bins=bins, color=BARREL_C, label=f"≥ 8 = β-barrel ({int((d['is_beta_barrel'] == True).sum()):,})")
    ax.axvline(7.5, ls="--", color="#555555", lw=1.2)
    ann = []
    for g in BARREL_MARKERS:
        r = d[d["gene"].str.lower() == g.lower()]
        if len(r) and pd.notna(r.iloc[0]["n_tm_strand"]):
            ann.append((g, int(r.iloc[0]["n_tm_strand"])))
    if ann:
        ax.text(0.97, 0.97, "\n".join(f"{g}  {n}" for g, n in sorted(ann, key=lambda t: t[1])),
                transform=ax.transAxes, fontsize=SS, va="top", ha="right", color="#2B2333",
                family="monospace", linespacing=1.35)
    ax.legend(fontsize=SS, frameon=False, loc="upper left")
    stylia.label(ax, xlabel="predicted transmembrane β-strands", ylabel="proteins",
                 title=f"β-barrel evidence — {orgname}")

    # ---- panel 3: TM helices per compartment ----
    ax = axs.next()
    present = [c for c in CLASSES if (d["localization"] == c).sum() > 0]
    data = [d.loc[d["localization"] == c, "n_tm_helix"].fillna(0).to_numpy() for c in present]
    bp = ax.boxplot(data, patch_artist=True, widths=0.55, showfliers=False,
                    medianprops=dict(color="#2B2333", lw=1.4))
    for patch, c in zip(bp["boxes"], present):
        patch.set_facecolor(LOC.LOC_CLASS_COLOR[c]); patch.set_edgecolor("#8A8A96"); patch.set_alpha(0.9)
    ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in present], fontsize=SS)
    stylia.label(ax, xlabel="", ylabel="TM α-helices (TMbed)",
                 title=f"Polytopic membrane proteins — {orgname}")

    # ---- panel 4: architecture composition per compartment ----
    ax = axs.next()
    feats = [
        ("β-barrel", BARREL_C, lambda s: s["is_beta_barrel"] == True),                       # noqa: E712
        ("TM helix", BURIED_C, lambda s: (s["is_beta_barrel"] != True) & (s["n_tm_helix"].fillna(0) > 0)),  # noqa: E712
        ("signal peptide only", NPG[1],
         lambda s: (s["is_beta_barrel"] != True) & (s["n_tm_helix"].fillna(0) == 0)          # noqa: E712
         & (s["has_signal_peptide_tmbed"] == True)),                                          # noqa: E712
        ("soluble", "#C9C9C7",
         lambda s: (s["is_beta_barrel"] != True) & (s["n_tm_helix"].fillna(0) == 0)          # noqa: E712
         & (s["has_signal_peptide_tmbed"] != True)),                                          # noqa: E712
    ]
    bottom = np.zeros(len(present))
    for name, colour, pred in feats:
        vals = np.array([100 * pred(d[d["localization"] == c]).sum()
                         / max((d["localization"] == c).sum(), 1) for c in present])
        ax.bar(range(len(present)), vals, bottom=bottom, color=colour, label=name, width=0.7)
        bottom += vals
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([LOC.LOC_ABBREV[c] for c in present], fontsize=SS)
    ax.set_ylim(0, 100)
    # Bars fill the axes edge to edge, so the legend needs its own ground to stay legible.
    ax.legend(fontsize=SS - 1, ncol=2, loc="lower center", frameon=True, facecolor="white",
              framealpha=0.93, edgecolor="none")
    stylia.label(ax, xlabel="", ylabel="% of compartment",
                 title=f"What each compartment is built from — {orgname}")

    # ---- panel 5: export signals and lipoprotein sorting ----
    ax = axs.next()
    n_sp_tmbed = int((d["has_signal_peptide_tmbed"] == True).sum())      # noqa: E712
    n_sp_uni = int((d["has_signal_peptide"] == True).sum())              # noqa: E712
    n_lipo = int((d["is_lipoprotein"] == True).sum())                    # noqa: E712
    lol = d["lipoprotein_sorting"].dropna().value_counts()
    cats = [("signal peptide\n(TMbed)", n_sp_tmbed, NPG[1]),
            ("signal peptide\n(UniProt)", n_sp_uni, "#C9C9C7"),
            ("Sec/SPII\nlipoprotein", n_lipo, "#B05CC8"),
            ("→ outer mem.\n(Lol +2)", int(lol.get("outer_membrane", 0)), BARREL_C),
            ("→ inner mem.\n(Lol +2)", int(lol.get("inner_membrane", 0)), BURIED_C)]
    bars = ax.bar(range(len(cats)), [n for _, n, _ in cats], color=[c for _, _, c in cats], width=0.66)
    ax.bar_label(bars, padding=2, fontsize=SS, fmt="%d")
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels([lab for lab, _, _ in cats], fontsize=SS - 1)
    ax.margins(y=0.2)
    src = d["signalp_source"].dropna()
    note = "SignalP 6.0" if (src == "signalp6").any() else "lipobox fallback (SignalP not installed)"
    ax.text(0.98, 0.97, note, transform=ax.transAxes, fontsize=SS, va="top", ha="right",
            color="#C98A1E" if "fallback" in note else "#555555")
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Export signals — {orgname}")

    # ---- panel 6: topology overrides ----
    ax = axs.next()
    ov = d.loc[d["localization_override"].fillna("") != "", "localization_override"].value_counts()
    if len(ov):
        labs = [OVERRIDE_LABEL.get(k, k) for k in ov.index]
        cols = [OVERRIDE_COLOR.get(k, "#C9C9C7") for k in ov.index]
        bars = ax.bar(range(len(ov)), ov.to_numpy(), color=cols, width=0.6)
        ax.bar_label(bars, padding=2, fontsize=SS, fmt="%d")
        ax.set_xticks(range(len(ov)))
        ax.set_xticklabels(labs, fontsize=SS - 1)
        ax.margins(y=0.30)
        ax.text(0.97, 0.60,
                "topology may overrule a\npredicted call only — never an\nexperimental or curated one",
                transform=ax.transAxes, fontsize=SS, ha="right", va="top", color="#555555",
                linespacing=1.4)
    else:
        ax.text(0.5, 0.5, "no overrides fired", ha="center", va="center", transform=ax.transAxes,
                fontsize=SS, color="#555555")
    stylia.label(ax, xlabel="", ylabel="proteins", title=f"Topology overrides — {orgname}")

    out = LOC.REPO_ROOT / "output" / "plots" / f"09k_topology_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)}")


if __name__ == "__main__":
    main()
