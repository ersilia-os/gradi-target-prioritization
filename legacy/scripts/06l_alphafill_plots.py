"""2x3 overview of the AlphaFill ligand-transfer annotation (06d), structure-based binding evidence.

AlphaFill (alphafill.eu) transplants ligands/cofactors from homologous PDB entries onto the
AlphaFold model: if a homolog co-crystallised a small molecule, that ligand is placed at the
structurally equivalent site — evidence that *something* binds there. Each transplant carries the
donor homolog's sequence identity and a global structural RMSD, so AlphaFill's distinctive signal
is transfer CONFIDENCE (high identity + low RMSD), which these panels foreground.

"Drug-like" is the STRICT tier from src/ligandability.py (a real bound molecule that is NOT a
promiscuous cofactor/nucleotide, free amino acid, sugar, lipid/detergent, buffer/cryo or solvent);
the broad "any transplant" count stays in the 06d CSV. Same house style as 06i/06j/06k: stylia
"slide" format, NPG palette, white/dark in-bar labels. Every panel is specific to the single
`--organism` of the slide.

  1  coverage & two-tier funnel    AlphaFill available / any-ligand transplant / drug-like transplant
  2  drug-like vs cofactor-only    of proteins with a transplant, how many are genuinely drug-like
  3  transplant confidence map     donor identity (%) vs global RMSD (Å) per drug-like protein
  4  top targets                   # distinct drug-like transplanted ligands
  5  most frequent transplants     which drug-like ligands recur (residual additives flagged)
  6  transplant depth              # distinct drug-like ligands per protein (proteins with >=1)

Reads output/results/<org>/<prefix>_alphafill.csv (06d). Output: output/plots/06l_alphafill_<prefix>.png
(one slide per --organism). Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)  # stylia rmtree's this on import; ensure it exists
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

REPO_ROOT = L.REPO_ROOT
TOP_N = 8
NPG = stylia.CategoricalPalette("npg").colors
ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}                 # NPG per organism (matches 06i/06j/06k)
DRUGLIKE_COLOR = "#00A087"                                     # NPG green
COFACTOR_COLOR = "#B7B7B7"
ID_MIN = 40.0    # donor identity (%) above which a transplant is "confident"
RMSD_MAX = 2.0   # global RMSD (Å) below which the structural fit is "confident"
ORGANISMS = [("kpneumoniae", "kp", "K. pneumoniae"), ("ecoli", "ec", "E. coli K-12")]


def druglike_set(s) -> list[str]:
    """Strict drug-like ligands from the broad `alphafill_ligand_ids` column."""
    if not isinstance(s, str) or not s:
        return []
    return [x for x in s.split(";") if x and L.is_druglike_ligand(x)]


def load(prefix: str, org: str) -> pd.DataFrame:
    d = pd.read_csv(L.results_dir(org) / f"{prefix}_alphafill.csv")
    d["available"] = d["alphafill_available"].fillna(False).astype(bool)
    d["n_broad"] = d["alphafill_n_ligand_any"].fillna(0)
    d["n_dl"] = d["alphafill_n_druglike"].fillna(0)
    d["ident_pct"] = d["alphafill_best_identity"].astype(float) * 100.0
    d["rmsd"] = d["alphafill_best_global_rmsd"].astype(float)
    d["dl_ligs"] = d["alphafill_ligand_ids"].map(druglike_set)
    return d


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    orgname = {p: n for _, p, n in ORGANISMS}[prefix]
    d = load(prefix, org)
    genes = L.load_genes(org)
    col = ORG_COLOR.get(prefix, NPG[0])

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)  # 16:9 slide, 2x3 panels

    # ---- panel 1: coverage & two-tier funnel (this organism) ----
    ax = axs.next()
    labels = ["AlphaFill\navailable", "Any-ligand\ntransplant", "Drug-like\ntransplant"]
    counts = [int(d["available"].sum()), int((d["n_broad"] > 0).sum()), int((d["n_dl"] > 0).sum())]
    bars = ax.bar(range(len(labels)), counts, color=[NPG[5], COFACTOR_COLOR, DRUGLIKE_COLOR])
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels)
    ax.set_ylabel("Number of proteins"); ax.set_xlabel("")
    ax.margins(y=0.18)
    ax.set_title(f"AlphaFill coverage — {orgname}")
    ax.text(0.5, -0.22, f"of {len(d):,} proteins (counts nested left→right)",
            transform=ax.transAxes, ha="center", va="top", color="#777777", fontsize="small")

    # ---- panel 2: drug-like vs cofactor/artifact-only, of proteins with a transplant ----
    ax = axs.next()
    n_dl = int((d["n_dl"] > 0).sum())
    n_cof = int(((d["n_broad"] > 0) & (d["n_dl"] == 0)).sum())
    b1 = ax.bar([0], [n_dl], color=DRUGLIKE_COLOR, label="Drug-like transplant")
    b2 = ax.bar([0], [n_cof], bottom=[n_dl], color=COFACTOR_COLOR, label="Cofactor / artifact-only")
    ax.bar_label(b1, label_type="center", color="white")
    ax.bar_label(b2, label_type="center", color="#444444")
    ax.bar_label(b2, labels=[f"{n_dl + n_cof}"], padding=2)
    ax.set_xticks([0]); ax.set_xticklabels(["Proteins with\na transplant"]); ax.set_xlim(-0.9, 0.9)
    ax.set_ylabel("Number of proteins"); ax.set_xlabel("")
    ax.margins(y=0.18)
    ax.set_title(f"Transplant quality split — {orgname}")
    ax.legend(loc="upper right")

    # ---- panel 3: transplant confidence map — donor identity vs global RMSD (drug-like) ----
    # ~2.5k transplants overplot badly as a scatter; show density as a hexbin (darker = more),
    # on a white→organism-colour ramp, and highlight the high-confidence quadrant.
    ax = axs.next()
    sub = d[(d["n_dl"] > 0) & d["ident_pct"].notna() & d["rmsd"].notna()]
    ymax = min(12.0, float(sub["rmsd"].max()) * 1.05) if len(sub) else 12.0
    cmap = LinearSegmentedColormap.from_list("org_density", ["#FFFFFF", col])
    ax.hexbin(sub["ident_pct"], sub["rmsd"], gridsize=26, bins="log", cmap=cmap, mincnt=1,
              extent=(20, 100, 0, ymax))
    # high-confidence quadrant (close donor + tight fit): shade via axes-fraction xmin
    ax.axhspan(0, RMSD_MAX, xmin=(ID_MIN - 20) / 80.0, xmax=1.0,
               color="#00A087", alpha=0.10, zorder=0)
    ax.axhline(RMSD_MAX, color="#777777", linestyle="--", linewidth=0.7)
    ax.axvline(ID_MIN, color="#777777", linestyle="--", linewidth=0.7)
    n_conf = int(((sub["ident_pct"] >= ID_MIN) & (sub["rmsd"] <= RMSD_MAX)).sum())
    ax.set_ylim(0, ymax); ax.set_xlim(20, 100)
    ax.set_xlabel("Best donor identity (%)"); ax.set_ylabel("Global RMSD (Å)")
    ax.set_title(f"Transplant confidence — {orgname}")
    ax.text(0.97, 0.95, f"high-confidence: {n_conf}\n(≥{ID_MIN:.0f}% id, ≤{RMSD_MAX:.0f} Å)",
            transform=ax.transAxes, ha="right", va="top", color="#2c6e63", fontsize="small")

    # ---- panel 4: top targets by # distinct drug-like transplanted ligands ----
    ax = axs.next()
    top = d[d["n_dl"] > 0].sort_values("n_dl", ascending=False).head(TOP_N).iloc[::-1]
    ys = range(len(top))
    bars = ax.barh(list(ys), top["n_dl"].to_numpy(), color=col)
    ax.set_yticks(list(ys)); ax.set_yticklabels([genes.get(a) or a for a in top["uniprot_accession"]])
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_xlim(0, top["n_dl"].max() * 1.12 if len(top) else 1)
    ax.set_ylabel(""); ax.set_xlabel("# distinct drug-like transplanted ligands")
    ax.set_title(f"Top targets — {orgname}")

    # ---- panel 5: most frequent drug-like transplanted ligands ----
    ax = axs.next()
    freq = Counter()
    for s in d["dl_ligs"]:
        freq.update(s)
    common = freq.most_common(TOP_N)[::-1]
    ys = range(len(common))
    colors = [NPG[i % len(NPG)] for i in range(len(common))][::-1]
    bars = ax.barh(list(ys), [c for _, c in common], color=colors)
    ax.set_yticks(list(ys)); ax.set_yticklabels([lig for lig, _ in common])
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_xlim(0, (common[-1][1] if common else 1) * 1.12)
    ax.set_ylabel(""); ax.set_xlabel("Proteins with this ligand")
    ax.set_title(f"Most frequent transplants — {orgname}")

    # ---- panel 6: transplant depth — # drug-like ligands per protein (>=1) ----
    ax = axs.next()
    n = d.loc[d["n_dl"] > 0, "n_dl"].to_numpy()
    bin_labels = ["1", "2–4", "5–9", "10–19", "20+"]
    binned = [int((n == 1).sum()),
              int(((n >= 2) & (n <= 4)).sum()),
              int(((n >= 5) & (n <= 9)).sum()),
              int(((n >= 10) & (n <= 19)).sum()),
              int((n >= 20).sum())]
    bars = ax.bar(range(len(bin_labels)), binned, color=col)
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(len(bin_labels))); ax.set_xticklabels(bin_labels)
    ax.set_xlabel("# distinct drug-like ligands"); ax.set_ylabel("Number of proteins")
    ax.margins(y=0.15)
    ax.set_title(f"Transplant depth — {orgname}")
    if len(n):
        ax.text(0.97, 0.92, f"median {int(np.median(n))} · max {int(n.max())}",
                transform=ax.transAxes, ha="right", va="top", color="#777777", fontsize="small")

    out = REPO_ROOT / "output" / "plots" / f"06l_alphafill_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
