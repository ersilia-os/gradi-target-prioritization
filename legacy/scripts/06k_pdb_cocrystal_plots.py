"""2x3 overview of the PDB co-crystal annotation (06c), direct vs homolog-transferred.

Co-crystallised drug-like ligands are empirical proof a small molecule binds the protein (DIRECT:
the protein's own structure, or a >=95%-identity PDB chain) or a homolog (TRANSFERRED: a <95%
chain by sequence, or an ortholog by accession). "Drug-like" here is the STRICT tier from
src/ligandability.py: a real bound molecule that is NOT a promiscuous cofactor/nucleotide (ATP,
NAD, FAD, heme, Fe-S clusters, …); the broad "any bound ligand" counts stay in the 06c CSV. Same
house style as 06i/06j and the agnostic deck: stylia "slide" format, NPG palette, white/dark in-bar
labels. Every panel is specific to the single `--organism` of the slide — no content is shared
between the kp and ec slides.

  1  co-crystal evidence by route        own / seq>=95% / seq<95% / ortholog-acc
  2  own experimental structures         of proteins with a PDB structure, liganded vs apo-only
  3  homolog donor organisms             which organisms supply the transferred co-crystals
  4  top targets                         # distinct drug-like co-crystal ligands, coloured by route
  5  most frequent co-crystal ligands    which ligands recur across the proteome
  6  ligand evidence depth               # distinct ligands per protein (proteins with >=1)

Reads output/results/<org>/<prefix>_pdb_cocrystals.csv (06c). Output: output/plots/06k_pdb_<prefix>.png
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
from matplotlib.patches import Patch  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

REPO_ROOT = L.REPO_ROOT
TOP_N = 8
NPG = stylia.CategoricalPalette("npg").colors
ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}                  # NPG per organism (matches 06i/06j)
ROUTE_COLOR = {"direct": "#00A087", "ortholog": "#F39B7F"}      # green (direct) / salmon (transferred)
ROUTE_LABEL = {"direct": "Direct (own / ≥95% id)", "ortholog": "Homolog-transferred"}
APO_COLOR = "#B7B7B7"
# 4-way provenance for the route-breakdown panel
PROV = [
    ("pdb_lig_direct_has_druglike", "Own\nstructure", NPG[3]),
    ("pdb_lig_seqdirect_has_druglike", "Seq\n≥95%", NPG[2]),
    ("pdb_lig_seqhom_has_druglike", "Seq\n<95%", NPG[4]),
    ("pdb_lig_ortho_has_druglike", "Ortholog\n(acc)", NPG[5]),
]
LIG_COLS = ["pdb_lig_direct_ligand_ids", "pdb_lig_seqdirect_ligand_ids",
            "pdb_lig_seqhom_ligand_ids", "pdb_lig_ortho_ligand_ids"]
ORGANISMS = [("kpneumoniae", "kp", "K. pneumoniae"), ("ecoli", "ec", "E. coli K-12")]


def ligand_set(row) -> set[str]:
    """Strict drug-like ligands across all routes. The `*_ligand_ids` columns are the BROAD set
    (cofactors/nucleotides included); these plots show the strict drug-like tier, so filter here."""
    out: set[str] = set()
    for c in LIG_COLS:
        v = row.get(c)
        if isinstance(v, str) and v:
            out.update(t for t in v.split(";") if t and L.is_druglike_ligand(t))
    return out


def short_org(name: str) -> str:
    """Donor-species label -> abbreviated binomial, strain suffix dropped.

    The 06c column smushes most names ('Ecoli_K12_MG1655', 'Bsubtilis_168') but leaves a few as a
    proper binomial ('Homo_sapiens'). A second token that is all-lowercase alpha marks the proper
    case (-> 'H. sapiens'); otherwise the first token is 'Initial+species' and we split it
    ('Ecoli' -> 'E. coli', 'Vcholerae' -> 'V. cholerae')."""
    parts = str(name).split("_")
    if not parts or not parts[0]:
        return str(name)
    head = parts[0]
    if len(parts) >= 2 and parts[1].isalpha() and parts[1].islower():
        return f"{head[0]}. {parts[1]}"          # proper 'Genus species'
    return f"{head[0]}. {head[1:]}" if len(head) > 1 else head  # smushed 'Initial+species'


def load(prefix: str, org: str) -> pd.DataFrame:
    d = pd.read_csv(L.results_dir(org) / f"{prefix}_pdb_cocrystals.csv")
    bool_cols = [p[0] for p in PROV] + ["pdb_lig_any_has_druglike", "pdb_lig_anydirect_has_druglike",
                                        "pdb_lig_direct_has_structure", "pdb_lig_any_has_ligand"]
    for c in bool_cols:
        d[c] = d.get(c, False)
        d[c] = d[c].fillna(False).astype(bool)
    d["ligs"] = d.apply(ligand_set, axis=1)
    d["n_ligs"] = d["ligs"].map(len)
    d["route"] = np.where(d["pdb_lig_anydirect_has_druglike"], "direct", "ortholog")
    return d


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

    # ---- panel 1: co-crystal evidence by route (this organism) ----
    ax = axs.next()
    labels = [p[1] for p in PROV]
    counts = [int(d[p[0]].sum()) for p in PROV]
    colors = [p[2] for p in PROV]
    bars = ax.bar(range(len(PROV)), counts, color=colors)
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(len(PROV))); ax.set_xticklabels(labels)
    ax.set_xlabel(""); ax.set_ylabel("Proteins with a drug-like co-crystal")
    ax.margins(y=0.18)
    ax.set_title(f"Co-crystal evidence by route — {orgname}")
    # two-tier note: how many proteins bind ONLY a cofactor/nucleotide (excluded from drug-like)
    n_cofactor_only = int((d["pdb_lig_any_has_ligand"] & ~d["pdb_lig_any_has_druglike"]).sum())
    ax.text(0.5, -0.22,
            f"\"Drug-like\" excludes cofactors/nucleotides (ATP, NAD, FAD, heme, Fe-S…).\n"
            f"A further {n_cofactor_only:,} proteins bind only such a cofactor.",
            transform=ax.transAxes, ha="center", va="top", color="#777777", fontsize="small")

    # ---- panel 2: own experimental structures, liganded vs apo-only (this organism) ----
    ax = axs.next()
    n_struct = int(d["pdb_lig_direct_has_structure"].sum())
    n_lig = int((d["pdb_lig_direct_has_structure"] & d["pdb_lig_direct_has_druglike"]).sum())
    n_apo = n_struct - n_lig
    b1 = ax.bar([0], [n_lig], color=ROUTE_COLOR["direct"], label="With drug-like co-crystal")
    b2 = ax.bar([0], [n_apo], bottom=[n_lig], color=APO_COLOR, label="Apo / cofactor-only")
    ax.bar_label(b1, label_type="center", color="white")
    ax.bar_label(b2, label_type="center", color="#444444")
    ax.bar_label(b2, labels=[f"{n_struct}"], padding=2)
    ax.set_xticks([0]); ax.set_xticklabels(["Own PDB\nstructure"])
    ax.set_xlim(-0.9, 0.9)
    ax.set_xlabel(""); ax.set_ylabel("Number of proteins")
    ax.margins(y=0.18)
    ax.set_title(f"Own experimental structures — {orgname}")
    ax.legend(loc="upper right")

    # ---- panel 3: homolog donor organisms (transferred evidence, this organism) ----
    ax = axs.next()
    sub = d.loc[d["pdb_lig_ortho_has_druglike"], "pdb_lig_ortho_best_species"].dropna()
    src = sub.map(short_org)
    sc = src.value_counts().sort_values(ascending=True).tail(TOP_N)
    ys = range(len(sc))
    colors = [NPG[i % len(NPG)] for i in range(len(sc))][::-1]  # distinct hue per organism
    bars = ax.barh(list(ys), sc.to_numpy(), color=colors)
    ax.set_yticks(list(ys)); ax.set_yticklabels(sc.index.tolist())
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_xlim(0, (sc.max() * 1.12) if len(sc) else 1)
    ax.set_ylabel(""); ax.set_xlabel("Proteins (ortholog route)")
    ax.set_title(f"Homolog donor organisms — {orgname}")

    # ---- panel 4: top targets by # distinct co-crystal ligands, coloured by route ----
    ax = axs.next()
    top = d[d["n_ligs"] > 0].sort_values("n_ligs", ascending=False).head(TOP_N).iloc[::-1]
    ys = range(len(top))
    colors = [ROUTE_COLOR[r] for r in top["route"]]
    bars = ax.barh(list(ys), top["n_ligs"].to_numpy(), color=colors)
    ax.set_yticks(list(ys)); ax.set_yticklabels([genes.get(a) or a for a in top["uniprot_accession"]])
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_xlim(0, top["n_ligs"].max() * 1.12 if len(top) else 1)
    ax.set_ylabel(""); ax.set_xlabel("# distinct drug-like co-crystal ligands")
    ax.set_title(f"Top targets — {orgname}")
    ax.legend(handles=[Patch(color=ROUTE_COLOR["direct"], label=ROUTE_LABEL["direct"]),
                       Patch(color=ROUTE_COLOR["ortholog"], label=ROUTE_LABEL["ortholog"])],
              loc="lower right")

    # ---- panel 5: most frequent co-crystallised ligands, one vivid colour each ----
    ax = axs.next()
    freq = Counter()
    for s in d["ligs"]:
        freq.update(s)
    common = freq.most_common(TOP_N)[::-1]
    ys = range(len(common))
    colors = [NPG[i % len(NPG)] for i in range(len(common))][::-1]
    bars = ax.barh(list(ys), [c for _, c in common], color=colors)
    ax.set_yticks(list(ys)); ax.set_yticklabels([lig for lig, _ in common])
    ax.bar_label(bars, label_type="edge", padding=2)
    ax.set_xlim(0, (common[-1][1] if common else 1) * 1.12)
    ax.set_ylabel(""); ax.set_xlabel("Proteins with this ligand")
    ax.set_title(f"Most frequent co-crystal ligands — {orgname}")

    # ---- panel 6: ligand evidence depth — # distinct ligands per protein (>=1) ----
    ax = axs.next()
    n = d.loc[d["n_ligs"] > 0, "n_ligs"].to_numpy()
    bin_labels = ["1", "2–4", "5–9", "10–19", "20+"]
    binned = [int((n == 1).sum()),
              int(((n >= 2) & (n <= 4)).sum()),
              int(((n >= 5) & (n <= 9)).sum()),
              int(((n >= 10) & (n <= 19)).sum()),
              int((n >= 20).sum())]
    bars = ax.bar(range(len(bin_labels)), binned, color=ORG_COLOR.get(prefix, NPG[0]))
    ax.bar_label(bars, padding=2)
    ax.set_xticks(range(len(bin_labels))); ax.set_xticklabels(bin_labels)
    ax.set_xlabel("# distinct drug-like ligands"); ax.set_ylabel("Number of proteins")
    ax.margins(y=0.15)
    ax.set_title(f"Ligand evidence depth — {orgname}")
    ax.text(0.97, 0.92, f"median {int(np.median(n))} · max {int(n.max())}",
            transform=ax.transAxes, ha="right", va="top", color="#777777", fontsize="small")

    out = REPO_ROOT / "output" / "plots" / f"06k_pdb_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(out))
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
