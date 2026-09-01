"""Structure snapshots of the top druggable AlphaFold targets, with the predicted pocket highlighted.

A companion to the 06m pocket/druggability poster: instead of statistics, show the actual models.
For the top druggable targets (highest pLDDT-weighted pocket consensus from 06e) we render the
AlphaFold model as a ray-traced cartoon coloured by the official per-residue pLDDT palette
(orange <50 → yellow 50–70 → cyan 70–90 → blue 90+), and overlay the residues of the top
P2Rank-predicted binding pocket as green sticks + a translucent green surface — i.e. where a small
molecule is predicted to bind. Six targets are tiled into one slide per organism.

Two conda envs (see install.sh): the cartoons are ray-traced by PyMOL in the **`gradi-pymol`** env
(`conda create -n gradi-pymol -c conda-forge pymol-open-source`); target selection and the montage
run in `gradi`. This script runs in `gradi` and shells out to `gradi-pymol` for the rendering via
scripts/_06n_pymol_render.py.

Reads output/results/<org>/<prefix>_pockets.csv (06e consensus ranking), the per-protein pocket
residues from data/processed/<org>/pockets/p2rank_run/<acc>.pdb_predictions.csv, and the AlphaFold
cif models. Output: output/plots/06n_structures_<prefix>.png (one slide per --organism).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)  # stylia rmtree's this on import; ensure it exists
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

REPO_ROOT = L.REPO_ROOT
N_TARGETS = 6
PYMOL_ENV = "gradi-pymol"
HELPER = Path(__file__).resolve().parent / "_06n_pymol_render.py"
POCKET_COLOR = "#00A087"
PLDDT_BANDS = [("#0053D6", "≥90 very high"), ("#65CBF3", "70–90 confident"),
               ("#FFDB13", "50–70 low"), ("#FF7D45", "<50 very low")]
ORGNAME = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli K-12"}


def pocket_resids(org: str, acc: str) -> list[int]:
    """Residue numbers of the top-ranked P2Rank pocket (from its residue_ids column)."""
    p = L.processed_dir(org, "pockets", "p2rank_run") / f"{acc}.pdb_predictions.csv"
    if not p.exists():
        return []
    df = pd.read_csv(p, skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]
    if df.empty or "residue_ids" not in df.columns:
        return []
    top = df.sort_values("rank").iloc[0]
    out = set()
    for tok in str(top.get("residue_ids") or "").split():
        if "_" in tok:
            try:
                out.add(int(tok.split("_")[1]))
            except ValueError:
                pass
    return sorted(out)


def pick_targets(org: str, prefix: str) -> list[dict]:
    """Top druggable targets (by pocket consensus) that have a cif model and a P2Rank pocket."""
    d = pd.read_csv(L.results_dir(org) / f"{prefix}_pockets.csv")
    d["consensus"] = d["pocket_consensus_score"].fillna(0.0)
    genes = L.load_genes(org)
    cand = (d[(d["p2rank_n_pockets"].fillna(0) > 0) & d["af_mean_plddt"].notna()]
            .sort_values("consensus", ascending=False))
    jobs = []
    for _, r in cand.iterrows():
        acc = r["uniprot_accession"]
        cif = L.af_cif_path(org, acc)
        res = pocket_resids(org, acc)
        if not cif.exists() or not res:
            continue
        jobs.append({"acc": acc, "label": genes.get(acc) or acc, "cif": str(cif),
                     "consensus": round(float(r["consensus"]), 2),
                     "plddt": round(float(r["af_mean_plddt"])), "pocket": res})
        if len(jobs) >= N_TARGETS:
            break
    return jobs


def montage(jobs: list[dict], img_dir: Path, prefix: str, orgname: str) -> Path:
    imgs = sorted(img_dir.glob("*.png"))
    fig = plt.figure(figsize=(15.5, 9))
    fig.suptitle(f"Top druggable AlphaFold targets — {orgname}", fontsize=16, fontweight="bold", y=0.975)
    for i, (img, j) in enumerate(zip(imgs, jobs)):
        ax = fig.add_subplot(2, 3, i + 1)
        ax.imshow(mpimg.imread(str(img)))
        ax.set_axis_off()
        ax.set_title(f"{j['label']}   ·   consensus {j['consensus']:.2f}   ·   pLDDT {j['plddt']}",
                     fontsize=12, pad=3)
        ax.text(0.5, -0.04, f"{len(j['pocket'])} predicted pocket residues",
                transform=ax.transAxes, ha="center", va="top", color="#5C6E75", fontsize=9.5)
    handles = [plt.Line2D([0], [0], color=c, lw=5, label=lab) for c, lab in PLDDT_BANDS]
    handles.append(plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=POCKET_COLOR,
                              markeredgecolor="white", markersize=11, label="predicted binding pocket"))
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=11,
               bbox_to_anchor=(0.5, 0.01))
    fig.text(0.5, 0.045,
             "AlphaFold cartoon coloured by per-residue pLDDT; pocket = top P2Rank site (fpocket+P2Rank consensus)",
             ha="center", color="#7A8890", fontsize=9.5)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.10, wspace=0.02, hspace=0.14)
    out = REPO_ROOT / "output" / "plots" / f"06n_structures_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out), dpi=150)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = L.ORGANISMS[org]

    jobs = pick_targets(org, prefix)
    if not jobs:
        raise SystemExit(f"[{org}] no targets with a cif model + P2Rank pocket found")

    with tempfile.TemporaryDirectory() as td:
        job_path = Path(td) / "job.json"
        img_dir = Path(td) / "img"
        job_path.write_text(json.dumps({"organism": org, "prefix": prefix, "jobs": jobs}))
        conda = os.environ.get("CONDA_EXE", "conda")
        subprocess.run([conda, "run", "-n", PYMOL_ENV, "pymol", "-cq", str(HELPER),
                        "--", str(job_path), str(img_dir)], check=True)
        out = montage(jobs, img_dir, prefix, ORGNAME[org])
    print(f"[{org}] wrote {out.relative_to(REPO_ROOT)} ({len(jobs)} structures: "
          f"{', '.join(j['label'] for j in jobs)})")


if __name__ == "__main__":
    main()
