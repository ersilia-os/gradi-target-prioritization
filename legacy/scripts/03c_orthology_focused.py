"""Focused 3-proteome orthology analysis: K. pneumoniae x E. coli x human (TABLES ONLY).

A self-contained, target-prioritization view of the project's three focal reference proteomes
(fetched by 00a). A good BacPROTAC target is BROAD-SPECTRUM (shared by both bacteria) and
SELECTIVE (no human ortholog), so we run OrthoFinder over exactly the three proteomes and read
the orthogroup membership structure through that lens. This script writes TABLES only; plotting
is done downstream from these tables.

Method: OrthoFinder (DIAMOND) on the 3 proteomes -> orthogroups; DIAMOND reciprocal-best-hits per
pair for %identity ("orthology depth"). Runs in the `gradi` env and shells out to the
`gradi-ortho` binaries (orthofinder/diamond) via a PATH-injected subprocess.

Outputs (TSV, keyed by UniProt accession), under data/processed/other/orthology/:
  three_way_orthogroups.tsv          orthogroup, n_kp, n_ec, n_human, region (the Venn data)
  three_way_protein_categories.tsv   per kp/ec protein: orthogroup, has_other_bacterium, has_human, selectivity
  broad_spectrum_selective_shortlist.tsv   the broad_selective subset (kp & ec, no human ortholog)
  three_way_rbh_identity.tsv         pair, query_accession, target_accession, pident (orthology depth)

    conda activate gradi
    python scripts/03c_orthology_focused.py
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PROC_DIR = REPO_ROOT / "data" / "processed" / "other" / "orthology"
RAW_DIR = REPO_ROOT / "data" / "raw" / "other" / "orthology"
OF_INPUT = PROC_DIR / "of_input_3way"
DEFAULT_ORTHO_BIN = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

# label -> (organism folder, proteome file stem, short key)
PROTEOMES = {
    "HS11286": ("kpneumoniae", "UP000007841_HS11286", "kp"),
    "EcoliK12": ("ecoli", "UP000000625_EcoliK12", "ec"),
    "Human": ("human", "UP000005640_Human", "hs"),
}


def _run(cmd: list[str], ortho_bin: Path) -> None:
    env = os.environ.copy()
    env["PATH"] = f"{ortho_bin}{os.pathsep}{env.get('PATH', '')}"
    subprocess.run(cmd, check=True, env=env)


def bare_accession_fasta(label: str) -> Path:
    """Copy a proteome FASTA into RAW_DIR with headers rewritten to the bare UniProt accession."""
    organism, stem, _ = PROTEOMES[label]
    src = REPO_ROOT / "data" / "raw" / organism / "proteome" / f"{stem}.fasta"
    out = RAW_DIR / f"{label}.fasta"
    if out.exists() and out.stat().st_size > 0:
        return out
    lines = []
    for line in src.read_text().splitlines():
        if line.startswith(">"):
            lines.append(">" + line[1:].split("|")[1])
        else:
            lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return out


def proteome_genes(label: str) -> pd.DataFrame:
    organism, stem, _ = PROTEOMES[label]
    tsv = REPO_ROOT / "data" / "raw" / organism / "proteome" / f"{stem}.tsv"
    return pd.read_csv(tsv, sep="\t")[["Entry", "Gene Names"]].rename(
        columns={"Entry": "uniprot_accession", "Gene Names": "gene"}
    )


# --------------------------- OrthoFinder ---------------------------
def run_orthofinder(threads: int, force: bool, ortho_bin: Path) -> Path:
    if OF_INPUT.exists() and force:
        shutil.rmtree(OF_INPUT)
    OF_INPUT.mkdir(parents=True, exist_ok=True)
    for label in PROTEOMES:
        shutil.copy(bare_accession_fasta(label), OF_INPUT / f"{label}.fasta")
    existing = sorted(OF_INPUT.glob("OrthoFinder/Results_*"))
    if existing and not force:
        print(f"  reusing OrthoFinder results {existing[-1].name}")
        return existing[-1]
    if not (ortho_bin / "orthofinder").exists():
        sys.exit(
            f"orthofinder not found in {ortho_bin} (use --ortho-bin / see install.sh)."
        )
    print("  running OrthoFinder (DIAMOND) on the 3 proteomes ...")
    _run(["orthofinder", "-f", str(OF_INPUT), "-t", str(threads)], ortho_bin)
    return sorted(OF_INPUT.glob("OrthoFinder/Results_*"))[-1]


def orthogroup_membership(results_dir: Path) -> pd.DataFrame:
    """One row per orthogroup with per-species gene lists + counts + Venn region.

    Uses Orthogroups.tsv (assigned) + Orthogroups_UnassignedGenes.tsv (species-specific singletons)
    so every protein is represented exactly once.
    """
    og_dir = results_dir / "Orthogroups"
    frames = [pd.read_csv(og_dir / "Orthogroups.tsv", sep="\t", dtype=str)]
    unassigned = og_dir / "Orthogroups_UnassignedGenes.tsv"
    if unassigned.exists():
        frames.append(pd.read_csv(unassigned, sep="\t", dtype=str))
    og = pd.concat(frames, ignore_index=True).fillna("")

    def genes(cell: str) -> list[str]:
        return [g.strip() for g in cell.split(",") if g.strip()]

    rows = []
    for _, r in og.iterrows():
        per = {lab: genes(r.get(lab, "")) for lab in PROTEOMES}
        counts = {lab: len(per[lab]) for lab in PROTEOMES}
        present = [
            tag
            for lab, tag in (("HS11286", "kp"), ("EcoliK12", "ec"), ("Human", "hs"))
            if counts[lab] > 0
        ]
        rows.append(
            {
                "orthogroup": r["Orthogroup"],
                "n_kp": counts["HS11286"],
                "n_ec": counts["EcoliK12"],
                "n_human": counts["Human"],
                "region": "+".join(present),
                "_kp": per["HS11286"],
                "_ec": per["EcoliK12"],
            }
        )
    return pd.DataFrame(rows)


def protein_categories(og: pd.DataFrame) -> pd.DataFrame:
    """Per bacterial protein: conservation flags + target-selectivity category."""
    rows = []
    for _, r in og.iterrows():
        in_kp, in_ec, in_hs = r.n_kp > 0, r.n_ec > 0, r.n_human > 0
        for organism, acc_list, other_present in (
            ("kpneumoniae", r["_kp"], in_ec),
            ("ecoli", r["_ec"], in_kp),
        ):
            cat = ("broad" if other_present else "narrow") + (
                "_human_homolog" if in_hs else "_selective"
            )
            for acc in acc_list:
                rows.append(
                    {
                        "organism": organism,
                        "uniprot_accession": acc,
                        "orthogroup": r.orthogroup,
                        "has_other_bacterium": other_present,
                        "has_human": in_hs,
                        "selectivity": cat,
                    }
                )
    df = pd.DataFrame(rows)
    genes = pd.concat(
        [proteome_genes("HS11286"), proteome_genes("EcoliK12")], ignore_index=True
    )
    return df.merge(genes, on="uniprot_accession", how="left")


# --------------------------- DIAMOND RBH (%identity) ---------------------------
def diamond_rbh(
    a: str, b: str, threads: int, evalue: float, ortho_bin: Path
) -> pd.DataFrame:
    """Reciprocal-best-hit %identity between two proteomes -> (query, target, pident) rows."""
    fa, fb = bare_accession_fasta(a), bare_accession_fasta(b)

    def best(q: Path, db_fa: Path, tag: str) -> pd.DataFrame:
        db = PROC_DIR / f"_db_{tag}.dmnd"
        out = PROC_DIR / f"_hits_{tag}.tsv"
        _run(
            ["diamond", "makedb", "--in", str(db_fa), "-d", str(db.with_suffix(""))],
            ortho_bin,
        )
        _run(
            [
                "diamond",
                "blastp",
                "-q",
                str(q),
                "-d",
                str(db),
                "-o",
                str(out),
                "-p",
                str(threads),
                "--evalue",
                str(evalue),
                "--max-target-seqs",
                "1",
                "--outfmt",
                "6",
                "qseqid",
                "sseqid",
                "pident",
                "bitscore",
            ],
            ortho_bin,
        )
        d = pd.read_csv(out, sep="\t", names=["q", "s", "pident", "bitscore"])
        return d.sort_values("bitscore", ascending=False).drop_duplicates("q")

    ab = best(fa, fb, f"{a}_{b}")  # a -> best in b
    ba = best(fb, fa, f"{b}_{a}")  # b -> best in a
    rev = dict(zip(ba.q, ba.s))
    keep = ab[[rev.get(s) == q for q, s in zip(ab.q, ab.s)]]  # reciprocal pairs
    return keep[["q", "s", "pident"]].rename(
        columns={"q": "query_accession", "s": "target_accession"}
    )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--evalue", type=float, default=1e-4)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--ortho-bin", type=Path, default=DEFAULT_ORTHO_BIN)
    args = ap.parse_args()

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROC_DIR.mkdir(parents=True, exist_ok=True)

    results_dir = run_orthofinder(args.threads, args.force, args.ortho_bin)
    og = orthogroup_membership(results_dir)
    cats = protein_categories(og)

    og.drop(columns=["_kp", "_ec"]).to_csv(
        PROC_DIR / "three_way_orthogroups.tsv", sep="\t", index=False
    )
    cats.to_csv(PROC_DIR / "three_way_protein_categories.tsv", sep="\t", index=False)
    shortlist = cats[cats.selectivity == "broad_selective"]
    shortlist.to_csv(
        PROC_DIR / "broad_spectrum_selective_shortlist.tsv", sep="\t", index=False
    )

    print(
        f"orthogroups: {len(og)} | region counts: {og['region'].value_counts().to_dict()}"
    )
    print(f"protein categories: {cats['selectivity'].value_counts().to_dict()}")
    print(f"broad-spectrum + selective shortlist: {len(shortlist)} proteins")

    print("DIAMOND RBH %identity per pair ...")
    pairs = {
        "kp-ec": ("HS11286", "EcoliK12"),
        "kp-human": ("HS11286", "Human"),
        "ec-human": ("EcoliK12", "Human"),
    }
    frames = []
    for name, (a, b) in pairs.items():
        d = diamond_rbh(a, b, args.threads, args.evalue, args.ortho_bin)
        d.insert(0, "pair", name)
        frames.append(d)
        print(
            f"  {name}: n={len(d)} median %id={d.pident.median():.1f}"
            if len(d)
            else f"  {name}: none"
        )
    rbh = pd.concat(frames, ignore_index=True)
    rbh.to_csv(PROC_DIR / "three_way_rbh_identity.tsv", sep="\t", index=False)

    print(f"Wrote 4 tables to {PROC_DIR.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
