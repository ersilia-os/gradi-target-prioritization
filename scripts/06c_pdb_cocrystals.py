"""PDB co-crystal ligand evidence, direct + ortholog-expanded (docs §2.2a).

Empirical proof that a small molecule binds the protein (or a homolog): a protein–ligand
co-crystal in the PDB. Two buckets:
  * DIRECT   — the anchor protein's own PDB structures. Read from the 04c coverage table
               (output/results/<org>/<prefix>_pdb_coverage.csv: pdb_ids / pdb_ligands), then
               filtered to drug-like ligands.
  * ORTHOLOG — for every ortholog (03a long table) we ask PDBe SIFTS which experimental PDB
               structures map to it, then which of those carry a drug-like bound ligand.
Ligand curation is two-tier (src/ligandability.py): `*_ligand_ids` / `*_n_ligand_any` are the BROAD
set (any real bound molecule, water/ions/buffer excluded), while `*_has_druglike` / `*_n_druglike`
are the STRICT drug-like tier that additionally drops promiscuous cofactors/nucleotides
(ATP, NAD, FAD, heme, Fe-S clusters, …). The strict tier is the headline signal.

PDBe endpoints (same as 04c):
  POST https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/{accs}
  GET  https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id}
Ortholog responses cached under data/processed/other/pdb_orthologs/ (shared across organisms).

Output (keyed by UniProt accession): output/results/<org>/<prefix>_pdb_cocrystals.csv
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

BEST_STRUCTURES = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/"
LIGAND_MONOMERS = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/"
BATCH = 100

PDB_SEQRES = L.REPO_ROOT / "data" / "raw" / "other" / "pdb" / "pdb_seqres.txt"
PDB_PROT_FAA = L.REPO_ROOT / "data" / "processed" / "other" / "pdb_orthologs" / "pdb_seqres_protein.faa"
SEQ_DIRECT_PIDENT = 95.0  # >= this %identity to a PDB chain == effectively the same protein


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def post_best_structures(accs: list[str]) -> dict:
    r = requests.post(
        BEST_STRUCTURES,
        data=",".join(accs),
        headers={"Content-Type": "application/x-www-form-urlencoded", "Connection": "close"},
        timeout=(10, 90),
    )
    if r.status_code == 404:
        return {}
    r.raise_for_status()
    return r.json()


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def get_ligands(pdb_id: str) -> list[str]:
    r = requests.get(LIGAND_MONOMERS + pdb_id, headers={"Connection": "close"}, timeout=(10, 60))
    if r.status_code == 404:
        return []
    r.raise_for_status()
    data = r.json().get(pdb_id, [])
    return sorted({m["chem_comp_id"] for m in data if m.get("chem_comp_id")})


def cache_best_structures(accs: list[str], cache_dir: Path, workers: int) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accs if not (cache_dir / f"{a}.json").exists()]
    if not todo:
        return
    chunks = [todo[i : i + BATCH] for i in range(0, len(todo), BATCH)]
    print(f"  best_structures: {len(todo)} accs in {len(chunks)} POSTs ...", flush=True)
    done = 0

    def run(chunk):
        res = post_best_structures(chunk)
        for acc in chunk:
            (cache_dir / f"{acc}.json").write_text(json.dumps(res.get(acc)))
        return len(chunk)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for k in pool.map(run, chunks):
            done += k
            if done % 2000 == 0:
                print(f"    {done}/{len(todo)}", flush=True)


def cache_ligands(pdb_ids: list[str], cache_dir: Path, workers: int) -> dict[str, list[str]]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    out: dict[str, list[str]] = {}
    todo = []
    for pid in pdb_ids:
        p = cache_dir / f"{pid}.json"
        if p.exists():
            out[pid] = json.loads(p.read_text())
        else:
            todo.append(pid)
    if todo:
        print(f"  ligands: {len(todo)} uncached PDB entries ...", flush=True)

        def run(pid):
            ligs = get_ligands(pid)
            (cache_dir / f"{pid}.json").write_text(json.dumps(ligs))
            return pid, ligs

        with ThreadPoolExecutor(max_workers=workers) as pool:
            for pid, ligs in pool.map(run, todo):
                out[pid] = ligs
    return out


def pdb_ids_for(acc: str, cache_dir: Path) -> list[str]:
    p = cache_dir / f"{acc}.json"
    if not p.exists():
        return []
    data = json.loads(p.read_text())
    if not data:
        return []
    return sorted({d["pdb_id"] for d in data if d.get("pdb_id")})


def direct_evidence(org: str, prefix: str) -> dict[str, dict]:
    """From the 04c coverage table: drug-like ligands on the anchor's own structures."""
    path = L.results_dir(org) / f"{prefix}_pdb_coverage.csv"
    out: dict[str, dict] = {}
    if not path.exists():
        print(f"  [warn] {path} missing; direct PDB evidence empty", flush=True)
        return out
    df = pd.read_csv(path)
    for _, r in df.iterrows():
        acc = r["uniprot_accession"]
        ligs_raw = str(r.get("pdb_ligands") or "")
        ligs = L.ligands(  # broad set (cofactors retained; strict tier derived downstream)
            [x for x in ligs_raw.replace(",", ";").split(";") if x and x != "nan"]
        )
        out[acc] = {
            "has_structure": bool(r.get("pdb_has_structure")),
            "ligs": ligs,
            "pdb_ids": str(r.get("pdb_ids") or "").replace(",", ";"),
        }
    return out


def build_pdb_protein_fasta() -> None:
    """One-time: extract protein chains from pdb_seqres into a DIAMOND-ready FASTA (id=<pdb>_<chain>)."""
    if PDB_PROT_FAA.exists() and PDB_PROT_FAA.stat().st_size > 0:
        return
    if not PDB_SEQRES.exists():
        raise FileNotFoundError(
            f"{PDB_SEQRES} missing — download pdb_seqres.txt.gz from RCSB derived_data."
        )
    PDB_PROT_FAA.parent.mkdir(parents=True, exist_ok=True)
    print("  building protein-only pdb_seqres FASTA ...", flush=True)
    keep = False
    n = 0
    with open(PDB_SEQRES) as fin, open(PDB_PROT_FAA, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                keep = "mol:protein" in line
                if keep:
                    fout.write(">" + line[1:].split()[0] + "\n")
                    n += 1
            elif keep:
                fout.write(line)
    print(f"  wrote {PDB_PROT_FAA.name} ({n} protein chains)", flush=True)


def seq_pdb_evidence(org: str, lig_dir: Path, workers: int) -> dict[str, dict]:
    """DIAMOND the proteome vs pdb_seqres; per protein, drug-like ligands on matched PDB entries,
    split into direct (>=95% id, same protein) vs homolog. Catches Kp structures that SIFTS keys
    under other-strain accessions."""
    build_pdb_protein_fasta()
    hits_tsv = L.processed_dir(org, "pdb") / "pdb_seqres_diamond.tsv"
    L.run_diamond_blastp(L.proteome_fasta(org), PDB_PROT_FAA, hits_tsv,
                         threads=8, min_id=30, query_cover=50, max_target_seqs=20)
    hits = L.load_diamond_hits(hits_tsv)
    if hits.empty:
        return {}
    hits["q"] = hits["qseqid"].map(L.acc_from_header)
    hits["pdbid"] = hits["sseqid"].str.split("_").str[0].str.lower()

    uniq_pdbs = sorted(set(hits["pdbid"]))
    print(f"  [seq] {len(hits)} chain hits across {len(uniq_pdbs)} PDB entries; fetching ligands ...", flush=True)
    lig_map = cache_ligands(uniq_pdbs, lig_dir, workers)
    pdb_ligands = {pid: L.ligands(ligs) for pid, ligs in lig_map.items()}  # broad

    out: dict[str, dict] = {}
    for q, sub in hits.groupby("q"):
        d_ligs, d_pdbs, d_pident = set(), set(), 0.0
        h_ligs, h_pdbs, h_pident = set(), set(), 0.0
        for r in sub.itertuples():
            dl = pdb_ligands.get(r.pdbid, [])
            if not dl:
                continue
            if r.pident >= SEQ_DIRECT_PIDENT:
                d_ligs.update(dl); d_pdbs.add(r.pdbid); d_pident = max(d_pident, r.pident)
            else:
                h_ligs.update(dl); h_pdbs.add(r.pdbid); h_pident = max(h_pident, r.pident)
        out[q] = {
            "seq_direct_ligs": sorted(d_ligs), "seq_direct_pdbs": sorted(d_pdbs), "seq_direct_pident": d_pident,
            "seq_hom_ligs": sorted(h_ligs), "seq_hom_pdbs": sorted(h_pdbs), "seq_hom_pident": h_pident,
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--skip-seq", action="store_true", help="skip the pdb_seqres sequence route")
    args = ap.parse_args()

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    accs = L.load_accessions(org)
    if args.limit:
        accs = accs[: args.limit]
    acc_set = set(accs)

    direct = direct_evidence(org, prefix)

    # ortholog accessions for the in-scope anchors
    orth = L.load_orthologs(org)
    orth = orth[orth["anchor_uniprot"].isin(acc_set)].copy()
    orth_accs = sorted({str(x) for x in orth["target_uniprot"] if str(x) and str(x) != "nan"})

    best_dir = L.REPO_ROOT / "data" / "processed" / "other" / "pdb_orthologs" / "best_structures"
    lig_dir = L.REPO_ROOT / "data" / "processed" / "other" / "pdb_orthologs" / "ligands"
    print(f"[{org}] {len(orth_accs)} ortholog accessions to resolve in PDBe", flush=True)
    cache_best_structures(orth_accs, best_dir, args.workers)

    # gather all PDB ids that appear for any ortholog, fetch their ligands
    acc_pdbs = {a: pdb_ids_for(a, best_dir) for a in orth_accs}
    all_pdbs = sorted({p for v in acc_pdbs.values() for p in v})
    print(f"[{org}] {len(all_pdbs)} unique ortholog PDB entries with structures", flush=True)
    lig_map = cache_ligands(all_pdbs, lig_dir, args.workers)
    pdb_ligands = {pid: L.ligands(ligs) for pid, ligs in lig_map.items()}  # broad

    # sequence route (catches Kp structures keyed under other-strain accessions in SIFTS)
    seq_ev = {} if args.skip_seq else seq_pdb_evidence(org, lig_dir, args.workers)

    # per-anchor aggregation
    orth_by_anchor = {a: sub for a, sub in orth.groupby("anchor_uniprot")}

    def _strict(ligs) -> list[str]:
        """Drug-like subset of a broad ligand list (drops cofactors/nucleotides)."""
        return [x for x in ligs if L.is_druglike_ligand(x)]

    rows = []
    for acc in accs:
        d = direct.get(acc, {"has_structure": False, "ligs": [], "pdb_ids": ""})
        d_broad = list(d["ligs"]); d_strict = _strict(d_broad)
        row = {
            "uniprot_accession": acc,
            "pdb_lig_direct_has_structure": d["has_structure"],
            "pdb_lig_direct_has_druglike": len(d_strict) > 0,
            "pdb_lig_direct_n_druglike": len(d_strict),
            "pdb_lig_direct_n_ligand_any": len(d_broad),
            "pdb_lig_direct_ligand_ids": ";".join(d_broad),  # broad (incl. cofactors)
            "pdb_lig_direct_pdb_ids": d["pdb_ids"],
        }
        # ortholog
        o_ligs: set[str] = set()
        o_pdbs: set[str] = set()
        best = None  # (pident, species)
        sub = orth_by_anchor.get(acc)
        if sub is not None:
            for _, r in sub.iterrows():
                tacc = str(r["target_uniprot"])
                hit_ligs: set[str] = set()
                for pid in acc_pdbs.get(tacc, []):
                    dl = pdb_ligands.get(pid, [])
                    if dl:
                        hit_ligs.update(dl)
                        o_pdbs.add(pid)
                if hit_ligs:
                    o_ligs.update(hit_ligs)
                    try:
                        pid_val = float(r.get("pident"))
                    except (TypeError, ValueError):
                        pid_val = 0.0
                    if best is None or pid_val > best[0]:
                        best = (pid_val, r.get("species"))
        o_broad = sorted(o_ligs); o_strict = _strict(o_broad)
        row["pdb_lig_ortho_has_druglike"] = len(o_strict) > 0
        row["pdb_lig_ortho_n_druglike"] = len(o_strict)
        row["pdb_lig_ortho_n_ligand_any"] = len(o_broad)
        row["pdb_lig_ortho_ligand_ids"] = ";".join(o_broad)  # broad (incl. cofactors)
        row["pdb_lig_ortho_pdb_ids"] = ";".join(sorted(o_pdbs))
        row["pdb_lig_ortho_best_species"] = best[1] if best else ""
        row["pdb_lig_ortho_best_pident"] = round(best[0], 1) if best else None

        # sequence route (se.*_ligs are broad sets from pdb_ligands)
        se = seq_ev.get(acc, {})
        sd_broad = list(se.get("seq_direct_ligs", [])); sd_strict = _strict(sd_broad)
        sh_broad = list(se.get("seq_hom_ligs", [])); sh_strict = _strict(sh_broad)
        row["pdb_lig_seqdirect_has_druglike"] = len(sd_strict) > 0
        row["pdb_lig_seqdirect_n_druglike"] = len(sd_strict)
        row["pdb_lig_seqdirect_n_ligand_any"] = len(sd_broad)
        row["pdb_lig_seqdirect_ligand_ids"] = ";".join(sd_broad)  # broad
        row["pdb_lig_seqdirect_pdb_ids"] = ";".join(se.get("seq_direct_pdbs", []))
        row["pdb_lig_seqdirect_best_pident"] = round(se.get("seq_direct_pident", 0.0), 1) if sd_broad else None
        row["pdb_lig_seqhom_has_druglike"] = len(sh_strict) > 0
        row["pdb_lig_seqhom_n_druglike"] = len(sh_strict)
        row["pdb_lig_seqhom_n_ligand_any"] = len(sh_broad)
        row["pdb_lig_seqhom_ligand_ids"] = ";".join(sh_broad)  # broad

        # strict (drug-like) aggregates — the headline signal
        row["pdb_lig_any_has_druglike"] = bool(
            row["pdb_lig_direct_has_druglike"]
            or row["pdb_lig_ortho_has_druglike"]
            or row["pdb_lig_seqdirect_has_druglike"]
            or row["pdb_lig_seqhom_has_druglike"]
        )
        # "direct" now means: own structure OR a >=95%-identity structure (same protein, any strain)
        row["pdb_lig_anydirect_has_druglike"] = bool(
            row["pdb_lig_direct_has_druglike"] or row["pdb_lig_seqdirect_has_druglike"]
        )
        # broad aggregate — any bound ligand incl. cofactors (kept for the two-tier comparison)
        row["pdb_lig_any_has_ligand"] = bool(
            d_broad or o_broad or sd_broad or sh_broad
        )
        rows.append(row)

    df = pd.DataFrame(rows)
    out = L.results_dir(org) / f"{prefix}_pdb_cocrystals.csv"
    df.to_csv(out, index=False)
    print(
        f"[{org}] wrote {out} ({len(df)} proteins; "
        f"own-struct druglike={int(df['pdb_lig_direct_has_druglike'].sum())}, "
        f"seq-direct(≥95%)={int(df['pdb_lig_seqdirect_has_druglike'].sum())}, "
        f"ortholog-acc={int(df['pdb_lig_ortho_has_druglike'].sum())}, "
        f"seq-homolog={int(df['pdb_lig_seqhom_has_druglike'].sum())}, "
        f"any={int(df['pdb_lig_any_has_druglike'].sum())})",
        flush=True,
    )


if __name__ == "__main__":
    main()
