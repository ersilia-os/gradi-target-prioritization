"""Per-protein experimental PDB sequence coverage for a reference proteome (docs §1.2a).

For every protein we ask PDBe (SIFTS) which experimental PDB structures map onto its UniProt
sequence, and summarise: how much of the sequence they cover (union across all structures, and
the best single structure), the best resolution/method, the covering PDB IDs, and whether any
covering structure is ligand-bound (holo).

Source — PDBe SIFTS "best_structures" (one batched endpoint gives coverage + residue ranges +
resolution + method per structure):
  POST/GET https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/{accession}
Apo/holo — per covering PDB entry:
  GET https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id}
  (404 / no-data == apo; otherwise the bound monomers, filtered against an ignore-list of
  water / ions / cryo-buffer additives).

Organism selected with --organism (kpneumoniae default, or ecoli). Raw PDBe responses are
cached under data/processed/<organism>/pdb/ for resumability. Output, keyed by UniProt accession:
  output/results/<organism>/<prefix>_pdb_coverage.csv
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

REPO_ROOT = Path(__file__).resolve().parents[1]
ORGANISMS = {
    "kpneumoniae": ("UP000007841_HS11286", "kp"),
    "ecoli": ("UP000000625_EcoliK12", "ec"),
}

BEST_STRUCTURES = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/"
LIGAND_MONOMERS = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/"
BATCH = 100  # accessions per best_structures POST

# Apo/holo heuristic: chem-comp IDs that do NOT count as a "real" ligand.
LIGAND_IGNORE = {
    "HOH",
    "DOD",  # water
    # monatomic ions
    "NA",
    "K",
    "LI",
    "CL",
    "BR",
    "IOD",
    "F",
    "MG",
    "CA",
    "ZN",
    "MN",
    "FE",
    "FE2",
    "CU",
    "CU1",
    "NI",
    "CD",
    "CO",
    "HG",
    "BA",
    "SR",
    "CS",
    "RB",
    "AL",
    "GA",
    # common buffer / cryo / crystallisation additives
    "SO4",
    "PO4",
    "PI",
    "NO3",
    "NH4",
    "ACT",
    "FMT",
    "GOL",
    "EDO",
    "PEG",
    "PGE",
    "PG4",
    "1PE",
    "P6G",
    "MPD",
    "DMS",
    "TRS",
    "EPE",
    "MES",
    "BME",
    "IMD",
    "BO3",
    "CAC",
    "ACY",
    "FLC",
    "CIT",
    "TLA",
    "MLI",
    "OXL",
    "SCN",
    "AZI",
    "PER",
    "DTT",
    "BTB",
    "MRD",
    "POL",
    "2PE",
    "12P",
    "15P",
}


def read_lengths(tsv_path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc, seq = row.get("Entry"), row.get("Sequence")
            if acc and seq:
                lengths[acc] = len(seq.strip())
    return lengths


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def post_best_structures(accs: list[str]) -> dict:
    r = requests.post(
        BEST_STRUCTURES,
        data=",".join(accs),
        headers={
            "Content-Type": "application/x-www-form-urlencoded",
            "Connection": "close",
        },
        timeout=(10, 90),
    )
    if r.status_code == 404:
        return {}
    r.raise_for_status()
    return r.json()


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def get_ligands(pdb_id: str) -> list[str]:
    r = requests.get(
        LIGAND_MONOMERS + pdb_id, headers={"Connection": "close"}, timeout=(10, 60)
    )
    if r.status_code == 404:
        return []  # apo: no bound monomers
    r.raise_for_status()
    data = r.json().get(pdb_id, [])
    return sorted({m["chem_comp_id"] for m in data if m.get("chem_comp_id")})


def union_residues(segments: list[tuple[int, int]]) -> int:
    """Total residues covered by a set of inclusive [start, end] intervals (merged)."""
    merged: list[list[int]] = []
    for s, e in sorted(segments):
        if s is None or e is None:
            continue
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return sum(e - s + 1 for s, e in merged)


def fetch_best_structures(accs: list[str], cache_dir: Path, workers: int) -> None:
    """Populate cache_dir/{acc}.json for every accession (JSON list, or 'null' if no PDB)."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accs if not (cache_dir / f"{a}.json").exists()]
    if not todo:
        return
    chunks = [todo[i : i + BATCH] for i in range(0, len(todo), BATCH)]
    print(
        f"  best_structures: {len(todo)} uncached accessions in {len(chunks)} POSTs ..."
    )
    done = 0

    def run(chunk):
        res = post_best_structures(chunk)
        for acc in chunk:
            (cache_dir / f"{acc}.json").write_text(
                json.dumps(res.get(acc))
            )  # list or null
        return len(chunk)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for k in pool.map(run, chunks):
            done += k
            print(f"    {done}/{len(todo)}", flush=True)


def fetch_ligands(
    pdb_ids: list[str], cache_dir: Path, workers: int
) -> dict[str, list[str]]:
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
        print(f"  ligands: {len(todo)} uncached PDB entries ...")

        def run(pid):
            ligs = get_ligands(pid)
            (cache_dir / f"{pid}.json").write_text(json.dumps(ligs))
            return pid, ligs

        with ThreadPoolExecutor(max_workers=workers) as pool:
            for pid, ligs in pool.map(run, todo):
                out[pid] = ligs
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    stem, prefix = ORGANISMS[args.organism]
    tsv = REPO_ROOT / "data" / "raw" / args.organism / "proteome" / f"{stem}.tsv"
    pdb_dir = REPO_ROOT / "data" / "processed" / args.organism / "pdb"
    bs_dir, lig_dir = pdb_dir / "best_structures", pdb_dir / "ligands"
    results_dir = REPO_ROOT / "output" / "results" / args.organism
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / f"{prefix}_pdb_coverage.csv"

    lengths = read_lengths(tsv)
    accs = list(lengths)
    print(f"{args.organism}: {len(accs)} proteins from {tsv.relative_to(REPO_ROOT)}")

    # 1) best_structures for all accessions (cached)
    fetch_best_structures(accs, bs_dir, args.workers)

    # 2) collect covering structures + the set of PDB ids needing ligand lookup
    per_acc: dict[str, list[dict]] = {}
    covering_pdbs: set[str] = set()
    for acc in accs:
        structs = json.loads((bs_dir / f"{acc}.json").read_text()) or []
        per_acc[acc] = structs
        covering_pdbs.update(s["pdb_id"] for s in structs)

    # 3) ligand info for covering PDB ids (cached) -> real (non-ignored) ligands per entry
    raw_ligs = fetch_ligands(sorted(covering_pdbs), lig_dir, args.workers)
    real_ligs = {
        pid: [c for c in ligs if c.upper() not in LIGAND_IGNORE]
        for pid, ligs in raw_ligs.items()
    }

    # 4) assemble per-protein rows
    rows = []
    for acc in accs:
        L = lengths[acc]
        structs = per_acc[acc]
        row = {
            "uniprot_accession": acc,
            "pdb_has_structure": bool(structs),
            "pdb_n_structures": 0,
            "pdb_n_chains": 0,
            "pdb_ids": "",
            "pdb_best_id": "",
            "pdb_best_chain": "",
            "pdb_coverage_union": 0.0,
            "pdb_coverage_best_single": 0.0,
            "pdb_best_resolution_A": None,
            "pdb_best_method": "",
            "pdb_has_holo": False,
            "pdb_holo_ids": "",
            "pdb_ligands": "",
            "pdb_apo_holo_per_structure": "",
            "source": "pdbe_sifts",
        }
        if structs:
            # group rows (one per chain mapping) by pdb entry
            by_pdb: dict[str, list[dict]] = {}
            for s in structs:
                by_pdb.setdefault(s["pdb_id"], []).append(s)
            # union coverage across ALL structures/chains
            all_segs = [(s["unp_start"], s["unp_end"]) for s in structs]
            row["pdb_coverage_union"] = round(min(union_residues(all_segs) / L, 1.0), 4)
            # best single structure = max per-entry union coverage (tie-break: lower resolution)
            entry_cov = {
                pid: union_residues([(s["unp_start"], s["unp_end"]) for s in rs]) / L
                for pid, rs in by_pdb.items()
            }

            def entry_res(pid):
                vals = [
                    s["resolution"]
                    for s in by_pdb[pid]
                    if s.get("resolution") is not None
                ]
                return min(vals) if vals else float("inf")

            best_pid = max(
                by_pdb, key=lambda p: (round(entry_cov[p], 6), -entry_res(p))
            )
            row["pdb_coverage_best_single"] = round(min(entry_cov[best_pid], 1.0), 4)
            row["pdb_best_id"] = best_pid
            row["pdb_best_chain"] = by_pdb[best_pid][0].get("chain_id", "")
            row["pdb_n_structures"] = len(by_pdb)
            row["pdb_n_chains"] = len(
                {(s["pdb_id"], s.get("chain_id")) for s in structs}
            )
            row["pdb_ids"] = ";".join(sorted(by_pdb))
            res_vals = [
                s["resolution"] for s in structs if s.get("resolution") is not None
            ]
            if res_vals:
                best = min(
                    structs,
                    key=lambda s: (
                        s["resolution"]
                        if s.get("resolution") is not None
                        else float("inf")
                    ),
                )
                row["pdb_best_resolution_A"] = best["resolution"]
                row["pdb_best_method"] = best.get("experimental_method", "")
            else:
                row["pdb_best_method"] = structs[0].get("experimental_method", "")
            # apo/holo
            holo_ids, all_real, per_struct = [], set(), []
            for pid in sorted(by_pdb):
                ligs = real_ligs.get(pid, [])
                if ligs:
                    holo_ids.append(pid)
                    all_real.update(ligs)
                    per_struct.append(f"{pid}:{('|'.join(ligs))}")
                else:
                    per_struct.append(f"{pid}:apo")
            row["pdb_has_holo"] = bool(holo_ids)
            row["pdb_holo_ids"] = ";".join(holo_ids)
            row["pdb_ligands"] = ";".join(sorted(all_real))
            row["pdb_apo_holo_per_structure"] = ";".join(per_struct)
        rows.append(row)

    df = pd.DataFrame(rows).sort_values("uniprot_accession").reset_index(drop=True)

    # ---- verification ----
    assert len(df) == len(accs), f"expected {len(accs)} rows, got {len(df)}"
    assert df["uniprot_accession"].is_unique, "duplicate accessions"
    assert df["pdb_coverage_union"].between(0, 1).all(), "union coverage out of [0,1]"
    bad = df[df["pdb_coverage_union"] + 1e-9 < df["pdb_coverage_best_single"]]
    assert bad.empty, f"union < best_single for {len(bad)} proteins"
    n_pdb = int(df["pdb_has_structure"].sum())
    n_holo = int(df["pdb_has_holo"].sum())
    print(
        f"With PDB: {n_pdb}/{len(df)} ({100 * n_pdb / len(df):.1f}%); holo: {n_holo}; "
        f"covering PDB entries: {len(covering_pdbs)}"
    )
    cov = df[df.pdb_has_structure]
    if n_pdb:
        print(
            f"  median union coverage (PDB subset): {cov.pdb_coverage_union.median():.2f}; "
            f"≥80% covered: {int((cov.pdb_coverage_union >= 0.8).sum())}"
        )

    df.to_csv(csv_path, index=False)
    print(f"Wrote {csv_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
