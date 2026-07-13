"""AlphaFill ligand-transfer evidence (docs §2.2b).

AlphaFill (alphafill.eu) transplants ligands and cofactors from homologous PDB entries onto the
AlphaFold model. For every protein we query the AlphaFill JSON metadata (keyed by UniProt
accession, same AFDB id as our cached models) and summarise the transplanted ligands as a
structure-based binding signal — *something* binds a homolog at a structurally equivalent site.

API:  GET https://alphafill.eu/v1/aff/<accession>/json  -> {"hits":[{alignment{identity},
      global_rmsd, pdb_id, transplants:[{analogue_id, clash{clash_count}, ...}]}]}
Responses cached under data/processed/<org>/alphafill/<acc>.json (resumable). Ligand curation is
two-tier (src/ligandability.py): `alphafill_n_ligand_any` / `alphafill_ligand_ids` are the BROAD set
(any real transplanted molecule), while `alphafill_n_druglike` is the STRICT drug-like tier that
drops promiscuous cofactors/nucleotides (ATP, NAD, FAD, heme, …). The strict tier is the headline.

Output (keyed by UniProt accession): output/results/<org>/<prefix>_alphafill.csv
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

API = "https://alphafill.eu/v1/aff/{}/json"


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=1, min=1, max=15))
def fetch_alphafill(acc: str) -> dict | None:
    r = requests.get(API.format(acc), headers={"Connection": "close"}, timeout=(10, 60))
    if r.status_code in (404, 400):
        return None  # no AlphaFill entry
    r.raise_for_status()
    try:
        return r.json()
    except ValueError:
        return None


def ensure_cache(accs: list[str], cache_dir: Path, workers: int) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    todo = [a for a in accs if not (cache_dir / f"{a}.json").exists()]
    if not todo:
        return
    print(f"  alphafill: {len(todo)} uncached accessions ...", flush=True)
    done = 0

    def run(acc):
        try:
            data = fetch_alphafill(acc)
        except Exception as e:  # noqa: BLE001
            data = {"_error": str(e)[:200]}
        (cache_dir / f"{acc}.json").write_text(json.dumps(data))
        return acc

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for _ in pool.map(run, todo):
            done += 1
            if done % 200 == 0:
                print(f"    {done}/{len(todo)}", flush=True)


def summarise(acc: str, data: dict | None) -> dict:
    row = {
        "uniprot_accession": acc,
        "alphafill_available": False,
        "alphafill_n_hits": 0,
        "alphafill_n_transplants": 0,
        "alphafill_n_druglike": 0,
        "alphafill_n_ligand_any": 0,
        "alphafill_ligand_ids": "",
        "alphafill_best_identity": None,
        "alphafill_best_global_rmsd": None,
        "alphafill_best_ligand": "",
    }
    if not data or "hits" not in data or not data["hits"]:
        return row
    hits = data["hits"]
    row["alphafill_available"] = True
    row["alphafill_n_hits"] = len(hits)
    broad: list[str] = []     # any real transplanted molecule (cofactors incl.)
    druglike: list[str] = []  # strict drug-like (cofactors/nucleotides excluded)
    best = None  # (key, ligand, identity, rmsd) — best among drug-like transplants
    for h in hits:
        ident = (h.get("alignment") or {}).get("identity")
        grmsd = h.get("global_rmsd")
        for t in h.get("transplants", []) or []:
            lig = (t.get("analogue_id") or "").strip().upper()
            if not lig:
                continue
            row["alphafill_n_transplants"] += 1
            if L.is_ligand(lig):
                broad.append(lig)
            if L.is_druglike_ligand(lig):
                druglike.append(lig)
                key = (
                    ident if ident is not None else 0.0,
                    -(grmsd if grmsd is not None else 99.0),
                )
                if best is None or key > best[0]:
                    best = (key, lig, ident, grmsd)
    row["alphafill_n_druglike"] = len(set(druglike))
    row["alphafill_n_ligand_any"] = len(set(broad))
    row["alphafill_ligand_ids"] = ";".join(sorted(set(broad)))  # broad (incl. cofactors)
    if best is not None:
        _, lig, ident, grmsd = best
        row["alphafill_best_ligand"] = lig
        row["alphafill_best_identity"] = round(ident, 4) if ident is not None else None
        row["alphafill_best_global_rmsd"] = round(grmsd, 3) if grmsd is not None else None
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    accs = L.load_accessions(org)
    if args.limit:
        accs = accs[: args.limit]

    cache_dir = L.processed_dir(org, "alphafill")
    ensure_cache(accs, cache_dir, args.workers)

    rows = []
    for acc in accs:
        p = cache_dir / f"{acc}.json"
        data = json.loads(p.read_text()) if p.exists() else None
        if isinstance(data, dict) and "_error" in data:
            data = None
        rows.append(summarise(acc, data))

    df = pd.DataFrame(rows)
    out = L.results_dir(org) / f"{prefix}_alphafill.csv"
    df.to_csv(out, index=False)
    n_avail = int(df["alphafill_available"].sum())
    n_drug = int((df["alphafill_n_druglike"] > 0).sum())
    print(
        f"[{org}] wrote {out} ({len(df)} proteins; {n_avail} with AlphaFill, "
        f"{n_drug} with a drug-like transplant)",
        flush=True,
    )


if __name__ == "__main__":
    main()
