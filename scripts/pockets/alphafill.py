"""AlphaFill transplants: ligands placed on our AlphaFold models by structural superposition.

    data/source/alphafill/<species>/<acc>.json            cached API responses
    data/processed/pockets/evidence/transplants_<sp>.tsv  LONG: one row per transplant, UNFILTERED

**This script does not decide what a ligand is.** It fetches and flattens; `holo.py` classifies and
counts. Keeping them apart avoids a circular dependency -- the drug-likeness vocabulary has to
cover the codes AlphaFill uses, so the raw transplant table must exist first.

WHY ALPHAFILL AND NOT OUR OWN HOMOLOGY TRANSFER
-------------------------------------------------
An earlier version of this axis transferred ligands itself: DIAMOND against BioLiP's holo chains,
a binding-site span test, and the best hit's % identity as the score. That is a hand-rolled
AlphaFill, and worse at the job -- AlphaFill (Hekkelman et al., *Nat Methods* 2023) superposes
structures rather than sequences and validates each transplant geometrically. Measured on the same
proteomes with the same drug-likeness rule: **AlphaFill 1,503 Kp / 1,162 Ec proteins, the
hand-rolled route 184 / 172.** The repo's own precedent is the same -- studiedness had a
hand-rolled k-NN transfer removed "on instruction to use well-established tools only".

The earlier `alphafill_check.py` concluded AlphaFill "adds little" and was deleted with this
rewrite: it imposed a 40% identity floor **AlphaFill does not use**, which threw away most of
what AlphaFill finds. Its donors sit at ~30% median identity by design.

WHAT ALPHAFILL ALREADY DECIDED, SO WE DO NOT
----------------------------------------------
Donors are PDB structures at **>= 25% identity over >= 85 aligned residues**; the ligand is
superposed using backbone atoms within 6 A and the fit is reported as `local_rmsd`. No RMSD
threshold is published -- only the metric -- so **nothing is filtered here** (owner's call,
2026-10-03): every transplant ships with its `local_rmsd`, donor identity and clash count, and
anyone who wants a stricter set can filter the evidence table.

**`analogue_id`, not `compound_id`, is the ligand.** When the donor holds a close analogue,
AlphaFill records the donor's compound as `compound_id` and the compound it stands in for as
`analogue_id` -- ANP->ATP, ACO->COA, AGS->ATP, TGG->GSH (2.7% of transplants). The analogue is
what AlphaFill asserts the protein binds, so it is the identity classified and counted; the
donor's own code ships beside it.

Run with the `gradi` env. Kp and Ec are cached from v1; S. aureus is fetched on first run.
  python scripts/pockets/alphafill.py
  python scripts/pockets/alphafill.py --species saureus --workers 6
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
CACHE = REPO_ROOT / "data" / "source" / "alphafill"
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "pockets" / "evidence"
API = "https://alphafill.eu/v1/aff/{acc}/json"

COLS = ["uniprot_ac", "ligand", "donor_ligand", "donor_pdb", "donor_chain", "identity",
        "global_rmsd", "local_rmsd", "clash_count"]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def fetch_one(species: str, acc: str, session: requests.Session) -> str:
    """Cache one response. Returns a status; never raises."""
    out = CACHE / species / f"{acc}.json"
    if out.exists():
        return "cached"
    for attempt in range(4):
        try:
            r = session.get(API.format(acc=acc), timeout=60)
            if r.status_code == 404:
                out.write_text("null")        # no AlphaFill entry: a real answer, cached as one
                return "no_entry"
            r.raise_for_status()
            # An HTTP 200 is not evidence of data -- it must parse, and be null or carry `hits`.
            payload = r.json()
            if payload is not None and "hits" not in payload:
                return "bad_payload"
            out.write_text(json.dumps(payload))
            return "fetched"
        except (requests.RequestException, ValueError):
            time.sleep(2 ** attempt)
    return "failed"


def flatten(species: str, accs: list[str]) -> tuple[pd.DataFrame, dict[str, int]]:
    rows, stats = [], {"no_file": 0, "null": 0, "no_hits": 0, "with_hits": 0}
    for acc in accs:
        p = CACHE / species / f"{acc}.json"
        if not p.exists():
            stats["no_file"] += 1
            continue
        d = json.loads(p.read_text())
        if d is None:
            stats["null"] += 1
            continue
        hits = d.get("hits")
        if not hits:
            stats["no_hits"] += 1
            continue
        stats["with_hits"] += 1
        for h in hits:
            ident = 100 * float(h["alignment"]["identity"])
            for t in h.get("transplants") or []:
                rows.append({
                    "uniprot_ac": acc,
                    "ligand": t["analogue_id"],
                    "donor_ligand": t["compound_id"],
                    "donor_pdb": h["pdb_id"],
                    "donor_chain": h.get("pdb_asym_id"),
                    "identity": round(ident, 2),
                    "global_rmsd": h.get("global_rmsd"),
                    "local_rmsd": t.get("local_rmsd"),
                    "clash_count": (t.get("clash") or {}).get("clash_count"),
                })
    return pd.DataFrame(rows, columns=COLS), stats


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say("pockets/alphafill.py -- AlphaFill transplants, fetched and flattened (nothing filtered)")
    say(f"  in : {API.format(acc='<acc>')}")
    say("  out: data/source/alphafill/<species>/  +  evidence/transplants_<species>.tsv")
    say("=" * 92)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    summary = []
    for sp in args.species:
        accs = P.load(sp)["uniprot_ac"].tolist()
        (CACHE / sp).mkdir(parents=True, exist_ok=True)
        todo = [a for a in accs if not (CACHE / sp / f"{a}.json").exists()]
        say(f"\n[{sp}] {len(accs):,} proteins, {len(accs) - len(todo):,} cached, "
            f"{len(todo):,} to fetch")
        if todo:
            t0, status = time.time(), {}
            with requests.Session() as s, ThreadPoolExecutor(max_workers=args.workers) as ex:
                futs = {ex.submit(fetch_one, sp, a, s): a for a in todo}
                for i, fut in enumerate(as_completed(futs), 1):
                    status[futs[fut]] = fut.result()
                    if i % 250 == 0 or i == len(todo):
                        say(f"    {i:,}/{len(todo):,}  ({time.time() - t0:.0f}s)")
            counts = pd.Series(status).value_counts()
            say("    " + ", ".join(f"{k} {v:,}" for k, v in counts.items()))
            if counts.get("failed", 0):
                sys.exit(f"FATAL {counts['failed']:,} responses failed after retries -- re-run; "
                         "cached work is kept.")

        long, stats = flatten(sp, accs)
        out = EVIDENCE_DIR / f"transplants_{sp}.tsv"
        long.to_csv(out, sep="\t", index=False)
        summary.append({"species": sp, "proteins": len(accs), **stats,
                        "transplants": len(long),
                        "distinct_ligands": long["ligand"].nunique(),
                        "analogue_substituted": int((long["ligand"]
                                                     != long["donor_ligand"]).sum())})
        say(f"  wrote {out.relative_to(REPO_ROOT)}  ({len(long):,} transplants)")

    say("\n" + "-" * 92)
    say("SUMMARY -- UNFILTERED; holo.py applies drug-likeness and counts")
    say(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
