"""Query ChEMBL for ligand activity records keyed by UniProt accession.

Input: the union of Kp UniProt accessions (from the proteome TSV) and
ortholog UniProt accessions (from ``data/processed/kp_orthodb_orthologs.tsv``).

For each accession:
1. GET ``/chembl/api/data/target.json?target_components__accession=<UNIPROT>``
   → collect every ``target_chembl_id`` returned.
2. For each TID: GET ``/chembl/api/data/activity.json?target_chembl_id=<TID>
   &standard_type__in=Ki,Kd,IC50&pchembl_value__isnull=false`` paginated until
   exhausted or ``--max-activities-per-uniprot`` reached.

Outputs ``data/processed/chembl_ligand_counts.tsv`` with columns:
    uniprot, n_chembl_targets, chembl_any, chembl_10um, chembl_1um,
    chembl_best_pchembl, target_chembl_ids, truncated

All HTTP responses are cached under ``tmp/chembl_cache/`` so reruns are
idempotent and respectful of the EBI API. Use ``--no-cache`` to force
re-fetches.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

from tenacity import retry, stop_after_attempt, wait_exponential


CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_UNIPROT_MAPPING_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
    "chembl_uniprot_mapping.txt"
)

PROCESSED_DIR = Path("data/processed")
CACHE_DIR = Path("tmp/chembl_cache")
KP_PROTEOME = Path("data/raw/klebsiella_pneumoniae_proteome.tsv")
ORTHOLOGS_TSV = PROCESSED_DIR / "klebsiella_pneumoniae_orthodb_orthologs.tsv"
MAPPING_CACHE = Path("tmp/chembl_uniprot_mapping.txt")

USER_AGENT = "ersilia-gradi-target-prioritization (mailto:miquel@ersilia.io)"


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, max=60))
def _http_get_json(url: str, timeout: int = 300) -> dict:
    req = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": USER_AGENT}
    )
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8"))


def _cache_load(path: Path) -> Optional[dict]:
    if path.exists() and path.stat().st_size > 0:
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            return None
    return None


def _cache_save(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def fetch_targets_for_uniprot(uniprot: str, use_cache: bool = True) -> list[str]:
    """Return all target_chembl_ids that list this UniProt as a component."""
    cache = CACHE_DIR / f"target_{uniprot}.json"
    payload: Optional[dict] = _cache_load(cache) if use_cache else None
    if payload is None:
        params = urllib.parse.urlencode(
            {"target_components__accession": uniprot, "limit": 50}
        )
        payload = _http_get_json(f"{CHEMBL_BASE}/target.json?{params}")
        _cache_save(cache, payload)
        time.sleep(0.2)
    tids: list[str] = []
    for tgt in payload.get("targets", []) or []:
        tid = tgt.get("target_chembl_id")
        if tid and tid not in tids:
            tids.append(tid)
    return tids


def fetch_activities_for_target(
    tid: str, max_activities: int, use_cache: bool = True
) -> tuple[list[dict], bool]:
    """Page through ``/activity.json`` for one target_chembl_id.

    Returns ``(activities, truncated)``.
    """
    base_params = {
        "target_chembl_id": tid,
        "standard_type__in": "Ki,Kd,IC50",
        "pchembl_value__isnull": "false",
        "limit": 1000,
    }
    initial_qs = urllib.parse.urlencode(base_params)
    url: Optional[str] = f"{CHEMBL_BASE}/activity.json?{initial_qs}"
    rows: list[dict] = []
    page = 0
    truncated = False
    while url:
        cache = CACHE_DIR / f"activity_{tid}_p{page}.json"
        payload: Optional[dict] = _cache_load(cache) if use_cache else None
        if payload is None:
            payload = _http_get_json(url)
            _cache_save(cache, payload)
            time.sleep(0.2)
        for act in payload.get("activities", []) or []:
            rows.append(act)
            if len(rows) >= max_activities:
                truncated = True
                break
        if truncated:
            break
        next_path = (payload.get("page_meta") or {}).get("next")
        if not next_path:
            break
        url = (
            next_path
            if next_path.startswith("http")
            else "https://www.ebi.ac.uk" + next_path
        )
        page += 1
    return rows, truncated


def summarize_activities(activities: list[dict]) -> dict[str, float]:
    """Collapse a list of activity records into tiered counts."""
    n_any = 0
    n_10um = 0
    n_1um = 0
    best_pchembl = float("nan")
    for a in activities:
        pv_raw = a.get("pchembl_value")
        if pv_raw in (None, "", "None"):
            continue
        try:
            pv = float(pv_raw)
        except (TypeError, ValueError):
            continue
        n_any += 1
        if pv >= 5.0:
            n_10um += 1
        if pv >= 6.0:
            n_1um += 1
        if math.isnan(best_pchembl) or pv > best_pchembl:
            best_pchembl = pv
    return {
        "chembl_any": n_any,
        "chembl_10um": n_10um,
        "chembl_1um": n_1um,
        "chembl_best_pchembl": best_pchembl,
    }


def load_uniprot_set(kp_proteome: Path, orthologs_tsv: Path) -> list[str]:
    """Return the deduplicated list of UniProt accessions to query."""
    accs: list[str] = []
    seen: set[str] = set()
    with open(kp_proteome, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index("Entry") if "Entry" in header else 0
        for line in f:
            cells = line.rstrip("\n").split("\t")
            if not cells:
                continue
            acc = cells[idx].strip()
            if acc and acc not in seen:
                seen.add(acc)
                accs.append(acc)
    if orthologs_tsv.exists():
        with open(orthologs_tsv, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = header.index("ortholog_uniprot")
            for line in f:
                cells = line.rstrip("\n").split("\t")
                if idx >= len(cells):
                    continue
                acc = cells[idx].strip()
                if acc and acc not in seen:
                    seen.add(acc)
                    accs.append(acc)
    return accs


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--kp-proteome",
        type=Path,
        default=KP_PROTEOME,
        help="Path to Kp proteome TSV (default: %(default)s).",
    )
    parser.add_argument(
        "--orthologs",
        type=Path,
        default=ORTHOLOGS_TSV,
        help="Path to ortholog TSV from script 02 (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROCESSED_DIR / "chembl_ligand_counts.tsv",
        help="Destination TSV (default: %(default)s).",
    )
    parser.add_argument(
        "--max-activities-per-uniprot",
        type=int,
        default=20000,
        help="Cap activities per UniProt (collapsed across TIDs) (default: %(default)s).",
    )
    parser.add_argument(
        "--limit-uniprots",
        type=int,
        default=None,
        help="Only process the first N UniProts (for smoke-testing).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Force re-fetch (ignore tmp/chembl_cache).",
    )
    args = parser.parse_args()

    accs = load_uniprot_set(args.kp_proteome, args.orthologs)
    if args.limit_uniprots:
        accs = accs[: args.limit_uniprots]
    print(f"Processing {len(accs)} UniProt accessions ...")

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)

    columns = [
        "uniprot",
        "n_chembl_targets",
        "chembl_any",
        "chembl_10um",
        "chembl_1um",
        "chembl_best_pchembl",
        "target_chembl_ids",
        "truncated",
    ]
    n_with_any = 0
    n_truncated = 0
    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for i, acc in enumerate(accs, start=1):
            try:
                tids = fetch_targets_for_uniprot(acc, use_cache=not args.no_cache)
            except Exception as e:
                print(f"  [{i}/{len(accs)}] {acc}: target ERROR {e!r}", file=sys.stderr)
                tids = []
            activities: list[dict] = []
            truncated = False
            for tid in tids:
                remaining = args.max_activities_per_uniprot - len(activities)
                if remaining <= 0:
                    truncated = True
                    break
                try:
                    rows, this_trunc = fetch_activities_for_target(
                        tid, remaining, use_cache=not args.no_cache
                    )
                except Exception as e:
                    print(
                        f"  [{i}/{len(accs)}] {acc} TID {tid}: activity ERROR {e!r}",
                        file=sys.stderr,
                    )
                    continue
                activities.extend(rows)
                if this_trunc:
                    truncated = True
            summary = summarize_activities(activities)
            best = summary["chembl_best_pchembl"]
            best_str = "" if math.isnan(best) else f"{best:.2f}"
            f.write(
                "\t".join(
                    [
                        acc,
                        str(len(tids)),
                        str(summary["chembl_any"]),
                        str(summary["chembl_10um"]),
                        str(summary["chembl_1um"]),
                        best_str,
                        ",".join(tids),
                        "True" if truncated else "False",
                    ]
                )
                + "\n"
            )
            if summary["chembl_any"] > 0:
                n_with_any += 1
            if truncated:
                n_truncated += 1
            if i % 50 == 0 or i == len(accs):
                print(
                    f"  [{i}/{len(accs)}] {n_with_any} UniProts with ChEMBL activity, "
                    f"{n_truncated} truncated"
                )

    print(f"Wrote {len(accs)} rows to {args.output}")
    print(
        f"Summary: {n_with_any}/{len(accs)} UniProts with ≥1 ChEMBL activity; "
        f"{n_truncated} truncated"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
