"""Compute per-protein structural coverage for a reference proteome.

For each UniProt entry in the chosen proteome, this script reports two
coverage metrics:

1. **PDB sequence coverage** from PDBe's SIFTS bulk mapping
   (``pdb_chain_uniprot.csv.gz``): the maximum fraction of the UniProt
   sequence covered by any single PDB chain, plus the union of all chain
   intervals (catches multi-domain proteins resolved piecewise).
2. **AlphaFold DB structural coverage at high confidence** from AFDB's
   per-prediction API: the fraction of residues with pLDDT > 70
   (= ``fractionPlddtVeryHigh`` + ``fractionPlddtConfident``), plus the
   stricter pLDDT > 90 fraction and the mean global pLDDT.

Output: ``data/processed/<slug>_structural_coverage.tsv``. Caches the
SIFTS file at ``tmp/pdb_chain_uniprot.csv.gz`` and one JSON per accession
under ``tmp/afdb/``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import statistics
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

from _common import DEFAULT_ORGANISM, ORGANISMS


UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
SIFTS_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/"
    "pdb_chain_uniprot.csv.gz"
)
AFDB_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_SIFTS_CACHE = Path("tmp/pdb_chain_uniprot.csv.gz")
DEFAULT_AFDB_CACHE_DIR = Path("tmp/afdb")
DEFAULT_WORKERS = 8


# ---------- HTTP helpers ----------

def http_get_text(url: str, timeout: int = 600) -> str:
    request = urllib.request.Request(url, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def http_download(url: str, output_path: Path, timeout: int = 1200) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url)
    with urllib.request.urlopen(request, timeout=timeout) as response, open(
        output_path, "wb"
    ) as out:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


# ---------- UniProt ----------

def fetch_proteome_metadata(proteome_id: str) -> list[dict[str, str]]:
    """Return one dict per protein with accession, gene_names, length, xref_pdb, xref_alphafolddb."""
    params = urllib.parse.urlencode(
        {
            "query": f"proteome:{proteome_id}",
            "format": "tsv",
            "fields": "accession,gene_names,length,xref_pdb,xref_alphafolddb",
        }
    )
    tsv = http_get_text(f"{UNIPROT_STREAM_URL}?{params}")
    rows = []
    lines = tsv.splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    idx = {col: i for i, col in enumerate(header)}
    for line in lines[1:]:
        if not line:
            continue
        cells = line.split("\t")
        # Pad to header length in case trailing fields are empty.
        while len(cells) < len(header):
            cells.append("")
        rows.append(
            {
                "accession": cells[idx["Entry"]],
                "gene_names": cells[idx.get("Gene Names", -1)] if "Gene Names" in idx else "",
                "length": cells[idx["Length"]],
                "xref_pdb": cells[idx["PDB"]] if "PDB" in idx else "",
                "xref_alphafolddb": cells[idx["AlphaFoldDB"]] if "AlphaFoldDB" in idx else "",
            }
        )
    return rows


# ---------- SIFTS ----------

def ensure_sifts_cache(cache_path: Path) -> None:
    if cache_path.exists():
        print(f"Using cached SIFTS mapping at {cache_path}")
        return
    print(f"Downloading SIFTS mapping -> {cache_path}")
    http_download(SIFTS_URL, cache_path)


def parse_sifts_for_accessions(
    cache_path: Path, accessions: set[str]
) -> tuple[dict[str, list[tuple[str, str, int, int]]], int]:
    """Stream-parse SIFTS, keeping only rows for ``accessions``.

    Returns ``(per_accession_rows, skipped_count)`` where
    ``per_accession_rows[acc]`` is a list of ``(pdb_id, chain_id, sp_beg, sp_end)``.
    """
    per_acc: dict[str, list[tuple[str, str, int, int]]] = defaultdict(list)
    skipped = 0
    with gzip.open(cache_path, "rt", encoding="utf-8", newline="") as f:
        first = f.readline()
        if not first.startswith("#"):
            # No comment line — rewind by re-opening.
            f.seek(0)
        reader = csv.DictReader(f)
        for row in reader:
            acc = row.get("SP_PRIMARY", "")
            if acc not in accessions:
                continue
            try:
                sp_beg = int(row["SP_BEG"])
                sp_end = int(row["SP_END"])
            except (KeyError, TypeError, ValueError):
                skipped += 1
                continue
            if sp_beg > sp_end:
                skipped += 1
                continue
            per_acc[acc].append((row["PDB"], row["CHAIN"], sp_beg, sp_end))
    return per_acc, skipped


def union_intervals(intervals: Iterable[tuple[int, int]]) -> int:
    """Return the total number of integer residues covered by [a, b] intervals (inclusive)."""
    sorted_iv = sorted(intervals)
    total = 0
    cur_lo: int | None = None
    cur_hi: int | None = None
    for lo, hi in sorted_iv:
        if cur_lo is None:
            cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi + 1:
            cur_hi = max(cur_hi, hi)
        else:
            total += cur_hi - cur_lo + 1
            cur_lo, cur_hi = lo, hi
    if cur_lo is not None:
        total += cur_hi - cur_lo + 1
    return total


def compute_pdb_coverage(
    sifts_rows: list[tuple[str, str, int, int]], length: int
) -> dict[str, object]:
    """Compute PDB coverage metrics for one accession.

    Returns dict with: pdb_chain_count, pdb_best_id, pdb_best_chain,
    pdb_max_coverage, pdb_total_coverage_union.
    """
    if not sifts_rows or length <= 0:
        return {
            "pdb_chain_count": 0,
            "pdb_best_id": "",
            "pdb_best_chain": "",
            "pdb_max_coverage": "",
            "pdb_total_coverage_union": "",
        }
    chains: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    all_intervals: list[tuple[int, int]] = []
    for pdb_id, chain_id, sp_beg, sp_end in sifts_rows:
        chains[(pdb_id, chain_id)].append((sp_beg, sp_end))
        all_intervals.append((sp_beg, sp_end))

    best_id, best_chain, best_cov = "", "", 0.0
    for (pdb_id, chain_id), intervals in chains.items():
        covered = union_intervals(intervals)
        cov = min(covered / length, 1.0)
        if cov > best_cov:
            best_id, best_chain, best_cov = pdb_id, chain_id, cov

    total_union = min(union_intervals(all_intervals) / length, 1.0)
    return {
        "pdb_chain_count": len(chains),
        "pdb_best_id": best_id,
        "pdb_best_chain": best_chain,
        "pdb_max_coverage": f"{best_cov:.3f}",
        "pdb_total_coverage_union": f"{total_union:.3f}",
    }


# ---------- AFDB ----------

def fetch_afdb_prediction(
    accession: str,
    cache_dir: Path,
    max_attempts: int = 3,
) -> dict:
    """Fetch and cache AFDB prediction summary for one accession.

    Returns a dict ``{"status": int, "data": list|dict|None}``. Caches 200,
    404, and 400 permanently; retries transient failures with backoff.
    """
    cache_path = cache_dir / f"{accession}.json"
    if cache_path.exists():
        with open(cache_path, "r", encoding="utf-8") as f:
            return json.load(f)

    url = AFDB_API_URL.format(accession=accession)
    last_exc: Exception | None = None
    for attempt in range(max_attempts):
        try:
            # Jitter inside the worker to avoid synchronized hammer.
            time.sleep(random.uniform(0.05, 0.20))
            request = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=30) as response:
                body = response.read().decode("utf-8")
                result = {"status": response.status, "data": json.loads(body) if body else None}
                break
        except urllib.error.HTTPError as e:
            if e.code in (404, 400):
                result = {"status": e.code, "data": None}
                break
            if e.code in (429, 500, 502, 503, 504) and attempt < max_attempts - 1:
                time.sleep(2**attempt + random.uniform(0, 0.5))
                last_exc = e
                continue
            return {"status": e.code, "data": None, "_transient": True}
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_exc = e
            if attempt < max_attempts - 1:
                time.sleep(2**attempt + random.uniform(0, 0.5))
                continue
            return {"status": -1, "data": None, "_transient": True, "_error": str(e)}
    else:  # noqa: PLW0120 (only entered if loop completes without break)
        return {"status": -1, "data": None, "_transient": True, "_error": str(last_exc)}

    cache_dir.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_suffix(".json.tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(result, f)
    tmp.replace(cache_path)
    return result


def extract_afdb_metrics(result: dict) -> dict[str, str]:
    blank = {
        "afdb_model_id": "",
        "afdb_global_plddt": "",
        "afdb_high_confidence_fraction": "",
        "afdb_very_high_confidence_fraction": "",
    }
    if not result or result.get("status") != 200:
        return blank
    data = result.get("data")
    if not isinstance(data, list) or not data:
        return blank
    entry = data[0]
    very_high = entry.get("fractionPlddtVeryHigh")
    confident = entry.get("fractionPlddtConfident")
    high = (very_high or 0) + (confident or 0) if (very_high is not None or confident is not None) else None
    return {
        "afdb_model_id": entry.get("modelEntityId", "") or "",
        "afdb_global_plddt": f"{entry['globalMetricValue']:.2f}" if entry.get("globalMetricValue") is not None else "",
        "afdb_high_confidence_fraction": f"{high:.3f}" if high is not None else "",
        "afdb_very_high_confidence_fraction": f"{very_high:.3f}" if very_high is not None else "",
    }


def fetch_afdb_for_accessions(
    accessions: list[str], cache_dir: Path, workers: int
) -> dict[str, dict]:
    results: dict[str, dict] = {}
    to_fetch = [
        a for a in accessions if not (cache_dir / f"{a}.json").exists()
    ]
    cached = len(accessions) - len(to_fetch)
    if cached:
        print(f"AFDB cache hits: {cached}/{len(accessions)}")
    # Load already-cached entries first.
    for acc in accessions:
        cp = cache_dir / f"{acc}.json"
        if cp.exists():
            with open(cp, "r", encoding="utf-8") as f:
                results[acc] = json.load(f)

    if not to_fetch:
        return results

    print(f"Fetching AFDB predictions for {len(to_fetch)} accessions "
          f"({workers} workers)")
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(fetch_afdb_prediction, acc, cache_dir): acc
            for acc in to_fetch
        }
        for fut in as_completed(futures):
            acc = futures[fut]
            try:
                results[acc] = fut.result()
            except Exception as e:  # noqa: BLE001
                results[acc] = {"status": -1, "data": None, "_error": str(e)}
            done += 1
            if done % 500 == 0 or done == len(to_fetch):
                print(f"  AFDB progress: {done}/{len(to_fetch)}")
    return results


# ---------- Output + summary ----------

OUTPUT_COLUMNS = [
    "accession",
    "gene_names",
    "length",
    "pdb_xref_count",
    "pdb_chain_count",
    "pdb_best_id",
    "pdb_best_chain",
    "pdb_max_coverage",
    "pdb_total_coverage_union",
    "afdb_model_id",
    "afdb_global_plddt",
    "afdb_high_confidence_fraction",
    "afdb_very_high_confidence_fraction",
]


def write_tsv(rows: list[dict[str, object]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(str(row.get(c, "")) for c in OUTPUT_COLUMNS) + "\n")


def print_summary(rows: list[dict[str, object]], sifts_skipped: int, plddt_threshold: float) -> None:
    n = len(rows)
    lengths = [int(r["length"]) for r in rows if str(r["length"]).isdigit()]
    print(f"\nTotal proteins: {n}")
    if lengths:
        print(
            f"  Length: mean={statistics.mean(lengths):.0f}  "
            f"median={statistics.median(lengths):.0f}  "
            f"min={min(lengths)}  max={max(lengths)}"
        )

    # PDB
    with_pdb = [r for r in rows if r["pdb_chain_count"] and int(r["pdb_chain_count"]) > 0]
    covs = [float(r["pdb_max_coverage"]) for r in with_pdb if r["pdb_max_coverage"]]
    print(f"\nPDB (via SIFTS):")
    print(f"  Proteins with any PDB chain: {len(with_pdb)} ({100 * len(with_pdb) / max(n, 1):.1f}%)")
    if covs:
        ge50 = sum(1 for c in covs if c >= 0.5)
        ge90 = sum(1 for c in covs if c >= 0.9)
        print(f"  max_coverage >= 0.50: {ge50}    >= 0.90: {ge90}")
        print(f"  max_coverage  mean={statistics.mean(covs):.3f}  median={statistics.median(covs):.3f}")

    # Sanity: UniProt xref_pdb says X but SIFTS gives 0
    diverge = sum(
        1
        for r in rows
        if int(r["pdb_xref_count"] or 0) > 0 and int(r["pdb_chain_count"] or 0) == 0
    )
    if diverge:
        print(
            f"  {diverge} proteins have a UniProt PDB xref but no SIFTS rows "
            f"(typical TrEMBL/SIFTS lag)."
        )
    if sifts_skipped:
        print(f"  SIFTS rows skipped (malformed/invalid SP_BEG/SP_END): {sifts_skipped}")

    # AFDB
    with_afdb = [r for r in rows if r["afdb_model_id"]]
    plddts = [float(r["afdb_global_plddt"]) for r in with_afdb if r["afdb_global_plddt"]]
    print(f"\nAlphaFold DB:")
    print(f"  Proteins with model: {len(with_afdb)} ({100 * len(with_afdb) / max(n, 1):.1f}%)")
    if plddts:
        print(f"  global pLDDT: mean={statistics.mean(plddts):.2f}  median={statistics.median(plddts):.2f}")
    high_fracs = [float(r["afdb_high_confidence_fraction"]) for r in with_afdb if r["afdb_high_confidence_fraction"]]
    if high_fracs:
        print(f"  fraction pLDDT>{plddt_threshold:.0f} (mean over modelled): {statistics.mean(high_fracs):.3f}")
        bins = [0] * 10
        for f in high_fracs:
            i = min(int(f * 10), 9)
            bins[i] += 1
        print(f"  histogram of high-confidence fraction (10% bins):")
        for i, c in enumerate(bins):
            lo = i * 10
            hi = (i + 1) * 10
            bar = "#" * min(c // max(1, max(bins) // 40 or 1), 40)
            print(f"    {lo:>3d}-{hi:<3d}%  {c:>5d}  {bar}")


# ---------- Main ----------

def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--organism",
        choices=sorted(ORGANISMS),
        default=DEFAULT_ORGANISM,
        help="Organism shortname (default: %(default)s).",
    )
    parser.add_argument(
        "--proteome-id",
        default=None,
        help="Override the UniProt proteome identifier for --organism.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Destination TSV. Defaults to "
            f"{DEFAULT_PROCESSED_DIR}/<slug>_structural_coverage.tsv."
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help="Concurrent workers for AFDB API (default: %(default)s).",
    )
    parser.add_argument(
        "--sifts-cache",
        type=Path,
        default=DEFAULT_SIFTS_CACHE,
        help="Path to cached SIFTS pdb_chain_uniprot.csv.gz (default: %(default)s).",
    )
    parser.add_argument(
        "--afdb-cache-dir",
        type=Path,
        default=DEFAULT_AFDB_CACHE_DIR,
        help="Directory for per-accession AFDB JSON cache (default: %(default)s).",
    )
    parser.add_argument(
        "--plddt-threshold",
        type=float,
        default=70.0,
        help=(
            "pLDDT cutoff for the labelled summary line (default: %(default)s). "
            "Both >70 and >90 fractions are always written to the TSV."
        ),
    )
    args = parser.parse_args()

    organism = ORGANISMS[args.organism]
    proteome_id = args.proteome_id or organism["proteome_id"]
    output_path = (
        args.output or DEFAULT_PROCESSED_DIR / f"{organism['slug']}_structural_coverage.tsv"
    )

    # 1. UniProt metadata.
    print(f"Fetching UniProt metadata for {organism['label']} ({proteome_id})")
    uniprot_rows = fetch_proteome_metadata(proteome_id)
    print(f"  Got {len(uniprot_rows)} proteins")
    accessions = {r["accession"] for r in uniprot_rows}

    # 2. SIFTS cache + parse.
    ensure_sifts_cache(args.sifts_cache)
    print("Parsing SIFTS (streaming, filtered to this proteome)")
    sifts_per_acc, sifts_skipped = parse_sifts_for_accessions(args.sifts_cache, accessions)
    n_with_sifts = sum(1 for v in sifts_per_acc.values() if v)
    print(f"  SIFTS rows for {n_with_sifts} of {len(accessions)} proteome accessions "
          f"({sifts_skipped} rows skipped)")

    # 3. AFDB fetch (only for accessions with an AFDB xref).
    afdb_candidates = [
        r["accession"] for r in uniprot_rows if r["xref_alphafolddb"].strip()
    ]
    print(f"AFDB candidates (have xref_alphafolddb): {len(afdb_candidates)}")
    afdb_results = fetch_afdb_for_accessions(
        afdb_candidates, args.afdb_cache_dir, args.workers
    )

    # 4. Build output rows.
    rows: list[dict[str, object]] = []
    for r in uniprot_rows:
        acc = r["accession"]
        try:
            length = int(r["length"])
        except (TypeError, ValueError):
            length = 0
        xref_count = sum(1 for t in r["xref_pdb"].split(";") if t.strip())
        pdb_metrics = compute_pdb_coverage(sifts_per_acc.get(acc, []), length)
        afdb_metrics = extract_afdb_metrics(afdb_results.get(acc, {}))
        rows.append(
            {
                "accession": acc,
                "gene_names": r["gene_names"],
                "length": length,
                "pdb_xref_count": xref_count,
                **pdb_metrics,
                **afdb_metrics,
            }
        )

    write_tsv(rows, output_path)
    print(f"\nWrote {len(rows)} rows to {output_path}")
    print_summary(rows, sifts_skipped, args.plddt_threshold)
    return 0


if __name__ == "__main__":
    sys.exit(main())
