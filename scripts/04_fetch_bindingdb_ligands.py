"""Index BindingDB affinity records by UniProt accession.

BindingDB has no clean per-target REST API. Their UniProt mapping file
contains only monomer→UniProt links (no affinities), so we download the
full bulk TSV and index it offline.

Workflow
--------
1. Download ``BindingDB_All_<release>_tsv.zip`` (~554 MB) into
   ``data/raw/bindingdb/`` if not present. The release filename is
   discovered by HEAD-checking the BindingDB downloads index, or pinned
   via ``--release YYYYMM``.
2. Stream-parse the unzipped TSV in chunks.
3. For each row whose target UniProt is in the Kp + ortholog union set,
   compute best affinity = min(Ki, Kd, IC50, EC50) in nM (drop ``>``/``<``
   qualified values — BindingDB convention for non-quantitative records).
4. Aggregate per UniProt to:
       bindingdb_any, bindingdb_10um, bindingdb_1um, bindingdb_best_pchembl

bindingdb_best_pchembl = ``9 - log10(best_nM)`` so it's directly
comparable to ChEMBL's ``pchembl_value``.

Output: ``data/processed/bindingdb_ligand_counts.tsv``.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Optional

import pandas as pd
from tenacity import retry, stop_after_attempt, wait_exponential


BINDINGDB_INDEX = "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/Download.jsp"
BINDINGDB_RAW_DIR = Path("data/raw/bindingdb")
PROCESSED_DIR = Path("data/processed")
KP_PROTEOME = Path("data/raw/klebsiella_pneumoniae_proteome.tsv")
ORTHOLOGS_TSV = PROCESSED_DIR / "klebsiella_pneumoniae_orthodb_orthologs.tsv"

USER_AGENT = "ersilia-gradi-target-prioritization (mailto:miquel@ersilia.io)"

# BindingDB stores up to N target chains per row, each with its own UniProt
# Primary ID columns (one for SwissProt, one for TrEMBL). We collect every
# "Primary ID of Target Chain <n>" column dynamically from the header so we
# don't miss multi-chain complexes or species split between SwissProt/TrEMBL.
COL_UNIPROT_RX = re.compile(
    r"^UniProt \((?:SwissProt|TrEMBL)\) Primary ID of Target Chain \d+$"
)
AFFINITY_COLS = ["Ki (nM)", "IC50 (nM)", "Kd (nM)", "EC50 (nM)"]


@retry(stop=stop_after_attempt(3), wait=wait_exponential(multiplier=2, max=60))
def _http_get_text(url: str, timeout: int = 60) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="ignore")


def resolve_release(explicit: Optional[str]) -> str:
    """Return a release id like '202506' (YYYYMM)."""
    if explicit:
        return explicit
    html = _http_get_text(BINDINGDB_INDEX)
    matches = re.findall(r"BindingDB_All_(\d{6})_tsv\.zip", html)
    if not matches:
        raise RuntimeError(
            f"Could not find BindingDB_All release on {BINDINGDB_INDEX}"
        )
    return sorted(matches)[-1]


def download_bulk(release: str, dest: Path) -> None:
    """Download BindingDB_All_<release>_tsv.zip into ``dest`` if missing.

    Uses the direct ``/rwd/bind/downloads/`` URL which serves the raw ZIP
    (content-type: application/zip). The SDFdownload.jsp wrapper returns
    an HTML acknowledgement page instead of the file.
    """
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return
    if dest.exists() and dest.stat().st_size <= 1_000_000:
        # Previous partial / HTML response — overwrite.
        dest.unlink()
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://www.bindingdb.org/rwd/bind/downloads/BindingDB_All_{release}_tsv.zip"
    print(f"Downloading {url} -> {dest}")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=1800) as resp, open(dest, "wb") as out:
        total = 0
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            total += len(chunk)
            if total % (50 * 1024 * 1024) < (1024 * 1024):
                print(f"  ... {total / (1024 * 1024):.0f} MB")
    print(f"Wrote {dest.stat().st_size / (1024 * 1024):.1f} MB to {dest}")


def open_bindingdb_tsv(zip_path: Path):
    """Return a file-like text stream for the single TSV inside the zip."""
    zf = zipfile.ZipFile(zip_path, "r")
    tsv_names = [n for n in zf.namelist() if n.lower().endswith(".tsv")]
    if not tsv_names:
        raise RuntimeError(f"No .tsv member found inside {zip_path}")
    return zf, zf.open(tsv_names[0], "r")


def load_uniprot_set(kp_proteome: Path, orthologs_tsv: Path) -> set[str]:
    accs: set[str] = set()
    with open(kp_proteome, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = header.index("Entry") if "Entry" in header else 0
        for line in f:
            cells = line.rstrip("\n").split("\t")
            if cells:
                accs.add(cells[idx].strip())
    if orthologs_tsv.exists():
        with open(orthologs_tsv, encoding="utf-8") as f:
            header = f.readline().rstrip("\n").split("\t")
            idx = header.index("ortholog_uniprot")
            for line in f:
                cells = line.rstrip("\n").split("\t")
                if idx < len(cells):
                    accs.add(cells[idx].strip())
    accs.discard("")
    return accs


def _to_nm(series: pd.Series) -> pd.Series:
    """Coerce a BindingDB affinity column to numeric nM, dropping qualified values."""
    return pd.to_numeric(series, errors="coerce")


def aggregate(zip_path: Path, allowed_uniprots: set[str], chunksize: int) -> pd.DataFrame:
    zf, fh = open_bindingdb_tsv(zip_path)
    try:
        reader = pd.read_csv(
            fh,
            sep="\t",
            dtype=str,
            chunksize=chunksize,
            on_bad_lines="skip",
            low_memory=True,
            quoting=3,  # csv.QUOTE_NONE — BindingDB doesn't quote, embedded quotes ok
        )
        per_uniprot: dict[str, dict[str, float]] = {}
        n_chunks = 0
        n_rows_total = 0
        n_rows_kept = 0
        uniprot_cols: Optional[list[str]] = None
        for chunk in reader:
            n_chunks += 1
            n_rows_total += len(chunk)
            if uniprot_cols is None:
                uniprot_cols = [c for c in chunk.columns if COL_UNIPROT_RX.match(c)]
                if not uniprot_cols:
                    raise RuntimeError(
                        "BindingDB TSV has no 'UniProt ... Primary ID of Target Chain N' "
                        f"columns; got {list(chunk.columns)[:10]}"
                    )
                print(f"  found {len(uniprot_cols)} UniProt-per-chain columns")
            # Match if ANY chain's UniProt is in the allowed set.
            row_mask = pd.Series(False, index=chunk.index)
            row_uniprots = pd.Series([""] * len(chunk), index=chunk.index)
            for col in uniprot_cols:
                in_set = chunk[col].isin(allowed_uniprots)
                # Keep the first matching uniprot per row for grouping below.
                first_match = (~row_mask) & in_set
                row_uniprots.loc[first_match] = chunk.loc[first_match, col]
                row_mask = row_mask | in_set
            if not row_mask.any():
                continue
            sub_chunk = chunk[row_mask]
            sub_uniprot = row_uniprots[row_mask]
            n_rows_kept += len(sub_chunk)
            present_cols = [c for c in AFFINITY_COLS if c in sub_chunk.columns]
            if not present_cols:
                continue
            affinity = pd.concat(
                [_to_nm(sub_chunk[c]) for c in present_cols], axis=1
            )
            affinity.columns = present_cols
            best_nm = affinity.min(axis=1)
            valid = best_nm.notna() & (best_nm > 0)
            if not valid.any():
                continue
            sub = pd.DataFrame(
                {
                    "uniprot": sub_uniprot.values[valid.values],
                    "best_nm": best_nm.values[valid.values],
                }
            )
            for uni, grp in sub.groupby("uniprot"):
                acc = per_uniprot.setdefault(
                    uni,
                    {"bindingdb_any": 0, "bindingdb_10um": 0, "bindingdb_1um": 0,
                     "best_nm": float("inf")},
                )
                acc["bindingdb_any"] += len(grp)
                acc["bindingdb_10um"] += int((grp["best_nm"] <= 10000).sum())
                acc["bindingdb_1um"] += int((grp["best_nm"] <= 1000).sum())
                cur_min = float(grp["best_nm"].min())
                if cur_min < acc["best_nm"]:
                    acc["best_nm"] = cur_min
            if n_chunks % 5 == 0:
                print(
                    f"  chunk {n_chunks}: scanned {n_rows_total} rows, "
                    f"matched {n_rows_kept} so far, {len(per_uniprot)} UniProts touched"
                )
        print(
            f"Done. Scanned {n_rows_total} rows; kept {n_rows_kept} matching the "
            f"ortholog set; {len(per_uniprot)} UniProts have ≥1 numeric affinity."
        )
    finally:
        fh.close()
        zf.close()

    rows = []
    for uni, acc in per_uniprot.items():
        best_nm = acc["best_nm"]
        best_pchembl = "" if math.isinf(best_nm) else f"{9.0 - math.log10(best_nm):.2f}"
        rows.append(
            {
                "uniprot": uni,
                "bindingdb_any": acc["bindingdb_any"],
                "bindingdb_10um": acc["bindingdb_10um"],
                "bindingdb_1um": acc["bindingdb_1um"],
                "bindingdb_best_pchembl": best_pchembl,
            }
        )
    return pd.DataFrame(rows).sort_values("uniprot").reset_index(drop=True)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--release",
        default=None,
        help="BindingDB release id YYYYMM (default: auto-detect latest).",
    )
    parser.add_argument(
        "--kp-proteome", type=Path, default=KP_PROTEOME, help="Kp proteome TSV."
    )
    parser.add_argument(
        "--orthologs",
        type=Path,
        default=ORTHOLOGS_TSV,
        help="Ortholog TSV from script 02.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=PROCESSED_DIR / "bindingdb_ligand_counts.tsv",
    )
    parser.add_argument(
        "--chunksize", type=int, default=200_000, help="pandas read_csv chunksize."
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Don't try to download; just use whatever zip is already in data/raw/bindingdb/.",
    )
    args = parser.parse_args()

    release = resolve_release(args.release)
    print(f"BindingDB release: {release}")
    zip_path = BINDINGDB_RAW_DIR / f"BindingDB_All_{release}_tsv.zip"
    if not args.skip_download:
        download_bulk(release, zip_path)
    else:
        existing = sorted(BINDINGDB_RAW_DIR.glob("BindingDB_All_*_tsv.zip"))
        if not existing:
            print(
                "No BindingDB zip found in data/raw/bindingdb/ and --skip-download set.",
                file=sys.stderr,
            )
            return 1
        zip_path = existing[-1]
        print(f"Using existing zip {zip_path}")

    accs = load_uniprot_set(args.kp_proteome, args.orthologs)
    print(f"UniProt union set: {len(accs)} accessions to match against BindingDB")

    print(f"Indexing {zip_path} ...")
    df = aggregate(zip_path, accs, args.chunksize)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(df)} rows to {args.output}")
    n_potent = int((pd.to_numeric(df["bindingdb_best_pchembl"], errors="coerce") >= 6).sum())
    print(
        f"Summary: {len(df)} UniProts with ≥1 affinity; "
        f"{n_potent} with a ≤1 µM measurement"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
