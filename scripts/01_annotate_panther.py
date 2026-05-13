"""Annotate a downloaded proteome with PANTHER protein family classifications.

Re-queries UniProt for the same reference proteome with the ``xref_panther``
field, then resolves PANTHER family/subfamily IDs to human-readable names
using PANTHER's HMM classification file (cached under ``tmp/``).

Writes a TSV to ``data/processed/<organism>_panther.tsv`` with columns:
``accession, gene_names, panther_family_id, panther_family_name,
panther_subfamily_id, panther_subfamily_name``. Multiple hits per protein
are joined with ``|``.
"""

from __future__ import annotations

import argparse
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

from _common import DEFAULT_ORGANISM, ORGANISMS


UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
PANTHER_DIR_URL = (
    "https://data.pantherdb.org/ftp/hmm_classifications/current_release/"
)

DEFAULT_PROCESSED_DIR = Path("data/processed")
PANTHER_CACHE_DIR = Path("tmp")


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


def fetch_proteome_panther_tsv(proteome_id: str) -> str:
    params = urllib.parse.urlencode(
        {
            "query": f"proteome:{proteome_id}",
            "format": "tsv",
            "fields": "accession,gene_names,xref_panther",
        }
    )
    return http_get_text(f"{UNIPROT_STREAM_URL}?{params}")


def resolve_panther_classifications_url() -> str:
    """Find the current PANTHER HMM classifications filename from the dir index."""
    html = http_get_text(PANTHER_DIR_URL)
    # Filenames look like PANTHER19.0_HMM_classifications
    import re

    matches = re.findall(r'href="(PANTHER[0-9.]+_HMM_classifications)"', html)
    if not matches:
        raise RuntimeError(
            f"Could not find PANTHER HMM classifications file under {PANTHER_DIR_URL}"
        )
    return PANTHER_DIR_URL + matches[0]


def load_panther_name_map(cache_path: Path) -> dict[str, str]:
    """Return ``{PTHR_ID: family_name}`` from PANTHER's HMM classifications file."""
    name_map: dict[str, str] = {}
    with open(cache_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            pthr_id, name = parts[0], parts[1]
            if pthr_id and name:
                name_map[pthr_id] = name
    return name_map


def split_panther_ids(raw: str) -> tuple[list[str], list[str]]:
    """Split a UniProt ``xref_panther`` value into (family_ids, subfamily_ids)."""
    families: list[str] = []
    subfamilies: list[str] = []
    for token in raw.split(";"):
        token = token.strip()
        if not token:
            continue
        if ":" in token:
            subfamilies.append(token)
        else:
            families.append(token)
    return families, subfamilies


def annotate(
    uniprot_tsv: str, name_map: dict[str, str]
) -> tuple[list[dict[str, str]], Counter]:
    rows: list[dict[str, str]] = []
    family_counter: Counter = Counter()
    lines = uniprot_tsv.splitlines()
    header = lines[0].split("\t")
    idx = {col: i for i, col in enumerate(header)}
    for line in lines[1:]:
        if not line:
            continue
        cells = line.split("\t")
        accession = cells[idx["Entry"]]
        gene_names = cells[idx["Gene Names"]] if "Gene Names" in idx else ""
        panther_raw = cells[idx["PANTHER"]] if "PANTHER" in idx else ""
        families, subfamilies = split_panther_ids(panther_raw)
        family_names = [name_map.get(fid, "") for fid in families]
        subfamily_names = [name_map.get(sid, "") for sid in subfamilies]
        for fid in families:
            family_counter[(fid, name_map.get(fid, ""))] += 1
        rows.append(
            {
                "accession": accession,
                "gene_names": gene_names,
                "panther_family_id": "|".join(families),
                "panther_family_name": "|".join(family_names),
                "panther_subfamily_id": "|".join(subfamilies),
                "panther_subfamily_name": "|".join(subfamily_names),
            }
        )
    return rows, family_counter


def write_tsv(rows: list[dict[str, str]], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "accession",
        "gene_names",
        "panther_family_id",
        "panther_family_name",
        "panther_subfamily_id",
        "panther_subfamily_name",
    ]
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for row in rows:
            f.write("\t".join(row[c] for c in columns) + "\n")


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
            f"{DEFAULT_PROCESSED_DIR}/<organism>_panther.tsv."
        ),
    )
    parser.add_argument(
        "--panther-cache",
        type=Path,
        default=None,
        help=(
            "Path to a cached PANTHER HMM classifications file. "
            f"Defaults to {PANTHER_CACHE_DIR}/<filename-from-PANTHER-ftp>."
        ),
    )
    parser.add_argument(
        "--top",
        type=int,
        default=10,
        help="Show this many top families in the summary (default: %(default)s).",
    )
    args = parser.parse_args()

    organism = ORGANISMS[args.organism]
    proteome_id = args.proteome_id or organism["proteome_id"]
    output_path = args.output or DEFAULT_PROCESSED_DIR / f"{organism['slug']}_panther.tsv"

    panther_url = resolve_panther_classifications_url()
    panther_filename = panther_url.rsplit("/", 1)[-1]
    panther_cache = args.panther_cache or PANTHER_CACHE_DIR / panther_filename
    if not panther_cache.exists():
        print(f"Downloading PANTHER classifications -> {panther_cache}")
        http_download(panther_url, panther_cache)
    else:
        print(f"Using cached PANTHER classifications at {panther_cache}")

    print(f"Loading PANTHER family names from {panther_cache}")
    name_map = load_panther_name_map(panther_cache)
    print(f"Loaded {len(name_map)} PANTHER IDs")

    print(
        f"Fetching UniProt PANTHER cross-references for "
        f"{organism['label']} ({proteome_id})"
    )
    uniprot_tsv = fetch_proteome_panther_tsv(proteome_id)

    rows, family_counter = annotate(uniprot_tsv, name_map)
    write_tsv(rows, output_path)

    annotated = sum(1 for r in rows if r["panther_family_id"])
    print(
        f"Wrote {len(rows)} proteins to {output_path} "
        f"({annotated} with PANTHER family, {len(rows) - annotated} without)"
    )
    print(f"\nTop {args.top} PANTHER families:")
    for (fid, fname), count in family_counter.most_common(args.top):
        label = fname or "(no name)"
        print(f"  {count:>4}  {fid}  {label}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
