"""Download a bacterial reference proteome from UniProt.

Streams a tab-separated file from the UniProt REST API into data/raw/.
Default columns: accession, gene_names, sequence.

Pick the species with --organism (e.g. kpneumoniae, ecoli); override the
strain with --proteome-id and the destination with --output if needed.
"""

from __future__ import annotations

import argparse
import sys
import urllib.parse
import urllib.request
from pathlib import Path

from _common import DEFAULT_ORGANISM, ORGANISMS


UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
DEFAULT_FIELDS = "accession,gene_names,sequence"
DEFAULT_RAW_DIR = Path("data/raw")


def build_url(proteome_id: str, fields: str) -> str:
    params = urllib.parse.urlencode(
        {
            "query": f"proteome:{proteome_id}",
            "format": "tsv",
            "fields": fields,
        }
    )
    return f"{UNIPROT_STREAM_URL}?{params}"


def download(url: str, output_path: Path, timeout: int = 600) -> int:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"Accept": "text/plain"})
    bytes_written = 0
    with urllib.request.urlopen(request, timeout=timeout) as response, open(
        output_path, "wb"
    ) as out_file:
        while True:
            chunk = response.read(1024 * 1024)
            if not chunk:
                break
            out_file.write(chunk)
            bytes_written += len(chunk)
    return bytes_written


def count_entries(path: Path) -> int:
    with open(path, "r", encoding="utf-8") as f:
        return max(sum(1 for _ in f) - 1, 0)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--organism",
        choices=sorted(ORGANISMS),
        default=DEFAULT_ORGANISM,
        help=(
            "Organism shortname (default: %(default)s). "
            + "; ".join(f"{k}: {v['label']}" for k, v in ORGANISMS.items())
            + "."
        ),
    )
    parser.add_argument(
        "--proteome-id",
        default=None,
        help="Override the UniProt proteome identifier for --organism.",
    )
    parser.add_argument(
        "--fields",
        default=DEFAULT_FIELDS,
        help="Comma-separated UniProt field names (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Destination TSV path. Defaults to "
            f"{DEFAULT_RAW_DIR}/<organism>_proteome.tsv."
        ),
    )
    args = parser.parse_args()

    organism = ORGANISMS[args.organism]
    proteome_id = args.proteome_id or organism["proteome_id"]
    output_path = args.output or DEFAULT_RAW_DIR / f"{organism['slug']}_proteome.tsv"

    url = build_url(proteome_id, args.fields)
    print(
        f"Downloading {organism['label']} proteome {proteome_id} "
        f"-> {output_path}"
    )
    n_bytes = download(url, output_path)
    n_entries = count_entries(output_path)
    print(
        f"Wrote {n_entries} entries ({n_bytes / 1024:.1f} KiB) to {output_path}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
