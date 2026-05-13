"""Shared constants for the scripts/ pipeline.

Each entry in ``ORGANISMS`` maps a CLI-friendly shortname to the canonical
UniProt reference proteome typically used for that species in the
literature, plus a ``slug`` used to build per-script output filenames as
``f"{slug}_{suffix}.tsv"``.

Imported as ``from _common import ORGANISMS`` by the sibling scripts;
Python sets ``sys.path[0]`` to the script's directory, so this works
regardless of the caller's working directory.
"""

from __future__ import annotations


ORGANISMS: dict[str, dict[str, str]] = {
    "kpneumoniae": {
        "proteome_id": "UP000007841",
        "label": "Klebsiella pneumoniae subsp. pneumoniae HS11286",
        "slug": "klebsiella_pneumoniae",
    },
    "ecoli": {
        "proteome_id": "UP000000625",
        "label": "Escherichia coli K-12 MG1655",
        "slug": "escherichia_coli",
    },
}

DEFAULT_ORGANISM = "kpneumoniae"
