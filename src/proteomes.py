"""Load the stage-00 proteome tables.

Stage 00 writes four tables, one per species, and nothing else at the top level:

    data/processed/00_proteomes/{kpneumoniae,ecoli,saureus,human}.tsv

`load_all()` is the stacked view. It exists so the pipeline does not need a `proteins.parquet`
artifact: the first run produced one and it turned out to be a byte-for-byte concat of these four
files -- 19 MB carrying no information the tables did not already have. Same for the `id_bridge.tsv`
that `id_bridge()` replaces here.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTEOME_DIR = REPO_ROOT / "data" / "processed" / "00_proteomes"
ACCESSORY_DIR = PROTEOME_DIR / "accessory"

SPECIES = ("kpneumoniae", "ecoli", "saureus", "human")


def load(species: str) -> pd.DataFrame:
    """One species table, 11 columns, keyed on `uniprot_ac`."""
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = PROTEOME_DIR / f"{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/00_download_proteomes.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column added back."""
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_locus_tags(species: str) -> pd.DataFrame:
    """`uniprot_ac`, `locus_tag`, `locus_tag_all` — kept out of the main table to keep it simple."""
    path = ACCESSORY_DIR / f"{species}_locus_tags.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/00_download_proteomes.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def with_locus_tags(species: str) -> pd.DataFrame:
    """The species table with its locus tags joined back on, for identifier-based joins."""
    return load(species).merge(load_locus_tags(species), on="uniprot_ac", how="left")


def load_annotation(species: str) -> pd.DataFrame:
    """The wide xref layer: kegg, string, embl, eggnog, biocyc, interpro, pfam, panther, go, ec, pdb.

    Note `eggnog` and `biocyc` are 0% for K. pneumoniae -- do not assume otherwise.
    """
    path = ACCESSORY_DIR / f"{species}_annotation.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/00_download_proteomes.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def id_bridge(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """Long form `uniprot_ac x (namespace, identifier)`, for joining external datasets.

    Pulls the locus tags in from accessory/ automatically, so this stays a one-call way to reach
    every identifier a protein is known by.
    """
    rows = []
    for sp in species:
        t = load(sp)
        t = t.merge(load_locus_tags(sp), on="uniprot_ac", how="left").fillna("")
        for ns in ("locus_tag", "locus_tag_all", "refseq", "geneid", "gene_name"):
            for ac, val in zip(t["uniprot_ac"], t[ns]):
                for tok in str(val or "").replace(";", " ").split():
                    rows.append((ac, sp, ns, tok))
    return pd.DataFrame(rows, columns=["uniprot_ac", "species", "namespace", "identifier"]) \
             .drop_duplicates()
