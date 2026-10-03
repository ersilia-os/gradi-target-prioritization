"""Load the stage-00 proteome tables.

Stage 00 writes four tables, one per species, and nothing else at the top level:

    data/processed/proteomes/proteome_{kpneumoniae,ecoli,saureus,human}.tsv

`load_all()` is the stacked view. It exists so the pipeline does not need a `proteins.parquet`
artifact: the first run produced one and it turned out to be a byte-for-byte concat of these four
files -- 19 MB carrying no information the tables did not already have. Same for the `id_bridge.tsv`
that `id_bridge()` replaces here.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTEOME_DIR = REPO_ROOT / "data" / "processed" / "proteomes"
EVIDENCE_DIR = PROTEOME_DIR / "evidence"
SCRATCH_DIR = PROTEOME_DIR / "scratch"
SPECIES = ("kpneumoniae", "ecoli", "saureus", "human")


DELIVERABLE_COLUMNS = ["uniprot_ac", "is_reviewed", "gene_name", "protein_name", "sequence"]


def load(species: str) -> pd.DataFrame:
    """`proteome_<species>.tsv` — 5 columns, keyed on `uniprot_ac`. **This file defines THE ROW
    ORDER** every other matrix in the project follows.

        uniprot_ac  is_reviewed  gene_name  protein_name  sequence

    **`sequence` stays here deliberately.** The project's standing rule is *map by sequence, not by
    accession* — HS11286 is a dark TrEMBL proteome whose accessions ChEMBL, BindingDB and the PDB
    never use — so the column every external join needs belongs in the table everything loads.

    **`gene_name` is NOT a join key.** It is 63.4% on Kp and 44.6% on Sa; join on `locus_tag`
    (`load_locus_tags()` / `with_locus_tags()`), which is 100% on Kp and 98.4% on Sa, and on
    `uniprot_ac` for human, which has no locus tags at all.

    The provenance columns — `gene_name_source`, `gene_synonyms`, `refseq`, `geneid` — moved to
    `evidence/proteome_full_<species>.tsv` on 2026-10-03, via `load_full()`.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = PROTEOME_DIR / f"proteome_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/proteomes/download.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_full(species: str) -> pd.DataFrame:
    """The 9-column table, in `evidence/` — the deliverable plus its provenance columns.

    Adds `gene_name_source` (which of the three filling tiers supplied the symbol, so an inferred
    name is never mistaken for a curated one), `gene_synonyms`, `refseq` and `geneid`. **`geneid`
    is what the literature axis keys on** — `gene2pubmed.py` and `pubtator.py` both join NCBI
    counts through it.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = EVIDENCE_DIR / f"proteome_full_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/proteomes/download.py first")
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
    path = EVIDENCE_DIR / f"locus_tags_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/proteomes/download.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def with_locus_tags(species: str) -> pd.DataFrame:
    """The species table with its locus tags joined back on, for identifier-based joins."""
    return load(species).merge(load_locus_tags(species), on="uniprot_ac", how="left")


def load_annotation(species: str) -> pd.DataFrame:
    """The wide xref layer: kegg, string, embl, eggnog, biocyc, interpro, pfam, panther, go, ec, pdb.

    Note `eggnog` and `biocyc` are 0% for K. pneumoniae -- do not assume otherwise.
    """
    path = EVIDENCE_DIR / f"annotation_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/proteomes/download.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def id_bridge(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """Long form `uniprot_ac x (namespace, identifier)`, for joining external datasets.

    Pulls the locus tags in from evidence/ + scratch/ automatically, so this stays a one-call way to reach
    every identifier a protein is known by.
    """
    rows = []
    for sp in species:
        t = load_full(sp)   # needs refseq + geneid, which the deliverable no longer carries
        t = t.merge(load_locus_tags(sp), on="uniprot_ac", how="left").fillna("")
        for ns in ("locus_tag", "locus_tag_all", "refseq", "geneid", "gene_name"):
            for ac, val in zip(t["uniprot_ac"], t[ns]):
                for tok in str(val or "").replace(";", " ").split():
                    rows.append((ac, sp, ns, tok))
    return pd.DataFrame(rows, columns=["uniprot_ac", "species", "namespace", "identifier"]) \
             .drop_duplicates()
