"""Load the stage-00 proteome tables.

Stage 00 writes four tables, one per species, and nothing else at the top level:

    data/processed/proteomes/proteome_{kpneumoniae,ecoli,saureus,human}.tsv

`load_all()` is the stacked view. It exists so the pipeline does not need a `proteins.parquet`
artifact: the first run produced one and it turned out to be a byte-for-byte concat of these four
files -- 19 MB carrying no information the tables did not already have. Same for the `id_bridge.tsv`
that `id_bridge()` replaces here.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
PROTEOME_DIR = REPO_ROOT / "data" / "processed" / "proteomes"
EVIDENCE_DIR = PROTEOME_DIR / "evidence"
SCRATCH_DIR = PROTEOME_DIR / "scratch"
SPECIES = ("kpneumoniae", "ecoli", "saureus", "human")


DELIVERABLE_COLUMNS = ["uniprot_ac", "gene_name", "protein_name", "sequence",
                       "proteomes_evidence"]

# A protein name that declares unknown function. `DUF` is literally "Domain of Unknown Function",
# so `DUF1176 domain-containing protein` belongs here -- while `Lipoprotein`, `Oxidoreductase` and
# `N-acetyltransferase domain-containing protein` name real functional classes and do not.
GENERIC_NAME = re.compile(r"uncharacteri[sz]ed|hypothetical|unknown function|\bDUF\d+", re.I)


def identity_evidence(full: pd.DataFrame) -> pd.Series:
    """`proteomes_evidence`, 1-3, from the 9-column table. The axis's half of the standard pair.

        3  the entry carries its OWN identity -- SwissProt-reviewed, or a gene symbol on the
           anchor entry itself -- AND a specific protein name
        2  one of those two
        1  neither: no symbol of its own, and a generic name

    Kp 2,205 / 2,469 / 1,054 · Ec 0 / 648 / 3,755 · Sa 1,510 / 485 / 894 ·
    human 0 / 490 / 19,926.

    **It replaces `is_reviewed`, which was nearly degenerate per species** -- E. coli and human are
    100% reviewed and K. pneumoniae is 7 of 5,728, so on three of four proteomes the boolean
    separated nothing. The raw flag stays in `evidence/proteome_full_<species>.tsv`.

    **"Own identity" is `is_reviewed OR gene_name_source == "anchor"`, and the second disjunct is
    what rescues Kp**: only 7 Kp entries are reviewed, but 1,055 carry a gene symbol on the entry
    itself. The two FILLED tiers (`species_exact`, `species_uniref90`) deliberately do NOT count --
    that is the naming-gap work this axis did (Kp 18.4% -> 63.4%) and it is inference, not the
    protein's own record. Mistaking the two would make an inferred symbol look like evidence.

    **E. coli and human have no level 1, and that is correct.** Both are 100% SwissProt, so every
    entry has been read by a curator and reaches at least 2.

    **There is no `proteome_consensus`** -- identity is not a magnitude. The convention is evidence
    always, consensus where the axis has one.
    """
    missing = [c for c in ("is_reviewed", "gene_name_source", "gene_name", "protein_name")
               if c not in full.columns]
    if missing:
        raise KeyError(f"identity_evidence needs the FULL table; absent: {missing}. "
                       f"Use load_full(), not load().")
    reviewed = full["is_reviewed"].astype(str).str.lower().isin(("true", "1")) \
        if full["is_reviewed"].dtype == object else full["is_reviewed"].astype(bool)
    own = reviewed | full["gene_name_source"].fillna("").eq("anchor")
    named = full["gene_name"].fillna("").astype(str).ne("")
    specific = ~full["protein_name"].fillna("").astype(str).str.contains(GENERIC_NAME)

    level = pd.Series(1, index=full.index, dtype="int64")
    level[(own | (named & specific)).to_numpy()] = 2
    level[(own & specific).to_numpy()] = 3
    return level.astype("Int64")


def load(species: str) -> pd.DataFrame:
    """`proteome_<species>.tsv` — 5 columns, keyed on `uniprot_ac`. **This file defines THE ROW
    ORDER** every other matrix in the project follows.

        uniprot_ac  gene_name  protein_name  sequence  proteomes_evidence

    **`proteomes_evidence` (1-3) replaced `is_reviewed` on 2026-10-04** — see `identity_evidence()`
    above for the ladder and why. The boolean was nearly degenerate per species (Ec and human 100%
    reviewed, Kp 7 of 5,728); it survives byte-identically in `evidence/proteome_full_<species>.tsv`
    via `load_full()`. There is no `proteome_consensus`: identity is not a magnitude.

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
    return _read_identity(path)


def _read_identity(path: Path) -> pd.DataFrame:
    """`dtype=str` everywhere EXCEPT the flag, which must come back a real bool.

    Everything else here is genuinely text -- accessions, names and sequences -- and
    `keep_default_na=False` keeps an empty gene name as `""` rather than NaN. But
    `is_reviewed` round-tripped as the STRINGS "True"/"False", which fails in two silent ways:
    `df[df.is_reviewed]` raises (pandas reads it as a column selection), and
    `df.is_reviewed == True` matches NOTHING because a string never equals a bool. Both of
    those read as "no reviewed proteins" rather than as an error.
    """
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if "is_reviewed" in df.columns:
        df["is_reviewed"] = df["is_reviewed"].map({"True": True, "False": False}).astype(bool)
    if "proteomes_evidence" in df.columns:
        df["proteomes_evidence"] = pd.to_numeric(df["proteomes_evidence"],
                                                errors="coerce").astype("Int64")
    return df


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
    return _read_identity(path)


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
