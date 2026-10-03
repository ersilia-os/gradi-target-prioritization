"""Load the stage-05 orthology tables: who is an ortholog of whom, and how similar they are.

Stage 05 answers the same question two ways, and keeps them apart:

    data/processed/orthology/orthologs.tsv     the DISCRETE matrix -- sparse, ortholog pairs only
    data/processed/orthology/neighbors.tsv     the CONTINUOUS matrix -- sparse, top-5 per species
    data/processed/orthology/orthology_<species>.tsv            uniprot_ac + 2 columns
    data/processed/orthology/evidence/orthology_<species>.tsv   DENSE, 29 columns

Why two matrices, and why a third table
---------------------------------------
v1 shipped orthology as membership alone. Its own retrospective calls this **trap 1**: the kp->ec
table's `pident`, `coverage` and `bitscore` columns were entirely empty, so *"any axis that
thresholds transfer on percent identity silently drops everything"*. The continuous matrix exists so
no consumer ever has to threshold on a column nobody populated.

The dense per-protein table exists for a subtler reason. A **sparse** matrix cannot express a zero --
absence of a row *is* the zero -- so "this protein has no E. coli ortholog" is not readable off it.
`load(species)` gives one row for every accession with `n_orthologs_<sp>` counts, and there a 0 is a
value you can read.

A zero here is a MEASURED zero
------------------------------
That claim rests on how the groups were built, not on convention. OrthoFinder was run **de novo on
our own four proteomes**, so every input protein is either assigned an orthogroup or named in
`Orthogroups_UnassignedGenes.tsv` -- both measured outcomes -- and the stage exits non-zero unless
assigned plus unassigned accounts for all 33,436 accessions exactly. Every accession is likewise a
query in all four DIAMOND searches; `searched` records it.

**No accession-LOOKUP database is used, and the distinction matters.** A lookup that returns nothing
because the accession is not cross-referenced is an *unknown*; a search that returns nothing because
the sequence matched nothing is a *measured* negative. On this anchor the two diverge violently:
UniProt's eggNOG xref covers **0.00% of Kp** while eggNOG-mapper, run over the same proteins in stage
02, assigned **91.5%** -- and v1's OrthoDB attempt pivoted through gene symbols and reached 18.3%,
which is just Kp's 18.4% native gene-name coverage.

So the rule is *search, never look up* (CLAUDE.md: "map by sequence, not by accession, when reaching
an external database"), and every group source here obeys it: OrthoFinder de novo, and OrthoDB by
DIAMOND against each species' own OrthoDB assembly. `orthodb_<species>.tsv` adds the ABSOLUTE,
panel-independent grouping that OrthoFinder's de novo orthogroups cannot be.

The two ortholog calls are never merged
---------------------------------------
`is_ortholog_orthofinder` and `is_rbh` sit side by side. They disagree: v1 computed both on the
Kp-Ec pair (3,179 vs 3,003 pairs) and never compared them. Merging would hide that, and v1's
retrospective already called an analogous max-pool over two activators "inflation".

`is_ortholog_orthofinder` is **cross-species by definition** -- orthology is a between-species
relation and OrthoFinder's `Orthologues/` emits nothing within a species. A `False` on a
within-species pair is a category error, not a measured zero. Within-species relatedness is
paralogy, which `same_orthogroup` and `n_paralogs` carry.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
ORTHOLOGY_DIR = REPO_ROOT / "data" / "processed" / "orthology"
EVIDENCE_DIR = ORTHOLOGY_DIR / "evidence"
SCRATCH_DIR = ORTHOLOGY_DIR / "scratch"
# All FOUR proteomes -- unlike `src/function.py`, human is in scope here. It is the whole point of
# the selectivity question: a target with a human ortholog is an off-target liability.
SPECIES = ("kpneumoniae", "ecoli", "saureus", "human")

# ---------------------------------------------------------------- the chosen configuration

TOP_K = 5  # nearest neighbours kept per (query, target species)
EVALUE = 1e-3
SENSITIVITY = "very-sensitive"

# Against v1's plain DIAMOND defaults, and the reason is the direction of the error. v1 found a human
# hit for only 732/5,728 Kp proteins (12.8%). Under-detecting human homology makes a target look MORE
# selective than it is, which is the expensive mistake here -- it advances a candidate that in fact
# has a human ortholog. Identity and coverage are kept as columns, never as filters, so a missing row
# means "no detectable homology at this e-value", which is a far stronger statement than v1's.

# ---------------------------------------------------------------- v1, as pre-registered expectations

# Measured by v1 and quoted across four of its documents. If a v2 run does not land near these,
# something is wrong with the run, not with v1. Written to `evidence/control.tsv` every run.
V1_KP_WITH_ECOLI_ORTHOLOG = 3179  # of 5,728 = 55.5%, the ceiling on Kp<->Ec transfer
V1_KP_EC_RBH_PAIRS = 3003
V1_RBH_MEDIAN_PIDENT = {"kpneumoniae__ecoli": 86.0, "kpneumoniae__human": 37.0,
                        "ecoli__human": 36.7}
V1_THREE_WAY_OG_ALL_SPECIES = 608
V1_THREE_WAY_SINGLE_COPY = 303
V1_KP_WITH_ANY_HUMAN_HIT = 732  # 12.8% at default sensitivity -- v2 should EXCEED this

TOLERANCE_FRACTION = 0.15  # how far an ortholog count may drift from v1 before the stage complains

# ---------------------------------------------------------------- OrthoDB: the ABSOLUTE groups

# OrthoFinder's orthogroups are inferred de novo and are therefore a function of the PANEL: recall
# rises with the number of species in the run (measured -- a 2-species run yielded 631 ortholog
# groups where a 4-species run yielded 2,568 Kp-Ec pairs). For a grouping that does not move when the
# query changes, the groups have to be defined elsewhere, over a fixed set of genomes. That is what
# OrthoDB is, and `orthodb_<species>.tsv` carries it alongside rather than instead.
ORTHODB_VERSION = "odb12v2"
ORTHODB_BASE = "https://data.orthodb.org/v12/download/odb_data_dump"

# v12 REPLACES v11 here; the two id spaces must never be joined. OrthoDB's own README says an
# "OG unique id (not stable and re-used between releases)", so a v11 id relabelled as v12 would be
# silently wrong. Every row carries `orthodb_version` for exactly this reason.
#
# The assignment is no longer "DIAMOND against this species' own OrthoDB assembly". That searched one
# organism (for Kp, a single 4,975-protein strain) and capped the anchor at 74.6%. v12 searches
# representatives of EVERY domain-level orthogroup and transfers a group only at levels our own
# lineage passes through. Measured: Kp 74.6 -> 91.1%, Ec 88.8 -> 95.1%, Sa 88.3 -> 91.3%,
# human -> 95.5%.

# Our species' OrthoDB organism code and the assembly OrthoDB built it from. Three of the four
# organisms are present; the anchor is not, and borrows a within-species strain.
ORTHODB_ORG = {
    "kpneumoniae": "72407_0",  # K. pneumoniae subsp. pneumoniae -- NOT HS11286, see below
    "ecoli": "83333_0",        # Escherichia coli K-12
    "saureus": "93061_0",      # S. aureus subsp. aureus NCTC 8325 -- our exact strain
    "human": "9606_0",         # Homo sapiens (T2T-CHM13v2.0)
}
ORTHODB_ASSEMBLY = {
    "72407_0": "GCF_001645745.1",
    "83333_0": "GCF_000974885.1",
    "93061_0": "GCF_000013425.1",
    "9606_0": "GCF_009914755.1",
}
ORTHODB_ORG_EXPECTED = {"83333_0": 4098, "93061_0": 2573, "9606_0": 20345, "72407_0": 4975}

# A group transfers only at a level OUR organism belongs to -- a Kp protein matching a Bacillus gene
# whose only group is Bacillus-level has learned nothing about Kp. From level2species column 4.
ORTHODB_LINEAGE = {
    "kpneumoniae": (2, 1224, 1236, 91347, 543, 570, 72407),
    "ecoli": (2, 1224, 1236, 91347, 543, 561, 83333),
    "saureus": (2, 1239, 91061, 1385, 90964, 1279, 93061),
    "human": (2759, 33208, 7742, 32523, 40674, 9347, 314146, 9443, 314295, 9604, 9606),
}

# Identity floor, DECOY-CALIBRATED rather than inherited. 33,436 composition-preserving shuffles of
# our own sequences produced 0 hits at 40%/50% and only 2 at 25%/50%, while the 40% floor cost ~10
# points of coverage on Kp and ~22 on S. aureus -- and E. coli correctness is FLAT across floors
# (73.3% at 40/50 vs 72.6% at 20/50). CLAUDE.md's >=40% rule is for annotation TRANSFER, a stricter
# task; do not conflate the two.
ORTHODB_MIN_PIDENT = 25.0
ORTHODB_MIN_COVERAGE = 50.0

# `orthodb_confidence` is the winning group's share of total bitscore at the domain level, and it is
# the ONLY good trust signal: measured on E. coli, >0.9 is 95.5% correct and <0.5 only 41.4%, while
# percent identity barely separates (95-100% identity is still just 79.0%). Filter on this, not on
# identity. Sequence-tier calls are ~65% correct overall, so for K. pneumoniae -- which has NO
# OrthoDB representation and is therefore 100% sequence-tier -- treat low-confidence rows as weak.
ORTHODB_HIGH_CONFIDENCE = 0.9

# ACCESSION JOINS INTO ORTHODB DO NOT WORK, measured three ways. OrthoDB's own `uniprot_id` field
# matches our accessions for only **62.1% of E. coli and 54.8% of human** (95.4% for S. aureus, whose
# exact strain OrthoDB happens to carry) because OrthoDB builds from RefSeq assemblies whose UniProt
# cross-references point at other substrains and at unreviewed entries our reviewed-only human
# proteome excludes. Joining on RefSeq instead is worse: 5.9% / 33.9%.
#
# So EVERY species is assigned by sequence, not just the anchor -- CLAUDE.md's standing rule, arrived
# at here for the third time. We DIAMOND each proteome against the protein set of its own OrthoDB
# assembly and join on `protein_id`, which OrthoDB populates for 100% of genes.
ORTHODB_SAME_PROTEIN_PIDENT = 99.5  # at/above this the OrthoDB entry IS our protein, not a transfer

# K. pneumoniae HS11286 is absent from OrthoDB entirely. Measured, via the REST API: `P0A6G7` (Ec
# ClpP) returns 6 groups and `Q16740` (human CLPP) returns 10, while `A0A0H3GJ69` and the *reviewed*
# `A0A0H3GGB5`/`A0A0H3GKK9` return 0. Same dark-TrEMBL problem that gives UniProt's eggNOG xref
# 0.00% on this proteome -- and confirmed here on OrthoDB's own side, not just UniProt's.
#
# So the anchor is assigned by SEQUENCE, per CLAUDE.md's standing rule ("map by sequence, not by
# accession, when reaching an external database"). OrthoDB indexes exactly one K. pneumoniae, so the
# transfer is WITHIN-species, which is the nearly-lossless case: stage 04 measured a within-species
# DIAMOND join at 96.4% success, median 100% identity.
MIN_TRANSFER_PIDENT = 95.0    # CLAUDE.md's ">=95% identity = direct" rule
MIN_TRANSFER_COVERAGE = 0.80

# **OrthoDB HAS NO ROOT LEVEL.** Its top levels are Bacteria (2), Archaea (2157) and Eukaryota
# (2759) -- 1,003 levels, none spanning domains. So OrthoDB CANNOT place a bacterial protein and a
# human protein in the same group, and no column here answers "does this target have a human
# ortholog". That question is answered twice over already: `neighbors.tsv` (continuous, complete, no
# coverage gap) and OrthoFinder's joint run. What OrthoDB adds is bacterial conservation breadth,
# over 17,551 bacterial species -- which is the question it is genuinely best at.
ORTHODB_DOMAIN_LEVEL = {"kpneumoniae": 2, "ecoli": 2, "saureus": 2, "human": 2759}
ORTHODB_LEVEL_NAME = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota"}

# The verdict vocabulary. A closed set that never contains "unknown": 100% of accessions leave with
# one of these, which is what makes a `no_group` a measured zero rather than a gap.
# `assigned_by_uniprot` means the protein IS an OrthoDB gene, so these are its own groups and no
# inference happened; `assigned_by_sequence` is a transfer from homologs, decided by a bitscore vote
# across all hits, with `orthodb_confidence` recording how contested that vote was so it can be
# judged. (v11's `assigned_same_protein` is gone -- the identity-threshold self-match it named no
# longer exists as a separate route.)
ORTHODB_VERDICTS = ("assigned_by_uniprot", "assigned_by_sequence", "no_group")

NEIGHBOR_COLUMNS = (
    "query_ac", "query_species", "target_ac", "target_species", "rank",
    "pident", "ppos", "alnlen", "qcov", "scov", "evalue", "bitscore", "bitscore_norm",
    "is_ortholog_orthofinder", "is_rbh",
)
ORTHOLOG_COLUMNS = (
    "query_ac", "query_species", "target_ac", "target_species",
    "is_ortholog_orthofinder", "is_rbh", "same_orthogroup", "orthogroup",
)

_BOOL = {"True": True, "False": False, "": False}


def _check(species: str) -> None:
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")


def _path(directory: Path, name: str) -> Path:
    path = directory / name
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/orthology/orthofinder.py first")
    return path


def _read(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for col in df.columns:
        if col.startswith(("is_", "same_", "in_", "has_", "searched")):
            df[col] = df[col].map(lambda v: _BOOL.get(v, bool(v)))
        elif (col in ("rank", "orthogroup_size", "orthodb_n_levels", "orthodb_n_candidate_ogs",
                      "level", "bacterial_panel_size")
              or col.startswith("n_")):
            # `n_` prefix, not a hand-kept list. The list was the bug: `n_bacterial_orthologs` was
            # added to the writer and came back as TEXT, so `<= panel_size` raised TypeError --
            # the same shape as the `best_pchembl` string-dtype trap in the ligands axis. Any new
            # count column is now numeric automatically; a writer and a loader that must be edited
            # in lockstep will eventually not be.
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
        elif col in ("orthodb_confidence", "orthodb_match_pident",
                     "bacterial_panel_orthologs"):
            # Float, not str: `orthodb_confidence > 0.9` is the documented way to filter these rows,
            # and leaving it as text makes that comparison raise instead of work.
            # `bacterial_panel_orthologs` is here because it is a FRACTION whose name reads like a
            # count -- it matches none of the prefixes above and fell through to `str`, so
            # `.mean()` concatenated 5,728 strings and any sort was lexical. Exactly the trap the
            # `n_` comment above describes, one column later.
            df[col] = pd.to_numeric(df[col], errors="coerce")
        elif col in ("pident", "ppos", "alnlen", "qcov", "scov", "evalue", "bitscore",
                     "bitscore_norm") or col.startswith(("best_identity_", "best_bitscore_")):
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


# The deliverable's dtypes, asserted on every load. `_read()` is deliberately permissive -- it
# also parses raw OrthoFinder output whose columns are named after proteomes -- so the strict
# check lives here, where the schema is fixed and small. It exists because
# `bacterial_panel_orthologs` shipped as TEXT: it matched none of `_read()`'s prefixes, so
# `.mean()` concatenated 5,728 strings and every sort was lexical.
DELIVERABLE_DTYPES = {"uniprot_ac": "object", "has_human_ortholog": "bool",
                      "bacterial_panel_orthologs": "float64"}


def load_orthologs() -> pd.DataFrame:
    """The discrete matrix: one row per ordered pair called an ortholog by *either* method.

    Sparse. A pair absent from this table was called an ortholog by neither method -- which is a
    measured negative, guaranteed by the coverage guard, not an absence of evidence.
    """
    return _read(_path(ORTHOLOGY_DIR, "orthologs.tsv"))


def load_neighbors() -> pd.DataFrame:
    """The continuous matrix: the `TOP_K` nearest neighbours in each target species.

    Self-hits are dropped, so the within-species block is a protein's nearest **paralogs**, not
    itself. `rank` is 1..k within each (query_ac, target_species).
    """
    return _read(_path(ORTHOLOGY_DIR, "neighbors.tsv"))


def load(species: str) -> pd.DataFrame:
    """`orthology_<species>.tsv` -- THE DELIVERABLE. Three columns, the three bacteria only.

        uniprot_ac  has_human_ortholog  bacterial_panel_orthologs

    **`has_human_ortholog`** is the selectivity liability. Under-detecting human homology would
    make a target look MORE selective than it is, which is why the search runs `--very-sensitive`.

    **`bacterial_panel_orthologs` IS A FRACTION, 0-1, not a count** -- the name reads like one, so
    check the scale before using it. It is OVER 28, not 26 -- the 26 tier-C comparator proteomes plus
    the three bacterial anchors, minus this protein's own species. It counts SPECIES sharing the
    orthogroup, never proteins, so a paralog pair does not inflate it. Median 0.536 Kp / 0.571 Ec
    / **0.179 Sa** -- that last is *S. aureus*'s Gram-positive isolation against a panel that is
    mostly Gram-negative, not a defect.

    **A 0 is MEASURED.** OrthoFinder runs de novo on our own FASTAs and the stage exits unless
    every protein is accounted for, so "no bacterial ortholog" is a finding, not a lookup miss.

    **Human has no deliverable** -- the panel columns are bacteria-only. Use `load_dense("human")`.

    Everything else -- orthogroup, paralogs, per-species ortholog counts and identities -- is in
    `load_dense()`.
    """
    _check(species)
    if species == "human":
        raise ValueError("human has no orthology deliverable (the panel is bacterial); "
                         'use load_dense("human")')
    df = _read(_path(ORTHOLOGY_DIR, f"orthology_{species}.tsv"))
    wrong = {c: str(df[c].dtype) for c, want in DELIVERABLE_DTYPES.items()
             if c in df.columns and str(df[c].dtype) != want}
    if wrong:
        raise TypeError(
            f"orthology_{species}.tsv has the wrong dtypes: {wrong} (expected "
            f"{ {c: DELIVERABLE_DTYPES[c] for c in wrong} }). A column that stays `object` is "
            f"text pretending to be a number -- add a rule to _read() above.")
    return df


def load_dense(species: str) -> pd.DataFrame:
    """The 29-column dense view, in `evidence/` -- one row per protein, every accession present.

    This is the table to read a zero off: `n_orthologs_<sp> == 0` means the protein was in the
    OrthoFinder run and in every DIAMOND search, and no ortholog was found. Carries `orthogroup`,
    `in_orthogroup`, `orthogroup_size`, `n_paralogs`, the five per-target-species columns for each
    of the four proteomes, `n_bacterial_orthologs` and `bacterial_panel_size`.
    """
    _check(species)
    return _read(_path(EVIDENCE_DIR, f"orthology_{species}.tsv"))


BACTERIA = ("kpneumoniae", "ecoli", "saureus")


def load_all(species: tuple[str, ...] = BACTERIA) -> pd.DataFrame:
    """The three bacteria's deliverables stacked, with a `species` column first. Human is
    excluded by default because the panel columns are bacterial and it has no deliverable."""
    frames = []
    for sp in species:
        df = load(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_dense_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """Every species' dense table, stacked, with a `species` column first. Human included."""
    frames = []
    for sp in species:
        df = load_dense(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def orthologs_of(accession: str, target_species: str | None = None) -> pd.DataFrame:
    """Every ortholog of one protein, optionally restricted to one target species."""
    df = load_orthologs()
    out = df[df["query_ac"] == accession]
    if target_species is not None:
        _check(target_species)
        out = out[out["target_species"] == target_species]
    return out


def neighbors_of(accession: str, target_species: str | None = None) -> pd.DataFrame:
    """One protein's nearest neighbours, best first."""
    df = load_neighbors()
    out = df[df["query_ac"] == accession]
    if target_species is not None:
        _check(target_species)
        out = out[out["target_species"] == target_species]
    return out.sort_values(["target_species", "rank"])


def load_orthogroups() -> pd.DataFrame:
    """Per-orthogroup composition: size and per-species member counts."""
    return _read(_path(EVIDENCE_DIR, "orthogroups.tsv"))


def load_disagreement() -> pd.DataFrame:
    """Pairs where OrthoFinder and RBH disagree -- the free validation set v1 never looked at."""
    return _read(_path(EVIDENCE_DIR, "method_disagreement.tsv"))


def control() -> pd.DataFrame:
    """This run's numbers beside v1's, for every pre-registered quantity above."""
    return _read(_path(EVIDENCE_DIR, "control.tsv"))


def load_orthodb(species: str) -> pd.DataFrame:
    """The ABSOLUTE orthogroup table for one species -- one row per protein, every accession present.

    `orthodb_verdict` is always one of `ORTHODB_VERDICTS`, never empty: `assigned_by_uniprot` (the protein
    is in OrthoDB, group read directly), `assigned_by_sequence` (the anchor, group inherited from a
    within-species match at >=95% identity, with `orthodb_transfer_pident` recording it), or
    `no_group` -- **a measured zero.**

    Unlike OrthoFinder's `orthogroup`, these ids do not depend on which species were in our run.

    **EMPTY STRING IS NOT A VALUE, and `notna()` alone will lie to you here.** `orthodb_og_domain`
    reads 5,728 non-null on Kp, but 416 of those are `""` -- a protein with no group -- so the real
    coverage is 5,312 (92.7%, the documented figure) and not 100%. `pd.NA` written to TSV returns
    as an empty string, and under `string` dtype that is a perfectly valid non-null entry. The same
    defect was found and fixed in `src/ligandability._coerce_precedents` for `exact_target`, so it
    is a property of this repo's TSV round-trips rather than a one-off. The group columns are
    nulled here; `orthodb_verdict` is the column to filter on in any case, since it says WHY a
    group is absent.
    """
    _check(species)
    d = _read(_path(ORTHOLOGY_DIR, f"orthodb_{species}.tsv"))
    for c in ("orthodb_og_domain", "orthodb_og_narrow", "orthodb_og_name"):
        if c in d.columns:
            d[c] = d[c].astype("string").replace("", pd.NA)
    return d


def load_orthodb_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """Every species' OrthoDB table, stacked, with a `species` column first."""
    frames = []
    for sp in species:
        df = load_orthodb(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_orthodb_long() -> pd.DataFrame:
    """Every (protein, orthogroup, level) assignment -- all 1,003 levels, not just the pivoted ones."""
    return _read(_path(SCRATCH_DIR, "orthodb_groups_long.tsv"))


def manifest() -> pd.DataFrame:
    """What the run did: per species counts, parameters, timings."""
    return _read(_path(EVIDENCE_DIR, "manifest.tsv"))
