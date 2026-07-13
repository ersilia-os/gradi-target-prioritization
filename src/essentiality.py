"""Shared helpers for the essentiality assessment stage (docs §4, scripts/07*).

Centralises what the 07a–07l scripts all need. Deliberately thin: it **reuses** the ligandability
helper module (`src/ligandability.py`) for everything generic — organism config, proteome loading,
the 03a ortholog table, the DIAMOND-by-sequence engine, and the output-path helpers — so the two
stages share exactly one implementation of those. On top of that it adds the essentiality-specific
pieces:

  * locus-tag / gene-symbol / b-number bridges built from the proteome TSV (HS11286 carries
    `KPHS_*` locus tags; E. coli K-12 carries Blattner `b####` and Keio `JW####` ids);
  * `map_strain_by_sequence()` — DIAMOND the anchor proteome against a locus-tag-keyed strain
    proteome (KPPR1/KPNIH1/NJST258/ECL8) and keep the best hit per anchor above an identity floor,
    so essentiality calls made on those strains can be transferred onto HS11286 by sequence;
  * the tunable scoring constants + the `essentiality_ladder` used by the 07h merge.

Keying is always by UniProt accession (project convention). Run with the `gradi` conda env.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

# Reuse the ligandability helpers verbatim (single source of truth for the generic machinery).
from src import ligandability as L  # noqa: E402

# Re-export the bits the 07* scripts lean on, so they can `from src import essentiality as E`
# and reach everything through one module.
ORGANISMS = L.ORGANISMS
run_diamond_blastp = L.run_diamond_blastp
load_diamond_hits = L.load_diamond_hits
acc_from_header = L.acc_from_header
load_accessions = L.load_accessions
load_genes = L.load_genes
load_orthologs = L.load_orthologs
ortholog_map = L.ortholog_map
orthologs_long_tsv = L.orthologs_long_tsv
proteome_fasta = L.proteome_fasta
proteome_tsv = L.proteome_tsv
results_dir = L.results_dir
processed_dir = L.processed_dir

# Human-readable organism names for plot titles (matches the ligandability plot convention).
ORG_DISPLAY: dict[str, str] = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli K-12"}

# Locus-tag prefixes seen in each anchor proteome's Gene Names column.
LOCUS_PREFIX: dict[str, str] = {"kpneumoniae": "KPHS_", "ecoli": "b"}


# --------------------------------------------------------------------------- paths
def essentiality_raw_dir(organism: str, *parts: str) -> Path:
    d = REPO_ROOT / "data" / "raw" / organism / "essentiality"
    for p in parts:
        d = d / p
    d.mkdir(parents=True, exist_ok=True)
    return d


def essentiality_processed_dir(organism: str, *parts: str) -> Path:
    return L.processed_dir(organism, "essentiality", *parts)


def legacy_essentiality_dir(*parts: str) -> Path:
    """The pre-reorganization essentiality material (ECL8 xlsx, sources.tsv, curated highlights)."""
    d = REPO_ROOT / "data" / "raw" / "legacy" / "essentiality"
    for p in parts:
        d = d / p
    return d


# --------------------------------------------------------------------------- id bridges
def _iter_gene_tokens(organism: str):
    """Yield (accession, token) for every whitespace token in the proteome Gene Names column."""
    tsv = L.proteome_tsv(organism)
    df = pd.read_csv(tsv, sep="\t", usecols=["Entry", "Gene Names"])
    for acc, gn in zip(df["Entry"], df["Gene Names"].fillna("")):
        acc = str(acc).strip()
        if not acc:
            continue
        for tok in str(gn).split():
            yield acc, tok


def locus_to_uniprot(organism: str) -> dict[str, str]:
    """Locus-tag -> UniProt accession, from the proteome TSV Gene Names column.

    HS11286 uses `KPHS_*`; E. coli K-12 uses Blattner `b####`. Case-insensitive keys.
    """
    prefix = LOCUS_PREFIX[organism]
    out: dict[str, str] = {}
    if organism == "ecoli":
        pat = re.compile(r"^b\d{4}$")
        for acc, tok in _iter_gene_tokens(organism):
            if pat.match(tok):
                out.setdefault(tok.upper(), acc)
    else:
        for acc, tok in _iter_gene_tokens(organism):
            if tok.startswith(prefix):
                out.setdefault(tok.upper(), acc)
    return out


def gene_to_uniprot(organism: str) -> dict[str, str]:
    """Lowercased gene symbol -> UniProt accession (first-seen wins).

    The gene symbol is the first Gene-Names token that is NOT a locus tag / ordered-locus id
    (i.e. not KPHS_*, not b####, not JW####). Used as the fallback bridge when a source table
    is keyed only by gene name (e.g. the ECL8 no-Prokka fallback, curated CRISPRi highlights).
    """
    locus_like = re.compile(r"^(KPHS_|b\d{4}$|JW)", re.IGNORECASE)
    out: dict[str, str] = {}
    tsv = L.proteome_tsv(organism)
    df = pd.read_csv(tsv, sep="\t", usecols=["Entry", "Gene Names"])
    for acc, gn in zip(df["Entry"], df["Gene Names"].fillna("")):
        acc = str(acc).strip()
        for tok in str(gn).split():
            if tok and not locus_like.match(tok):
                out.setdefault(tok.lower(), acc)
                break
    return out


def gene_aliases_to_uniprot(organism: str) -> dict[str, str]:
    """Lowercased gene symbol (EVERY alias, not just the first) -> UniProt accession.

    Published screens key on the primary gene name OR a synonym; this maps all Gene-Names symbol
    tokens (excluding locus ids KPHS_/b####/JW####) so screen tables land robustly. First-seen wins.
    """
    locus_like = re.compile(r"^(KPHS_|KPN_|VK055|b\d{4}$|JW)", re.IGNORECASE)
    out: dict[str, str] = {}
    df = pd.read_csv(L.proteome_tsv(organism), sep="\t", usecols=["Entry", "Gene Names"])
    for acc, gn in zip(df["Entry"], df["Gene Names"].fillna("")):
        acc = str(acc).strip()
        for tok in str(gn).split():
            if tok and not locus_like.match(tok):
                out.setdefault(tok.lower(), acc)
    return out


def jw_to_uniprot(organism: str = "ecoli") -> dict[str, str]:
    """Keio JW clone id (e.g. JW0001) -> UniProt accession, from the proteome Gene Names column.

    JW ids carry a `.N` suffix in some UniProt entries (e.g. `JW1527.1`); both the exact token and
    the suffix-stripped form are indexed so a bare `JW1527` from a screen still resolves.
    """
    jw = re.compile(r"^JW\d+", re.IGNORECASE)
    out: dict[str, str] = {}
    for acc, tok in _iter_gene_tokens(organism):
        if jw.match(tok):
            out.setdefault(tok.upper(), acc)
            out.setdefault(tok.split(".")[0].upper(), acc)
    return out


# --------------------------------------------------------------------------- E. coli -> Kp transfer
def transfer_ecoli_to_kp(ec_value_by_acc: dict[str, float] | pd.Series,
                         reduce: str = "max") -> dict[str, float]:
    """Lift an E.-coli-accession-keyed signal onto K. pneumoniae HS11286 anchors via orthology.

    Uses the 03a `kp_orthologs_long` table restricted to species `Ecoli_K12_MG1655`
    (anchor_uniprot = HS11286, target_uniprot = E. coli K-12). For each Kp anchor with one or more
    E. coli orthologs, `reduce` ("max" | "mean") aggregates the ortholog values. Returns
    {kp_accession: value}. ~3,100+ Kp proteins have an E. coli ortholog.
    """
    ec = dict(ec_value_by_acc) if not isinstance(ec_value_by_acc, pd.Series) else ec_value_by_acc.to_dict()
    orth = load_orthologs("kpneumoniae")
    sub = orth[orth["species"] == "Ecoli_K12_MG1655"][["anchor_uniprot", "target_uniprot"]].copy()
    sub["val"] = sub["target_uniprot"].map(ec)
    sub = sub.dropna(subset=["val"])
    if sub.empty:
        return {}
    agg = sub.groupby("anchor_uniprot")["val"].max() if reduce == "max" else sub.groupby("anchor_uniprot")["val"].mean()
    return agg.to_dict()


# --------------------------------------------------------------------------- sequence transfer
# Identity/coverage floor for transferring a per-gene essentiality call from another strain's
# protein onto the HS11286 anchor by DIAMOND best hit. Same-species strains sit far above this;
# the floor guards against spurious low-identity hits (mirrors the ligandability TRANSFER_FLOOR).
STRAIN_MIN_PIDENT = 50.0
STRAIN_MIN_QCOV = 70.0


def map_strain_by_sequence(
    organism: str,
    strain_faa: Path,
    out_tsv: Path,
    min_pident: float = STRAIN_MIN_PIDENT,
    min_qcov: float = STRAIN_MIN_QCOV,
) -> pd.DataFrame:
    """DIAMOND the anchor proteome against a locus-tag-keyed strain proteome; best hit per anchor.

    `strain_faa` must be a FASTA whose headers are the strain's locus tags (e.g. `>VK055_01234`).
    Returns a DataFrame with columns: uniprot_accession, strain_locus, pident, qcov, bitscore —
    one row per anchor protein that has a hit above the floor. The raw DIAMOND hits are cached at
    `out_tsv` (idempotent via `run_diamond_blastp`).
    """
    L.run_diamond_blastp(L.proteome_fasta(organism), strain_faa, out_tsv)
    hits = L.load_diamond_hits(out_tsv)
    if hits.empty:
        return pd.DataFrame(columns=["uniprot_accession", "strain_locus", "pident", "qcov", "bitscore"])
    hits["uniprot_accession"] = hits["qseqid"].map(L.acc_from_header)
    hits = hits[(hits["pident"] >= min_pident) & (hits["qcovhsp"] >= min_qcov)]
    hits = hits.sort_values("bitscore", ascending=False).drop_duplicates("uniprot_accession")
    out = hits[["uniprot_accession", "sseqid", "pident", "qcovhsp", "bitscore"]].rename(
        columns={"sseqid": "strain_locus", "qcovhsp": "qcov"}
    )
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------- scoring (07h)
# All weights/thresholds are module-level so the composite blend is auditable + tunable, exactly
# like src/ligandability is consumed by 06g. Unavailable tracks are dropped with renormalized
# weights in 07h (not zero-filled), so a gated/deferred track lowers confidence, not the score.
W_EXPERIMENTAL = 0.40   # direct Kp Tn-seq / CRISPRi essential evidence (strongest)
W_ECOLI_TRANSFER = 0.20  # E. coli essential set lifted via ortholog
W_PREDICTOR = 0.40       # consensus of the computational predictors
# within the predictor consensus, ProteomeLM (the LM SOTA) is weighted above the classical methods
PRED_WEIGHTS = {"proteomelm": 0.5, "geptop": 0.3, "fba": 0.2}

TIER_ESSENTIAL = 0.60    # composite at/above -> essential (unless overridden by a hard call)
TIER_LIKELY = 0.35       # composite at/above -> likely_essential
GEPTOP_CUTOFF = 0.24     # Geptop 2.0 default essential-score cutoff
FBA_GROWTH_CUTOFF = 0.01  # single-gene KO growth ratio below this -> in-silico essential


def essentiality_ladder(call: str | None, score: float | None = None) -> float:
    """Map a categorical experimental essentiality call to a [0,1] sub-score.

    essential -> 1.0 ; conditional/vulnerable -> 0.7 ; unclear/intermediate -> 0.4 ;
    non-essential -> 0.0 ; unknown/missing -> 0.0. An explicit numeric `score` (already in [0,1])
    overrides the categorical mapping when provided.
    """
    if score is not None and not pd.isna(score):
        return float(min(1.0, max(0.0, score)))
    c = (call or "").strip().lower()
    if c in ("essential", "ess"):
        return 1.0
    if c in ("conditional", "conditionally_essential", "vulnerable", "vulnerability"):
        return 0.7
    if c in ("unclear", "intermediate", "reduced_fitness"):
        return 0.4
    return 0.0
