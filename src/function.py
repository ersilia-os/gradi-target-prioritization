"""Load the stage-02 functional annotation.

Stage 02 assigns each bacterial protein broad functional classes, under two schemes:

    data/processed/function/cog_{kpneumoniae,ecoli,saureus}.tsv      COG2024 category letter
    data/processed/function/goslim_{...}.tsv                         goslim_prokaryote terms
    data/processed/function/eggnog_{...}.tsv                         eggNOG-mapper, full output

COG says *which orthologous group*; GO slim says *what it does*. Both come from established tools
(COGclassifier / RPS-BLAST vs NCBI's CDD, and eggNOG-mapper v2), and **neither reaches 100%**.
That is deliberate: a protein no established tool can annotate gets an empty row rather than a
guess. Roughly a quarter of each proteome is short, uncharacterised and carries no Pfam or
InterPro either -- nothing annotates it.

Use `load(species)` for the two classifications side by side, `load_eggnog()` for the rest of what
eggNOG-mapper produced (KEGG KO, EC, BRITE, CAZy, PFAMs, preferred names).

One row per protein, keyed on `uniprot_ac`, every protein present. A protein RPS-BLAST found no COG
for keeps its row with empty strings -- absence of a category is a fact about the protein, not a
missing record, and dropping it would silently shrink every downstream join.

Human is absent by construction, not by omission: COG2024 is 2,296 genomes of bacteria and archaea
with no eukaryotes, so there is no honest COG letter to give it.

Two columns are easy to confuse:

    cog_category      exactly one letter -- the broad class. Use this one.
    cog_category_all  the whole string. COG4862 is `KTN`; COG orders those letters by importance
                      and `cog_category` is the first. 9-16% of classified proteins are multi-letter,
                      so `cog_category` is a real choice being made, not a formality.

`R` (general function prediction only) and `S` (function unknown) are categories, but they are not
answers -- `informative()` drops them.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
FUNCTION_DIR = REPO_ROOT / "data" / "processed" / "function"
EVIDENCE_DIR = FUNCTION_DIR / "evidence"
SCRATCH_DIR = FUNCTION_DIR / "scratch"
# The three bacteria. Human is not classifiable here -- see the module docstring.
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# Categories that are classified but uninformative.
POORLY_CHARACTERIZED = ("R", "S")


def load_cog(species: str) -> pd.DataFrame:
    """One species' COG table: uniprot_ac, cog_id, cog_category, cog_category_all, cog_group, ..."""
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = EVIDENCE_DIR / f"cog_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/cog.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_cog_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column added back."""
    frames = []
    for sp in species:
        df = load_cog(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_cog_counts(species: str) -> pd.DataFrame:
    """Per-letter counts for one species: all 26 categories, zeros included."""
    path = EVIDENCE_DIR / f"cog_counts_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/cog.py first")
    return pd.read_csv(path, sep="\t", dtype={"n": int}, keep_default_na=False)


def load_goslim(species: str) -> pd.DataFrame:
    """One species' GO slim table, keyed on `uniprot_ac`, 100% covered.

    Three aspects, each with one chosen term (`goslim_mf`), its label (`goslim_mf_name`) and the
    full multi-label set (`goslim_mf_all`). A protein can have MF and no BP -- absence of an aspect
    is a fact, not a missing value.

    `goslim_source` is `curated` (UniProt's own GO) or `eggnog` (eggNOG-mapper v2, orthology-based
    transfer). Rows with neither are empty -- see the module docstring.

    `eggnog` rows are a tool's inference, not curation: `eggnog_evalue` and `eggnog_score` carry
    the evidence, and `confident()` filters on them.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = EVIDENCE_DIR / f"goslim_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/goslim.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def load_goslim_terms() -> pd.DataFrame:
    """The 97-term goslim_prokaryote vocabulary: go_id, name, aspect, depth."""
    path = EVIDENCE_DIR / "goslim_terms.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/goslim.py first")
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def load(species: str) -> pd.DataFrame:
    """`function_<species>.tsv` — the axis deliverable. One row per protein, four columns:

        uniprot_ac  cog_categories  goslim_terms

    Both term columns are `;`-joined lists in vocabulary order; **an empty string means the
    protein carries no term**, which is a complete row, not a missing one. Split with
    `.str.split(";")` after dropping empties.

    **A zero-length list means NOT ANNOTATED, not "the function was ruled out."** On
    K. pneumoniae 1,506 proteins (26.3%) carry no GO-slim term because nothing is known about
    them.

    **No evidence column, for either scheme** (owner's call, 2026-10-03). COG never had one to
    carry — `cogclassifier` is 1:1 with "has a category", measured on all three species. GO-slim's
    separated only **322 proteins out of 13,020** and lives on byte-identically as `evidence` in
    `evidence/goslim_matrix_<species>.tsv` and as `goslim_source` in `evidence/goslim_<species>.tsv`.
    **So this table does not say whether a GO term is UniProt-curated or inferred from an eggNOG
    orthogroup** — read one of those two files before treating a term as curated.

    **What this form CANNOT say, and the matrices can.** A packed list cannot distinguish a term
    that is merely unannotated from one the organism structurally cannot reach: 8 GO-slim terms
    are eukaryote/plant concepts, *S. aureus* has 19 unreachable because it is Gram-positive, and
    COG `Y` is nuclear structure. Those are **kept columns** in
    `evidence/{goslim,cog}_matrix_<species>.tsv` and that is where structural zeros live. Use
    `load_goslim_matrix()` / `load_cog_matrix()` for anything that needs a feature matrix, or that
    needs to tell "impossible" from "unknown".

    **Multi-label, deliberately** — built from the `*_all` columns: 34–52% of annotated proteins
    carry more than one slim term (max 11 on Kp), 12.4–12.6% of classified proteins more than one
    COG letter.

    The packed table and the matrices are **provably interchangeable**: `function/matrix.py`
    asserts the matrices round-trip to the long-form source AND that the packed columns re-expand
    to the matrices exactly.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = FUNCTION_DIR / f"function_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/matrix.py first")
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def load_long(species: str) -> pd.DataFrame:
    """Both schemes' long-form SOURCE tables, joined on `uniprot_ac` — names, e-values, OGs."""
    return load_cog(species).merge(load_goslim(species), on="uniprot_ac", how="outer")


def load_eggnog(species: str) -> pd.DataFrame:
    """Everything eggNOG-mapper produced for one species, keyed on `uniprot_ac`.

    Beyond `gos` (which feeds the GO-slim table) this carries `kegg_ko`, `kegg_pathway`,
    `kegg_module`, `brite`, `ec`, `cazy`, `pfams`, `preferred_name`, `description` and
    `cog_category` -- a second, independent opinion on the COG letter. Empty strings where
    eggNOG-mapper found no orthologous group.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = EVIDENCE_DIR / f"eggnog_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/eggnog.py first")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def confident(df: pd.DataFrame, max_evalue: float = 1e-10) -> pd.DataFrame:
    """Curated rows, plus eggNOG rows whose seed-ortholog hit was at least this good.

    eggNOG-mapper already applies its own thresholds, so this is a second, stricter cut for
    analyses that want curation-grade evidence only. Lower it knowingly.
    """
    ev = pd.to_numeric(df.get("eggnog_evalue"), errors="coerce")
    return df[(df["goslim_source"] == "curated") | (ev <= max_evalue)]


def informative(df: pd.DataFrame) -> pd.DataFrame:
    """Drop the unclassified rows and the R/S ones -- what is left is an actual functional call."""
    return df[(df["cog_category"] != "") & (~df["cog_category"].isin(POORLY_CHARACTERIZED))]


# ---------------------------------------------------------------- the matrices (the deliverable)

GOSLIM_MATRIX_TERMS = 97
COG_MATRIX_LETTERS = 26


def load_goslim_matrix(species: str) -> pd.DataFrame:
    """`goslim_matrix_<species>.tsv` — one row per protein, 97 binary GO-slim columns, + `evidence`.

    **Complete by construction**: every protein in the proteome has a row. A protein with no
    annotation is an **all-zero row with `evidence == "none"`**, never a missing row.

    Three things to hold:

    1. **A zero means NOT ANNOTATED, not absent.** On K. pneumoniae 1,506 proteins (26.3%) are
       all-zero because nothing is known about them — not because the function was ruled out. Read
       `evidence` (`curated` | `eggnog` | `none`) before treating a zero as a measured negative.
    2. **Multi-label, deliberately.** Built from the `*_all` columns, not the single chosen term:
       34–41% of MF-annotated and 32–52% of BP-annotated proteins carry more than one slim term
       (max 7). `goslim_<aspect>` in the source table keeps only one and is not what this is.
    3. **Some columns are always zero, and that is information.** All 97 terms are present in every
       species file so the three stack; 8 never occur anywhere (thylakoid, photosynthesis,
       extracellular matrix, protein tag, histone binding, ECM organization, nutrient reservoir,
       cell adhesion mediator — all eukaryote/plant), and *S. aureus* has 19 always-zero columns
       because it is Gram-positive. Dropping them would give the species different shapes.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = EVIDENCE_DIR / f"goslim_matrix_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/matrix.py first")
    return pd.read_csv(path, sep="\t")


def load_cog_matrix(species: str) -> pd.DataFrame:
    """`cog_matrix_<species>.tsv` — one row per protein, 26 binary COG-letter columns, + `evidence`.

    Same contract as `load_goslim_matrix`: complete, multi-label, zeros mean *not annotated*
    (`evidence` is `cogclassifier` | `none`; 1,200 Kp / 684 Ec / 777 Sa are unclassified).

    **Multi-label matters here too**: 12.4–12.6% of classified proteins carry more than one letter
    (max 4). Reconciled against the independently-computed `evidence/cog_counts_<species>.tsv`,
    which counts only the CHOSEN letter — chosen + extra-from-multi-label == matrix total, exactly
    (Kp 4,528 + 591 = 5,119).

    **`Y` (nuclear structure) is always zero** — no bacterium has one. Column order follows the
    vendored `data/source/cdd/cog_func_category.tsv`, which groups related categories
    (information storage → cellular processes → metabolism → poorly characterised) rather than
    sorting alphabetically. It is vendored precisely so the schema is not defined by a pip install.
    """
    if species not in SPECIES:
        raise ValueError(f"unknown species {species!r}; expected one of {SPECIES}")
    path = EVIDENCE_DIR / f"cog_matrix_{species}.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/matrix.py first")
    return pd.read_csv(path, sep="\t")


def load_goslim_matrix_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column. Safe because every file has the same 97 columns."""
    return _stack(load_goslim_matrix, species)


def load_cog_matrix_all(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """All species stacked, with a `species` column. Safe because every file has the same 26 columns."""
    return _stack(load_cog_matrix, species)


def _stack(loader, species: tuple[str, ...]) -> pd.DataFrame:
    frames = []
    for sp in species:
        df = loader(sp)
        df.insert(0, "species", sp)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    return out


def matrix_manifest() -> pd.DataFrame:
    """Per species and matrix: rows, columns, annotated, all-zero rows, dead columns, evidence mix."""
    path = EVIDENCE_DIR / "matrix_manifest.tsv"
    if not path.exists():
        raise FileNotFoundError(f"{path} -- run scripts/function/matrix.py first")
    return pd.read_csv(path, sep="\t")
