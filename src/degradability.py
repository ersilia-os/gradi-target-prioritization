"""Loader for the Clp-protease degradability annotation produced by
``scripts/03_annotate_clp_degradability.py``.

The annotation TSV is keyed by UniProt accession but ``src/assemble.py`` joins
the v1 target table by gene symbol (cross-strain locus tags do not share an ID
space). We collapse to one row per gene symbol here; on the rare collisions
(paralogs sharing a symbol) the maximum score wins, which matches how the
assemble layer treats its other gene-symbol-keyed evidence frames.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

CLP_TSV = Path("data/processed/klebsiella_pneumoniae_clp_degradability.tsv")

HEADLINE_COLUMNS = [
    "clp_degradability_score",
    "clp_degradability_tier",
    "degron_feature_score",
    "cterm_ssra_like",
    "nterm_destabilizing",
    "ecoli_clp_trapped",
    "ecoli_halflife_class",
]


def _empty() -> pd.DataFrame:
    return pd.DataFrame(columns=["gene_symbol", *HEADLINE_COLUMNS])


def load_clp_degradability(path: str | Path = CLP_TSV) -> pd.DataFrame:
    """Return one-row-per-gene_symbol DataFrame with the headline columns.

    Score columns are coerced to float; boolean columns are kept as strings
    (matching the TSV) and converted to native bool. Missing TSV → empty
    frame (so a fresh checkout without `data/processed/` still assembles)."""
    p = Path(path)
    if not p.exists():
        return _empty()

    df = pd.read_csv(p, sep="\t")
    if df.empty:
        return _empty()

    df = df[df["gene_symbol"].fillna("").astype(str) != ""].copy()
    df["gene_symbol"] = df["gene_symbol"].astype(str).str.lower()

    for col in ["clp_degradability_score", "degron_feature_score"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["cterm_ssra_like", "nterm_destabilizing", "ecoli_clp_trapped"]:
        df[col] = df[col].astype(str).str.lower() == "true"

    # Collapse paralogs sharing a gene symbol: max-score row wins.
    df = (
        df.sort_values("clp_degradability_score", ascending=False)
        .drop_duplicates(subset=["gene_symbol"], keep="first")
        .reset_index(drop=True)
    )
    return df[["gene_symbol", *HEADLINE_COLUMNS]]
