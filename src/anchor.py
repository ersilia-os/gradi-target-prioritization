"""Load the HS11286 anchor proteome and the E. coli K-12 reference proteome
from the pre-staged UniProt TSV dumps in `data/raw/`.

The TSV format is `Entry<tab>Gene Names<tab>Sequence`. `Gene Names` is a
space-separated list that always contains the locus tag (KPHS_NNNNN for Kp,
b#### / JW#### for E. coli) plus optionally a gene symbol and synonyms.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

KP_LOCUS_RX = re.compile(r"(KPHS_[0-9p]+)")
EC_BNUM_RX = re.compile(r"\b(b\d{4})\b")
EC_JW_RX = re.compile(r"\b(JW\d+(?:\.\d+)?)\b")
GENE_SYMBOL_RX = re.compile(r"^[a-z][a-zA-Z0-9_]{2,5}$")


@dataclass
class ProteomeRow:
    uniprot: str
    locus_tag: str | None
    gene_symbol: str | None
    synonyms: list[str]
    sequence: str


def _split_names(names: str) -> list[str]:
    return [tok for tok in (names or "").split() if tok]


def _first_gene_symbol(tokens: list[str], locus_rx: re.Pattern[str]) -> str | None:
    """Pick the first token that looks like a canonical gene symbol — lowercase
    initial letter, 3–6 chars, and not a locus tag / Blattner / Keio ID."""
    for tok in tokens:
        if locus_rx.search(tok):
            continue
        if EC_BNUM_RX.search(tok) or EC_JW_RX.search(tok):
            continue
        if GENE_SYMBOL_RX.match(tok):
            return tok
    return None


def load_kp_anchor(path: str | Path) -> pd.DataFrame:
    """Return a one-row-per-KPHS-protein DataFrame.

    Columns: uniprot, kp_locus_tag, kp_gene_symbol, synonyms, sequence.
    Plasmid-encoded proteins (KPHS_pNNNNNN) are kept and flagged via the locus
    tag itself; the pipeline downstream can filter chromosomal-only if desired.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"Entry": "uniprot", "Gene Names": "_names", "Sequence": "sequence"})
    df["_tokens"] = df["_names"].fillna("").apply(_split_names)
    df["kp_locus_tag"] = df["_tokens"].apply(
        lambda toks: next((m.group(1) for tok in toks if (m := KP_LOCUS_RX.search(tok))), None)
    )
    df["kp_gene_symbol"] = df["_tokens"].apply(lambda toks: _first_gene_symbol(toks, KP_LOCUS_RX))
    df["synonyms"] = df["_tokens"].apply(
        lambda toks: [t for t in toks if not KP_LOCUS_RX.search(t) and not GENE_SYMBOL_RX.match(t)]
    )
    df["chromosomal"] = df["kp_locus_tag"].fillna("").str.match(r"^KPHS_\d+$")
    return df[
        ["uniprot", "kp_locus_tag", "kp_gene_symbol", "synonyms", "chromosomal", "sequence"]
    ]


def load_ec_reference(path: str | Path) -> pd.DataFrame:
    """Return one-row-per-MG1655-protein DataFrame.

    Columns: uniprot, ec_b_number, ec_jw_number, ec_gene_symbol, synonyms, sequence.
    """
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"Entry": "uniprot", "Gene Names": "_names", "Sequence": "sequence"})
    df["_tokens"] = df["_names"].fillna("").apply(_split_names)
    df["ec_b_number"] = df["_tokens"].apply(
        lambda toks: next((m.group(1) for tok in toks if (m := EC_BNUM_RX.search(tok))), None)
    )
    df["ec_jw_number"] = df["_tokens"].apply(
        lambda toks: next((m.group(1) for tok in toks if (m := EC_JW_RX.search(tok))), None)
    )
    df["ec_gene_symbol"] = df["_tokens"].apply(
        lambda toks: _first_gene_symbol(toks, EC_BNUM_RX)
    )
    df["synonyms"] = df["_tokens"].apply(
        lambda toks: [
            t for t in toks
            if not EC_BNUM_RX.search(t) and not EC_JW_RX.search(t) and not GENE_SYMBOL_RX.match(t)
        ]
    )
    return df[
        ["uniprot", "ec_b_number", "ec_jw_number", "ec_gene_symbol", "synonyms", "sequence"]
    ]
