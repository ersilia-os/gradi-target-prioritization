"""Parsers for the per-paper essentiality supplementary tables.

Each parser returns a long-form DataFrame with the columns
    [gene_symbol, locus_tag, call, score, source, condition, flavor]
where `flavor` is one of:
    - in_vitro_essential
    - in_vivo_lung / in_vivo_urine / in_vivo_serum
    - conditional_<stressor>
    - vulnerability_crispri
`call` is in {essential, fitness_defect, non_essential, unclear}.
`score` carries the effect size where available (e.g. logFC for Eichelberger
conditional screens) and is None otherwise. `source` is the short paper key
matching `data/raw/essentiality/literature/sources.tsv`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

EICH_DIR = Path("data/raw/essentiality/literature/eichelberger2024_ECL8")
ZHU_HIGHLIGHTS = Path(
    "data/raw/essentiality/literature/zhu2023_kp_crispri/curated_highlights.tsv"
)


# ---------- Eichelberger 2024 (ECL8) -------------------------------------------------

def _eich_call_from_flags(row: pd.Series) -> str:
    if int(row.get("Essential", 0) or 0) == 1:
        return "essential"
    if int(row.get("Unclear", 0) or 0) == 1:
        return "unclear"
    if int(row.get("Non-essential", 0) or 0) == 1:
        return "non_essential"
    return "unclear"


def load_eichelberger_in_vitro() -> pd.DataFrame:
    p = EICH_DIR / "elife-88971-fig1-data1-v1.xlsx"
    df = pd.read_excel(p, sheet_name=0, skiprows=1)
    df["call"] = df.apply(_eich_call_from_flags, axis=1)
    return pd.DataFrame(
        dict(
            gene_symbol=df["Gene_name"].astype(str).str.split("_").str[0],
            locus_tag=df["locus_tag"],
            call=df["call"],
            score=df["Insertion_index_score"],
            source="eichelberger2024",
            condition="LB outgrowth",
            flavor="in_vitro_essential",
        )
    )


def _eich_conditional(file: str, condition: str, flavor: str, q_threshold: float = 0.05) -> pd.DataFrame:
    df = pd.read_excel(EICH_DIR / file, sheet_name=0, skiprows=1)
    df = df.dropna(subset=["locus_tag", "logFC", "q.value"])
    df["call"] = "non_essential"
    sig = df["q.value"] < q_threshold
    df.loc[sig & (df["logFC"] < 0), "call"] = "fitness_defect"
    df.loc[sig & (df["logFC"] > 0), "call"] = "fitness_advantage"
    return pd.DataFrame(
        dict(
            gene_symbol=df["gene_name"].astype(str).str.split("_").str[0],
            locus_tag=df["locus_tag"],
            call=df["call"],
            score=df["logFC"],
            source="eichelberger2024",
            condition=condition,
            flavor=flavor,
        )
    )


def load_eichelberger_urine() -> pd.DataFrame:
    return _eich_conditional(
        "elife-88971-fig4-data1-v1.xlsx",
        condition="pooled human urine vs LB",
        flavor="in_vivo_urine",
    )


def load_eichelberger_serum() -> pd.DataFrame:
    return _eich_conditional(
        "elife-88971-fig6-data1-v1.xlsx",
        condition="pooled human serum vs heat-inactivated control",
        flavor="in_vivo_serum",
    )


# ---------- Zhu 2023 (Kp Mobile-CRISPRi-seq, paywalled — highlights only) -----------

def load_zhu_highlights() -> pd.DataFrame:
    if not ZHU_HIGHLIGHTS.exists():
        return pd.DataFrame(columns=["gene_symbol", "locus_tag", "call", "score", "source", "condition", "flavor"])
    df = pd.read_csv(ZHU_HIGHLIGHTS, sep="\t")
    return pd.DataFrame(
        dict(
            gene_symbol=df["gene_symbol"],
            locus_tag=None,
            call="fitness_defect",
            score=None,
            source="zhu2023",
            condition=df["condition"],
            flavor="vulnerability_crispri",
        )
    )


# ---------- Manual-stage parsers (fire when the user drops the file) ----------------
# These are intentionally minimal stubs; flesh out once the gated supplementary
# tables are downloaded into the expected_local_path locations from sources.tsv.

def load_bachman_lung_if_present() -> pd.DataFrame:
    p = Path("data/raw/essentiality/literature/bachman2015_KPPR1")
    candidates = list(p.glob("mbo003152358sd*.xls*"))
    if not candidates:
        return _empty()
    # Bachman 2015 Data Set S1 carries the per-gene insertion site fitness; column
    # layout is (locus_tag VK055_*, gene, fold-change, p-value, q-value, etc.).
    frames: list[pd.DataFrame] = []
    for f in candidates:
        df = pd.read_excel(f, sheet_name=0)
        cols = {c.lower(): c for c in df.columns}
        lt = cols.get("locus_tag") or cols.get("locus tag") or next(iter(df.columns))
        gn = cols.get("gene") or cols.get("gene_symbol") or lt
        fc = cols.get("log2fc") or cols.get("fold_change") or cols.get("log2_fc") or None
        q = cols.get("q-value") or cols.get("q_value") or cols.get("qvalue") or None
        df = df.rename(columns={lt: "locus_tag", gn: "gene_symbol"})
        if fc:
            df = df.rename(columns={fc: "score"})
        else:
            df["score"] = None
        if q is not None:
            df["call"] = "non_essential"
            df.loc[df[q] < 0.05, "call"] = "fitness_defect"
        else:
            df["call"] = "fitness_defect"  # whole table is the hit list
        frames.append(
            df.assign(
                source="bachman2015",
                condition="mouse pneumonia, 24h",
                flavor="in_vivo_lung",
            )[["gene_symbol", "locus_tag", "call", "score", "source", "condition", "flavor"]]
        )
    return pd.concat(frames, ignore_index=True)


def load_ramage_in_vitro_if_present() -> pd.DataFrame:
    p = Path("data/raw/essentiality/literature/ramage2017_KPNIH1")
    candidates = list(p.glob("JB.00352-17_zjb999094540sd2.xlsx"))
    if not candidates:
        return _empty()
    df = pd.read_excel(candidates[0], sheet_name=0)
    # Ramage 2017 DSet S2 is a sheet of essential gene rows; KPNIH1_* locus tags.
    cols = {c.lower(): c for c in df.columns}
    lt = cols.get("locus_tag") or cols.get("locus tag") or next(iter(df.columns))
    gn = cols.get("gene") or cols.get("gene_symbol") or lt
    df = df.rename(columns={lt: "locus_tag", gn: "gene_symbol"})
    df["call"] = "essential"
    df["score"] = None
    return df.assign(
        source="ramage2017",
        condition="LB agar, MKP103, consensus of Tn-seq + arrayed library",
        flavor="in_vitro_essential",
    )[["gene_symbol", "locus_tag", "call", "score", "source", "condition", "flavor"]]


def load_goodall_Ec_if_present() -> pd.DataFrame:
    p = Path("data/raw/essentiality/literature/goodall2018_Ec_BW25113")
    candidates = list(p.glob("*.xlsx"))
    if not candidates:
        return _empty()
    df = pd.read_excel(candidates[0], sheet_name=0)
    cols = {c.lower(): c for c in df.columns}
    bnum = cols.get("b_number") or cols.get("b#") or cols.get("locus_tag") or next(iter(df.columns))
    gn = cols.get("gene") or cols.get("gene_symbol") or bnum
    df = df.rename(columns={bnum: "locus_tag", gn: "gene_symbol"})
    df["call"] = "essential"
    df["score"] = None
    return df.assign(
        source="goodall2018",
        condition="LB outgrowth, BW25113, TraDIS",
        flavor="in_vitro_essential_Ec",  # tagged Ec — feeds the orthology layer, not the Kp anchor directly
    )[["gene_symbol", "locus_tag", "call", "score", "source", "condition", "flavor"]]


def _empty() -> pd.DataFrame:
    return pd.DataFrame(
        columns=["gene_symbol", "locus_tag", "call", "score", "source", "condition", "flavor"]
    )


def load_all() -> pd.DataFrame:
    frames = [
        load_eichelberger_in_vitro(),
        load_eichelberger_urine(),
        load_eichelberger_serum(),
        load_zhu_highlights(),
        load_bachman_lung_if_present(),
        load_ramage_in_vitro_if_present(),
        load_goodall_Ec_if_present(),
    ]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return _empty()
    return pd.concat(frames, ignore_index=True)
