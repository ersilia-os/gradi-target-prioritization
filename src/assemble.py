"""Assemble the v1 K. pneumoniae target annotation table.

Output schema (one row per HS11286 KPHS_* locus tag):

    kp_locus_tag, kp_gene_symbol, uniprot, product, chromosomal,
    plfam_id, pgfam_id, has_plfam,

    ess_in_vitro_call,        ess_in_vitro_score,        ess_in_vitro_sources,
    ess_in_vivo_lung_call,    ess_in_vivo_lung_score,    ess_in_vivo_lung_sources,
    ess_in_vivo_urine_call,   ess_in_vivo_urine_score,   ess_in_vivo_urine_sources,
    ess_in_vivo_serum_call,   ess_in_vivo_serum_score,   ess_in_vivo_serum_sources,
    ess_vulnerability_call,   ess_vulnerability_score,   ess_vulnerability_sources,
    ess_Ec_inferred_call,     ess_Ec_inferred_via,       ess_Ec_inferred_sources,

    clp_degradability_score,  clp_degradability_tier,    degron_feature_score,
    cterm_ssra_like,          nterm_destabilizing,
    ecoli_clp_trapped,        ecoli_halflife_class,

    notes

Joins are made by gene symbol because cross-strain locus tags (ecl8_*, VK055_*,
KPNIH1_*) do not share a common ID space. Symbol-based join is imperfect; the
notes column records ambiguity (same symbol matching multiple anchor rows).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from .anchor import load_kp_anchor, load_ec_reference
from .conservation import load_bvbrc_features
from .degradability import HEADLINE_COLUMNS as CLP_HEADLINE_COLUMNS
from .degradability import load_clp_degradability
from .essentiality import load_all

KP_PROTEOME = Path("data/raw/klebsiella_pneumoniae_proteome.tsv")
EC_PROTEOME = Path("data/raw/escherichia_coli_proteome.tsv")

FLAVORS_DIRECT = {
    "in_vitro_essential": "ess_in_vitro",
    "in_vivo_lung": "ess_in_vivo_lung",
    "in_vivo_urine": "ess_in_vivo_urine",
    "in_vivo_serum": "ess_in_vivo_serum",
    "vulnerability_crispri": "ess_vulnerability",
}

CALL_PRIORITY = ["essential", "fitness_defect", "unclear", "fitness_advantage", "non_essential"]


def _consensus_call(calls: list[str]) -> str:
    seen = [c for c in calls if c]
    if not seen:
        return ""
    for c in CALL_PRIORITY:
        if c in seen:
            return c
    return seen[0]


def _join_flavor_block(anchor: pd.DataFrame, evidence: pd.DataFrame, prefix: str) -> pd.DataFrame:
    """Collapse multi-source evidence for one flavor onto the anchor by gene symbol."""
    if evidence.empty:
        anchor[f"{prefix}_call"] = ""
        anchor[f"{prefix}_score"] = pd.NA
        anchor[f"{prefix}_sources"] = ""
        return anchor

    ev = evidence.dropna(subset=["gene_symbol"]).copy()
    ev["gene_symbol"] = ev["gene_symbol"].str.lower()
    grouped = (
        ev.groupby("gene_symbol")
        .agg(
            call=("call", lambda s: _consensus_call(list(s))),
            score=("score", lambda s: ",".join(str(round(float(x), 3)) for x in s if pd.notna(x))),
            sources=("source", lambda s: ",".join(sorted(set(s)))),
        )
        .reset_index()
    )
    a = anchor.copy()
    a["_sym"] = a["kp_gene_symbol"].fillna("").str.lower()
    a = a.merge(grouped, how="left", left_on="_sym", right_on="gene_symbol")
    a = a.drop(columns=["_sym", "gene_symbol"])
    a = a.rename(
        columns={
            "call": f"{prefix}_call",
            "score": f"{prefix}_score",
            "sources": f"{prefix}_sources",
        }
    )
    a[f"{prefix}_call"] = a[f"{prefix}_call"].fillna("")
    a[f"{prefix}_sources"] = a[f"{prefix}_sources"].fillna("")
    return a


def _join_ec_inferred(
    anchor: pd.DataFrame, ec_ref: pd.DataFrame, ec_evidence: pd.DataFrame
) -> pd.DataFrame:
    """For each KPHS row with a gene symbol that matches an E. coli reference protein
    flagged essential, set ess_Ec_inferred_call='essential'."""
    if ec_evidence.empty:
        anchor["ess_Ec_inferred_call"] = ""
        anchor["ess_Ec_inferred_via"] = ""
        anchor["ess_Ec_inferred_sources"] = ""
        return anchor

    ev = ec_evidence[ec_evidence["flavor"] == "in_vitro_essential_Ec"].copy()
    # gene_symbol in evidence is E. coli's; match against E. coli reference symbol set
    ec_essential_symbols = {
        s.lower() for s in ev["gene_symbol"].dropna().astype(str) if s
    } & {s.lower() for s in ec_ref["ec_gene_symbol"].dropna().astype(str) if s}

    sources = ",".join(sorted(ev["source"].dropna().unique()))
    a = anchor.copy()
    sym = a["kp_gene_symbol"].fillna("").str.lower()
    a["ess_Ec_inferred_call"] = sym.isin(ec_essential_symbols).map(
        {True: "essential", False: ""}
    )
    a["ess_Ec_inferred_via"] = a["ess_Ec_inferred_call"].apply(
        lambda v: "gene-symbol match (Ec ortholog inferred)" if v else ""
    )
    a["ess_Ec_inferred_sources"] = a["ess_Ec_inferred_call"].apply(
        lambda v: sources if v else ""
    )
    return a


def build_annotation() -> pd.DataFrame:
    anchor = load_kp_anchor(KP_PROTEOME)
    # Keep one row per locus tag (collapse plasmid duplicates that lack locus_tag)
    anchor = anchor.dropna(subset=["kp_locus_tag"]).drop_duplicates(subset=["kp_locus_tag"])
    ec_ref = load_ec_reference(EC_PROTEOME)
    bvbrc = load_bvbrc_features()
    evidence = load_all()

    anchor = anchor.merge(bvbrc, how="left", on="kp_locus_tag")

    for flavor, prefix in FLAVORS_DIRECT.items():
        sub = evidence[evidence["flavor"] == flavor]
        anchor = _join_flavor_block(anchor, sub, prefix)

    anchor = _join_ec_inferred(anchor, ec_ref, evidence)

    clp = load_clp_degradability()
    if not clp.empty:
        sym = anchor["kp_gene_symbol"].fillna("").str.lower()
        anchor = anchor.assign(_sym=sym).merge(
            clp, how="left", left_on="_sym", right_on="gene_symbol"
        ).drop(columns=["_sym", "gene_symbol"])
    else:
        for c in CLP_HEADLINE_COLUMNS:
            anchor[c] = pd.NA

    cols = [
        "kp_locus_tag", "kp_gene_symbol", "uniprot", "product", "chromosomal",
        "plfam_id", "pgfam_id", "has_plfam",
        "ess_in_vitro_call", "ess_in_vitro_score", "ess_in_vitro_sources",
        "ess_in_vivo_lung_call", "ess_in_vivo_lung_score", "ess_in_vivo_lung_sources",
        "ess_in_vivo_urine_call", "ess_in_vivo_urine_score", "ess_in_vivo_urine_sources",
        "ess_in_vivo_serum_call", "ess_in_vivo_serum_score", "ess_in_vivo_serum_sources",
        "ess_vulnerability_call", "ess_vulnerability_score", "ess_vulnerability_sources",
        "ess_Ec_inferred_call", "ess_Ec_inferred_via", "ess_Ec_inferred_sources",
        *CLP_HEADLINE_COLUMNS,
    ]
    for c in cols:
        if c not in anchor.columns:
            anchor[c] = ""
    return anchor[cols].sort_values("kp_locus_tag").reset_index(drop=True)


def main() -> None:
    out_dir = Path("output/results")
    out_dir.mkdir(parents=True, exist_ok=True)
    df = build_annotation()
    df.to_csv(out_dir / "kp_target_annotation_v1.csv", index=False)
    try:
        df.to_parquet(out_dir / "kp_target_annotation_v1.parquet", index=False)
    except Exception as e:
        print(f"parquet write skipped: {e}")
    print(f"wrote {len(df)} rows to {out_dir / 'kp_target_annotation_v1.csv'}")
    # Tiny summary for log readability
    for col in [
        "ess_in_vitro_call",
        "ess_in_vivo_urine_call",
        "ess_in_vivo_serum_call",
        "ess_vulnerability_call",
        "ess_Ec_inferred_call",
    ]:
        nonblank = (df[col] != "").sum()
        print(f"  {col}: {nonblank} non-blank rows")


if __name__ == "__main__":
    main()
