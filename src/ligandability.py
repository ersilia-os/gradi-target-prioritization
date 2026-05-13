"""Ligandability evidence loader — joins ChEMBL and BindingDB per-UniProt
ligand counts to the K. pneumoniae anchor table via OrthoDB ortholog groups.

For each Kp protein we emit two parallel blocks:

- ``lig_kp_*`` — counts that come ONLY from the Kp protein's own UniProt
  accession. Almost always zero (Kp TrEMBL entries are rarely targets of
  published assays directly).
- ``lig_ortho_*`` — counts aggregated across all UniProts in the same
  OrthoDB group (i.e. across all bacterial species OrthoDB tracks for that
  protein family). This is where chemical-matter signal actually surfaces.

Best-pchembl values are aggregated with ``max`` so the ortholog block
reports the most potent reported binder anywhere in the family.

Provenance columns let downstream users see which species contributed:
``ortholog_uniprots``, ``ortholog_species``, ``orthodb_group_ids``. The
``ortholog_uniprots`` and ``ortholog_species`` strings are truncated at 50
items with a ``+N more`` suffix to keep the CSV readable.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


ORTHOLOGS_TSV = Path("data/processed/klebsiella_pneumoniae_orthodb_orthologs.tsv")
CHEMBL_TSV = Path("data/processed/chembl_ligand_counts.tsv")
BINDINGDB_TSV = Path("data/processed/bindingdb_ligand_counts.tsv")

CHEMBL_COUNT_COLS = ["chembl_any", "chembl_10um", "chembl_1um"]
CHEMBL_BEST_COL = "chembl_best_pchembl"
BINDINGDB_COUNT_COLS = ["bindingdb_any", "bindingdb_10um", "bindingdb_1um"]
BINDINGDB_BEST_COL = "bindingdb_best_pchembl"


def load_orthologs(path: Path = ORTHOLOGS_TSV) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=[
                "kp_uniprot",
                "kp_locus_tag",
                "ortholog_uniprot",
                "ortholog_organism_name",
                "orthodb_group_id",
            ]
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return df


def load_chembl(path: Path = CHEMBL_TSV) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=["uniprot"] + CHEMBL_COUNT_COLS + [CHEMBL_BEST_COL]
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for c in CHEMBL_COUNT_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df[CHEMBL_BEST_COL] = pd.to_numeric(df[CHEMBL_BEST_COL], errors="coerce")
    return df[["uniprot"] + CHEMBL_COUNT_COLS + [CHEMBL_BEST_COL]]


def load_bindingdb(path: Path = BINDINGDB_TSV) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(
            columns=["uniprot"] + BINDINGDB_COUNT_COLS + [BINDINGDB_BEST_COL]
        )
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for c in BINDINGDB_COUNT_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df[BINDINGDB_BEST_COL] = pd.to_numeric(df[BINDINGDB_BEST_COL], errors="coerce")
    return df[["uniprot"] + BINDINGDB_COUNT_COLS + [BINDINGDB_BEST_COL]]


def _truncate_list(items: list[str], cap: int = 50) -> str:
    items = [x for x in items if x]
    if len(items) <= cap:
        return ",".join(items)
    return ",".join(items[:cap]) + f",+{len(items) - cap} more"


def _aggregate_per_kp_locus(
    ortho: pd.DataFrame,
    chembl: pd.DataFrame,
    bindingdb: pd.DataFrame,
) -> pd.DataFrame:
    """Compute the ortholog-aggregate block keyed on kp_locus_tag."""
    ligand_by_uniprot = chembl.merge(bindingdb, on="uniprot", how="outer")
    count_cols = CHEMBL_COUNT_COLS + BINDINGDB_COUNT_COLS
    for c in count_cols:
        if c in ligand_by_uniprot.columns:
            ligand_by_uniprot[c] = (
                pd.to_numeric(ligand_by_uniprot[c], errors="coerce").fillna(0).astype(int)
            )
        else:
            ligand_by_uniprot[c] = 0
    for c in [CHEMBL_BEST_COL, BINDINGDB_BEST_COL]:
        if c not in ligand_by_uniprot.columns:
            ligand_by_uniprot[c] = pd.NA

    merged = ortho.merge(
        ligand_by_uniprot,
        how="left",
        left_on="ortholog_uniprot",
        right_on="uniprot",
    )
    for c in count_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0).astype(int)

    grouped = merged.groupby("kp_locus_tag", dropna=False)
    rows = []
    for kp_locus, sub in grouped:
        agg = {f"lig_ortho_{c}": int(sub[c].sum()) for c in count_cols}
        for src, col in [
            ("lig_ortho_chembl_best_pchembl", CHEMBL_BEST_COL),
            ("lig_ortho_bindingdb_best_pchembl", BINDINGDB_BEST_COL),
        ]:
            vals = pd.to_numeric(sub[col], errors="coerce").dropna()
            agg[src] = float(vals.max()) if not vals.empty else math.nan

        with_lig = sub[
            (sub["chembl_any"] > 0) | (sub["bindingdb_any"] > 0)
        ]
        uniprots = with_lig["ortholog_uniprot"].dropna().unique().tolist()
        species = with_lig["ortholog_organism_name"].dropna().unique().tolist()
        groups = sub["orthodb_group_id"].dropna().unique().tolist()
        groups = [g for g in groups if g]
        agg["ortholog_uniprots"] = _truncate_list(uniprots)
        agg["ortholog_species"] = _truncate_list(species)
        agg["orthodb_group_ids"] = ",".join(sorted(set(groups)))
        agg["kp_locus_tag"] = kp_locus
        rows.append(agg)
    return pd.DataFrame(rows)


def _direct_per_kp_locus(
    ortho: pd.DataFrame,
    chembl: pd.DataFrame,
    bindingdb: pd.DataFrame,
) -> pd.DataFrame:
    """Compute the Kp-direct block keyed on kp_locus_tag.

    Direct = the Kp protein's own UniProt accession only (no ortholog rollup).
    """
    direct = ortho.drop_duplicates(subset=["kp_locus_tag", "kp_uniprot"])[
        ["kp_locus_tag", "kp_uniprot"]
    ]
    direct = direct.merge(
        chembl, how="left", left_on="kp_uniprot", right_on="uniprot"
    ).drop(columns=["uniprot"], errors="ignore")
    direct = direct.merge(
        bindingdb, how="left", left_on="kp_uniprot", right_on="uniprot"
    ).drop(columns=["uniprot"], errors="ignore")
    for c in CHEMBL_COUNT_COLS + BINDINGDB_COUNT_COLS:
        direct[c] = pd.to_numeric(direct[c], errors="coerce").fillna(0).astype(int)
    direct[CHEMBL_BEST_COL] = pd.to_numeric(direct[CHEMBL_BEST_COL], errors="coerce")
    direct[BINDINGDB_BEST_COL] = pd.to_numeric(direct[BINDINGDB_BEST_COL], errors="coerce")
    rename = {
        "chembl_any": "lig_kp_chembl_any",
        "chembl_10um": "lig_kp_chembl_10um",
        "chembl_1um": "lig_kp_chembl_1um",
        CHEMBL_BEST_COL: "lig_kp_chembl_best_pchembl",
        "bindingdb_any": "lig_kp_bindingdb_any",
        "bindingdb_10um": "lig_kp_bindingdb_10um",
        "bindingdb_1um": "lig_kp_bindingdb_1um",
        BINDINGDB_BEST_COL: "lig_kp_bindingdb_best_pchembl",
    }
    direct = direct.rename(columns=rename).drop(columns=["kp_uniprot"])
    return direct


def load_ligandability(
    orthologs_path: Path = ORTHOLOGS_TSV,
    chembl_path: Path = CHEMBL_TSV,
    bindingdb_path: Path = BINDINGDB_TSV,
) -> pd.DataFrame:
    """Return one row per ``kp_locus_tag`` with direct + ortholog-aggregate
    ligandability columns + provenance."""
    ortho = load_orthologs(orthologs_path)
    chembl = load_chembl(chembl_path)
    bindingdb = load_bindingdb(bindingdb_path)
    if ortho.empty:
        return pd.DataFrame(columns=["kp_locus_tag"])

    ortho = ortho[ortho["kp_locus_tag"] != ""].copy()
    direct = _direct_per_kp_locus(ortho, chembl, bindingdb)
    ortho_block = _aggregate_per_kp_locus(ortho, chembl, bindingdb)

    merged = direct.merge(ortho_block, how="outer", on="kp_locus_tag")
    out_cols = [
        "kp_locus_tag",
        "lig_kp_chembl_any",
        "lig_kp_chembl_10um",
        "lig_kp_chembl_1um",
        "lig_kp_chembl_best_pchembl",
        "lig_kp_bindingdb_any",
        "lig_kp_bindingdb_10um",
        "lig_kp_bindingdb_1um",
        "lig_kp_bindingdb_best_pchembl",
        "lig_ortho_chembl_any",
        "lig_ortho_chembl_10um",
        "lig_ortho_chembl_1um",
        "lig_ortho_chembl_best_pchembl",
        "lig_ortho_bindingdb_any",
        "lig_ortho_bindingdb_10um",
        "lig_ortho_bindingdb_1um",
        "lig_ortho_bindingdb_best_pchembl",
        "ortholog_uniprots",
        "ortholog_species",
        "orthodb_group_ids",
    ]
    for c in out_cols:
        if c not in merged.columns:
            merged[c] = pd.NA
    int_cols = [c for c in out_cols if c.endswith(("_any", "_10um", "_1um"))]
    for c in int_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").fillna(0).astype(int)
    str_cols = ["ortholog_uniprots", "ortholog_species", "orthodb_group_ids"]
    for c in str_cols:
        merged[c] = merged[c].fillna("").astype(str)
    return merged[out_cols]
