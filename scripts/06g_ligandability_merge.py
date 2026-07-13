"""Merge ligandability tracks + disorder filter + composite score (docs §2.4 + result).

Joins every per-track table (06a ChEMBL, 06b BindingDB, 06c PDB co-crystals, 06d AlphaFill,
06e pockets, 06f AF2Bind placeholder) plus AlphaFold confidence (04a) on UniProt accession,
applies the disorder filter, and computes a transparent composite ligandability score + tier.

Composite = weighted blend of three normalised, independently-kept sub-scores:
  * evidence_binding    — ChEMBL + BindingDB potency (direct + bacterial orthologs; human is a
                          discounted homolog signal). pChEMBL/pAff >= 6  <=> <= 1 uM.
  * evidence_structural — PDB co-crystal drug-like ligand (direct > ortholog) or AlphaFill
                          drug-like transplant (scaled by donor identity).
  * evidence_pocket     — pLDDT-weighted fpocket/P2Rank consensus, further x disorder_penalty.
Disorder filter (§2.4): disorder_frac = af_frac_low + af_frac_very_low; penalty = 1 - frac.

All sub-scores and weights are kept as columns/constants so the blend is auditable + tunable.
Output: output/results/<org>/<prefix>_ligandability.csv   Run with the `gradi` env.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

# composite weights (tunable) — experimental binding strongest, then structural, then predicted
W_BINDING = 0.45
W_STRUCTURAL = 0.30
W_POCKET = 0.25
HUMAN_DISCOUNT = 0.6  # human-ortholog binding counts, but discounted (homolog, not selective)
TIER_TRACTABLE = 0.60
TIER_PARTIAL = 0.35
POCKET_STRONG = 0.50  # pLDDT-weighted pocket consensus that alone implies a druggable site
POCKET_WEAK = 0.30
STRUCT_STRONG = 0.70  # AlphaFill / homolog co-crystal strength implying a real ligand site
DISORDER_INTRACTABLE = 0.50  # mostly-disordered + no hard evidence -> intractable


def _read(org: str, suffix: str, cols: list[str] | None = None) -> pd.DataFrame | None:
    _, prefix = L.ORGANISMS[org]
    p = L.results_dir(org) / f"{prefix}_{suffix}.csv"
    if not p.exists():
        print(f"  [warn] missing track table {p.name}; its columns will default", flush=True)
        return None
    df = pd.read_csv(p, usecols=cols) if cols else pd.read_csv(p)
    return df.drop_duplicates("uniprot_accession")


def _potency_ladder(best_p, n_potent, n_tested) -> float:
    """tested-not-potent -> 0.35; potent -> 0.7..1.0 by best potency (pAff/pChEMBL 6->9)."""
    if not n_tested or n_tested <= 0:
        return 0.0
    if not n_potent or n_potent <= 0:
        return 0.35
    bp = best_p if best_p is not None and not pd.isna(best_p) else 6.0
    return float(min(1.0, 0.7 + 0.3 * np.clip((bp - 6.0) / 3.0, 0, 1)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = L.ORGANISMS[org]

    base = pd.DataFrame({"uniprot_accession": L.load_accessions(org)})
    genes = L.load_genes(org)
    base["gene"] = base["uniprot_accession"].map(genes)

    tracks = {
        "chembl": _read(org, "chembl"),
        "bindingdb": _read(org, "bindingdb"),
        "pdb": _read(org, "pdb_cocrystals"),
        "alphafill": _read(org, "alphafill"),
        "pockets": _read(org, "pockets"),
        "af2bind": _read(org, "af2bind"),
        "af": _read(org, "alphafold_structure",
                    ["uniprot_accession", "af_mean_plddt", "af_frac_low_plddt", "af_frac_very_low_plddt"]),
    }
    df = base
    for t in tracks.values():
        if t is not None:
            df = df.merge(t, on="uniprot_accession", how="left")

    # ---- disorder filter (§2.4)
    low = df.get("af_frac_low_plddt", pd.Series(0, index=df.index)).fillna(0.0)
    vlow = df.get("af_frac_very_low_plddt", pd.Series(0, index=df.index)).fillna(0.0)
    df["disorder_frac"] = (low + vlow).clip(0, 1)
    df["disorder_penalty"] = (1.0 - df["disorder_frac"]).clip(0, 1)

    # ---- evidence_binding (direct + bacterial primary; human discounted)
    def binding_row(r) -> float:
        prim = max(
            _potency_ladder(r.get("chembl_any_best_pchembl"), r.get("chembl_any_n_potent"), r.get("chembl_any_n_compounds")),
            _potency_ladder(r.get("bindingdb_any_best_paff"), r.get("bindingdb_any_n_potent"), r.get("bindingdb_any_n_compounds")),
        )
        human = max(
            _potency_ladder(r.get("chembl_human_best_pchembl"), r.get("chembl_human_n_potent"), r.get("chembl_human_n_compounds")),
            _potency_ladder(r.get("bindingdb_human_best_paff"), r.get("bindingdb_human_n_potent"), r.get("bindingdb_human_n_compounds")),
        )
        return float(max(prim, HUMAN_DISCOUNT * human))

    df["evidence_binding"] = df.apply(binding_row, axis=1)

    # ---- evidence_structural (PDB co-crystal > ortholog; AlphaFill scaled by identity)
    def structural_row(r) -> float:
        s_pdb = 0.0
        # own structure OR a >=95%-identity PDB chain (same protein, any strain) == strong
        if bool(r.get("pdb_lig_direct_has_druglike")) or bool(r.get("pdb_lig_seqdirect_has_druglike")):
            s_pdb = 1.0
        elif bool(r.get("pdb_lig_ortho_has_druglike")) or bool(r.get("pdb_lig_seqhom_has_druglike")):
            s_pdb = 0.7
        s_af = 0.0
        n_af = r.get("alphafill_n_druglike")
        if n_af is not None and not pd.isna(n_af) and n_af > 0:
            ident = r.get("alphafill_best_identity")
            ident = float(ident) if ident is not None and not pd.isna(ident) else 0.3
            s_af = float(0.4 + 0.5 * np.clip(ident, 0, 1))
        return float(max(s_pdb, s_af))

    df["evidence_structural"] = df.apply(structural_row, axis=1)

    # ---- evidence_pocket (pLDDT-weighted consensus x disorder penalty)
    pcs = df.get("pocket_consensus_score", pd.Series(0.0, index=df.index)).fillna(0.0)
    df["evidence_pocket"] = (pcs * df["disorder_penalty"]).clip(0, 1)

    # ---- composite
    df["ligandability_score"] = (
        W_BINDING * df["evidence_binding"]
        + W_STRUCTURAL * df["evidence_structural"]
        + W_POCKET * df["evidence_pocket"]
    ).round(4)

    # hard-evidence flags
    any_potent = (
        (df.get("chembl_any_n_potent", 0).fillna(0) > 0)
        | (df.get("bindingdb_any_n_potent", 0).fillna(0) > 0)
    )
    # hard structural evidence = own PDB co-crystal OR a >=95%-id PDB co-crystal (same protein)
    direct_cocrystal = df.get("pdb_lig_direct_has_druglike", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    seqdirect_cocrystal = df.get("pdb_lig_seqdirect_has_druglike", pd.Series(False, index=df.index)).fillna(False).astype(bool)
    df["has_hard_evidence"] = any_potent | direct_cocrystal | seqdirect_cocrystal
    df["human_ligandable_family"] = (
        (df.get("chembl_human_n_potent", 0).fillna(0) > 0)
        | (df.get("bindingdb_human_n_potent", 0).fillna(0) > 0)
    )

    def tier_row(r) -> str:
        # mostly-disordered with no hard binding evidence -> intractable regardless of pockets
        if r["disorder_frac"] >= DISORDER_INTRACTABLE and not r["has_hard_evidence"]:
            return "intractable"
        # tractable: concrete binding evidence, a strong homolog ligand site, a strong druggable
        # pocket, or a high composite
        if (
            r["has_hard_evidence"]
            or r["evidence_structural"] >= STRUCT_STRONG
            or r["evidence_pocket"] >= POCKET_STRONG
            or r["ligandability_score"] >= TIER_TRACTABLE
        ):
            return "tractable"
        # partial: a plausible pocket, some structural transplant, any tested compounds, or mid score
        if (
            r["evidence_pocket"] >= POCKET_WEAK
            or r["evidence_structural"] > 0
            or r["evidence_binding"] > 0
            or r["ligandability_score"] >= TIER_PARTIAL
        ):
            return "partial"
        return "intractable"

    df["ligandability_tier"] = df.apply(tier_row, axis=1)

    # selectivity from the 03c three-way categories (broad_selective = prime: in other bacteria,
    # no human ortholog) — lets us emit the prioritized shortlist for presentations
    cats_path = L.REPO_ROOT / "data" / "processed" / "other" / "orthology" / "three_way_protein_categories.tsv"
    if cats_path.exists():
        cats = pd.read_csv(cats_path, sep="\t")
        sel = cats[cats["organism"] == org].set_index("uniprot_accession")["selectivity"]
        df["selectivity"] = df["uniprot_accession"].map(sel)
    else:
        df["selectivity"] = pd.NA

    out = L.results_dir(org) / f"{prefix}_ligandability.csv"
    df.to_csv(out, index=False)

    # prime shortlist: broad-spectrum + human-selective + tractable, best evidence first
    short_cols = [
        "uniprot_accession", "gene", "selectivity", "ligandability_score", "ligandability_tier",
        "has_hard_evidence", "evidence_binding", "evidence_structural", "evidence_pocket",
        "chembl_any_n_potent", "chembl_any_best_pchembl", "chembl_direct_acc", "chembl_direct_organism",
        "bindingdb_any_n_potent", "pdb_lig_seqdirect_has_druglike", "pdb_lig_seqdirect_pdb_ids",
        "alphafill_best_ligand", "fpocket_max_drug_score", "human_ligandable_family", "disorder_frac",
    ]
    short_cols = [c for c in short_cols if c in df.columns]
    shortlist = (
        df[(df["selectivity"] == "broad_selective") & (df["ligandability_tier"] == "tractable")]
        .sort_values(["has_hard_evidence", "ligandability_score"], ascending=False)[short_cols]
    )
    sl_path = L.results_dir(org) / f"{prefix}_ligandability_shortlist.csv"
    shortlist.to_csv(sl_path, index=False)
    print(f"  prime shortlist: {len(shortlist)} broad-selective tractable targets -> {sl_path.name}", flush=True)

    vc = df["ligandability_tier"].value_counts().to_dict()
    print(
        f"[{org}] wrote {out} ({len(df)} proteins, {df.shape[1]} cols)\n"
        f"  tiers: {vc}\n"
        f"  hard evidence: {int(df['has_hard_evidence'].sum())}; "
        f"mean score: {df['ligandability_score'].mean():.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
