"""Stage 08a — export per-protein prioritization data for the browser web app.

Joins the finished per-axis result CSVs (essentiality + ligandability + the ESM-C
projection and structural annotation) on the canonical ``uniprot_accession`` key,
curates a compact column set, derives 0–1 composite *component* values, and writes
one JSON per organism to ``app/data/{kp,ec}.json`` for the static SPA under ``app/``.

Design notes:
- Defensive: a desired column is emitted only if it exists on disk. New axes/columns
  can be added later without breaking older data (the front-end greys out what's absent).
- Composite components are the *available* 0–1 signals only — Essentiality, Ligandability,
  Broad-spectrum breadth, Human-selectivity. Novelty/Degradability/Expression have no data
  yet and are intentionally not emitted (the front-end shows them as "coming soon").
- The cross-axis weighting itself is done live in the browser; this script only exposes
  the normalized per-component values so the weights can be tuned interactively.

Run with the ``gradi`` conda env interpreter.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.ligandability import ORGANISMS, REPO_ROOT  # noqa: E402

RESULTS = REPO_ROOT / "output" / "results"
PROCESSED = REPO_ROOT / "data" / "processed"
OUT_DIR = REPO_ROOT / "app" / "data"

# Desired columns per source table. Intersected with what actually exists (defensive).
ESS_COLS = [
    "uniprot_accession", "gene", "selectivity",
    "essentiality_score", "essentiality_tier",
    "evidence_experimental", "evidence_transfer", "evidence_predictor",
    "experimentally_essential", "evidence_sources",
    "entero_pct_essential", "bacteria_pct_essential",
    "proteomelm_ess_score", "geptop_score", "fba_essential",
    "ec_transfer_essential", "n_ecoli_orthologs",
]
LIG_COLS = [
    "uniprot_accession",
    "ligandability_score", "ligandability_tier",
    "evidence_binding", "evidence_structural", "evidence_pocket",
    "has_hard_evidence", "disorder_frac", "human_ligandable_family",
    "chembl_any_n_compounds", "chembl_any_n_potent", "chembl_any_best_pchembl",
    "chembl_direct_organism",
    "bindingdb_any_n_potent",
    "pdb_lig_any_has_druglike", "pdb_lig_direct_pdb_ids", "pdb_lig_seqdirect_pdb_ids",
    "alphafill_best_ligand", "alphafill_n_druglike",
    "fpocket_max_drug_score", "p2rank_top_score", "pocket_consensus_score",
]
PROJ_COLS = ["uniprot_accession", "tsne_x", "tsne_y", "family"]
STRUCT_COLS = [
    "uniprot_accession", "af_available", "af_mean_plddt",
    "af_n_domains", "af_is_multidomain", "af_cif_url",
]
PDB_COLS = [
    "uniprot_accession", "pdb_has_structure", "pdb_n_structures",
    "pdb_ids", "pdb_best_resolution_A", "pdb_best_method", "pdb_has_holo",
]
POP_COLS = [
    "uniprot_accession", "popularity_score", "popularity_tier",
    "own_n_pubs", "own_reviewed", "n_orthologs",
    "best_homolog_gene", "best_homolog_organism", "best_homolog_n_pubs",
]
FAM_IP_COLS = [
    "uniprot_accession", "interpro_family_names", "interpro_superfamily_names",
    "interpro_domain_names", "pfam_ids", "n_interpro_entries",
]
FAM_PANTHER_COLS = [
    "uniprot_accession", "panther_family_names", "panther_subfamily_names",
]

# --- task-agnostic functional class (heuristic keyword map) -----------------
# Ordered: first matching class wins. Keys are matched (substring, lowercased)
# against the protein's combined InterPro/PANTHER family text + gene symbol.
FUNCTIONAL_CLASS_RULES = [
    ("ribosomal_translation", ["ribosom", "aminoacyl-trna", "trna synth", "trna ligase",
                               "elongation factor", "translation initiation", "release factor",
                               "trna-", "rrna", "peptidyl-trna"]),
    ("dna_replication_repair", ["dna polymerase", "dna-directed dna", "helicase", "topoisomerase",
                                "gyrase", "primase", "recombinase", "dna repair", "excinuclease",
                                "mismatch repair", "chromosome partition", "reca", "single-strand"]),
    ("signaling", ["two-component", "histidine kinase", "response regulator", "chemotaxis",
                   "sensor", "diguanylate", "cyclic-di", "signal transduction", "ggdef", "protein kinase"]),
    ("transcription_regulation", ["transcription", "sigma factor", "sigma-", "rna polymerase",
                                  "dna-directed rna", "repressor", "activator", "transcriptional regul",
                                  "helix-turn-helix", "winged helix", "hth", "lysr", "tetr", "arac",
                                  "regulator"]),
    ("transport", ["transport", "permease", "major facilitator", "mfs", "abc transporter",
                   "porin", "channel", "efflux", "symporter", "antiporter", "tonb", "secretion",
                   "translocase", "importer", "exporter"]),
    ("cell_envelope", ["outer membrane", "cell wall", "peptidoglycan", "murein", "lipopolysaccharide",
                       "lipoprotein", "cell division", "septum", "ftsz", "pilus", "fimbria", "flagell",
                       "capsule", "lps", "membrane"]),
    ("oxidoreductase", ["oxidoreductase", "dehydrogenase", "reductase", "oxidase", "cytochrome",
                        "ferredoxin", "peroxidase", "catalase", "nad(p)", "nadh", "flavo", "monooxygenase",
                        "dioxygenase"]),
    ("transferase", ["transferase", "synthase", "synthetase", "methyltransferase", "acyltransferase",
                     "glycosyltransferase", "aminotransferase", "phosphoribosyltransferase", "kinase"]),
    ("hydrolase_protease", ["hydrolase", "protease", "peptidase", "lipase", "esterase", "nuclease",
                            "amidase", "glycosidase", "phosphatase", "phosphodiesterase", "atpase",
                            "thioesterase", "deacetylase"]),
    ("lyase_isomerase_ligase", ["lyase", "isomerase", "ligase", "carboxylase", "aldolase", "dehydratase",
                                "mutase", "epimerase", "racemase", "hydratase", "decarboxylase", "cyclase"]),
    ("uncharacterized", ["duf", "uncharacter", "hypothetical", "not named", "unknown function",
                         "domain of unknown"]),
]
FUNCTIONAL_CLASS_TEXT_COLS = [
    "interpro_superfamily_names", "interpro_family_names", "panther_family_names",
    "interpro_domain_names", "gene",
]


def conservation_scores(organism: str) -> pd.Series:
    """0–1 cross-species ortholog spread: distinct non-human species with an
    ortholog / panel size. Presence-based (not essentiality). Index = anchor accession."""
    _pid, prefix = ORGANISMS[organism]
    path = REPO_ROOT / "data" / "processed" / "other" / "orthology" / f"{prefix}_orthologs_long.tsv"
    if not path.exists():
        print(f"  ! missing (conservation skipped): {path.name}")
        return pd.Series(dtype=float)
    d = pd.read_csv(path, sep="\t", usecols=["anchor_uniprot", "tier", "species"])
    d = d[(d["tier"] != "human") & (d["species"] != "Homo_sapiens")]
    panel = d["species"].nunique()
    if panel == 0:
        return pd.Series(dtype=float)
    spread = d.groupby("anchor_uniprot")["species"].nunique() / panel
    print(f"  conservation panel size: {panel} species")
    return spread.clip(0.0, 1.0)


def human_closeness_scores(organism: str) -> pd.Series:
    """0–1 closeness to the nearest human ortholog = max(pident)/100 over human-tier
    ortholog hits. 0 (via fillna downstream) = no human ortholog = maximally selective."""
    _pid, prefix = ORGANISMS[organism]
    path = REPO_ROOT / "data" / "processed" / "other" / "orthology" / f"{prefix}_orthologs_long.tsv"
    if not path.exists():
        print(f"  ! missing (human_closeness skipped): {path.name}")
        return pd.Series(dtype=float)
    d = pd.read_csv(path, sep="\t", usecols=["anchor_uniprot", "tier", "species", "pident"])
    d = d[(d["tier"] == "human") | (d["species"] == "Homo_sapiens")]
    d["pident"] = pd.to_numeric(d["pident"], errors="coerce")
    d = d.dropna(subset=["pident"])
    if d.empty:
        return pd.Series(dtype=float)
    return (d.groupby("anchor_uniprot")["pident"].max() / 100.0).clip(0.0, 1.0)


def _functional_class(row) -> str:
    parts = []
    for c in FUNCTIONAL_CLASS_TEXT_COLS:
        v = row.get(c)
        if isinstance(v, str) and v:
            parts.append(v.lower())
    text = " ; ".join(parts)
    if not text.strip():
        return "uncharacterized"
    for cls, kws in FUNCTIONAL_CLASS_RULES:
        if any(k in text for k in kws):
            return cls
    return "other"

# Columns to coerce to real booleans (CSV stores them as "True"/"False"/blank).
BOOL_COLS = {
    "experimentally_essential", "has_hard_evidence", "human_ligandable_family",
    "fba_essential", "ec_transfer_essential", "af_available", "af_is_multidomain",
    "pdb_has_structure", "pdb_has_holo",
}

# Composite components exposed to the weight sliders. Each maps to a 0–1 value.
COMPONENT_KEYS = ["comp_essentiality", "comp_ligandability", "comp_breadth",
                  "comp_human_selective", "comp_novelty"]


def _read_subset(path: Path, cols: list[str]) -> pd.DataFrame:
    """Read a CSV keeping only the desired columns that actually exist."""
    if not path.exists():
        print(f"  ! missing (skipped): {path.name}")
        return pd.DataFrame(columns=["uniprot_accession"])
    df = pd.read_csv(path, low_memory=False)
    keep = [c for c in cols if c in df.columns]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        print(f"  · {path.name}: absent cols {missing}")
    return df[keep].copy()


def _to_bool(series: pd.Series) -> pd.Series:
    def conv(v):
        if pd.isna(v):
            return None
        if isinstance(v, bool):
            return v
        s = str(v).strip().lower()
        if s in ("true", "1", "1.0", "yes", "t"):
            return True
        if s in ("false", "0", "0.0", "no", "f", ""):
            return False
        return None
    return series.map(conv)


def _human_selective(sel) -> float | None:
    """1.0 if the target has no meaningful human ortholog (safe), 0.0 if it does."""
    if pd.isna(sel):
        return None
    s = str(sel)
    if s.endswith("_selective"):
        return 1.0
    if s.endswith("_human_homolog"):
        return 0.0
    return None


def _clip01(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").clip(0.0, 1.0)


def _clean(v):
    """JSON-safe scalar: NaN/NaT -> None, numpy -> python, round floats."""
    if v is None:
        return None
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return None
        return round(v, 4)
    if isinstance(v, (pd.Timestamp,)):
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(v, "item"):  # numpy scalar
        return v.item()
    return v


def build_organism(organism: str) -> dict:
    pid, prefix = ORGANISMS[organism]
    rdir = RESULTS / organism
    print(f"\n[{organism}] prefix={prefix}")

    ess = _read_subset(rdir / f"{prefix}_essentiality.csv", ESS_COLS)
    lig = _read_subset(rdir / f"{prefix}_ligandability.csv", LIG_COLS)
    proj = _read_subset(rdir / f"{prefix}_esmc600m_projection.csv", PROJ_COLS)
    struct = _read_subset(rdir / f"{prefix}_alphafold_structure.csv", STRUCT_COLS)
    pdb = _read_subset(rdir / f"{prefix}_pdb_coverage.csv", PDB_COLS)
    pop = _read_subset(PROCESSED / organism / "bibliometric" / f"{prefix}_popularity.csv", POP_COLS)
    famdir = PROCESSED / organism / "families"
    ipro = _read_subset(famdir / f"{prefix}_interpro_annotation.csv", FAM_IP_COLS)
    panther = _read_subset(famdir / f"{prefix}_panther_annotation.csv", FAM_PANTHER_COLS)

    if ess.empty or "essentiality_score" not in ess.columns:
        raise SystemExit(f"FATAL: essentiality table for {organism} is missing/empty")

    df = ess
    for other in (lig, proj, struct, pdb, pop, ipro, panther):
        if not other.empty:
            df = df.merge(other, on="uniprot_accession", how="left")

    # task-agnostic functional class (heuristic)
    df["functional_class"] = df.apply(_functional_class, axis=1)

    # Overview display scores (not composite components):
    # Structure = pocket druggability (alias of evidence_pocket for independent labelling)
    if "evidence_pocket" in df.columns:
        df["structure_score"] = _clip01(df["evidence_pocket"])
    # Conservation = cross-species ortholog spread (presence-based, essentiality-independent)
    cons = conservation_scores(organism)
    df["conservation_score"] = df["uniprot_accession"].map(cons).fillna(0.0).clip(0.0, 1.0)
    # Human closeness (0–1): identity to nearest human ortholog; 0 = no human ortholog.
    hc = human_closeness_scores(organism)
    df["human_closeness"] = df["uniprot_accession"].map(hc).fillna(0.0).clip(0.0, 1.0)

    # Coerce booleans.
    for c in BOOL_COLS & set(df.columns):
        df[c] = _to_bool(df[c])

    # ---- composite components (0–1) --------------------------------------
    df["comp_essentiality"] = _clip01(df["essentiality_score"])
    if "ligandability_score" in df.columns:
        df["comp_ligandability"] = _clip01(df["ligandability_score"])
    breadth = pd.to_numeric(df.get("entero_pct_essential"), errors="coerce")
    if "bacteria_pct_essential" in df.columns:
        breadth = breadth.fillna(pd.to_numeric(df["bacteria_pct_essential"], errors="coerce"))
    df["comp_breadth"] = breadth.clip(0.0, 1.0)
    if "selectivity" in df.columns:
        df["comp_human_selective"] = df["selectivity"].map(_human_selective)
    # Novelty / neglectedness = inverse of bibliometric studiedness (dark = novel = 1).
    if "popularity_score" in df.columns:
        df["comp_novelty"] = (1.0 - _clip01(df["popularity_score"])).clip(0.0, 1.0)

    # Fallback display name.
    df["name"] = df["gene"].where(df["gene"].notna() & (df["gene"].astype(str) != ""),
                                  df["uniprot_accession"])

    # ---- coverage report -------------------------------------------------
    n = len(df)
    print(f"  rows: {n}")
    for c in COMPONENT_KEYS:
        if c in df.columns:
            cov = int(df[c].notna().sum())
            print(f"  coverage {c}: {cov}/{n} ({100*cov/n:.0f}%)")
        else:
            print(f"  coverage {c}: ABSENT")
    for c in ("structure_score", "conservation_score", "human_closeness"):
        if c in df.columns:
            nz = int((pd.to_numeric(df[c], errors="coerce") > 0).sum())
            print(f"  {c}: mean={df[c].mean():.2f}, >0 for {nz}/{n} ({100*nz/n:.0f}%)")
    if "interpro_family_names" in df.columns:
        fam_cov = int(df["interpro_family_names"].notna().sum())
        print(f"  coverage interpro_family_names: {fam_cov}/{n} ({100*fam_cov/n:.0f}%)")
    print("  functional_class distribution:")
    for cls, cnt in df["functional_class"].value_counts().items():
        print(f"    {cls:24s} {cnt:5d} ({100*cnt/n:.0f}%)")

    # ---- emit ------------------------------------------------------------
    # Drop columns that are dead weight in the payload (unused by the web app and
    # trivially reconstructable): the AlphaFold CIF URL is a template over the accession.
    df = df.drop(columns=[c for c in ("af_cif_url",) if c in df.columns])

    # Columnar wire format: emit rows as arrays of cell values in `columns` order
    # (not per-row objects), which removes the repeated key strings from the JSON.
    # The client rehydrates back to row-objects on load, so nothing else changes.
    columns = list(df.columns)
    rows = [[_clean(v) for v in rec] for rec in df.itertuples(index=False, name=None)]
    return {
        "organism": organism,
        "prefix": prefix,
        "proteome_id": pid,
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "n": n,
        "format": "columnar",
        "columns": columns,
        "components": COMPONENT_KEYS,
        "rows": rows,
    }


def pocket_resids(organism: str) -> dict:
    """Top (rank-1) P2Rank pocket residue numbers per protein, for the 3D viewer.

    Parses the raw P2Rank predictions (data/processed/<org>/pockets/p2rank_run/
    <ACC>.pdb_predictions.csv); the first data row is the highest-scoring pocket.
    residue_ids look like "A_324 A_399 ..." → we keep the integer residue numbers,
    which match the AlphaFold model's numbering (UniProt positions)."""
    base = PROCESSED / organism / "pockets" / "p2rank_run"
    out: dict[str, list[int]] = {}
    if not base.exists():
        return out
    for csvf in base.glob("*.pdb_predictions.csv"):
        acc = csvf.name[: -len(".pdb_predictions.csv")]
        try:
            with open(csvf) as fh:
                rdr = csv.reader(fh)
                hdr = [h.strip() for h in next(rdr)]
                ri = hdr.index("residue_ids")
                first = next(rdr, None)
                if not first or len(first) <= ri:
                    continue
                nums = []
                for tok in first[ri].split():
                    t = tok.strip().split("_")[-1]
                    if t.lstrip("-").isdigit():
                        nums.append(int(t))
                if nums:
                    out[acc] = nums
        except Exception:
            continue
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for organism, (_pid, prefix) in ORGANISMS.items():
        payload = build_organism(organism)
        out = OUT_DIR / f"{prefix}.json"
        out.write_text(json.dumps(payload, ensure_ascii=False, separators=(",", ":")))
        size_mb = out.stat().st_size / 1e6
        print(f"  -> wrote {out.relative_to(REPO_ROOT)} ({size_mb:.1f} MB)")
        # compact top-pocket residue map (separate lazy-loaded file for the 3D viewer)
        pk = pocket_resids(organism)
        pkout = OUT_DIR / f"pockets_{prefix}.json"
        pkout.write_text(json.dumps(pk, separators=(",", ":")))
        print(f"  -> wrote {pkout.relative_to(REPO_ROOT)} "
              f"({pkout.stat().st_size / 1e6:.1f} MB, {len(pk)} pockets)")
    print("\nDone.")


if __name__ == "__main__":
    main()
