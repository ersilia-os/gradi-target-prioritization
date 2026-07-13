"""Popularity / studiedness annotation (docs/01 §1.4) by orthology transfer.

The HS11286 proteome is bibliometrically dark (0.1% reviewed; UniProt pub counts flat at ~1 =
the genome paper), so direct literature signals are uninformative. Instead we transfer
"studiedness" from each protein's **best-characterized ortholog**: the orthology mapping from
03a links ~89% of Kp proteins to orthologs across 24 species, including the curation hubs
(E. coli K-12 is 94% reviewed). The best-studied ortholog (most curated literature) carries the
knowledge about the protein's function.

UniProt-only signals (no Europe PMC). Proteins with no ortholog -> `dark` (own annotation only).

Per anchor protein the score blends three [0,1] components (weights tunable):
  0.25 * own_annotation_norm
  0.45 * log1p(best_ortholog_pubs) / log1p(P95)
  0.30 * best_ortholog_annotation_norm        (norm = (score-1)/4)

Organism selected with --organism (kpneumoniae default, or ecoli). Output (one CSV):
  data/processed/<organism>/bibliometric/<prefix>_popularity.csv
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import io
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
META_FIELDS = "accession,reviewed,annotation_score,lit_pubmed_id,gene_primary,organism_name"
IDMAP = "https://rest.uniprot.org/idmapping"      # bulk metadata for ortholog targets
W_OWN, W_PUBS, W_ANN = 0.25, 0.45, 0.30          # score weights (sum 1)
TIER_CUTS = (0.33, 0.66)                          # dark < .33 <= studied <= .66 < well_studied

REPO_ROOT = Path(__file__).resolve().parents[1]
ORTHO_DIR = REPO_ROOT / "data" / "processed" / "other" / "orthology"
ORGANISMS = {
    "kpneumoniae": {"proteome": "UP000007841", "prefix": "kp"},
    "ecoli":       {"proteome": "UP000000625", "prefix": "ec"},
}


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def stream_tsv(query: str) -> pd.DataFrame:
    r = requests.get(STREAM, params={"query": query, "format": "tsv", "fields": META_FIELDS}, timeout=300)
    r.raise_for_status()
    return pd.read_csv(io.StringIO(r.text), sep="\t")


def n_pubs(cell) -> int:
    if not isinstance(cell, str) or not cell.strip():
        return 0
    return len([x for x in cell.split(";") if x.strip()])


def normalize_meta(df: pd.DataFrame, key: str = "Entry") -> pd.DataFrame:
    """Standardise a UniProt metadata TSV to: accession, reviewed, annotation, n_pubs, gene, organism.

    `key` is the column holding the accession to index by ("Entry" for a normal query,
    "From" for id-mapping results so it matches the queried ortholog accession).
    """
    out = pd.DataFrame({
        "accession": df[key],
        "reviewed": df["Reviewed"].eq("reviewed"),
        "annotation": pd.to_numeric(df["Annotation"], errors="coerce").fillna(0.0),
        "n_pubs": df["PubMed ID"].apply(n_pubs),
        "gene": df.get("Gene Names (primary)", pd.Series(dtype=str)).fillna(""),
        "organism": df.get("Organism", pd.Series(dtype=str)).fillna(""),
    })
    return out.dropna(subset=["accession"]).drop_duplicates("accession").set_index("accession")


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def submit_idmapping(ids: list[str]) -> str:
    r = requests.post(f"{IDMAP}/run",
                      data={"ids": ",".join(ids), "from": "UniProtKB_AC-ID", "to": "UniProtKB"},
                      timeout=300)
    r.raise_for_status()
    return r.json()["jobId"]


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def _idmap_status(job: str) -> dict:
    r = requests.get(f"{IDMAP}/status/{job}", timeout=120)
    r.raise_for_status()
    return r.json()


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def _idmap_results(job: str) -> str:
    r = requests.get(f"{IDMAP}/uniprotkb/results/stream/{job}",
                     params={"format": "tsv", "fields": META_FIELDS}, timeout=600)
    r.raise_for_status()
    return r.text


def fetch_targets_meta(accs: list[str], cache: Path) -> pd.DataFrame:
    """Bulk UniProt metadata for ortholog targets via id-mapping (one job), cached."""
    if cache.exists():
        print(f"Using cached ortholog metadata: {cache.relative_to(REPO_ROOT)}")
        return pd.read_csv(cache).set_index("accession")
    print(f"Submitting {len(accs)} ortholog targets to UniProt id-mapping ...")
    job = submit_idmapping(accs)
    while _idmap_status(job).get("jobStatus") in ("RUNNING", "NEW"):
        time.sleep(3)
    meta = normalize_meta(pd.read_csv(io.StringIO(_idmap_results(job)), sep="\t"), key="From")
    cache.parent.mkdir(parents=True, exist_ok=True)
    meta.reset_index().to_csv(cache, index=False)
    print(f"  cached {len(meta)} rows to {cache.relative_to(REPO_ROOT)}")
    return meta


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    spec = ORGANISMS[args.organism]
    proteome_id, prefix = spec["proteome"], spec["prefix"]

    out_dir = REPO_ROOT / "data" / "processed" / args.organism / "bibliometric"
    out_csv = out_dir / f"{prefix}_popularity.csv"
    meta_cache = REPO_ROOT / "data" / "raw" / args.organism / "bibliometric" / "ortholog_metadata.csv"

    # --- anchor proteome metadata (own studiedness backbone, 100% coverage) ---
    print(f"Fetching {args.organism} proteome metadata ...")
    own = normalize_meta(stream_tsv(f"proteome:{proteome_id}"))
    print(f"  {len(own)} anchor proteins")

    # --- orthologs + their metadata ---
    long = pd.read_csv(ORTHO_DIR / f"{prefix}_orthologs_long.tsv", sep="\t", low_memory=False)
    long = long[long["target_uniprot"].astype(str) != long["anchor_uniprot"].astype(str)]
    targets = sorted(long["target_uniprot"].dropna().astype(str).unique())
    tmeta = fetch_targets_meta(targets, meta_cache)

    # best-studied ortholog per anchor (max curated pubs, tie-break annotation score)
    j = long.join(tmeta, on="target_uniprot")
    j = j.dropna(subset=["n_pubs"])
    j = j.sort_values(["n_pubs", "annotation"], ascending=False)
    best = j.groupby("anchor_uniprot").first()
    n_orth = long.groupby("anchor_uniprot").size()

    # --- assemble per-protein table over the full proteome ---
    rows = []
    for acc, m in own.iterrows():
        b = best.loc[acc] if acc in best.index else None
        rows.append({
            "uniprot_accession": acc,
            "own_annotation_score": m["annotation"],
            "own_reviewed": bool(m["reviewed"]),
            "own_n_pubs": int(m["n_pubs"]),
            "n_orthologs": int(n_orth.get(acc, 0)),
            "best_homolog_uniprot": b["target_uniprot"] if b is not None else "",
            "best_homolog_organism": b["organism"] if b is not None else "",
            "best_homolog_gene": b["gene"] if b is not None else "",
            "best_homolog_annotation_score": float(b["annotation"]) if b is not None else 0.0,
            "best_homolog_n_pubs": int(b["n_pubs"]) if b is not None else 0,
        })
    df = pd.DataFrame(rows)

    # --- score + tier ---
    p95 = max(np.percentile(df["best_homolog_n_pubs"], 95), 1)
    own_norm = ((df["own_annotation_score"] - 1) / 4).clip(0, 1)
    pubs_norm = (np.log1p(df["best_homolog_n_pubs"]) / np.log1p(p95)).clip(0, 1)
    ann_norm = ((df["best_homolog_annotation_score"] - 1) / 4).clip(0, 1)
    df["popularity_score"] = (W_OWN * own_norm + W_PUBS * pubs_norm + W_ANN * ann_norm).round(4)
    lo, hi = TIER_CUTS
    df["popularity_tier"] = np.where(df["popularity_score"] < lo, "dark",
                                     np.where(df["popularity_score"] <= hi, "studied", "well_studied"))

    cols = ["uniprot_accession", "popularity_score", "popularity_tier",
            "own_annotation_score", "own_reviewed", "own_n_pubs", "n_orthologs",
            "best_homolog_uniprot", "best_homolog_organism", "best_homolog_gene",
            "best_homolog_annotation_score", "best_homolog_n_pubs"]
    df = df[cols]

    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_csv, index=False)

    # --- summary ---
    n = len(df)
    print(f"\n{args.organism}: {n} proteins | with a best ortholog: {(df.n_orthologs>0).sum()} "
          f"({100*(df.n_orthologs>0).mean():.1f}%)")
    print("tiers:", df["popularity_tier"].value_counts().to_dict())
    print("top organisms supplying the best homolog:")
    print(df.loc[df.best_homolog_uniprot != "", "best_homolog_organism"].value_counts().head(6).to_string())
    print(f"\nWrote {out_csv.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
