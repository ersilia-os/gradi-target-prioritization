"""AlphaFold structural annotation of a reference proteome (docs §1.2b).

For every protein we summarise its AlphaFold (AFDB) monomer model: availability, per-residue
confidence (pLDDT), and PAE-derived domain organisation. AFDB holds only single-chain MONOMER
models and does not predict oligomeric state, so monomer-vs-multimer is deliberately out of
scope here (a separate track). Two confidence axes are captured:

  - pLDDT  -> local/per-residue confidence (mean + confidence-bin fractions), read straight
             from the AFDB prediction endpoint (AFDB pre-computes these fractions).
  - PAE    -> global/relative-position confidence; we cluster the Predicted Aligned Error
             matrix into structural domains (the standard `pae_to_domains` graph approach) to
             report domain count, largest-domain fraction, and inter-domain PAE (how uncertain
             the relative orientation of domains is).

Organism selected with --organism (kpneumoniae default, or ecoli). Source: EBI AlphaFold DB
API (https://alphafold.ebi.ac.uk/api/prediction/{accession}). The model (.cif) and PAE
(.json.gz) are cached under data/processed/<organism>/alphafold/ for reuse by later structural
tracks (degron exposure §3, pocket detection §2). Re-running is resumable.

Outputs (keyed by UniProt accession, per CLAUDE.md):
  output/results/<organism>/<prefix>_alphafold_structure.csv
  output/plots/<organism>/<prefix>_alphafold_structure.png

Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
import networkx as nx
import numpy as np
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
# organism -> (proteome file stem, output prefix, display name)
ORGANISMS = {
    "kpneumoniae": ("UP000007841_HS11286", "kp", "K. pneumoniae HS11286"),
    "ecoli": ("UP000000625_EcoliK12", "ec", "E. coli K-12 MG1655"),
}

API = "https://alphafold.ebi.ac.uk/api/prediction"


class Paths:
    """Per-organism input/cache/output locations (derived from --organism)."""

    def __init__(self, organism: str):
        stem, prefix, name = ORGANISMS[organism]
        self.organism, self.prefix, self.name = organism, prefix, name
        self.tsv = REPO_ROOT / "data" / "raw" / organism / "proteome" / f"{stem}.tsv"
        af_dir = REPO_ROOT / "data" / "processed" / organism / "alphafold"
        self.cif = af_dir / "cif"
        self.pae = af_dir / "pae"
        self.pred = af_dir / "pred"
        self.results = REPO_ROOT / "output" / "results" / organism
        self.plots = REPO_ROOT / "output" / "plots" / organism
        self.csv = self.results / f"{prefix}_alphafold_structure.csv"
        self.plot = self.plots / f"{prefix}_alphafold_structure.png"

    def mkdirs(self):
        for d in (self.cif, self.pae, self.pred, self.results, self.plots):
            d.mkdir(parents=True, exist_ok=True)


# PAE -> domains (standard pae_to_domains parameters)
PAE_POWER = 1.0
PAE_CUTOFF = (
    12.0  # Å: ignore residue pairs less certain than this (sparsifies the graph)
)
RESOLUTION = (
    0.3  # Louvain resolution; calibrated so compact low-PAE proteins -> 1 domain
)
# (the 0.05-0.5 range is a stable plateau on this proteome; 1.0 over-fragments single domains)
MIN_DOMAIN_RESIDUES = 10

COLUMNS = [
    "uniprot_accession",
    "af_available",
    "af_model_id",
    "af_model_version",
    "af_modeled_len",
    "af_is_complex",
    "af_mean_plddt",
    "af_frac_very_high_plddt",
    "af_frac_confident_plddt",
    "af_frac_low_plddt",
    "af_frac_very_low_plddt",
    "af_frac_high_plddt",
    "af_pae_max",
    "af_pae_mean",
    "af_n_domains",
    "af_largest_domain_frac",
    "af_interdomain_pae_mean",
    "af_is_multidomain",
    "af_cif_url",
    "af_pae_url",
    "source",
]


def read_proteome(tsv_path: Path) -> list[str]:
    accs: list[str] = []
    with open(tsv_path) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("Entry"):
                accs.append(row["Entry"])
    return accs


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def _get_json(url: str):
    r = requests.get(url, headers={"Connection": "close"}, timeout=(10, 60))
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def _get_bytes(url: str) -> bytes:
    r = requests.get(url, headers={"Connection": "close"}, timeout=(10, 120))
    r.raise_for_status()
    return r.content


def pae_to_domains(pae: np.ndarray) -> list[set[int]]:
    """Partition residues into structural domains by community-detecting the PAE graph.

    Edge weight between residues i,j is 1/pae**power for pairs more certain than the cutoff
    (the symmetric mean PAE is used since PAE is asymmetric). Communities are found with
    Louvain; clusters smaller than MIN_DOMAIN_RESIDUES are dropped.
    """
    n = pae.shape[0]
    sym = (pae + pae.T) / 2.0
    g = nx.Graph()
    g.add_nodes_from(range(n))
    iu, ju = np.triu_indices(n, k=1)
    mask = sym[iu, ju] <= PAE_CUTOFF
    iu, ju = iu[mask], ju[mask]
    w = 1.0 / np.clip(sym[iu, ju], 0.25, None) ** PAE_POWER
    g.add_weighted_edges_from(zip(iu.tolist(), ju.tolist(), w.tolist()))
    communities = nx.community.louvain_communities(
        g, weight="weight", resolution=RESOLUTION, seed=0
    )
    return [c for c in communities if len(c) >= MIN_DOMAIN_RESIDUES]


def domain_metrics(pae: np.ndarray) -> dict:
    domains = pae_to_domains(pae)
    n = pae.shape[0]
    if not domains:
        return {
            "af_n_domains": 0,
            "af_largest_domain_frac": float("nan"),
            "af_interdomain_pae_mean": float("nan"),
            "af_is_multidomain": False,
        }
    largest = max(len(d) for d in domains)
    inter = float("nan")
    if len(domains) >= 2:
        sym = (pae + pae.T) / 2.0
        labels = np.full(n, -1)
        for k, d in enumerate(domains):
            labels[list(d)] = k
        assigned = labels >= 0
        sub = sym[np.ix_(assigned, assigned)]
        lab = labels[assigned]
        cross = lab[:, None] != lab[None, :]
        inter = float(sub[cross].mean()) if cross.any() else float("nan")
    return {
        "af_n_domains": len(domains),
        "af_largest_domain_frac": largest / n,
        "af_interdomain_pae_mean": inter,
        "af_is_multidomain": len(domains) >= 2,
    }


def process(acc: str, paths: Paths) -> dict:
    row = {c: None for c in COLUMNS}
    row["uniprot_accession"] = acc
    row["af_available"] = False

    # 1) prediction metadata (cached)
    pred_cache = paths.pred / f"{acc}.json"
    if pred_cache.exists():
        pred = json.loads(pred_cache.read_text()) or None
    else:
        pred = _get_json(f"{API}/{acc}")
        pred_cache.write_text(json.dumps(pred) if pred is not None else "null")
    if not pred:
        return row  # no AFDB model (e.g. > 2700 aa)
    e = pred[0]

    ver = e["latestVersion"]
    row.update(
        af_available=True,
        af_model_id=e.get("entryId"),
        af_model_version=ver,
        af_modeled_len=int(e["uniprotEnd"]) - int(e["uniprotStart"]) + 1,
        af_is_complex=bool(e.get("isComplex", False)),
        af_mean_plddt=e.get("globalMetricValue"),
        af_frac_very_high_plddt=e.get("fractionPlddtVeryHigh"),
        af_frac_confident_plddt=e.get("fractionPlddtConfident"),
        af_frac_low_plddt=e.get("fractionPlddtLow"),
        af_frac_very_low_plddt=e.get("fractionPlddtVeryLow"),
        af_cif_url=e.get("cifUrl"),
        af_pae_url=e.get("paeDocUrl"),
        source=f"alphafold_db_v{ver}",
    )
    vh = row["af_frac_very_high_plddt"] or 0.0
    cf = row["af_frac_confident_plddt"] or 0.0
    row["af_frac_high_plddt"] = vh + cf  # docs §1.2b: pLDDT > 70

    # 2) cache the model (.cif, plain so downstream tools can read it directly)
    cif_path = paths.cif / f"{e['entryId']}-model_v{ver}.cif"
    if not cif_path.exists() and e.get("cifUrl"):
        cif_path.write_bytes(_get_bytes(e["cifUrl"]))

    # 3) cache PAE (.json.gz) and compute domain metrics
    pae_path = paths.pae / f"{e['entryId']}-pae_v{ver}.json.gz"
    if pae_path.exists():
        with gzip.open(pae_path, "rt") as fh:
            pae_doc = json.load(fh)
    else:
        pae_doc = _get_json(e["paeDocUrl"])
        with gzip.open(pae_path, "wt") as fh:
            json.dump(pae_doc, fh)
    doc = pae_doc[0] if isinstance(pae_doc, list) else pae_doc
    pae = np.asarray(doc["predicted_aligned_error"], dtype=float)
    row["af_pae_max"] = float(doc.get("max_predicted_aligned_error", pae.max()))
    row["af_pae_mean"] = float(pae.mean())
    row.update(domain_metrics(pae))
    return row


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="process only first N (debug)")
    args = ap.parse_args()

    paths = Paths(args.organism)
    paths.mkdirs()

    accs = read_proteome(paths.tsv)
    if args.limit:
        accs = accs[: args.limit]
    n = len(accs)
    print(
        f"Annotating {n} {paths.name} proteins from AlphaFold DB ({args.workers} workers) ..."
    )

    rows: list[dict] = []
    done = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(lambda acc: process(acc, paths), accs):
            rows.append(row)
            done += 1
            if done % 250 == 0 or done == n:
                print(f"  {done}/{n}", flush=True)

    df = (
        pd.DataFrame(rows, columns=COLUMNS)
        .sort_values("uniprot_accession")
        .reset_index(drop=True)
    )

    # ---- verification ----
    n_avail = int(df["af_available"].sum())
    missing = df.loc[~df["af_available"], "uniprot_accession"].tolist()
    assert len(df) == n, f"expected {n} rows, got {len(df)}"
    assert df["uniprot_accession"].is_unique, "duplicate accessions"
    avail = df[df["af_available"]]
    assert avail["af_mean_plddt"].between(0, 100).all(), "pLDDT out of range"
    assert (avail["af_n_domains"] >= 1).all(), "modeled protein with 0 domains"
    print(f"AFDB models: {n_avail}/{n} ({100 * n_avail / n:.1f}%).")
    if missing:
        print(f"  NO MODEL ({len(missing)}): {', '.join(missing)}")
    print(
        f"  mean pLDDT: median={avail.af_mean_plddt.median():.1f}  "
        f"multidomain={int(avail.af_is_multidomain.sum())} "
        f"({100 * avail.af_is_multidomain.mean():.1f}%)"
    )

    df.to_csv(paths.csv, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].hist(avail["af_mean_plddt"], bins=50, color="#2c7fb8")
    axes[0].set_xlabel("mean pLDDT")
    axes[0].set_ylabel("proteins")
    axes[0].set_title("AlphaFold per-protein confidence")
    sc = axes[1].scatter(
        avail["af_modeled_len"],
        avail["af_n_domains"],
        c=avail["af_mean_plddt"],
        cmap="viridis",
        s=8,
        alpha=0.6,
    )
    axes[1].set_xlabel("modeled length (residues)")
    axes[1].set_ylabel("PAE domains")
    axes[1].set_title("Domain organisation vs length")
    fig.colorbar(sc, ax=axes[1], label="mean pLDDT")
    fig.suptitle(
        f"{paths.name} — AlphaFold structural annotation (n={n_avail}/{n} modeled)",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(paths.plot, dpi=200, facecolor="white")
    plt.close(fig)

    print(
        f"Wrote:\n  {paths.csv.relative_to(REPO_ROOT)}\n  {paths.plot.relative_to(REPO_ROOT)}\n"
        f"  cache: {paths.cif.relative_to(REPO_ROOT)}/ , {paths.pae.relative_to(REPO_ROOT)}/"
    )


if __name__ == "__main__":
    main()
