"""In-silico single-gene-deletion essentiality by FBA (docs §4.3d).

Flux-balance analysis on a genome-scale metabolic reconstruction: knock out each modelled gene and
ask whether biomass flux collapses. This is an orthogonal, mechanistic essentiality signal covering
the metabolic subset of the proteome.

Per organism we use the matching curated BiGG model:
  * kpneumoniae -> **iYL1228** (Liao 2011; K. pneumoniae MGH 78578; 1,229 genes as `KPN_` locus tags).
    Model genes are mapped onto HS11286 by SEQUENCE: we build a `KPN_`-keyed MGH78578 proteome FASTA
    from UniProt (proteome UP000000265, ordered-locus field) and DIAMOND the HS11286 proteome against
    it (the project's dark-proteome mapping convention).
  * ecoli -> **iML1515** (Monk 2017; E. coli K-12 MG1655; genes as Blattner `b####`), mapped directly
    via the proteome b-number bridge — no DIAMOND needed.

A gene is `fba_essential` if the single-KO growth ratio (KO biomass / wild-type biomass) < 0.01.
Proteins absent from the model are `fba_status="not_in_model"` (FBA only covers metabolism).

Output: output/results/<org>/<prefix>_ess_fba.csv
  uniprot_accession, fba_essential (bool/NA), fba_growth_ratio, mapped_model_gene,
  fba_mapping_pident, fba_status
Run with the `gradi` conda env interpreter (needs `cobra`). Model + KO results are cached.
"""

from __future__ import annotations

import argparse
import io
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

HEADERS = {"User-Agent": "Mozilla/5.0", "Connection": "close"}
MODELS = {
    "kpneumoniae": {"bigg": "iYL1228", "gene_type": "KPN", "mgh_proteome": "UP000000265"},
    "ecoli":       {"bigg": "iML1515", "gene_type": "bnum"},
}
GROWTH_CUTOFF = E.FBA_GROWTH_CUTOFF  # KO/WT biomass ratio below this -> essential


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=2, max=30))
def _download(url: str) -> bytes:
    r = requests.get(url, headers=HEADERS, timeout=(10, 120))
    r.raise_for_status()
    return r.content


def load_model(bigg: str):
    import cobra

    models_dir = E.REPO_ROOT / "data" / "raw" / "other" / "essentiality" / "models"
    models_dir.mkdir(parents=True, exist_ok=True)
    path = models_dir / f"{bigg}.json"
    if not path.exists():
        print(f"downloading BiGG model {bigg} ...", flush=True)
        path.write_bytes(_download(f"http://bigg.ucsd.edu/static/models/{bigg}.json"))
    return cobra.io.load_json_model(str(path))


def single_gene_deletion(model) -> pd.DataFrame:
    """Return per model-gene: growth_ratio, fba_essential. Cached by the caller."""
    import cobra
    from cobra.flux_analysis import single_gene_deletion as sgd

    wt = model.slim_optimize()
    print(f"  wild-type biomass = {wt:.4f}; {len(model.genes)} genes", flush=True)
    res = sgd(model, gene_list=model.genes, processes=1)  # processes=1: robust & fast enough (~1k LPs)
    # cobra returns a df indexed by frozenset({gene_id}); normalise to a single gene id
    res = res.reset_index(drop=True)
    res["gene"] = res["ids"].map(lambda s: next(iter(s)) if len(s) else "")
    res["growth"] = res["growth"].fillna(0.0)
    res["growth_ratio"] = (res["growth"] / wt).clip(lower=0).round(4)
    res["fba_essential"] = res["growth_ratio"] < GROWTH_CUTOFF
    return res[["gene", "growth_ratio", "fba_essential"]]


def build_mgh78578_fasta(proteome_id: str, out_faa: Path) -> Path:
    """KPN_-locus-keyed MGH78578 protein FASTA from UniProt (for the sequence mapping)."""
    if out_faa.exists() and out_faa.stat().st_size > 0:
        return out_faa
    import re

    print(f"fetching MGH78578 proteome {proteome_id} from UniProt ...", flush=True)
    url = (f"https://rest.uniprot.org/uniprotkb/stream?query=proteome:{proteome_id}"
           "&format=tsv&fields=accession,gene_names,sequence")
    txt = _download(url).decode()
    df = pd.read_csv(io.StringIO(txt), sep="\t")
    df.columns = ["accession", "gene_names", "sequence"]
    kpn_re = re.compile(r"\bKPN_\d+\b")  # the iYL1228 locus-tag namespace, embedded in Gene Names
    n = 0
    with open(out_faa, "w") as fh:
        for _, r in df.iterrows():
            m = kpn_re.search(str(r["gene_names"]))
            seq = str(r["sequence"]).strip()
            if m and seq and seq != "nan":
                fh.write(f">{m.group(0)}\n{seq}\n")
                n += 1
    print(f"  wrote {n} KPN_-keyed sequences -> {out_faa.relative_to(E.REPO_ROOT)}", flush=True)
    return out_faa


def map_model_genes(org: str, spec: dict, ko: pd.DataFrame) -> pd.DataFrame:
    """Attach uniprot_accession + mapping confidence to each model gene's KO result."""
    if spec["gene_type"] == "bnum":
        # E. coli: model gene ids are b-numbers -> direct proteome bridge
        b2acc = E.locus_to_uniprot("ecoli")  # {B####: acc}
        ko["uniprot_accession"] = ko["gene"].str.upper().map(b2acc)
        ko["fba_mapping_pident"] = np.where(ko["uniprot_accession"].notna(), 100.0, np.nan)
        return ko.dropna(subset=["uniprot_accession"])
    # K. pneumoniae: sequence-map HS11286 -> MGH78578 (KPN_) via DIAMOND
    faa = build_mgh78578_fasta(spec["mgh_proteome"],
                               E.essentiality_raw_dir(org, "strains") / "mgh78578_KPN.faa")
    hit_tsv = E.essentiality_processed_dir(org, "fba") / "hs11286_vs_mgh78578.tsv"
    m = E.map_strain_by_sequence(org, faa, hit_tsv)  # uniprot_accession, strain_locus(KPN_), pident, ...
    m = m.rename(columns={"strain_locus": "gene", "pident": "fba_mapping_pident"})
    out = ko.merge(m[["uniprot_accession", "gene", "fba_mapping_pident"]], on="gene", how="inner")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    spec = MODELS[org]

    ko_cache = E.essentiality_processed_dir(org, "fba") / f"{prefix}_{spec['bigg']}_ko.csv"
    if ko_cache.exists():
        ko = pd.read_csv(ko_cache)
        print(f"[{org}] using cached KO table {ko_cache.relative_to(E.REPO_ROOT)}", flush=True)
    else:
        t = time.time()
        model = load_model(spec["bigg"])
        ko = single_gene_deletion(model)
        ko.to_csv(ko_cache, index=False)
        print(f"[{org}] single-gene deletion done in {time.time()-t:.1f}s "
              f"({int(ko.fba_essential.sum())} in-silico essential of {len(ko)})", flush=True)

    mapped = map_model_genes(org, spec, ko)
    # best (most-essential, then highest identity) model gene per anchor protein
    mapped["_rank"] = mapped["fba_essential"].astype(int) * 1000 + mapped["fba_mapping_pident"].fillna(0)
    mapped = mapped.sort_values("_rank", ascending=False).drop_duplicates("uniprot_accession")

    accs = E.load_accessions(org)
    base = pd.DataFrame({"uniprot_accession": accs})
    df = base.merge(
        mapped[["uniprot_accession", "fba_essential", "growth_ratio", "gene", "fba_mapping_pident"]],
        on="uniprot_accession", how="left",
    ).rename(columns={"growth_ratio": "fba_growth_ratio", "gene": "mapped_model_gene"})
    df["fba_status"] = np.where(df["mapped_model_gene"].notna(), "in_model", "not_in_model")
    # proteins not in the model: essentiality unknown by FBA -> NA (not False)
    df["fba_essential"] = df["fba_essential"].astype("object").where(df["mapped_model_gene"].notna(), pd.NA)

    out = E.results_dir(org) / f"{prefix}_ess_fba.csv"
    df.to_csv(out, index=False)
    n_model = int((df["fba_status"] == "in_model").sum())
    n_ess = int((df["fba_essential"] == True).sum())  # noqa: E712
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)} ({len(df)} proteins; "
          f"{n_model} in {spec['bigg']}; {n_ess} FBA-essential)", flush=True)


if __name__ == "__main__":
    main()
