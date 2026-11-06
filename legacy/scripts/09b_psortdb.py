"""Stage 09b — PSORTb 3.0 calls, taken precomputed from PSORTdb (docs §5.1b).

PSORTb 3.0 remains the reference rule-based predictor for bacterial localization, and it is useful
here as an **independent, orthogonal opinion** to cross-check DeepLocPro (09c) against — its
SCL-BLAST / motif / signal / SVM modules share no machinery with a protein language model, so
agreement between the two is real corroboration rather than two views of the same embedding.

We do not *run* PSORTb. The previous version of this script drove `brinkmanlab/psortb_commandline`
under Docker; on Apple Silicon that container hung for ~19 h under Rosetta emulation and never
emitted a single prediction. The Brinkman lab already publishes PSORTb 3.0 results for every RefSeq
bacterial genome through PSORTdb (cPSORTdb), including both of our anchors, so we download the
per-genome table instead: same predictor, same version, zero local compute.

PSORTdb is keyed by RefSeq protein accession (`ref|YP_005221101.1`), so the table is joined onto
UniProt through the UniProt ID-mapping service. That currently resolves cleanly — every mapped
accession lands inside the reference proteome — and the script asserts on the coverage so a silent
degradation shows up as an error rather than as quietly missing rows. If RefSeq ever drifts far
enough that coverage collapses, the fallback is the house pattern for this repo: DIAMOND the anchor
proteome against the genome's RefSeq FASTA and transfer by best hit (`L.run_diamond_blastp`).

Output: output/results/<organism>/<prefix>_loc_psortdb.csv
Run with the `gradi` conda env. Network fetch (PSORTdb + UniProt ID mapping), both cached on disk.
"""
from __future__ import annotations

import argparse
import csv
import io
import sys
import time
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

# PSORTdb internal assembly ids, resolved once from db.psort.org/search (organism name -> id).
# HS11286 = GCF_000240185.1, MG1655 = GCF_000005845.2.
ASSEMBLY_ID = {"kpneumoniae": 446671, "ecoli": 449203}

PSORTDB_DOWNLOAD = "http://db.psort.org/search/results/download"
IDMAP_RUN = "https://rest.uniprot.org/idmapping/run"
IDMAP_STATUS = "https://rest.uniprot.org/idmapping/status/{job}"
IDMAP_RESULTS = "https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job}"

IDMAP_CHUNK = 5000          # UniProt's documented per-job id ceiling
MIN_COVERAGE = 0.50         # below this the RefSeq join is no longer trustworthy


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=3, max=60))
def fetch_psortdb(assembly: int, dest: Path) -> Path:
    """Download the per-genome PSORTb 3.0 table. Needs a (possibly empty) POST body or it 411s."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  using cached {dest.name}", flush=True)
        return dest
    r = requests.post(
        f"{PSORTDB_DOWNLOAD}?assembly={assembly}&id=",
        data={"id": ""},
        headers={"Content-Type": "application/x-www-form-urlencoded"},
        timeout=600,
    )
    r.raise_for_status()
    dest.write_text(r.text)
    return dest


@retry(stop=stop_after_attempt(4), wait=wait_exponential(multiplier=2, min=3, max=60))
def _idmap_chunk(ids: list[str]) -> dict[str, str]:
    job = requests.post(
        IDMAP_RUN,
        data={"from": "RefSeq_Protein", "to": "UniProtKB", "ids": ",".join(ids)},
        timeout=120,
    )
    job.raise_for_status()
    job_id = job.json()["jobId"]

    for _ in range(120):
        st = requests.get(IDMAP_STATUS.format(job=job_id), timeout=60)
        st.raise_for_status()
        body = st.json()
        if body.get("jobStatus") in ("RUNNING", "NEW"):
            time.sleep(3)
            continue
        break

    res = requests.get(IDMAP_RESULTS.format(job=job_id), params={"format": "tsv"}, timeout=300)
    res.raise_for_status()
    rows = csv.DictReader(io.StringIO(res.text), delimiter="\t")
    return {r["From"]: r["Entry"] for r in rows if r.get("Entry")}


def refseq_to_uniprot(ids: list[str], cache: Path) -> dict[str, str]:
    if cache.exists() and cache.stat().st_size > 0:
        df = pd.read_csv(cache, sep="\t")
        print(f"  using cached id mapping ({len(df)} ids)", flush=True)
        return dict(zip(df["refseq"], df["uniprot"]))

    mapping: dict[str, str] = {}
    for i in range(0, len(ids), IDMAP_CHUNK):
        chunk = ids[i:i + IDMAP_CHUNK]
        mapping.update(_idmap_chunk(chunk))
        print(f"  mapped {len(mapping)}/{len(ids)} RefSeq ids", flush=True)
    pd.DataFrame({"refseq": list(mapping), "uniprot": list(mapping.values())}).to_csv(
        cache, sep="\t", index=False)
    return mapping


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ASSEMBLY_ID), default="kpneumoniae")
    ap.add_argument("--refresh", action="store_true", help="ignore the on-disk caches")
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    assembly = ASSEMBLY_ID[org]
    raw = LOC.localization_raw_dir(org, "psortdb")

    tab = raw / f"assembly_{assembly}.tab"
    idcache = raw / f"assembly_{assembly}_refseq2uniprot.tsv"
    if args.refresh:
        tab.unlink(missing_ok=True)
        idcache.unlink(missing_ok=True)

    print(f"[{org}] PSORTdb assembly {assembly}", flush=True)
    fetch_psortdb(assembly, tab)
    rows = list(csv.DictReader(open(tab), delimiter="\t"))
    print(f"[{org}] {len(rows)} PSORTb records", flush=True)

    refseq_ids = [r["SeqID"].split("|")[-1] for r in rows]
    mapping = refseq_to_uniprot(refseq_ids, idcache)

    canon = set(LOC.load_accessions(org))
    recs = []
    for r, rid in zip(rows, refseq_ids):
        acc = mapping.get(rid)
        if not acc or acc not in canon:
            continue
        final = (r.get("Final_Localization") or "").strip()
        score = r.get("Final_Score") or ""
        recs.append({
            "uniprot_accession": acc,
            "psortb_localization": LOC.PSORTB_MAP.get(final.lower().replace(" ", ""), "unknown"),
            "psortb_score": float(score) if score.strip() else None,
            "psortb_raw": final,
            "refseq_accession": rid,
        })

    df = pd.DataFrame(recs).drop_duplicates(subset="uniprot_accession")
    coverage = len(df) / len(canon)
    print(f"[{org}] joined onto {len(df)}/{len(canon)} reference proteins ({coverage:.1%})",
          flush=True)
    if coverage < MIN_COVERAGE:
        raise SystemExit(
            f"[{org}] RefSeq->UniProt coverage {coverage:.1%} is below {MIN_COVERAGE:.0%}. "
            "The accession join is no longer reliable — fall back to DIAMOND best-hit transfer "
            "against the genome's RefSeq FASTA (src.ligandability.run_diamond_blastp)."
        )

    out = LOC.results_dir(org) / f"{prefix}_loc_psortdb.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(df)} proteins)", flush=True)
    print(f"[{org}] classes: {df['psortb_localization'].value_counts().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
