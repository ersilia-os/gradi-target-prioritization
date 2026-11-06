"""Geptop 2.0 essentiality prediction (docs §4.3b, orthology+phylogeny BASELINE).

Geptop 2.0 (Wen et al. 2019) scores gene essentiality by transferring DEG essential-gene calls from
37 reference bacterial genomes via best-hit orthology, weighting each reference by its evolutionary
closeness to the query. It is the classical, interpretable baseline that cross-checks the ProteomeLM
predictor: a gene whose ortholog is essential across many diverse bacteria is a confident call.

The upstream standalone is a Python-2 + NCBI-BLAST `.rar`; rather than fight that headless, we
reimplement the algorithm with DIAMOND (already installed) over the reference data bundled in the
`.rar` (staged by hand under data/raw/other/essentiality/geptop/Geptop2/: `datasets2/*.faa` = the 37
reference proteomes; `DEG2` = the pooled essential-gene ID list). Per Geptop:

  score_raw[gene] = Σ_ref  w_ref · [ best DIAMOND hit of `gene` in reference `ref` is DEG-essential ]

with the reference weight `w_ref` set to the median %identity of that reference's best hits to the
query proteome (a data-derived evolutionary-closeness proxy for Geptop's composition-vector distance;
documented simplification — see docs/essentiality_log.md). The raw score is min-max normalised to
[0,1] and thresholded at Geptop's default cutoff 0.24.

Output: output/results/<org>/<prefix>_ess_geptop.csv
  uniprot_accession, geptop_score, geptop_essential, geptop_status
Run with the `gradi` env (shells out to the gradi-ortho DIAMOND). Per-reference hits are cached.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

GEPTOP_DIR = E.REPO_ROOT / "data" / "raw" / "other" / "essentiality" / "geptop" / "Geptop2"
DATASETS = GEPTOP_DIR / "datasets2"
DEG2 = GEPTOP_DIR / "DEG2"
CUTOFF = E.GEPTOP_CUTOFF  # 0.24


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]

    if not DEG2.exists() or not DATASETS.exists():
        # graceful skip: emit a deferred placeholder so the merge schema stays uniform
        print(f"[{org}] Geptop reference data missing under {GEPTOP_DIR}; emitting placeholder. "
              "Extract Geptop_v2.0.rar with `unar` (see install.sh).", flush=True)
        accs = E.load_accessions(org)
        pd.DataFrame({"uniprot_accession": accs, "geptop_score": pd.NA,
                      "geptop_essential": pd.NA, "geptop_status": "deferred: no reference data"}
                     ).to_csv(E.results_dir(org) / f"{prefix}_ess_geptop.csv", index=False)
        return

    essential_ids = {ln.strip() for ln in open(DEG2) if ln.strip()}
    print(f"[{org}] DEG2 essential ids: {len(essential_ids)}", flush=True)

    accs = E.load_accessions(org)
    query_faa = E.proteome_fasta(org)
    refs = sorted(DATASETS.glob("*.faa"))
    hit_dir = E.essentiality_processed_dir(org, "geptop")

    # weighted vote accumulators
    raw = pd.Series(0.0, index=accs)
    wsum = 0.0
    for i, ref in enumerate(refs, 1):
        out_tsv = hit_dir / f"{prefix}_vs_{ref.stem}.tsv"
        E.run_diamond_blastp(query_faa, ref, out_tsv)  # cached
        hits = E.load_diamond_hits(out_tsv)
        if hits.empty:
            continue
        hits["acc"] = hits["qseqid"].map(E.acc_from_header)
        best = hits.sort_values("bitscore", ascending=False).drop_duplicates("acc")
        w_ref = float(np.median(best["pident"])) / 100.0  # evolutionary-closeness proxy
        best["ess"] = best["sseqid"].isin(essential_ids)
        contrib = best.set_index("acc")["ess"].astype(float) * w_ref
        raw = raw.add(contrib.reindex(accs).fillna(0.0), fill_value=0.0)
        wsum += w_ref
        if i % 10 == 0 or i == len(refs):
            print(f"  {i}/{len(refs)} references (w_ref[{ref.stem}]={w_ref:.2f})", flush=True)

    # normalise to [0,1] (Geptop min-max); guard degenerate range
    rng = raw.max() - raw.min()
    score = (raw - raw.min()) / rng if rng > 0 else raw * 0.0
    df = pd.DataFrame({
        "uniprot_accession": accs,
        "geptop_score": score.round(4).to_numpy(),
        "geptop_essential": (score >= CUTOFF).to_numpy(),
        "geptop_status": f"ok ({len(refs)} refs)",
    })
    out = E.results_dir(org) / f"{prefix}_ess_geptop.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)} ({len(df)} proteins; "
          f"{int(df.geptop_essential.sum())} essential at cutoff {CUTOFF})", flush=True)


if __name__ == "__main__":
    main()
