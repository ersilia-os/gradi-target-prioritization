"""Paralog clusters for a proteome, so cross-validation folds cannot leak.

Train and test for a per-organism endpoint are the SAME organism, so the twins to separate are
paralogs, not orthologs. A random split puts a protein's near-duplicate in the training fold and
reports a score the model did not earn -- stage 04's `_splitters` docstring records what that cost
v1's essentiality head, which used an ungrouped `cv=5`.

The anchor species get their groups from stage 05's OrthoFinder run. A tier-D SCREEN STRAIN
(ECL8, KPNIH1) is in no OrthoFinder run, so its groups are built here: MMseqs2 `easy-cluster` at
30% identity / 80% coverage -- the same threshold stage 04 inherited for its S. aureus clusters, so
the two axes group on a comparable definition of "too similar to split".

MMseqs2 comes from `gradi-ortho` (`GRADI_MMSEQS_BIN` overrides); it has no osx-arm64 build and must
not be installed into `gradi`.

Run with the `gradi` env:
    python scripts/essentiality/paralog_clusters.py --strain kpneumoniae__ecl8__GCA_000315385.1
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

NCBI_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "ncbi"
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality" / "scratch" / "paralog_clusters"
DEFAULT_MMSEQS = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin" / "mmseqs"

MIN_SEQ_ID = 0.30        # stage 04's MMseqs2 clustering threshold, reused deliberately
MIN_COV = 0.80

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def mmseqs_bin() -> str:
    p = Path(os.environ.get("GRADI_MMSEQS_BIN", DEFAULT_MMSEQS))
    if not p.exists():
        sys.exit(f"FATAL no mmseqs at {p}. It lives in `gradi-ortho` (osx-64, Rosetta); do NOT "
                 "install it into `gradi`.")
    return str(p)


def cluster(label: str, min_id: float, min_cov: float, refresh: bool) -> pd.DataFrame:
    out = OUT_DIR / f"{label}_id{int(min_id * 100)}_cov{int(min_cov * 100)}.tsv"
    if out.exists() and not refresh:
        say(f"  cached {out.relative_to(REPO_ROOT)}")
        return pd.read_csv(out, sep="\t")

    faa = NCBI_DIR / f"{label}.faa"
    if not faa.exists():
        sys.exit(f"FATAL {faa} missing -- run scripts/proteomes/download.py --tier D --only {label}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        pre = Path(td) / "clu"
        cmd = [mmseqs_bin(), "easy-cluster", str(faa), str(pre), str(Path(td) / "tmp"),
               "--min-seq-id", str(min_id), "-c", str(min_cov), "--cov-mode", "0",
               "--threads", "4", "-v", "1"]
        say(f"  mmseqs easy-cluster --min-seq-id {min_id} -c {min_cov}")
        r = subprocess.run(cmd, capture_output=True, text=True)
        tsv = Path(f"{pre}_cluster.tsv")
        if r.returncode != 0 or not tsv.exists():
            sys.exit(f"FATAL mmseqs failed (exit {r.returncode}):\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        d = pd.read_csv(tsv, sep="\t", names=["cluster", "member"])
    d.to_csv(out, sep="\t", index=False)
    say(f"  wrote {out.relative_to(REPO_ROOT)}")
    return d


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--strain", nargs="+", required=True)
    ap.add_argument("--min-seq-id", type=float, default=MIN_SEQ_ID)
    ap.add_argument("--min-cov", type=float, default=MIN_COV)
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    say("=" * 92)
    say("paralog clusters -- the CV grouping for screen-strain endpoints")
    say("=" * 92)
    for label in a.strain:
        say(f"\n{label}")
        d = cluster(label, a.min_seq_id, a.min_cov, a.refresh)
        n_prot, n_clu = d["member"].nunique(), d["cluster"].nunique()
        sizes = d.groupby("cluster").size()
        multi = int((sizes > 1).sum())
        say(f"  {n_prot} proteins -> {n_clu} clusters "
            f"({n_prot - n_clu} collapsed into a paralog family)")
        say(f"  families with >1 member: {multi}   largest: {int(sizes.max())}")
        say(f"  singletons: {int((sizes == 1).sum())} ({100 * (sizes == 1).sum() / n_clu:.1f}%)")
    say("\n" + "=" * 92)


if __name__ == "__main__":
    main()
