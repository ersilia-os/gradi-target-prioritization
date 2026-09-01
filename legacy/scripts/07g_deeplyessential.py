"""DeeplyEssential essentiality prediction (docs §4.3c) — SCAFFOLD / DEFERRED.

DeeplyEssential (Hasan & Lonardi 2020, github.com/ucrbioinfo/DeeplyEssential) is a DNA+protein deep
neural network trained on DEG. It is DEFERRED this round because: (1) no trained weights are released
(only code + training data), (2) the stack is Python-2.7 + TensorFlow-1.6 + Keras-2.1.5, effectively
unbuildable in the py3.11 `gradi` env, (3) it needs per-gene DNA which we do not currently stage, and
(4) the repository ships no license. The essentiality signal is already covered by three independent
predictors (ProteomeLM 4.3a, Geptop 4.3b, FBA 4.3d), so this is a low-value addition.

This script only emits a placeholder column so the 07h merge keeps a uniform schema (mirrors the
scripts/06f_af2bind.py deferred pattern). The node is marked `:::planned` in docs/04_essentiality.md.

Output (placeholder): output/results/<org>/<prefix>_ess_deeplyessential.csv
  uniprot_accession, deeplyessential_score (NaN), deeplyessential_essential (NA), deeplyessential_status
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--run", action="store_true",
                    help="attempt a real run (NOT IMPLEMENTED — no released weights / py2+TF1.6)")
    args = ap.parse_args()

    if args.run:
        raise SystemExit(
            "DeeplyEssential execution is deferred (track 4.3c): no released weights, license-less, "
            "and a Python-2.7 + TensorFlow-1.6 stack. See docs/04_essentiality.md §4.3c. "
            "Essentiality is covered by 07d (ProteomeLM), 07e (Geptop) and 07f (FBA)."
        )

    org = args.organism
    _, prefix = E.ORGANISMS[org]
    df = pd.DataFrame({
        "uniprot_accession": E.load_accessions(org),
        "deeplyessential_score": pd.NA,
        "deeplyessential_essential": pd.NA,
        "deeplyessential_status": "deferred",
    })
    out = E.results_dir(org) / f"{prefix}_ess_deeplyessential.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote placeholder {out.relative_to(E.REPO_ROOT)} "
          f"({len(df)} proteins; DeeplyEssential deferred)", flush=True)


if __name__ == "__main__":
    main()
