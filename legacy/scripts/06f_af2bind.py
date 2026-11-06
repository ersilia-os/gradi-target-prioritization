"""AF2Bind binding-site prediction (docs §2.3a) — SCAFFOLD / DEFERRED.

AF2Bind (https://github.com/sokrypton/af2bind) predicts small-molecule binding residues from
AlphaFold2's pair representation. A full-proteome run needs AF2 MSAs + a model forward pass per
protein (GPU-bound, effectively days on CPU), so this round it is DEFERRED: the node is marked
`:::planned` in docs/02_ligandability.md and this script only emits a placeholder column so the
06g merge has a uniform schema. The real implementation outline is kept below for when a GPU is
available; it can also be run on a small high-priority shortlist via --shortlist.

Output (placeholder): output/results/<org>/<prefix>_af2bind.csv
  uniprot_accession, af2bind_max_score (NaN), af2bind_n_binding_residues (NaN), af2bind_status

Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    ap.add_argument(
        "--run",
        action="store_true",
        help="actually run AF2Bind (NOT IMPLEMENTED here — needs the af2bind/ColabFold env + GPU)",
    )
    args = ap.parse_args()

    if args.run:
        raise SystemExit(
            "AF2Bind execution is deferred (track 2.3a). Wire up af2bind + AF2 params on a GPU "
            "host, predict per-residue binding scores on the AlphaFold models in "
            f"data/processed/{args.organism}/alphafold/cif/, then fill af2bind_max_score / "
            "af2bind_n_binding_residues. See docs/02_ligandability.md §2.3a."
        )

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    accs = L.load_accessions(org)
    df = pd.DataFrame(
        {
            "uniprot_accession": accs,
            "af2bind_max_score": pd.NA,
            "af2bind_n_binding_residues": pd.NA,
            "af2bind_status": "deferred",
        }
    )
    out = L.results_dir(org) / f"{prefix}_af2bind.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote placeholder {out} ({len(df)} proteins; AF2Bind deferred)", flush=True)


if __name__ == "__main__":
    main()
