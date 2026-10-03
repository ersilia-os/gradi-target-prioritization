"""ProteomeLM-Ess worker. Runs in the `gradi-plm-ess` env, NEVER in `gradi`.

`scripts/essentiality/proteomelm_ess.py` shells out to this file with the `gradi-plm-ess`
interpreter; nothing imports it. It reads a protein FASTA and writes one TSV:

    protein_id <TAB> p_essential <TAB> rank

plus a small JSON manifest recording the head revision, the backbone, the layer and the timings --
the things that determine the numbers and would otherwise be unrecoverable from the TSV alone.

Why a separate process at all
-----------------------------
The essentiality head ships in the authors' package from GitHub, which is a DIFFERENT build from
the `proteomelm` 1.0.0 installed in `gradi`. Measured before the split: the two backbones' forward
passes are numerically identical for our use (the only difference is upstream's added support for
pre-expanded 3D/4D attention masks, a branch a standard 2D mask never takes, plus a `polarize()`
helper that is never called). So the split is not about the model code -- it is about pip: letting
the git package resolve its own `torch`/`esm`/`transformers` inside `gradi` is the failure this
project has hit repeatedly, and stage 01's ESM-C is what breaks.

**This file must not be fed our own `proteomelm_<species>.npz`.** Those matrices are z-scored
genome-wide and window sequences above 4,096 aa; the head was fitted on raw `hidden_states[8]` with
sequences TRUNCATED at 4,096. Same model, different numbers -- the *External models* rule in
CLAUDE.md, and the reason this worker recomputes ESM-C and the backbone pass itself rather than
reusing anything.

CLASS 0 IS ESSENTIAL in the head's output (`ESSENTIAL_LABEL = 0`), so `p_essential` is the softmax
of class 0, not class 1. The authors' `essentiality_table` already applies that via
`config.essential_index`; this worker does not re-derive it, precisely so the convention has one
owner. The `gradi` side spot-checks the polarity against the ribosome anyway, because an inverted
essentiality column is well-formed, confident and exactly wrong.

Stdlib + the authors' package only. No pandas beyond what their `predict` returns, no `src/`:
importing this from `gradi` would drag the git build's dependency tree into that process, which is
the whole thing the env split prevents. Communication is by file, one way.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser(description="Score one proteome with ProteomeLM-Ess.")
    ap.add_argument("--fasta", required=True, type=Path, help="protein FASTA, one whole proteome")
    ap.add_argument("--out", required=True, type=Path, help="TSV to write")
    ap.add_argument("--manifest", type=Path, help="JSON to write (config + timings)")
    ap.add_argument("--head", default="Bitbol-Lab/ProteomeLM-ess")
    ap.add_argument("--device", default="cpu", help="cpu | mps | cuda")
    ap.add_argument("--esm-device", default=None, help="defaults to --device")
    args = ap.parse_args()

    # Imported here so `--help` works without the heavy stack.
    from proteomelm.essentiality import EssentialityPredictor, read_fasta

    ids, sequences = read_fasta(args.fasta)
    print(f"[worker] {len(ids):,} proteins from {args.fasta.name}", flush=True)

    started = time.time()
    predictor = EssentialityPredictor.from_pretrained(
        args.head, device=args.device, esm_device=args.esm_device or args.device)
    table = predictor.predict(ids, sequences)
    elapsed = time.time() - started

    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)
    print(f"[worker] wrote {len(table):,} rows to {args.out} in {elapsed:.1f}s", flush=True)

    if args.manifest:
        cfg = predictor.config
        args.manifest.write_text(json.dumps({
            "head": args.head,
            "backbone": cfg.backbone,
            "layer": cfg.layer,
            "esm_model": cfg.esm_model,
            "backbone_dtype": cfg.backbone_dtype,
            "essential_index": cfg.essential_index,
            "n_proteins": len(ids),
            "device": args.device,
            "seconds": round(elapsed, 1),
            "timings": {k: round(v, 1) for k, v in predictor.timings.items()},
        }, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
