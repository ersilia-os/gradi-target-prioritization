"""Stage 09c — DeepLocPro 1.0 subcellular-localization prediction (docs §5.1c).

The primary predictor for the localization axis. DeepLocPro is an ESM-2 650M backbone with an
attention-pooled classifier head, ensembled over 20 nested-CV checkpoints, emitting a probability
over six prokaryotic compartments. On the post-2010 Gram-negative benchmark it beats PSORTb 3.0
across the board (accuracy 0.74 vs 0.34, macro-F1 0.75 vs 0.35, multiclass MCC 0.69 vs 0.30) and,
unlike PSORTb, always returns a call — which is what lets us close the ~60% `unknown` tail that
the UniProt-only track (09a) leaves behind.

We drive `DeepLocPro.model.EnsembleModel` directly instead of shelling out to the `deeplocpro`
CLI, for three reasons:
  1. the CLI's `embed_batch()` gates GPU placement on `torch.cuda.is_available()`, so its own
     `-d mps` flag crashes on Apple Silicon ("Placeholder storage has not been allocated on MPS
     device") — driving the model ourselves puts the tokens and the weights on the same device;
  2. it writes a timestamped CSV with no resume, whereas the run is ~10k single-sequence forward
     passes and needs to survive an interruption;
  3. we want the full probability vector and a top1-top2 margin as a confidence, not just argmax.

Caches one small JSON per protein under
data/processed/<organism>/localization/deeplocpro/, the same resumable pattern as 06e/07*.

Output: output/results/<organism>/<prefix>_loc_deeplocpro.csv
Run with the **`gradi-loc`** conda env (DeepLocPro needs `fair-esm`, which collides with the
EvolutionaryScale `esm` package that `gradi` uses for the 01a ESM-C embeddings). See install.sh.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

# DeepLocPro's class order, verbatim from DeepLocPro/deeplocpro.py.
DLP_LABELS = [
    "Cell wall & surface",
    "Extracellular",
    "Cytoplasmic",
    "Cytoplasmic Membrane",
    "Outer Membrane",
    "Periplasmic",
]
CANON = [LOC.DEEPLOCPRO_MAP[lbl.lower()] for lbl in DLP_LABELS]

# ESM-2 cost grows quadratically with length and the longest Kp proteins are ~2-4k residues.
# Localization signal is overwhelmingly N-terminal (signal peptides, TM topology), so truncating
# the C-terminus is the standard compromise; we record which proteins it affected.
MAX_LENGTH = 2000


def load_model(device: torch.device):
    from DeepLocPro.model import EnsembleModel

    model = EnsembleModel()
    model.eval()
    model.to(device)
    return model


@torch.no_grad()
def predict_one(model, sequence: str, device: torch.device) -> list[float]:
    """Reimplements EnsembleModel.embed_batch for a single sequence, device-correct."""
    batch_converter = model.esm_alphabet.get_batch_converter()
    _labels, _strs, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)

    out = model.esm_model(toks, repr_layers=[33], return_contacts=False)["representations"][33]
    out = out.cpu()
    out[out != out] = 0.0  # ESM-2 can emit NaN on odd residues; DeepLocPro zeroes them too

    embedding = out.transpose(0, 1)[1:-1][:, 0]          # strip BOS/EOS -> (L, 1280)
    mask = torch.ones(embedding.shape[0])
    probs, _attn = model(embedding.unsqueeze(0).to(device), mask.unsqueeze(0).to(device))
    return [float(p) for p in probs[0]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--limit", type=int, default=None, help="first N proteins (debug)")
    ap.add_argument("--device", default="mps", choices=["mps", "cuda", "cpu"])
    ap.add_argument("--max-length", type=int, default=MAX_LENGTH)
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]

    device = torch.device(args.device)
    if args.device == "mps" and not torch.backends.mps.is_available():
        print("[warn] MPS unavailable, falling back to CPU", flush=True)
        device = torch.device("cpu")

    accs = LOC.load_accessions(org)
    if args.limit:
        accs = accs[: args.limit]
    seqs = LOC.read_fasta(LOC.proteome_fasta(org))
    cache = LOC.localization_processed_dir(org, "deeplocpro")

    todo = [a for a in accs if not (cache / f"{a}.json").exists()]
    print(f"[{org}] {len(accs)} proteins, {len(accs) - len(todo)} cached, {len(todo)} to run "
          f"on {device.type}", flush=True)

    if todo:
        model = load_model(device)
        t0 = time.time()
        for n, acc in enumerate(todo, 1):
            seq = seqs.get(acc)
            if not seq:
                continue
            truncated = len(seq) > args.max_length
            probs = predict_one(model, seq[: args.max_length], device)
            (cache / f"{acc}.json").write_text(json.dumps({"probs": probs, "truncated": truncated}))
            if n % 200 == 0 or n == len(todo):
                rate = n / (time.time() - t0)
                eta = (len(todo) - n) / rate / 60
                print(f"  {n}/{len(todo)}  {rate:.1f} prot/s  eta {eta:.0f} min", flush=True)

    rows = []
    for acc in accs:
        f = cache / f"{acc}.json"
        if not f.exists():
            continue
        rec = json.loads(f.read_text())
        probs = rec["probs"]
        ranked = sorted(range(len(probs)), key=lambda i: -probs[i])
        row = {
            "uniprot_accession": acc,
            "dlp_localization": CANON[ranked[0]],
            "dlp_confidence": round(probs[ranked[0]], 4),
            "dlp_margin": round(probs[ranked[0]] - probs[ranked[1]], 4),
            "dlp_truncated": rec["truncated"],
        }
        for cls, p in zip(CANON, probs):
            row[f"dlp_p_{cls}"] = round(p, 4)
        rows.append(row)

    df = pd.DataFrame(rows)
    out = LOC.results_dir(org) / f"{prefix}_loc_deeplocpro.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(df)} proteins)", flush=True)
    if len(df):
        print(f"[{org}] classes: {df['dlp_localization'].value_counts().to_dict()}", flush=True)


if __name__ == "__main__":
    main()
