"""DeepLocPro worker for stage 03. Runs in the `gradi-loc` env, NEVER in `gradi`.

Not a stage, which is why it lives in `scripts/workers/` rather than beside the numbered stages.
`scripts/localization/predict.py` shells out to this file with the `gradi-loc` interpreter; nothing
imports it. It reads a FASTA and writes one small JSON per protein into a cache directory:

    {"probs": [6 floats, DeepLocPro's own class order], "truncated": bool}

and nothing else. The label vocabulary, the Gram remap and every table live on the `gradi` side in
`src/localization.py` -- this file deliberately knows no project conventions.

Why a separate process at all
-----------------------------
DeepLocPro depends on `fair-esm`, which installs the same top-level `esm` package as the
EvolutionaryScale `esm` that stage 01 uses for ESM-C. They cannot coexist, so the predictor lives in
`gradi-loc` (see install.sh) and is reached across a process boundary.

`gradi-loc` also carries much newer pandas and numpy than `gradi`, so this file imports **stdlib and
torch only**. No pandas, no `src/`. Communication is by file, one way.

**Why this is not in `src/`.** In this repo `src/` means "imported in-process by `gradi`" -- it is a
package, and `proteomes`, `embeddings`, `function` and `localization` are all reached with
`from src import X`. Importing *this* file from `gradi` would pull `torch` + `DeepLocPro` +
`fair-esm` into that process, which is exactly the `esm` collision the env split exists to prevent.
Filing it under `src/` would put it next to modules that are safe to import, inside the package
`gradi` imports from. The reusable half of the predictor -- the vocabulary, the Gram remap, the
TMbed parsing, the loaders -- IS in `src/localization.py`; what is left here is only what cannot run
under the `gradi` interpreter at all.

TMbed needs no equivalent worker: its `tmbed` console script already has `gradi-loc`'s interpreter in
its shebang, so stage 03 can shell the binary directly.

Why we drive the model instead of the `deeplocpro` CLI
-----------------------------------------------------
1. **The CLI's `-d mps` flag is broken.** `EnsembleModel.embed_batch()` gates device placement on
   `torch.cuda.is_available()` (`DeepLocPro/model.py:45`), so with `-d mps` the weights move to MPS
   while the tokens stay on CPU: *"Placeholder storage has not been allocated on MPS device"*.
   `predict_one` below is that method rewritten to put both on the same device.
2. **No resume.** The CLI writes one timestamped `results_<YYYYmmdd-HHMMSS>.csv` at the end; this is
   ~13k single-sequence forward passes and must survive an interruption.
3. **It discards the probability vector**, and also writes a per-protein attention plot we do not
   want.

Invoked as:
    <gradi-loc python> scripts/localization/workers/deeplocpro.py --fasta Q.faa --cache DIR [--device auto]
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import torch

# ESM-2 attention is quadratic in length; the C-terminus is what gets cut because localization
# signal is overwhelmingly N-terminal. Kept in step with DLP_MAX_LENGTH in src/localization.py.
MAX_LENGTH = 2000


def read_fasta(path: Path) -> dict[str, str]:
    """Bare-accession FASTA -> {accession: sequence}. Headers are written bare by stage 03."""
    seqs: dict[str, str] = {}
    acc, chunks = None, []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if acc is not None:
                seqs[acc] = "".join(chunks)
            acc, chunks = line[1:].strip(), []
        elif acc is not None:
            chunks.append(line.strip())
    if acc is not None:
        seqs[acc] = "".join(chunks)
    return seqs


def pick_device(name: str) -> torch.device:
    if name == "auto":
        if torch.cuda.is_available():
            return torch.device("cuda")
        if torch.backends.mps.is_available():
            return torch.device("mps")
        return torch.device("cpu")
    if name == "mps" and not torch.backends.mps.is_available():
        print("[warn] MPS unavailable, falling back to CPU", flush=True)
        return torch.device("cpu")
    return torch.device(name)


@torch.no_grad()
def predict_one(model, sequence: str, device: torch.device) -> list[float]:
    """`EnsembleModel.embed_batch` for a single sequence, device-correct.

    The only change from DeepLocPro's own code is that `toks` follows the weights onto `device`
    instead of being moved only when CUDA is present.
    """
    batch_converter = model.esm_alphabet.get_batch_converter()
    _labels, _strs, toks = batch_converter([("seq", sequence)])
    toks = toks.to(device)

    out = model.esm_model(toks, repr_layers=[33], return_contacts=False)["representations"][33]
    out = out.cpu()
    out[out != out] = 0.0                      # ESM-2 can emit NaN; DeepLocPro zeroes them too

    embedding = out.transpose(0, 1)[1:-1][:, 0]        # strip BOS/EOS -> (L, 1280)
    mask = torch.ones(embedding.shape[0])
    probs, _attn = model(embedding.unsqueeze(0).to(device), mask.unsqueeze(0).to(device))
    return [float(p) for p in probs[0]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--cache", type=Path, required=True, help="one <accession>.json per protein")
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    ap.add_argument("--max-length", type=int, default=MAX_LENGTH)
    ap.add_argument("--report-every", type=int, default=200)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    args.cache.mkdir(parents=True, exist_ok=True)
    todo = [a for a, s in seqs.items() if s and not (args.cache / f"{a}.json").exists()]

    device = pick_device(args.device)
    print(f"    {len(seqs)} proteins, {len(seqs) - len(todo)} cached, {len(todo)} to run "
          f"on {device.type}", flush=True)
    if not todo:
        return

    from DeepLocPro.model import EnsembleModel
    model = EnsembleModel()
    model.eval()
    model.to(device)

    failures: list[str] = []
    t0 = time.time()
    for n, acc in enumerate(todo, 1):
        seq = seqs[acc]
        try:
            probs = predict_one(model, seq[: args.max_length], device)
        except Exception as exc:                          # noqa: BLE001 - name it, never swallow it
            failures.append(f"{acc} ({len(seq)} aa): {type(exc).__name__}: {exc}")
            continue
        # Atomic: a partial JSON must never be cached as a completed protein.
        tmp = args.cache / f"{acc}.json.tmp"
        tmp.write_text(json.dumps({"probs": probs, "truncated": len(seq) > args.max_length}))
        tmp.rename(args.cache / f"{acc}.json")
        if n % args.report_every == 0 or n == len(todo):
            rate = n / max(time.time() - t0, 1e-9)
            print(f"    {n}/{len(todo)}  {rate:.1f} prot/s  "
                  f"~{(len(todo) - n) / max(rate, 1e-9) / 60:.0f} min left", flush=True)

    if failures:
        sys.exit("FAILED on " + str(len(failures)) + " proteins:\n  " + "\n  ".join(failures[:20]))


if __name__ == "__main__":
    main()
