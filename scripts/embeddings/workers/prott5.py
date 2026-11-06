"""ProtT5-XL-U50 per-protein embeddings. A worker, not a stage.

Runs under the `gradi-loc` interpreter and is never imported. Two reasons, both load-bearing:

  * ProtT5 is loaded through `T5Tokenizer`, and transformers 5.x routes that through the tiktoken
    converter and dies with a spurious "`tiktoken` is required to read a `tiktoken` file". `gradi-loc`
    pins `transformers==4.44.2` for exactly this; `gradi` runs 4.48.1.
  * The weights are already on disk there -- TMbed ships `Rostlab/prot_t5_xl_half_uniref50-enc`
    (1024-d, fp16, encoder-only, 2.4 GB) inside its package. Nothing is downloaded.

POOLING
-------
Mean over residue positions, **excluding padding and excluding the terminal `</s>`**. ProtT5's
tokenizer appends EOS and no BOS, and the vector UniProt publishes is the mean over the sequence's
own residues. Getting this wrong yields an embedding that looks perfectly reasonable and is subtly
not the published one, which is why the caller checks a proteome UniProt publishes before trusting
one it does not. `--pooling` exposes the alternatives so the control can discriminate between them
rather than assume.

RESIDUE ENCODING
----------------
ProtT5 wants residues space-separated, and U/Z/O/B mapped to `X` -- its vocabulary has no token for
them, so skipping this silently maps rare residues to `<unk>`. Selenocysteine is real in these
proteomes (E. coli `fdnG`, `fdhF`), so it is not hypothetical.

LONG SEQUENCES
--------------
T5 self-attention is O(L^2), and *S. aureus* carries a 9,535-residue protein (36x the median) that
would ask for ~6 GB of attention in one pass. Sequences above `--max-len` are therefore embedded in
non-overlapping windows and pooled as the **length-weighted mean of window means**, which is exactly
the mean over all residues -- the windows differ from a single pass only in that each sees its own
context, not the whole chain. Every such protein is FLAGGED in the output (`chunked`), never
silently approximated. The default 4096 leaves every *labeled* S. aureus protein (max 2,478 aa) in a
single pass; ProtT5 uses relative position embeddings, so length carries no positional cliff, and
the caller's control validates single-pass agreement out to E. coli's longest (2,339 aa).

Sharded and resumable: the encoder on MPS has stalled before in this repo (an uninterruptible kernel
wait), so a killed run must not lose the batch. Shard size is in the filename, so changing it cannot
silently reuse mismatched shards.

IN   a .tsv with columns `uniprot_ac` and `sequence`
OUT  a .npz: accessions, embeddings (n, 1024) float32, model, pooling, dim

  ~/miniconda3/envs/gradi-loc/bin/python scripts/embeddings/workers/prott5.py --in p.tsv --out e.npz
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_MODEL = (Path.home() / "miniconda3" / "envs" / "gradi-loc" / "lib" / "python3.11"
                 / "site-packages" / "tmbed" / "models" / "t5")
EMBED_DIM = 1024
NONSTANDARD = re.compile(r"[UZOB]")


def encode(seq: str) -> str:
    """ProtT5's input convention: space-separated residues, U/Z/O/B -> X."""
    return " ".join(NONSTANDARD.sub("X", seq.strip().upper()))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", dest="out", required=True)
    ap.add_argument("--model", default=str(DEFAULT_MODEL))
    ap.add_argument("--shard-dir", default=None, help="resume cache; default alongside --out")
    ap.add_argument("--shard-size", type=int, default=250)
    ap.add_argument("--batch-tokens", type=int, default=6000,
                    help="max residues per batch; batches are built by length to limit padding")
    ap.add_argument("--device", default="mps", choices=["mps", "cpu"])
    ap.add_argument("--pooling", default="mean_no_eos",
                    choices=["mean_no_eos", "mean_with_eos"])
    ap.add_argument("--max-len", type=int, default=4096,
                    help="sequences longer than this are embedded in windows and pooled by length")
    args = ap.parse_args()

    import torch
    from transformers import T5EncoderModel, T5Tokenizer

    df = pd.read_csv(args.inp, sep="\t")
    for col in ("uniprot_ac", "sequence"):
        if col not in df.columns:
            sys.exit(f"FATAL {args.inp} has no `{col}` column (has {list(df.columns)})")
    # Sort by length so each batch is near-uniform and padding is small. The output is re-ordered
    # back to the input order at the end -- row alignment is a promise the caller relies on.
    df = df.reset_index(drop=True)
    order = df["sequence"].str.len().sort_values(ascending=False).index.to_numpy()

    shard_dir = Path(args.shard_dir or (Path(args.out).parent / "shards_prott5"))
    shard_dir.mkdir(parents=True, exist_ok=True)

    device = torch.device(args.device)
    print(f"loading ProtT5 from {args.model} on {device} ...", flush=True)
    tok = T5Tokenizer.from_pretrained(args.model, do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained(args.model).to(device).eval()
    # fp16 weights on CPU are slow and on MPS are fine; float() on CPU avoids a silent crawl.
    model = model.float() if args.device == "cpu" else model.half()
    print(f"loaded: {sum(p.numel() for p in model.parameters()) / 1e6:.0f}M params, "
          f"pooling={args.pooling}", flush=True)

    out = np.zeros((len(df), EMBED_DIM), dtype=np.float32)
    done = np.zeros(len(df), dtype=bool)
    chunked = np.zeros(len(df), dtype=bool)

    def embed_one_long(seq: str) -> np.ndarray:
        """Length-weighted mean of window means == the mean over every residue."""
        acc, total = np.zeros(EMBED_DIM, dtype=np.float64), 0
        for w in range(0, len(seq), args.max_len):
            piece = seq[w:w + args.max_len]
            e = tok.batch_encode_plus([encode(piece)], add_special_tokens=True,
                                      padding="longest", return_tensors="pt")
            i_, m_ = e["input_ids"].to(device), e["attention_mask"].to(device)
            with torch.no_grad():
                h = model(input_ids=i_, attention_mask=m_).last_hidden_state
            mm = m_.clone()
            if args.pooling == "mean_no_eos":
                mm[0, m_.sum(1) - 1] = 0
            mm = mm.unsqueeze(-1).to(h.dtype)
            v = ((h * mm).sum(1) / mm.sum(1).clamp(min=1)).float().cpu().numpy()[0]
            acc += v.astype(np.float64) * len(piece)
            total += len(piece)
        return (acc / max(total, 1)).astype(np.float32)

    n_shards = (len(df) + args.shard_size - 1) // args.shard_size
    t_start = time.time()
    for s in range(n_shards):
        rows = order[s * args.shard_size:(s + 1) * args.shard_size]
        shard = shard_dir / f"shard_{args.shard_size}_{s:05d}.npz"
        if shard.exists():
            z = np.load(shard)
            out[rows] = z["emb"]
            done[rows] = True
            continue

        emb = np.zeros((len(rows), EMBED_DIM), dtype=np.float32)
        pos = {r: idx for idx, r in enumerate(rows)}     # row id -> slot in this shard
        i = 0
        while i < len(rows):
            # Pack a batch under the residue budget; always take at least one sequence, so a single
            # protein longer than the budget still runs rather than looping forever.
            j, budget = i, 0
            while j < len(rows):
                L = len(df.at[rows[j], "sequence"])
                if j > i and budget + L > args.batch_tokens:
                    break
                budget += L
                j += 1
            chunk = rows[i:j]
            # Pull anything over the window budget out of the batch and window it on its own.
            long_rows = [r for r in chunk if len(df.at[r, "sequence"]) > args.max_len]
            if long_rows:
                for r in long_rows:
                    emb[pos[r]] = embed_one_long(df.at[r, "sequence"])
                    chunked[r] = True
                    print(f"  windowed {df.at[r, 'uniprot_ac']} "
                          f"({len(df.at[r, 'sequence'])} aa, {args.max_len} per window)", flush=True)
                chunk = np.array([r for r in chunk if r not in set(long_rows)])
                if len(chunk) == 0:
                    i = j
                    continue
            seqs = [encode(df.at[r, "sequence"]) for r in chunk]
            enc = tok.batch_encode_plus(seqs, add_special_tokens=True, padding="longest",
                                        return_tensors="pt")
            ids = enc["input_ids"].to(device)
            mask = enc["attention_mask"].to(device)
            with torch.no_grad():
                hs = model(input_ids=ids, attention_mask=mask).last_hidden_state
            m = mask.clone()
            if args.pooling == "mean_no_eos":
                # Drop the terminal </s>: it is attended (mask==1) but is not a residue.
                lens = mask.sum(1)
                m[torch.arange(len(chunk), device=device), lens - 1] = 0
            m = m.unsqueeze(-1).to(hs.dtype)
            pooled = (hs * m).sum(1) / m.sum(1).clamp(min=1)
            pooled = pooled.float().cpu().numpy()
            for n_, r in enumerate(chunk):
                emb[pos[r]] = pooled[n_]
            i = j
        tmp = shard.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, emb=emb)
        tmp.rename(tmp.with_name(shard.name))          # atomic: a half-written shard is never read
        out[rows] = emb
        done[rows] = True
        el = time.time() - t_start
        print(f"  shard {s + 1}/{n_shards}: {int(done.sum())}/{len(df)} proteins, "
              f"{el / 60:.1f} min elapsed, ~{el / (s + 1) * (n_shards - s - 1) / 60:.1f} min left",
              flush=True)

    if not done.all():
        sys.exit(f"FATAL {int((~done).sum())} proteins were never embedded")
    norms = np.linalg.norm(out, axis=1)
    if not np.isfinite(out).all() or (norms == 0).any():
        sys.exit(f"FATAL {int((~np.isfinite(out)).all(1).sum())} non-finite / "
                 f"{int((norms == 0).sum())} all-zero embeddings -- a zero row is a real position "
                 "in embedding space and must never be shipped as one")

    np.savez_compressed(args.out, accessions=df["uniprot_ac"].to_numpy().astype(object),
                        embeddings=out, model="Rostlab/prot_t5_xl_half_uniref50-enc",
                        pooling=args.pooling, dim=EMBED_DIM, chunked=chunked,
                        max_len=args.max_len)
    print(f"wrote {args.out}: {out.shape}, {int(chunked.sum())} windowed, "
          f"{(time.time() - t_start) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
