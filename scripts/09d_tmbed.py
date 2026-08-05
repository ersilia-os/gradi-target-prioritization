"""Stage 09d — TMbed transmembrane topology (docs §5.1d).

TMbed (Rostlab, Apache-2.0) predicts, per residue, one of `S` signal peptide, `H` transmembrane
α-helix, `B` transmembrane β-strand, `i` non-transmembrane inside, `o` non-transmembrane outside,
using ProtT5-XL-U50 embeddings. That single pass gives us the two things the localization axis
cannot get from a compartment classifier alone:

  * **the β-barrel call** — a β-barrel mis-assigned as cytoplasmic is the worst-case error for this
    axis (docs §5), so TMbed's `B` segments act as a hard override onto `outer_membrane` in 09g;
  * **the cytoplasmic-domain fraction** — the `i` residue share, which is what turns a flat
    "inner membrane" label into the graded Clp-accessibility score (an IM protein with a large
    cytoplasm-facing domain is reachable by ClpXP; a mostly-buried one is not).

Chosen over DeepTMHMM, which is distributed only through BioLib (Docker or cloud round-trip);
TMbed is a plain `pip install git+…` and needs no daemon.

**Runtime note:** TMbed gates GPU use on `torch.cuda.is_available()` (`tmbed/embed.py`), so on
Apple Silicon it runs ProtT5 on CPU. That is slow, hence the sharding below: each shard's
prediction file is cached under data/processed/<organism>/localization/tmbed/, so an interrupted
run resumes where it stopped rather than starting over.

Output: output/results/<organism>/<prefix>_loc_topology.csv
Run with the **`gradi-loc`** conda env. See install.sh.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

TMBED_BIN = "tmbed"
SHARD_SIZE = 250


def run_shard(fasta: Path, pred: Path, threads: int) -> None:
    cmd = [
        TMBED_BIN, "predict",
        "-f", str(fasta),
        "-p", str(pred),
        "--out-format", "1",          # 3-line, explicit inside/outside
        "--threads", str(threads),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--limit", type=int, default=None, help="first N proteins (debug)")
    ap.add_argument("--shard-size", type=int, default=SHARD_SIZE)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]

    accs = LOC.load_accessions(org)
    if args.limit:
        accs = accs[: args.limit]
    seqs = LOC.read_fasta(LOC.proteome_fasta(org))
    accs = [a for a in accs if seqs.get(a)]

    cache = LOC.localization_processed_dir(org, "tmbed")
    shards = [accs[i:i + args.shard_size] for i in range(0, len(accs), args.shard_size)]
    pending = [(n, s) for n, s in enumerate(shards) if not (cache / f"shard_{args.shard_size}_{n:04d}.pred").exists()]
    print(f"[{org}] {len(accs)} proteins in {len(shards)} shards, {len(pending)} to run", flush=True)

    t0 = time.time()
    for done, (n, shard) in enumerate(pending, 1):
        fa = cache / f"shard_{args.shard_size}_{n:04d}.fasta"
        LOC.write_fasta({a: seqs[a] for a in shard}, fa)
        tmp = cache / f"shard_{args.shard_size}_{n:04d}.pred.tmp"
        run_shard(fa, tmp, args.threads)
        tmp.rename(cache / f"shard_{args.shard_size}_{n:04d}.pred")   # atomic: a partial file is never cached
        fa.unlink(missing_ok=True)
        rate = done / (time.time() - t0)
        print(f"  shard {done}/{len(pending)}  eta {(len(pending) - done) / rate / 60:.0f} min",
              flush=True)

    feats: dict[str, dict] = {}
    for n in range(len(shards)):
        f = cache / f"shard_{args.shard_size}_{n:04d}.pred"
        if f.exists():
            feats.update(LOC.parse_tmbed(f))

    df = pd.DataFrame([{"uniprot_accession": a, **feats[a]} for a in accs if a in feats])
    out = LOC.results_dir(org) / f"{prefix}_loc_topology.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(df)} proteins)", flush=True)
    if len(df):
        print(f"[{org}] beta-barrels: {int(df['is_beta_barrel'].sum())} | "
              f"with TM helix: {int((df['n_tm_helix'] > 0).sum())} | "
              f"signal peptide: {int(df['has_signal_peptide_tmbed'].sum())}", flush=True)


if __name__ == "__main__":
    main()
