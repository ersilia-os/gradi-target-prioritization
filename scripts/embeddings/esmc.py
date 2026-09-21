"""ESM-C 600M mean-pooled embeddings, one vector per protein, one NPZ per species.

Stage 01. Reads the stage-00 tables and writes a 1152-dimensional embedding for every protein.

The model
---------
`esmc_600m` (ESMC_600M_202412), 1152-dim, mean-pooled over residues. This is **the largest ESM-C with
downloadable weights**: the installed `esm` package registers only ESMC_300M and ESMC_600M in
`pretrained.py`, while `esmc-6b` exists solely as a name check for the hosted Forge API, which needs a
token, ships sequences to a third party and is not reproducible offline. Weights come from the local
HuggingFace cache.

Scope
-----
The three bacteria only -- 13,020 proteins, 3,798,100 residues. Human is deliberately out: it held
75% of the residues and every pathological length (titin 34,350 aa, MUC16 14,507, MUC3B 13,477).
Add it back with `--species human` if that ever changes; nothing in the code excludes it.

Measured on this machine (MPS, torch 2.12) before any of this was written:

    100 aa 0.30 s | 257 aa 0.15 s | 500 aa 0.26 s | 994 aa 0.59 s
    2035 aa 1.98 s | 3163 aa 3.24 s | 9535 aa 22.0 s  <- Q2FYJ6, the longest, works fine

~1,700 aa/s in the typical 250-1,000 aa band, falling to ~975 aa/s by 3k as O(L^2) attention bites.
A full run is 42-63 min. No length cap is needed: only 9 proteins exceed 2,000 aa.

ProtT5 is deferred, and the reason is worth knowing before anyone tries again: UniProt publishes
precomputed ProtT5 per-protein embeddings, but for **only 8 proteomes**. E. coli K-12 and human are
among them; **Kp HS11286 and Sa NCTC 8325 both 404**. The all-Swiss-Prot file would cover 7 of 5,728
Kp proteins. Covering the anchor organism means running ProtT5-XL-U50 locally, which is a separate
decision.

Sharding
--------
v1's `01a` accumulated every vector in memory and wrote once at the end, so a killed run lost
everything. Here the work is cut into shards of 250 written under `scratch/shards_esmc/` via atomic
rename, so a partial file is never cached and a restart resumes. **The shard size is part of the cache
filename** -- changing `--shard-size` must not silently reuse mismatched shards (a trap v1 hit).

Output
------
    data/processed/embeddings/embeddings_<species>.npz
        accessions   object   (n,)       row-aligned to `embeddings`
        embeddings   float32  (n, 1152)  mean-pooled, BOS/EOS stripped
        model        scalar   "esmc_600m"
        pooling      scalar   "mean"
        dim          scalar   1152

Load it with `src.embeddings.load()` rather than by hand -- `accessions` is an object array, so a raw
`np.load` needs `allow_pickle=True`.

Run with the `gradi` env (NOT `gradi-loc`, whose `fair-esm` claims the same top-level `esm` package):
    python scripts/embeddings/esmc.py                          # the three bacteria
    python scripts/embeddings/esmc.py --limit 20 --species saureus   # smoke test
    python scripts/embeddings/esmc.py --species kpneumoniae --refresh
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

# Let any op unsupported on Apple MPS fall back to CPU instead of crashing.
# MUST be set before torch is imported -- order-sensitive, and silent if you get it wrong.
os.environ.setdefault("PYTORCH_ENABLE_MPS_FALLBACK", "1")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import torch  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402
from src import matrices as M  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
NCBI_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "ncbi"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
SHARD_DIR = SCRATCH_DIR / "shards_esmc"

MODEL_ID = "esmc_600m"
EMBED_DIM = 1152
POOLING = "mean"
DEFAULT_SPECIES = ("kpneumoniae", "ecoli", "saureus")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def pick_device(choice: str) -> str:
    if choice != "auto":
        return choice
    if torch.cuda.is_available():
        return "cuda"
    mps = getattr(torch.backends, "mps", None)
    if mps is not None and mps.is_available():
        return "mps"
    return "cpu"


def load_client(device: str):
    """Load ESM-C, with actionable messages for the two failure modes v1 documented."""
    try:
        from esm.models.esmc import ESMC
    except ImportError as exc:
        sys.exit(
            f"could not import the EvolutionaryScale `esm` package ({exc}).\n"
            "Run this with the `gradi` env: ~/miniconda3/envs/gradi/bin/python\n"
            "NOTE `gradi-loc` will NOT work -- its `fair-esm` claims the same top-level `esm` name."
        )
    try:
        return ESMC.from_pretrained(MODEL_ID).to(device).eval()
    except Exception as exc:  # noqa: BLE001 - we want to classify, then re-raise usefully
        text = str(exc).lower()
        if any(k in text for k in ("401", "403", "gated", "token", "authenticate", "login")):
            sys.exit(
                f"ESM-C weights look gated ({exc}).\n"
                "Authenticate once with `huggingface-cli login`, or export HF_TOKEN=..."
            )
        raise


def embed_one(client, seq: str) -> np.ndarray:
    """Mean-pooled embedding for one sequence, BOS/EOS stripped."""
    from esm.sdk.api import ESMProtein, LogitsConfig

    with torch.no_grad():
        pt = client.encode(ESMProtein(sequence=seq))
        out = client.logits(pt, LogitsConfig(sequence=True, return_embeddings=True))
        # embeddings[0] is (L+2, 1152): BOS, residues..., EOS. Pool over residues only.
        vec = out.embeddings[0][1:-1].mean(dim=0)
    return vec.float().cpu().numpy()


def shard_path(species: str, shard_size: int, idx: int) -> Path:
    # shard_size is in the name on purpose: changing it must not reuse mismatched shards.
    return SHARD_DIR / f"embeddings_{species}_{shard_size}_{idx:04d}.npz"


def load_strain_frame(faa: Path) -> pd.DataFrame:
    """(uniprot_ac, sequence) from a plain protein FASTA, keyed on its record id.

    For TRAINING proteomes that are not one of the four registry species -- the tier-D screen
    strains. `uniprot_ac` keeps its name so everything downstream is identical, but the ids here
    are RefSeq `WP_*` protein accessions, not UniProt. That is deliberate: the column is the row
    key, and renaming it per source would fork every consumer.
    """
    ids, seqs, cur, buf = [], [], None, []
    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                ids.append(cur); seqs.append("".join(buf))
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if cur:
        ids.append(cur); seqs.append("".join(buf))
    return pd.DataFrame({"uniprot_ac": ids, "sequence": seqs})


def run_species(client, species: str, shard_size: int, limit: int | None,
                refresh: bool, frame: pd.DataFrame | None = None,
                out_dir: Path | None = None, canonical: bool = True) -> dict:
    t0 = time.time()
    df = frame if frame is not None else P.load(species)[["uniprot_ac", "sequence"]]
    df = df[df["sequence"].str.strip().ne("")].sort_values("uniprot_ac").reset_index(drop=True)
    if limit:
        df = df.head(limit)
    n_total = len(df)
    residues = int(df["sequence"].str.len().sum())
    say(f"  {species}: {n_total} proteins, {residues:,} residues")

    shards = [df.iloc[i:i + shard_size] for i in range(0, n_total, shard_size)]
    if refresh:
        for i in range(len(shards)):
            shard_path(species, shard_size, i).unlink(missing_ok=True)
    pending = [i for i in range(len(shards))
               if not shard_path(species, shard_size, i).exists()]
    say(f"    {len(shards)} shards of {shard_size}; {len(pending)} to compute, "
        f"{len(shards) - len(pending)} cached")

    skipped: list[str] = []
    done_res = 0
    for pos, i in enumerate(pending, 1):
        part = shards[i]
        s0 = time.time()
        accs, vecs = [], []
        for acc, seq in zip(part["uniprot_ac"], part["sequence"]):
            try:
                vecs.append(embed_one(client, seq))
                accs.append(acc)
            except Exception as exc:  # noqa: BLE001 - name it, never swallow it
                skipped.append(f"{acc} ({len(seq)} aa): {type(exc).__name__}: {exc}")
        tmp = shard_path(species, shard_size, i).with_suffix(".tmp.npz")
        np.savez_compressed(tmp,
                            accessions=np.array(accs, dtype=object),
                            embeddings=np.vstack(vecs).astype(np.float32) if vecs
                            else np.zeros((0, EMBED_DIM), dtype=np.float32))
        tmp.rename(shard_path(species, shard_size, i))   # atomic: no partial file is ever cached
        done_res += int(part["sequence"].str.len().sum())
        dt = time.time() - s0
        rate = done_res / max(time.time() - t0, 1e-9)
        left = residues * (len(pending) - pos) / max(len(pending), 1) / max(rate, 1e-9)
        say(f"    shard {pos}/{len(pending)}  {len(accs)} proteins in {dt:5.1f}s  "
            f"({rate:,.0f} aa/s, ~{left / 60:.0f} min left)")

    # ---- assemble from the shard cache
    all_accs, all_vecs = [], []
    for i in range(len(shards)):
        z = np.load(shard_path(species, shard_size, i), allow_pickle=True)
        all_accs.extend(z["accessions"].tolist())
        if z["embeddings"].shape[0]:
            all_vecs.append(z["embeddings"])
    mat = np.vstack(all_vecs).astype(np.float32) if all_vecs \
        else np.zeros((0, EMBED_DIM), dtype=np.float32)

    # Canonical row order: every matrix in this project has the same rows in the same order, so
    # an embedding can be hstacked onto any other axis with no join. Shards are assembled in shard
    # order, which is whatever order the FASTA happened to be in -- reorder once, here at the end.
    if canonical:
        all_accs, mat = M.reindex_arrays(all_accs, mat, species)
    else:
        # A screen strain is NOT one of the four anchor proteomes, so it has no canonical order to
        # match and `src.matrices` must not be asked for one. These are training features, not a
        # deliverable matrix -- which is also why they are written under scratch/.
        order = np.argsort(np.array(all_accs, dtype=object))
        all_accs = [all_accs[i] for i in order]
        mat = mat[order]

    out = (out_dir or OUT_DIR) / f"embeddings_{species}.npz"
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out,
                        accessions=np.array(all_accs, dtype=object),
                        embeddings=mat,
                        model=np.array(MODEL_ID),
                        pooling=np.array(POOLING),
                        dim=np.array(EMBED_DIM))
    elapsed = time.time() - t0
    say(f"    wrote {out.name}  {mat.shape[0]} x {mat.shape[1]}  "
        f"({out.stat().st_size / 1e6:.1f} MB, {elapsed / 60:.1f} min)")
    if skipped:
        say(f"    !! {len(skipped)} proteins FAILED:")
        for s in skipped[:10]:
            say(f"       {s}")
    return dict(species=species, n=mat.shape[0], n_expected=n_total, dim=mat.shape[1],
                residues=residues, model=MODEL_ID, pooling=POOLING,
                seconds=round(elapsed, 1), skipped=len(skipped),
                sha256=hashlib.sha256(mat.tobytes()).hexdigest()[:16],
                fetched_at=datetime.now(timezone.utc).isoformat(timespec="seconds"))


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES), choices=list(P.SPECIES),
                    help="species to embed (default: the three bacteria; human is out of scope)")
    ap.add_argument("--strain", nargs="+", metavar="LABEL",
                    help="embed a tier-D screen strain by registry label instead of a species. "
                         "Reads data/source/uniprot/proteomes/ncbi/<LABEL>.faa and writes to "
                         "scratch/strains/ -- these are TRAINING FEATURES, not a deliverable "
                         "matrix, so they carry no canonical row order")
    ap.add_argument("--shard-size", type=int, default=250)
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    ap.add_argument("--limit", type=int, help="only the first N proteins per species (smoke test)")
    ap.add_argument("--refresh", action="store_true", help="discard cached shards and recompute")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    device = pick_device(args.device)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SHARD_DIR.mkdir(parents=True, exist_ok=True)

    rule("=")
    say("STAGE 01 - ESM-C 600M mean-pooled embeddings")
    rule("=")
    say(f"  in       : data/processed/proteomes/<species>.tsv")
    say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/  (+ scratch/shards_esmc/)")
    say(f"  model    : {MODEL_ID}  ({EMBED_DIM}-dim, {POOLING}-pooled)")
    say(f"  {'strains ' if args.strain else 'species '} : "
        f"{', '.join(args.strain or args.species)}")
    say(f"  device   : {device}   torch {torch.__version__}")
    say(f"  shards   : {args.shard_size} proteins each")
    if args.limit:
        say(f"  limit    : {args.limit} proteins per species (SMOKE TEST)")
    say()

    t0 = time.time()
    say("loading ESM-C ...")
    client = load_client(device)
    say(f"  ready in {time.time() - t0:.1f}s\n")

    rule()
    say("EMBED")
    rule()
    if args.strain:
        strain_dir = OUT_DIR / "scratch" / "strains"
        rows = []
        for lbl in args.strain:
            faa = NCBI_DIR / f"{lbl}.faa"
            if not faa.exists():
                sys.exit(f"FATAL {faa} missing -- run scripts/proteomes/download.py --tier D "
                         f"--only {lbl}")
            rows.append(run_species(client, lbl, args.shard_size, args.limit, args.refresh,
                                    frame=load_strain_frame(faa), out_dir=strain_dir,
                                    canonical=False))
    else:
        rows = [run_species(client, sp, args.shard_size, args.limit, args.refresh)
                for sp in args.species]
    say()

    man = pd.DataFrame(rows)
    man.insert(1, "device", device)
    man.to_csv(EVIDENCE_DIR / "esmc_manifest.tsv", sep="\t", index=False)

    rule()
    say("SUMMARY")
    rule()
    say(f"  {'species':<14} {'n':>7} {'expected':>9} {'dim':>5} {'residues':>11} "
        f"{'min':>7} {'skipped':>8}")
    for r in rows:
        flag = "" if r["n"] == r["n_expected"] else "   <- MISMATCH"
        say(f"  {r['species']:<14} {r['n']:>7} {r['n_expected']:>9} {r['dim']:>5} "
            f"{r['residues']:>11,} {r['seconds'] / 60:>7.1f} {r['skipped']:>8}{flag}")
    total_min = sum(r["seconds"] for r in rows) / 60
    say(f"\n  total {total_min:.1f} min   manifest -> evidence/esmc_manifest.tsv")
    bad = [r for r in rows if r["n"] != r["n_expected"] or r["skipped"]]
    if bad:
        sys.exit(f"\nFAILED: {len(bad)} species did not embed completely -- see above.")
    rule("=")
    say("stage 01 complete.")
    rule("=")


if __name__ == "__main__":
    main()
