"""ProtT5-XL-U50 per-protein embeddings -- a stage-01 sibling, and its validation control.

Stage 01 ships ESM-C 600M and **deliberately defers ProtT5**, for a reason this script does not
overturn: UniProt publishes precomputed ProtT5 vectors for only **7 reference proteomes**, and
*K. pneumoniae* HS11286 and *S. aureus* NCTC 8325 are both absent. Covering them means running
ProtT5-XL-U50 locally, which is what this does.

THE CONTROL IS THE POINT, NOT AN EXTRA
--------------------------------------
CLAUDE.md's standing instruction on ProtT5 is: *"If that is ever done, validate it against UniProt's
E. coli h5 first (same model, same sequences, so cosine must be ~ 1)."* **E. coli MG1655
(`UP000000625_83333`) IS one of the 7**, and it is our exact anchor -- so the same model over the
same sequences must reproduce the same vectors. That is a real, falsifiable check on a local
inference pipeline, and it pins down the one thing UniProt's README leaves unstated: how residues
are pooled to one vector. A wrong pooling produces a perfectly plausible embedding that is not the
published one; only this catches it.

So the control runs first and **nothing is written if it fails**. It is scored on a SAMPLE rather
than all 4,403 E. coli proteins -- deliberately weighted to the LONGEST ones, because that is where
a pooling, truncation or windowing bug actually shows; agreement on a 250-residue protein is nearly
free. `--pooling` variants are scored against the reference so the choice is measured, not assumed.

WHAT IT COSTS, AND WHY NO DOWNLOAD
----------------------------------
The weights are already on disk: TMbed (stage 03, `gradi-loc`) ships
`Rostlab/prot_t5_xl_half_uniref50-enc` -- 1024-d, fp16, encoder-only, 2.4 GB -- inside its package.
Nothing is fetched but UniProt's 10 MB reference h5. Inference runs in `gradi-loc` across a process
boundary (`GRADI_LOC_BIN` overrides), because ProtT5 needs `transformers==4.44.2`: 5.x routes
`T5Tokenizer` through the tiktoken converter and dies with a spurious tiktoken error.

Writes `data/processed/embeddings/prott5_<species>.npz` beside the ESM-C files, plus
`evidence/prott5_control.tsv` -- the per-protein cosine against UniProt, which is the evidence.

Run with the `gradi` env.
  python scripts/embeddings/prott5.py                      # ecoli (control) + saureus
  python scripts/embeddings/prott5.py --species ecoli --compare-pooling
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
RAW_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "prott5_precomputed"
WORKER = REPO_ROOT / "scripts" / "embeddings" / "workers" / "prott5.py"
DEFAULT_LOC_BIN = Path.home() / "miniconda3" / "envs" / "gradi-loc" / "bin" / "python"

EMBED_DIM = 1024
# The 7 proteomes UniProt publishes ProtT5 for. Ours is E. coli; Kp and Sa are NOT here, which is
# the whole reason this script exists.
UNIPROT_REFERENCE = {"ecoli": "UP000000625_83333"}
REFERENCE_URL = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/"
                 "embeddings/{d}/per-protein.h5")

# Calibrated in this script's own run, not inherited. Same model, same sequences: agreement should be
# essentially exact, so anything below this means the local pipeline differs from UniProt's.
COSINE_FLOOR = 0.99

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def loc_bin() -> str:
    p = Path(os.environ.get("GRADI_LOC_BIN", DEFAULT_LOC_BIN))
    if not p.exists():
        sys.exit(f"FATAL no gradi-loc interpreter at {p}. See install.sh -- and do NOT install "
                 "ProtT5's transformers pin into `gradi`.")
    return str(p)


def fetch_reference(species: str, refresh: bool) -> Path:
    """UniProt's own ProtT5 vectors for a proteome it publishes. ~10 MB."""
    import urllib.request

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / f"{UNIPROT_REFERENCE[species]}_per-protein.h5"
    if path.exists() and not refresh:
        say(f"  reference cached: {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
        return path
    url = REFERENCE_URL.format(d=UNIPROT_REFERENCE[species])
    say(f"  fetching {url}")
    tmp = path.with_suffix(".tmp")
    with urllib.request.urlopen(url, timeout=300) as r:
        declared = int(r.headers.get("Content-Length", 0))
        tmp.write_bytes(r.read())
    got = tmp.stat().st_size
    # An HTTP 200 is not evidence of data. Check the bytes against what the server declared.
    if declared and got != declared:
        sys.exit(f"FATAL short read: {got} bytes vs Content-Length {declared}")
    tmp.rename(path)
    say(f"  fetched {got / 1e6:.1f} MB (matches Content-Length)")
    return path


def control_sample(species: str, n: int, seed: int = 0) -> pd.DataFrame:
    """Proteins to validate on: the longest half, plus a random half.

    Length-weighted on purpose. Single-pass agreement is only *established* up to the longest
    protein actually checked, so the sample's maximum length is a reported property of the control,
    not an afterthought.
    """
    prot = P.load(species)[["uniprot_ac", "sequence"]].copy()
    prot["_len"] = prot["sequence"].str.len()
    longest = prot.nlargest(n // 2, "_len")
    rest = prot.drop(longest.index).sample(min(n - len(longest), len(prot) - len(longest)),
                                           random_state=seed)
    return pd.concat([longest, rest]).drop(columns="_len").reset_index(drop=True)


def embed(species: str, binary: str, pooling: str, device: str, refresh: bool,
          subset: pd.DataFrame | None = None, tag: str = "") -> Path:
    """Dispatch a proteome (or a subset, for the control) to the gradi-loc worker."""
    stem = f"prott5_{species}{tag}" + ("" if pooling == "mean_no_eos" else f"_{pooling}")
    out = (EVIDENCE_DIR if tag else OUT_DIR) / f"{stem}.npz"
    if out.exists() and not refresh:
        z = np.load(out, allow_pickle=True)
        say(f"  cached {out.name}: {z['embeddings'].shape}, pooling={z['pooling']}")
        return out

    prot = subset if subset is not None else P.load(species)[["uniprot_ac", "sequence"]]
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    tsv = SCRATCH_DIR / f"prott5_input_{stem}.tsv"
    prot.to_csv(tsv, sep="\t", index=False)
    say(f"  embedding {len(prot):,} proteins ({prot['sequence'].str.len().sum() / 1e6:.2f}M "
        f"residues), pooling={pooling}, device={device}")
    proc = subprocess.run([binary, str(WORKER), "--in", str(tsv), "--out", str(out),
                           "--pooling", pooling, "--device", device,
                           "--shard-dir", str(SCRATCH_DIR / f"shards_{stem}")],
                          text=True)
    if proc.returncode != 0 or not out.exists():
        sys.exit(f"FATAL the ProtT5 worker failed for {species} (exit {proc.returncode})")
    return out


def control(species: str, npz: Path, reference: Path) -> pd.DataFrame:
    """Per-protein cosine between our vectors and UniProt's. The evidence, not a formality."""
    import h5py

    z = np.load(npz, allow_pickle=True)
    ours = {a: v for a, v in zip(z["accessions"].astype(str), z["embeddings"].astype(np.float32))}
    rows = []
    with h5py.File(reference, "r") as f:
        for key in f.keys():
            if key not in ours:
                continue
            ref = np.asarray(f[key], dtype=np.float32)
            mine = ours[key]
            if ref.shape != mine.shape:
                sys.exit(f"FATAL dimension mismatch on {key}: ours {mine.shape} vs "
                         f"UniProt {ref.shape}")
            cos = float(ref @ mine / (np.linalg.norm(ref) * np.linalg.norm(mine)))
            rows.append({"uniprot_ac": key, "cosine": cos})
    if not rows:
        sys.exit("FATAL no accession in UniProt's h5 matched our proteome -- the control cannot "
                 "run, so nothing here is validated")
    return pd.DataFrame(rows).sort_values("cosine").reset_index(drop=True)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=["saureus"])
    ap.add_argument("--control-n", type=int, default=300,
                    help="E. coli proteins to validate against UniProt (half the longest)")
    ap.add_argument("--pooling", default="mean_no_eos",
                    choices=["mean_no_eos", "mean_with_eos"])
    ap.add_argument("--compare-pooling", action="store_true",
                    help="embed E. coli both ways and score each against UniProt")
    ap.add_argument("--device", default="mps", choices=["mps", "cpu"])
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    species = list(dict.fromkeys(args.species))

    rule("=")
    say("ProtT5-XL-U50 per-protein embeddings (stage-01 sibling)")
    rule("=")
    say(f"  model      Rostlab/prot_t5_xl_half_uniref50-enc, {EMBED_DIM}-d, fp16, encoder-only")
    say("  weights    already on disk in gradi-loc (TMbed ships them) -- nothing downloaded")
    say(f"  species    {', '.join(species)}   (E. coli always carries the control, separately)")
    say(f"  control    cosine vs UniProt {UNIPROT_REFERENCE['ecoli']}, floor {COSINE_FLOOR}")
    say("  NOTE       S. aureus and K. pneumoniae are NOT among UniProt's 7 published proteomes,")
    say("             which is why they must be computed and why E. coli must be checked")
    rule("=")

    binary = loc_bin()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rule()
    say("REFERENCE")
    rule()
    ref = fetch_reference("ecoli", args.refresh)

    poolings = ["mean_no_eos", "mean_with_eos"] if args.compare_pooling else [args.pooling]
    rule()
    say("CONTROL - E. coli, ours vs UniProt's own ProtT5")
    rule()
    sample = control_sample("ecoli", args.control_n)
    say(f"  sample     {len(sample)} E. coli proteins, lengths "
        f"{sample.sequence.str.len().min()}-{sample.sequence.str.len().max()} aa "
        f"(half are the longest in the proteome)")
    best, results = None, []
    for pooling in poolings:
        npz = embed("ecoli", binary, pooling, args.device, args.refresh,
                    subset=sample, tag="_control")
        c = control("ecoli", npz, ref)
        med, lo, frac = c.cosine.median(), c.cosine.min(), float((c.cosine >= COSINE_FLOOR).mean())
        longest = int(sample.sequence.str.len().max())
        say(f"  {pooling:<14} n={len(c):>5}  median cosine {med:.6f}  worst {lo:.6f}  "
            f"{100 * frac:.2f}% above {COSINE_FLOOR}  (single-pass validated to {longest} aa)")
        results.append({"pooling": pooling, "n": len(c), "median_cosine": round(med, 6),
                        "min_cosine": round(lo, 6), "frac_above_floor": round(frac, 4)})
        if best is None or med > best[1]:
            best = (pooling, med, c)
    pd.DataFrame(results).to_csv(EVIDENCE_DIR / "prott5_pooling.tsv", sep="\t", index=False)

    pooling, med, c = best
    c.to_csv(EVIDENCE_DIR / "prott5_control.tsv", sep="\t", index=False)
    if med < COSINE_FLOOR:
        sys.exit(f"\nFAILED: best pooling ({pooling}) reaches median cosine {med:.6f}, below the "
                 f"{COSINE_FLOOR} floor. The local ProtT5 pipeline does NOT reproduce UniProt's "
                 "published vectors, so no species is written -- a plausible-looking embedding that "
                 "differs from the published one is worse than none.")
    say(f"  PASS: {pooling} reproduces UniProt at median cosine {med:.6f}")

    rule()
    say("EMBED")
    rule()
    for sp in species:
        embed(sp, binary, pooling, args.device, args.refresh)

    rule()
    say("SUMMARY")
    rule()
    for sp in species:
        p = OUT_DIR / f"prott5_{sp}.npz"
        if p.exists():
            z = np.load(p, allow_pickle=True)
            say(f"  {sp:<12} {z['embeddings'].shape}  pooling={z['pooling']}  "
                f"{p.stat().st_size / 1e6:.1f} MB")
    rule("=")
    say("ProtT5 embeddings complete.")
    rule("=")


if __name__ == "__main__":
    main()
