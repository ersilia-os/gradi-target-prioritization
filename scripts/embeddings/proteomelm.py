"""ProteomeLM contextualised embeddings -- proteome in, embedding matrix out (stage 07).

ProteomeLM (Malbranke, Zalaffi & Bitbol, PNAS 2026, 10.1073/pnas.2524201123; papers in
docs/papers/) is a set-transformer over a WHOLE proteome: one token per protein, each token being
that protein's mean-pooled ESM-C 600M embedding. Its output is what ProteomeLM-Ess consumes, so this
script is the front half of the essentiality axis. The label corpus is the other half, already built
by essentiality/labels.py and essentiality/deg_proteomes.py.

SELF-CONTAINED FROM THE SEQUENCES ONWARD
----------------------------------------
Sequences come from stage 00 (`src.proteomes.load`), which is the repo's canonical proteome table.
Everything after that is computed here: this does NOT import src/embeddings.py, does NOT read stage
01's .npz, and does NOT touch stage 05's OrthoDB tables. So the ESM-C pooling convention and the
ProteomeLM forward are ours alone, and the comparison against stage 01 stays a CONTROL rather than a
dependency -- agreement validates both pipelines precisely because neither feeds the other.

An earlier version downloaded the proteomes from UniProt itself. That was dropped because it was
measured to be pointless: the download was byte-identical to stage 00 (4,403/4,403 E. coli sequences,
no extras either way), while adding a real failure mode -- UniProt's stream endpoint is chunked, has
no Content-Length to verify against, and dropped mid-transfer on the 5,728-protein K. pneumoniae
fetch (`http.client.IncompleteRead`). Reading the local table cannot fail that way.

Why there is no OrthoDB here
----------------------------
ProteomeLM takes a per-protein "functional encoding" (`group_embeds`) which, during TRAINING, is a
vector sampled from the protein's OrthoDB ancestral path. At INFERENCE the authors' own released code
does not use OrthoDB at all: `prepare_ppi(..., use_odb: bool = False, ...)` is the default (annotated
`# TODO: use odb on the fly`), `ppi/feature_extraction.py` hard-codes `use_odb=False` twice, and
`utils/embedding.py:163` returns the ESM-C vectors AS `group_embeds` with the comment "to avoid
relying on ODB". Their reference notebook does the same. So each protein's own ESM-C embedding is its
functional encoding -- the documented fallback -- and `group_embeds_mode` records that in every
output file so it can never be mistaken for the OrthoDB variant.

The recipe, and where each number comes from
--------------------------------------------
  1. proteome FASTA from UniProt, UNFILTERED (see the reviewed-only trap below)
  2. ESM-C 600M, mean over NON-PAD tokens -- so BOS/EOS are INSIDE the mean, which is ProteomeLM's
     own convention (`average_representation`), not stage 01's. Measured difference against stripping
     them: cosine 0.999970 median, 0.998291 worst. Negligible, but this tool is faithful.
  3. ProteomeLM-L forward over the whole proteome at once, `group_embeds=None`
  4. hidden_states[8], z-scored with the genome-wide mean/SD -- "the best performing version of
     ProteomeLM-Ess is the one trained on the embeddings of layer 8 of ProteomeLM-L"

Measured facts, so nobody re-derives them
-----------------------------------------
  * the forward is CHEAP: 1.3 s for all 4,403 E. coli proteins with ProteomeLM-M on CPU fp32; L is
    ~3-4x that (18 layers at width 1152 vs 12 at 768). Still seconds.
  * protein ORDER IS IRRELEVANT: permuting the input and un-permuting the output agrees to
    max|diff| 1.1e-05. There are no positional embeddings; `max_position_embeddings: 512` in the HF
    config is an inert inherited DistilBert field, NOT a proteome-size cap.
  * proteome CONTEXT MATTERS: the same 2,000 proteins embedded inside the full E. coli proteome
    versus alone differ at cosine 0.9676. That difference is the entire reason this stage exists,
    and it is why `--limit` is not a free smoke test.

Traps
-----
  * **Do not fetch with `reviewed:true` for bacteria.** The paper's own `download_proteome` defaults
    to it; for K. pneumoniae HS11286 that returns a handful of entries instead of 5,728. Human is the
    opposite -- its unfiltered proteome is 147,506 TrEMBL-bloated entries -- so human alone is
    fetched reviewed-only, and that asymmetry is asserted, not assumed.
  * **The ProteomeLM forward is NOT shardable.** It runs over the whole proteome at once, so a shard
    boundary would change the values. Only the ESM-C step shards. Do not "fix" this by copying stage
    01's cache.
  * `output_attentions=False` matters: it keeps SDPA in play so the N x N attention matrix is never
    materialised. N = 5,728 for Kp.

Output
  data/source/uniprot/proteome_fasta/<label>/proteome.fasta        as downloaded, + SOURCE.md
  data/processed/embeddings/
      proteomelm_<label>.npz     accessions - embeddings (n,dim) - model - layer - n_layers - dim -
                                 group_embeds_mode - esmc_pooling - proteome_id
      scratch/esmc_<label>.npz the expensive intermediate, kept
      scratch/shards_proteomelm/          ESM-C resume cache
      evidence/proteomelm_manifest.tsv
  Load with src.proteomelm.load() rather than by hand -- the accessions array needs allow_pickle.

Run with the `gradi` env (NOT `gradi-loc`: its fair-esm claims the same top-level `esm` name).
  python scripts/embeddings/proteomelm.py --dry-run
  python scripts/embeddings/proteomelm.py --proteome UP000000625 --label ecoli --limit 400
  python scripts/embeddings/proteomelm.py
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
sys.path.insert(0, str(REPO_ROOT))   # only so `src.proteomes` can supply the input sequences
from src import matrices as M  # noqa: E402
OUT_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
SHARD_DIR = SCRATCH_DIR / "shards_proteomelm"

ESMC_MODEL = "esmc_600m"
EMBED_DIM = 1152
# ProteomeLM's own pooling: mean over non-pad tokens, so BOS/EOS are included. Named in the output.
ESMC_POOLING = "mean_with_bos_eos"
PROTEOMELM_REPO = "Bitbol-Lab/ProteomeLM-{size}"
DEFAULT_SIZE = "L"     # the paper's best for essentiality; M being cached is not a reason to use it
DEFAULT_LAYER = 8      # layer 8 of L's 18 -- normalised depth 0.44
GROUP_EMBEDS_MODE = "self"

# label -> (UniProt proteome id, expected protein count, reviewed-only?)
# The counts are the assertion targets; a silent short download is the failure mode this guards.
PROTEOMES = {
    "kpneumoniae": ("UP000007841", 5728, False),
    "ecoli": ("UP000000625", 4403, False),
    "saureus": ("UP000008816", 2889, False),
    "human": ("UP000005640", 20416, True),   # reviewed-only: unfiltered is 147,506 TrEMBL entries
}
DEFAULT_LABELS = ("kpneumoniae", "ecoli", "saureus")

MIN_DIM_SD = 1e-4          # below this the matrix has collapsed and is not an embedding
PERM_TOLERANCE = 1e-3      # measured 1.1e-05; a real break would be orders of magnitude worse

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
    mps = getattr(torch.backends, "mps", None)   # older torch has no .mps attribute at all
    if mps is not None and mps.is_available():
        return "mps"
    return "cpu"


def read_proteome(label: str, limit: int | None) -> tuple[list[str], list[str]]:
    """Sequences from stage 00, sorted by accession. Returns (accessions, sequences).

    Sorting is deliberate: it makes shard boundaries deterministic across runs and independent of
    whatever order the source table happens to be in.
    """
    from src import proteomes as P  # noqa: PLC0415 - the one pipeline dependency, and only for input

    pid, n_expected, _ = PROTEOMES[label]
    df = P.load(label)[["uniprot_ac", "sequence"]]
    df = df[df["sequence"].str.strip().ne("")].sort_values("uniprot_ac").reset_index(drop=True)
    if len(df) != n_expected:
        sys.exit(f"FATAL {label}: stage 00 has {len(df):,} proteins with a sequence, expected "
                 f"{n_expected:,}. Re-run scripts/proteomes/download.py.")
    if df["uniprot_ac"].duplicated().any():
        sys.exit(f"FATAL {label}: duplicate accessions in the stage-00 table")
    accs = df["uniprot_ac"].tolist()
    seqs = df["sequence"].tolist()
    say(f"  {label:<13} {len(accs):>6,} proteins, {sum(len(x) for x in seqs):>10,} residues  "
        f"({pid}, from stage 00)")
    if limit:
        accs, seqs = accs[:limit], seqs[:limit]
    return accs, seqs


def load_esmc(device: str):
    """Load ESM-C, classifying the two failure modes v1 documented."""
    try:
        from esm.models.esmc import ESMC
    except ImportError as exc:
        sys.exit(
            f"could not import the EvolutionaryScale `esm` package ({exc}).\n"
            "Run this with the `gradi` env: ~/miniconda3/envs/gradi/bin/python\n"
            "NOTE `gradi-loc` will NOT work -- its `fair-esm` claims the same top-level `esm` name."
        )
    try:
        return ESMC.from_pretrained(ESMC_MODEL).to(device).eval()
    except Exception as exc:  # noqa: BLE001 - classify, then re-raise usefully
        text = str(exc).lower()
        if any(k in text for k in ("401", "403", "gated", "token", "authenticate", "login")):
            sys.exit(f"ESM-C weights look gated ({exc}).\n"
                     "Authenticate once with `huggingface-cli login`, or export HF_TOKEN=...")
        raise


def shard_path(label: str, shard_size: int, idx: int) -> Path:
    # shard_size is in the name on purpose: changing it must not reuse mismatched shards.
    return SHARD_DIR / f"esmc_{label}_{shard_size}_{idx:04d}.npz"


def esmc_embeddings(label: str, accs: list[str], seqs: list[str], device: str,
                    shard_size: int, refresh: bool) -> np.ndarray:
    """ESM-C 600M, mean over NON-PAD tokens -- ProteomeLM's pooling, BOS/EOS included.

    The only slow step (~1,800 aa/s measured), so it caches. The cache is keyed BY ACCESSION, not by
    position: shards store the accessions they contain and assembly indexes by name. That means
    re-ordering the input, or adding proteins to it, reuses everything already computed instead of
    invalidating it -- which a position-keyed cache would, and which cost 15 minutes of recomputation
    once already.

    An accession-keyed cache has one dangerous failure mode: if the SEQUENCE behind an accession ever
    changes (a UniProt revision, a different proteome release, a switch of sequence source), the
    stale vector still matches by name and is served silently. So every shard also stores a sha256 of
    each sequence, verified on load; a mismatch is treated as not-cached and recomputed, and it is
    named in the log rather than absorbed. Shards written before this guard existed carry no hash and
    are reported as unverified -- the ones in this repo were checked by recomputation (45 proteins
    across the three species, max|diff| exactly 0.0) before that was accepted.
    """
    SHARD_DIR.mkdir(parents=True, exist_ok=True)
    if refresh:
        for f in SHARD_DIR.glob(f"esmc_{label}_*.npz"):
            f.unlink(missing_ok=True)

    want_sha = {a: hashlib.sha256(seq.encode()).hexdigest()[:16]
                for a, seq in zip(accs, seqs)}
    lut: dict[str, np.ndarray] = {}
    stale, unverified = [], 0
    for f in sorted(SHARD_DIR.glob(f"esmc_{label}_*.npz")):
        try:
            z = np.load(f, allow_pickle=True)
            shas = z["seq_sha"].astype(str) if "seq_sha" in z.files else None
            if shas is None:
                unverified += len(z["accessions"])
            for k, (a, v) in enumerate(zip(z["accessions"].astype(str), z["embeddings"])):
                if shas is not None and a in want_sha and shas[k] != want_sha[a]:
                    stale.append(a)          # sequence changed: the cached vector is for the old one
                    continue
                lut[a] = v
        except Exception as exc:  # noqa: BLE001 - a corrupt shard is recomputed, never trusted
            say(f"    unreadable shard {f.name} ({type(exc).__name__}); ignoring it")
    if unverified:
        say(f"    {unverified:,} cached embeddings predate the sequence-hash guard "
            "(verified by recomputation; use --refresh to rebuild from scratch)")
    if stale:
        say(f"    !! {len(stale):,} cached embeddings are for a DIFFERENT sequence than stage 00 "
            f"now has, e.g. {stale[:3]} -- recomputing those")
    todo = [i for i, a in enumerate(accs) if a not in lut]
    say(f"  cache: {len(lut):,} embeddings on disk; {len(todo):,} of {len(accs):,} to compute")

    if todo:
        model = load_esmc(device)
        pad = model.tokenizer.pad_token_id
        residues_total = sum(len(seqs[i]) for i in todo)
        batches = [todo[i:i + shard_size] for i in range(0, len(todo), shard_size)]
        existing = len(list(SHARD_DIR.glob(f"esmc_{label}_*.npz")))
        t0, done_res = time.time(), 0
        for pos, batch in enumerate(batches, 1):
            vecs, s0 = [], time.time()
            with torch.no_grad():
                for i in batch:
                    ids = model._tokenize([seqs[i]]).long().to(device)
                    out = model(ids).embeddings              # (1, L+2, 1152)
                    mask = (ids != pad).unsqueeze(-1)
                    summed = (out * mask).sum(dim=1)
                    count = mask.sum(dim=1).clamp(min=1)
                    vecs.append((summed / count).squeeze(0).float().cpu().numpy())
            mat = np.vstack(vecs).astype(np.float32)
            dest = SHARD_DIR / f"esmc_{label}_{shard_size}_{existing + pos - 1:04d}.npz"
            tmp = dest.with_suffix(".tmp.npz")
            np.savez_compressed(
                tmp, accessions=np.array([accs[i] for i in batch], dtype=object),
                embeddings=mat,
                seq_sha=np.array([want_sha[accs[i]] for i in batch], dtype=object))
            tmp.rename(dest)      # atomic: a partial shard is never cached
            for a, v in zip((accs[i] for i in batch), mat):
                lut[a] = v
            done_res += sum(len(seqs[i]) for i in batch)
            rate = done_res / max(time.time() - t0, 1e-9)
            left = (residues_total - done_res) / max(rate, 1e-9)
            say(f"    batch {pos}/{len(batches)}  {len(batch)} proteins in {time.time() - s0:5.1f}s "
                f"({rate:,.0f} aa/s, ~{left / 60:.0f} min left)")

    missing = [a for a in accs if a not in lut]
    if missing:
        sys.exit(f"FATAL {label}: {len(missing):,} proteins never embedded, e.g. {missing[:3]}")
    return np.vstack([lut[a] for a in accs]).astype(np.float32)


def load_proteomelm(size: str, device: str):
    from proteomelm import ProteomeLMForMaskedLM
    repo = PROTEOMELM_REPO.format(size=size)
    model = ProteomeLMForMaskedLM.from_pretrained(repo).to(device).float().eval()
    return model, repo


def proteomelm_forward(model, esmc: np.ndarray, layer: int, device: str) -> np.ndarray:
    """One forward over the WHOLE proteome. Not shardable -- a shard boundary changes the values.

    `group_embeds=None` makes each protein its own functional encoding, which is what the model does
    internally (`if group_embeds is None: group_embeds = inputs_embeds.clone()`) and what the
    authors' released inference path does.
    `output_attentions=False` keeps SDPA in play, so the N x N attention is never materialised.
    """
    x = torch.tensor(esmc, dtype=torch.float32, device=device).unsqueeze(0)
    with torch.no_grad():
        out = model(inputs_embeds=x, group_embeds=None,
                    output_hidden_states=True, output_attentions=False)
    hs = out.hidden_states
    if not (0 <= layer < len(hs)):
        sys.exit(f"FATAL layer {layer} out of range: this model exposes {len(hs)} hidden states "
                 f"(0..{len(hs) - 1})")
    return hs[layer].squeeze(0).float().cpu().numpy()


def zscore_genome(mat: np.ndarray) -> np.ndarray:
    """Normalise with the GENOME-WIDE mean and SD -- 'Protein embeddings are normalized using the
    genome-wide mean and SD before being given as input to ProteomeLM-Ess'."""
    mu, sd = mat.mean(axis=0, keepdims=True), mat.std(axis=0, keepdims=True)
    return (mat - mu) / np.maximum(sd, 1e-8)


def check_permutation(model, esmc: np.ndarray, layer: int, device: str, base: np.ndarray) -> float:
    """Order must not matter: there are no positional embeddings. Measured 1.1e-05."""
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(esmc))
    shuffled = proteomelm_forward(model, esmc[perm], layer, device)
    back = np.empty_like(shuffled)
    back[perm] = shuffled
    return float(np.abs(base - back).max())


def check_context(model, esmc: np.ndarray, layer: int, device: str, base: np.ndarray) -> float:
    """Context must matter, or ProteomeLM has collapsed to a per-protein encoder. Measured 0.9676."""
    half = max(2, len(esmc) // 2)
    sub = proteomelm_forward(model, esmc[:half], layer, device)
    a, b = base[:half], sub
    cos = (a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-12)
    return float(np.median(cos))


def crosscheck_stage01(label: str, accs: list[str], esmc: np.ndarray) -> str:
    """REPORTED, never depended on. Stage 01 strips BOS/EOS; we keep them. Expect cosine ~0.99997.

    Agreement validates both pipelines exactly because neither feeds the other. A missing stage-01
    file is not a failure -- this tool is self-contained.
    """
    p = REPO_ROOT / "data" / "processed" / "embeddings" / f"embeddings_{label}.npz"
    if not p.exists():
        return "stage 01 absent (fine -- this tool is self-contained)"
    z = np.load(p, allow_pickle=True)
    lut = {a: i for i, a in enumerate(z["accessions"].astype(str))}
    idx = [(i, lut[a]) for i, a in enumerate(accs) if a in lut]
    if not idx:
        return "no shared accessions with stage 01"
    ours = esmc[[i for i, _ in idx]]
    theirs = z["embeddings"][[j for _, j in idx]].astype(np.float32)
    cos = (ours * theirs).sum(1) / (np.linalg.norm(ours, axis=1) * np.linalg.norm(theirs, axis=1)
                                    + 1e-12)
    return (f"{len(idx):,} shared: cosine median {np.median(cos):.6f}, min {cos.min():.6f} "
            f"(expected ~0.99997 -- the BOS/EOS convention difference)")


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--label", nargs="+", default=list(DEFAULT_LABELS), choices=list(PROTEOMES),
                   help="which proteomes (default: the three bacteria; human needs reviewed-only)")
    ap.add_argument("--size", default=DEFAULT_SIZE, choices=["XS", "S", "M", "L"],
                   help=f"ProteomeLM size (default {DEFAULT_SIZE} -- the paper's best for essentiality)")
    ap.add_argument("--layer", type=int, default=DEFAULT_LAYER,
                   help=f"hidden state to keep (default {DEFAULT_LAYER}, i.e. layer 8 of L's 18)")
    ap.add_argument("--shard-size", type=int, default=250)
    ap.add_argument("--device", choices=["auto", "cuda", "mps", "cpu"], default="auto")
    ap.add_argument("--limit", type=int,
                   help="only the first N proteins (SMOKE TEST -- changes the values, see below)")
    ap.add_argument("--refresh", action="store_true", help="re-download and recompute")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and write nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    device = pick_device(args.device)
    pre = "smoke_" if args.limit else ""

    rule("=")
    say("STAGE 07 - ProteomeLM contextualised embeddings")
    rule("=")
    say(f"  proteomes : {', '.join(args.label)}")
    say(f"  esm-c     : {ESMC_MODEL}  ({EMBED_DIM}-dim, pooling={ESMC_POOLING})")
    say(f"  proteomelm: {PROTEOMELM_REPO.format(size=args.size)}  layer {args.layer}, "
        f"z-scored genome-wide")
    say(f"  group_embeds: {GROUP_EMBEDS_MODE}  (the authors' own inference default -- no OrthoDB)")
    say(f"  device    : {device}   torch {torch.__version__}")
    say("  sequences : data/processed/proteomes/proteome_<label>.tsv  (stage 00)")
    say(f"  out       : {OUT_DIR.relative_to(REPO_ROOT)}/{pre}proteomelm_<label>.npz")
    if args.limit:
        say("")
        say(f"  !! --limit {args.limit}: this is NOT a free smoke test. ProteomeLM embeds a protein")
        say("     IN THE CONTEXT of the whole proteome, so truncating changes the values of the")
        say("     proteins that remain (measured: cosine 0.9676 full vs half). Output goes to")
        say("     scratch/smoke_* and is NOT comparable to a full run.")
    rule("=")
    if args.dry_run:
        say("dry run -- nothing fetched, nothing written.")
        for lab in args.label:
            pid, n, rev = PROTEOMES[lab]
            say(f"    {lab:<13} {pid}  {n:>6} proteins  reviewed={rev}")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    rule()
    say("SEQUENCES  (stage 00)")
    rule()
    data = {lab: read_proteome(lab, args.limit) for lab in args.label}

    rule()
    say(f"ESM-C  ({ESMC_POOLING})")
    rule()
    esmc = {}
    for lab, (accs, seqs) in data.items():
        say(f"  {lab}: {len(accs):,} proteins")
        # No limit-specific cache key: the cache is accession-keyed, so a --limit run warms the
        # same store the full run reads.
        esmc[lab] = esmc_embeddings(lab, accs, seqs, device, args.shard_size, args.refresh)

    rule()
    say(f"PROTEOMELM-{args.size}  (whole proteome per forward; NOT shardable)")
    rule()
    t0 = time.time()
    model, repo = load_proteomelm(args.size, device)
    n_layers = int(model.config.n_layers)
    say(f"  loaded {repo}: dim {model.config.dim}, {n_layers} layers, "
        f"{sum(p.numel() for p in model.parameters()):,} params  ({time.time() - t0:.0f}s)")

    rows = []
    for lab, (accs, _) in data.items():
        t1 = time.time()
        raw = proteomelm_forward(model, esmc[lab], args.layer, device)
        mat = zscore_genome(raw)
        secs = time.time() - t1
        say(f"  {lab:<13} {mat.shape}  in {secs:.1f}s")

        # --- controls
        if not np.isfinite(mat).all():
            sys.exit(f"FATAL {lab}: non-finite values in the embedding")
        sd_min = float(mat.std(axis=0).min())
        if sd_min < MIN_DIM_SD:
            sys.exit(f"FATAL {lab}: a dimension has SD {sd_min:.2e} -- the matrix has collapsed")
        perm = check_permutation(model, esmc[lab], args.layer, device, raw)
        if perm > PERM_TOLERANCE:
            sys.exit(f"FATAL {lab}: permutation changed the output by {perm:.2e} "
                     f"(> {PERM_TOLERANCE}); the wiring is wrong")
        ctx = check_context(model, esmc[lab], args.layer, device, raw)
        say(f"                permutation max|diff| {perm:.2e} (ok)   "
            f"context cosine full-vs-half {ctx:.4f}")
        if ctx > 0.999:
            say("                !! context cosine ~1.0: ProteomeLM is behaving as a per-protein "
                "encoder, which would remove the reason for this stage")
        say(f"                stage-01 cross-check: {crosscheck_stage01(lab, accs, esmc[lab])}")

        out = (SCRATCH_DIR if args.limit else OUT_DIR) / f"{pre}proteomelm_{lab}.npz"
        # Canonical row order -- see src/matrices.py. Skipped under --limit, whose smoke output is
        # a subset by construction and so cannot be canonical.
        if not args.limit:
            accs, mat = M.reindex_arrays(accs, mat, lab)
        np.savez_compressed(
            out,
            accessions=np.array(accs, dtype=object), embeddings=mat,
            model=np.array(repo), layer=np.array(args.layer), n_layers=np.array(n_layers),
            dim=np.array(mat.shape[1]), group_embeds_mode=np.array(GROUP_EMBEDS_MODE),
            esmc_pooling=np.array(ESMC_POOLING), proteome_id=np.array(PROTEOMES[lab][0]))
        rows.append(dict(label=lab, proteome_id=PROTEOMES[lab][0], n=mat.shape[0],
                         n_expected=len(accs), dim=mat.shape[1], model=repo, layer=args.layer,
                         group_embeds_mode=GROUP_EMBEDS_MODE, esmc_pooling=ESMC_POOLING,
                         perm_maxdiff=f"{perm:.2e}", context_cosine=round(ctx, 4),
                         seconds=round(secs, 1), device=device,
                         sha256=hashlib.sha256(mat.tobytes()).hexdigest()[:16],
                         path=str(out.relative_to(REPO_ROOT)),
                         run_at=datetime.now(timezone.utc).isoformat(timespec="seconds")))

    rule()
    say("SUMMARY")
    rule()
    say(f"  {'label':<13} {'n':>7} {'dim':>5} {'layer':>6} {'perm':>10} {'ctx cos':>8} {'sec':>6}")
    for r in rows:
        flag = "" if r["n"] == r["n_expected"] else "   <- MISMATCH"
        say(f"  {r['label']:<13} {r['n']:>7,} {r['dim']:>5} {r['layer']:>6} "
            f"{r['perm_maxdiff']:>10} {r['context_cosine']:>8.4f} {r['seconds']:>6.1f}{flag}")
    man = pd.DataFrame(rows)
    man_path = (SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}proteomelm_manifest.tsv"
    man.to_csv(man_path, sep="\t", index=False)

    rule()
    say("OUTPUTS")
    rule()
    for r in rows:
        p = REPO_ROOT / r["path"]
        say(f"  {r['path']}  ({p.stat().st_size / 1e6:.1f} MB)")
    say(f"  {man_path.relative_to(REPO_ROOT)}")

    # Integrity re-read: the file on disk is what downstream will consume.
    for r in rows:
        z = np.load(REPO_ROOT / r["path"], allow_pickle=True)
        if z["embeddings"].shape[0] != r["n"] or str(z["group_embeds_mode"]) != GROUP_EMBEDS_MODE:
            sys.exit(f"FATAL {r['label']}: {r['path']} did not round-trip")
    say(f"  re-read check: {len(rows)} file(s), shapes and group_embeds_mode intact")

    bad = [r for r in rows if r["n"] != r["n_expected"]]
    if bad:
        sys.exit(f"\nFAILED: {len(bad)} proteome(s) incomplete -- see above.")
    rule("=")
    say("stage 07 proteomelm embeddings complete.")
    rule("=")


if __name__ == "__main__":
    main()
