"""Localization -- one subcellular compartment per protein, from sequence alone.

Stage 03. Reads the stage-00 tables and gives every bacterial protein a compartment call plus a
topology score, using nothing but the amino-acid sequence. No accession lookup, no annotation
transfer, no curated database: the same method for all three species, so the label means the same
thing in each.

Why this axis exists
--------------------
The cytoplasmic Clp machinery can only degrade what it can physically reach, so localization is the
gate on BacPROTAC-style targeting. A protein with no assignment is a protein silently dropped from
every downstream shortlist -- which is exactly what v1 did. Its axis shipped UniProt-curated labels
only, leaving 59.6% of Kp and 48.6% of Ec at `localization = unknown`, and those proteins vanished
from the webapp's `degrader` preset without a word.

Why sequence only
-----------------
Because annotation depth is not evenly distributed and pretending otherwise makes the deliverable
mean different things in different species. UniProt has **1** experimentally-evidenced localization
for Kp HS11286 against 803 for E. coli K-12. v1 patched around that with four more tracks; all four
are unusable here. PSORTb (via PSORTdb) is keyed on RefSeq assembly ids and in v1 never won a single
call, because DeepLocPro always outranked it. SignalP 6.0 is licensed and was never installed.
STEPdb 2.0 is E. coli only. Ortholog transfer needs orthology this pipeline does not have yet -- and
none of them help *S. aureus*, which is new in v2.

So the two tracks that were actually load-bearing in v1 are the whole of this stage.

DeepLocPro -- the compartment
-----------------------------
DeepLocPro 1.0 (Moreno et al. 2024, Bioinformatics 40:btae677): ESM-2 650M with an attention-pooled
head, ensembled over 20 nested-CV checkpoints, emitting a softmax over six prokaryotic compartments.
On the post-2010 Gram-negative benchmark it beats PSORTb 3.0 across the board (accuracy 0.74 vs
0.34, macro-F1 0.75 vs 0.35, MCC 0.69 vs 0.30) and, unlike PSORTb, **always returns a call**. That
is what makes 100% coverage a property of the method rather than something to chase.

    cytoplasm -> cytoplasmic_membrane -> periplasm -> outer_membrane -> cell_wall_surface
    -> extracellular

The Gram-positive trap
----------------------
DeepLocPro takes a Gram group, and `positive` does NOT merely mask the two Gram-negative-only
classes -- `DeepLocPro/utils.py:remap_probabilities` **adds their probability into Extracellular**.
So *S. aureus* has four reachable classes and its `extracellular` count absorbs whatever the model
wanted to call periplasmic. v1 never met this; it had no Gram-positive.

The remap is arithmetic on the output vector, not a different forward pass. So the model runs ONCE,
the raw un-remapped six-vector is persisted to evidence/ + scratch/, and the remap is applied only when
choosing a label. Costs nothing and keeps the inflation measurable -- the SPOT CHECKS section
prints how much mass it moved.

TMbed -- the topology
---------------------
TMbed (Rostlab, Apache-2.0) is not a second classifier. It labels every residue with one of
`H`/`h` transmembrane helix, `B`/`b` transmembrane strand, `S` signal peptide, `i` non-membrane
inside, `o` non-membrane outside, from ProtT5-XL-U50 embeddings. Two numbers come out of that
string and into the deliverable:

    cytoplasmic_fraction   count('i') / length -- the share of the protein on the cytoplasmic side
    has_signal_peptide     'S' in labels -- exported, and therefore not reachable

`cytoplasmic_fraction` is the score: 1.0 wholly cytoplasmic, ~0.5 a membrane protein with real
cytoplasm-facing domains, ~0.0 exported or a barrel. `has_signal_peptide` says *why* a fraction is
near zero -- exported rather than membrane-buried -- and makes one contradiction visible: a
signal-peptide-bearing protein called `cytoplasm` is almost certainly a mis-call.

`--out-format 4`, not v1's format 1, which collapsed `h` into `H` and `b` into `B` and threw away
segment orientation for no saving. Everything else derivable from the label string (helix count,
strand count, longest cytoplasmic run) goes to evidence/ + scratch/ rather than being discarded: the run is
CPU-only hours and cannot be partially redone, so dropping a free number is the expensive choice.
The raw label string is kept for the same reason.

Two environments
----------------
Both predictors need `gradi-loc`, never `gradi`: DeepLocPro depends on `fair-esm`, which claims the
same top-level `esm` package as the EvolutionaryScale `esm` that stage 01 uses for ESM-C. This
script runs in `gradi` and reaches across a process boundary -- `scripts/localization/workers/deeplocpro.py` under
the `gradi-loc` interpreter, and the `tmbed` console script, whose shebang already points there.
Override with GRADI_LOC_BIN. Weights are ~4.9 GB (ESM-2 650M + ProtT5-XL-U50) and already cached.

Scope
-----
The three bacteria only -- 13,020 proteins. Human is out and cannot be otherwise: DeepLocPro is
prokaryote-only, and a eukaryotic call would need DeepLoc 2.x and a disjoint vocabulary.

Output
------
    evidence/deeplocpro_<species>.tsv     uniprot_ac  localization  confidence
    evidence/tmbed_<species>.tsv          uniprot_ac  cytoplasmic_fraction  has_signal_peptide

**Neither is the deliverable.** `scripts/localization/merge.py` stacks the two into
`data/processed/localization/localization_<species>.tsv`, the axis's single table; these two are
the per-predictor record it is built from. They stay separate ON DISK because the two tracks run
independently and are separately resumable -- TMbed alone is CPU-only hours -- so merging inside
this script would mean a `--only` run could not write a deliverable at all.

    evidence/deeplocpro_probabilities_<species>.tsv   raw un-remapped 6-vector, margin, truncated
    evidence/tmbed_topology_<species>.tsv             n_tm_helix, n_tm_strand, cyto_longest, length
    evidence/tmbed_labels_<species>.tsv               the raw per-residue label string
    evidence/deeplocpro_counts_<species>.tsv          per-compartment counts, zeros included
    scratch/deeplocpro_cache/<species>_L<max>/       one JSON per protein (resumable)
    scratch/tmbed_shards/<species>_<size>_<NNNN>.pred   shard cache (resumable)
    evidence/{deeplocpro,tmbed}_manifest.tsv

Two predictors, one table. They answer different questions and neither derives from the other, so
the merged table keeps both side by side and never reduces them to one call: `merge.py` stacks the
columns, it does not arbitrate between them. No accessibility score is computed here -- that is a
modelling decision for whichever stage consumes this.

Run with the `gradi` env:
    python scripts/localization/predict.py                                 # both predictors, 3 species
    python scripts/localization/predict.py --species ecoli --limit 20      # smoke test, ~2 min
    python scripts/localization/predict.py --only deeplocpro               # ~95 min, skip TMbed's hours
    python scripts/localization/predict.py --only tmbed --refresh          # rebuild every shard
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import localization as LOC  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import matrices as M  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "localization"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
WORKER = REPO_ROOT / "scripts" / "localization" / "workers" / "deeplocpro.py"

DEFAULT_SPECIES = ("kpneumoniae", "ecoli", "saureus")
DEFAULT_SHARD_SIZE = 250
TMBED_OUT_FORMAT = "4"

# Neither predictor is installable into `gradi` -- see the module docstring.
DEFAULT_GRADI_LOC_BIN = Path.home() / "miniconda3" / "envs" / "gradi-loc" / "bin"

DLP_OUT_COLUMNS = ["uniprot_ac", "localization", "confidence"]
TMBED_OUT_COLUMNS = ["uniprot_ac", "cytoplasmic_fraction", "has_signal_peptide"]
TMBED_ACCESSORY_COLUMNS = ["uniprot_ac", "n_tm_helix", "n_tm_strand",
                           "cyto_longest_segment", "tmbed_length"]

# Proteins whose compartment is textbook. This is v1's 20-marker E. coli panel, which is what caught
# its single worst classification bug (40 peripheral inner-membrane proteins, DnaK among them, binned
# as inner_membrane instead of cytoplasm). Reported per gene, and enforced as a RATE, not per gene:
# DeepLocPro's own benchmark accuracy is 0.74, so demanding 20/20 from a predictor would be a
# spurious failure waiting to happen. The measured rate belongs in the run log.
GRAM_NEGATIVE_MARKERS: dict[str, str] = {
    "ompA": "outer_membrane", "ompC": "outer_membrane", "ompF": "outer_membrane",
    "lamB": "outer_membrane", "btuB": "outer_membrane", "lptD": "outer_membrane",
    "fhuA": "outer_membrane", "tolC": "outer_membrane",
    "malE": "periplasm", "dsbA": "periplasm", "oppA": "periplasm", "mglB": "periplasm",
    "secY": "cytoplasmic_membrane", "acrB": "cytoplasmic_membrane",
    "atpB": "cytoplasmic_membrane", "ftsW": "cytoplasmic_membrane",
    "groEL": "cytoplasm", "dnaK": "cytoplasm", "rpoB": "cytoplasm", "clpP": "cytoplasm",
}

MARKERS: dict[str, dict[str, str]] = {
    "ecoli": GRAM_NEGATIVE_MARKERS,
    # Kp takes the same panel. `ompC`/`ompF` are absent by biology, not by annotation gap: the
    # Klebsiella orthologs are OmpK36 and OmpK35 and are not named ompC/ompF here. `phoE` is a
    # named Kp outer-membrane porin, so it stands in. Absent genes are skipped and named.
    "kpneumoniae": {**GRAM_NEGATIVE_MARKERS, "phoE": "outer_membrane"},
    "saureus": {
        "clpP": "cytoplasm", "clpC": "cytoplasm", "clpX": "cytoplasm", "rpoB": "cytoplasm",
        "secY": "cytoplasmic_membrane", "atpB": "cytoplasmic_membrane",
        "spa": "cell_wall_surface",       # LPXTG-anchored; tests a class DeepLocPro rarely calls
        "hly": "extracellular",           # alpha-hemolysin (`hla` is the synonym, not the name)
    },
}
MARKER_MIN_AGREEMENT = 0.80

# TMbed's strand counts are the axis's best-validated signal: v1's came out textbook. The invariant
# enforced is the *call*, not the count -- a barrel has >= 8 strands, and TolC does not, because it
# is a trimer contributing 4 strands per monomer. The counts themselves are printed for comparison.
#
# GRAM-NEGATIVE ONLY, and not as a convenience: a Gram-positive organism has no outer membrane, so
# it has no beta-barrel to find. Worse, the same gene symbol can name a different protein -- Sa's
# `fhuA` is "Iron compound ABC transporter, ATP-binding protein" (cytoplasmic), not Ec's
# "Ferrichrome outer membrane transporter". CLAUDE.md's "join on locus_tag, not gene name" applies
# to hand-written marker panels too.
BARREL_MARKERS: dict[str, int] = {
    "ompA": 8, "ompC": 16, "ompF": 16, "lamB": 18, "btuB": 22, "fhuA": 22, "lptD": 26,
}
NOT_A_BARREL = {"tolC": 4}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------- the gradi-loc env

def ensure_gradi_loc() -> tuple[Path, Path, str]:
    """Locate the `gradi-loc` python and `tmbed`, or exit with the fix. Returns versions too."""
    bin_dir = Path(os.environ.get("GRADI_LOC_BIN", DEFAULT_GRADI_LOC_BIN))
    python, tmbed = bin_dir / "python", bin_dir / "tmbed"
    missing = [str(p) for p in (python, tmbed) if not p.exists()]
    if missing:
        sys.exit(
            "the `gradi-loc` env is not where expected.\n"
            f"  missing: {', '.join(missing)}\n"
            "  it must stay SEPARATE from `gradi`: DeepLocPro needs `fair-esm`, which claims the\n"
            "  same top-level `esm` package as the EvolutionaryScale `esm` stage 01 uses, so\n"
            "  installing it into `gradi` silently breaks embedding generation.\n"
            "  See install.sh for the recipe, or point GRADI_LOC_BIN at the right bin/ directory."
        )
    probe = subprocess.run(
        [str(python), "-c",
         "import DeepLocPro, tmbed, torch, transformers; "
         "print(DeepLocPro.__version__ if hasattr(DeepLocPro,'__version__') else '1.0.0', "
         "tmbed.__version__ if hasattr(tmbed,'__version__') else '?', "
         "torch.__version__, transformers.__version__)"],
        capture_output=True, text=True,
    )
    if probe.returncode != 0:
        sys.exit(f"`gradi-loc` python cannot import the predictors:\n{probe.stderr.strip()}")
    dlp_v, tmbed_v, torch_v, tf_v = probe.stdout.split()
    return python, tmbed, f"DeepLocPro {dlp_v} · tmbed {tmbed_v} · torch {torch_v} · transformers {tf_v}"


# ---------------- input

def species_frame(species: str, limit: int | None) -> pd.DataFrame:
    """The (uniprot_ac, sequence) rows this stage predicts on, in a stable order."""
    df = P.load(species)[["uniprot_ac", "sequence"]]
    df = df[df["sequence"].str.strip().ne("")].sort_values("uniprot_ac").reset_index(drop=True)
    return df.head(limit) if limit else df


def write_query_fasta(df: pd.DataFrame, path: Path) -> int:
    """Write the query FASTA with headers = bare accessions.

    Deliberately NOT data/source/uniprot/proteomes/*.fasta: those carry `sp|A5A616|MGTS_ECOLI` headers, so
    every prediction would need parsing back to an accession. Building it here makes the FASTA
    header *be* the join key, and guarantees the query set is exactly the rows in the deliverable.
    """
    with path.open("w", encoding="utf-8") as f:
        for ac, seq in zip(df["uniprot_ac"], df["sequence"]):
            f.write(f">{ac}\n{seq}\n")
    return len(df)


def cache_key(accessions: list[str], **params) -> str:
    """Identity of a cache: which accessions were asked, under which parameters.

    A prediction file cannot tell you what it covers, so the coverage is recorded beside it. A
    20-protein smoke-test cache silently serving a 5,728-protein run is the v1 short-payload
    failure mode, and it is what this guards.
    """
    tag = "|".join(f"{k}={v}" for k, v in sorted(params.items()))
    return hashlib.md5(("\n".join(sorted(accessions)) + "|" + tag).encode()).hexdigest()


# ---------------- DeepLocPro

def run_deeplocpro(species, df, python, device, refresh, limit) -> dict:
    """Predict, then assemble the deliverable and the raw-probability evidence table."""
    t0 = time.time()
    pre = "smoke_" if limit else ""
    cache = SCRATCH_DIR / "deeplocpro_cache" / f"{pre}{species}_L{LOC.DLP_MAX_LENGTH}"
    gram = LOC.GRAM_GROUP[species]
    say(f"  {species:<14} {len(df):>6} proteins   gram={gram}")

    if refresh and cache.exists():
        shutil.rmtree(cache)
    cache.mkdir(parents=True, exist_ok=True)
    # Counted before the worker runs, so `seconds` in the manifest is interpretable: after a fully
    # cached re-run it is assembly time only, not prediction time.
    n_cached_before = sum(1 for a in df["uniprot_ac"] if (cache / f"{a}.json").exists())

    with tempfile.TemporaryDirectory() as tmp:
        query = Path(tmp) / f"{species}.faa"
        write_query_fasta(df, query)
        proc = subprocess.run(
            [str(python), str(WORKER), "--fasta", str(query), "--cache", str(cache),
             "--device", device],
            text=True,
        )
    if proc.returncode != 0:
        sys.exit(f"FAILED: the gradi-loc DeepLocPro worker exited {proc.returncode} for {species}.")

    rows, prob_rows, missing = [], [], []
    for acc in df["uniprot_ac"]:
        f = cache / f"{acc}.json"
        if not f.exists():
            missing.append(acc)
            continue
        rec = json.loads(f.read_text())
        probs = rec["probs"]
        label, conf, margin = LOC.label_from_probs(probs, gram)
        rows.append({"uniprot_ac": acc, "localization": label, "confidence": round(conf, 4)})
        prob_rows.append({"uniprot_ac": acc, "margin": round(margin, 4),
                          "truncated": bool(rec["truncated"]),
                          **{f"p_{cls}": round(p, 4) for cls, p in zip(LOC.DLP_CANON, probs)}})
    if missing:
        sys.exit(f"FAILED: {len(missing)} {species} proteins have no DeepLocPro result "
                 f"(first few: {', '.join(missing[:5])}).")

    out = pd.DataFrame(rows)[DLP_OUT_COLUMNS]
    probs_df = pd.DataFrame(prob_rows)[
        ["uniprot_ac"] + [f"p_{c}" for c in LOC.DLP_CANON] + ["margin", "truncated"]]

    bad = sorted(set(out["localization"]) - set(LOC.LOC_CLASSES))
    if bad:
        sys.exit(f"FAILED: {species} produced labels outside LOC_CLASSES: {bad}")
    if gram == "positive":
        illegal = sorted(set(out["localization"]) - set(LOC.GRAM_POSITIVE_CLASSES))
        if illegal:
            sys.exit(f"FAILED: {species} is Gram-positive but was given {illegal} -- the remap "
                     "did not run.")
    sums = probs_df[[f"p_{c}" for c in LOC.DLP_CANON]].sum(axis=1)
    if not sums.between(0.999, 1.001).all():
        sys.exit(f"FAILED: {species} has probability vectors that do not sum to 1.0 "
                 f"(min {sums.min():.4f}, max {sums.max():.4f}).")

    out_path = (SCRATCH_DIR / f"smoke_deeplocpro_{species}.tsv") if limit \
        else (EVIDENCE_DIR / f"deeplocpro_{species}.tsv")
    if not limit:
        out = M.reindex(out, species)   # canonical row order -- see src/matrices.py
    out.to_csv(out_path, sep="\t", index=False)
    probs_df.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}deeplocpro_probabilities_{species}.tsv", sep="\t", index=False)

    counts = (pd.DataFrame({"localization": list(LOC.LOC_CLASSES)})
                .merge(out["localization"].value_counts().rename("n"),
                       left_on="localization", right_index=True, how="left")
                .fillna({"n": 0}))
    counts["n"] = counts["n"].astype(int)
    counts["pct"] = (100 * counts["n"] / max(len(out), 1)).round(2)
    counts.insert(0, "species", species)
    counts.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}deeplocpro_counts_{species}.tsv", sep="\t", index=False)

    n_trunc = int(probs_df["truncated"].sum())
    moved = float((probs_df["p_periplasm"] + probs_df["p_outer_membrane"]).mean()) \
        if gram == "positive" else 0.0
    return {
        "species": species, "gram_group": gram, "n": len(out), "n_expected": len(df),
        "truncated": n_trunc, "mean_confidence": round(float(out["confidence"].mean()), 4),
        "mean_margin": round(float(probs_df["margin"].mean()), 4),
        "remapped_mass": round(moved, 4), "n_computed": len(df) - n_cached_before,
        "model": "DeepLocPro 1.0 (ESM-2 650M, 20-checkpoint ensemble)",
        "max_length": LOC.DLP_MAX_LENGTH, "device": device,
        "seconds": round(time.time() - t0, 1),
        "path": str(out_path.relative_to(REPO_ROOT)),
        "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


# ---------------- TMbed

def shard_path(species: str, shard_size: int, idx: int, pre: str = "") -> Path:
    # shard_size is in the name on purpose: changing it must not reuse mismatched shards.
    return SCRATCH_DIR / "tmbed_shards" / f"{pre}{species}_{shard_size}_{idx:04d}.pred"


def run_tmbed(species, df, tmbed, shard_size, threads, refresh, limit) -> dict:
    """Shard, predict per shard with resume, then parse every label string into features."""
    t0 = time.time()
    pre = "smoke_" if limit else ""
    (SCRATCH_DIR / "tmbed_shards").mkdir(parents=True, exist_ok=True)
    accs = list(df["uniprot_ac"])
    seqs = dict(zip(df["uniprot_ac"], df["sequence"]))
    shards = [accs[i:i + shard_size] for i in range(0, len(accs), shard_size)]
    say(f"  {species:<14} {len(accs):>6} proteins in {len(shards)} shards of {shard_size}")

    key_path = SCRATCH_DIR / f".{pre}tmbed_{species}.key"
    key = cache_key(accs, shard_size=shard_size, out_format=TMBED_OUT_FORMAT)
    stale = key_path.exists() and key_path.read_text().strip() != key
    if refresh or stale:
        if stale:
            say(f"  {'':<14} [rerun] cached shards were built for a different query set")
        for i in range(len(shards) + 64):
            shard_path(species, shard_size, i, pre).unlink(missing_ok=True)

    pending = [i for i in range(len(shards)) if not shard_path(species, shard_size, i, pre).exists()]
    say(f"  {'':<14} {len(shards) - len(pending)} cached, {len(pending)} to run on CPU "
        f"(threads={threads})")

    for pos, i in enumerate(pending, 1):
        s0 = time.time()
        with tempfile.TemporaryDirectory() as tmp:
            fa = Path(tmp) / f"{species}_{i:04d}.faa"
            write_query_fasta(df[df["uniprot_ac"].isin(shards[i])], fa)
            dest = shard_path(species, shard_size, i, pre)
            out = subprocess.run(
                [str(tmbed), "predict", "-f", str(fa), "-p", str(dest.with_suffix(".tmp")),
                 "--out-format", TMBED_OUT_FORMAT, "--threads", str(threads)],
                capture_output=True, text=True,
            )
            if out.returncode != 0:
                sys.exit(f"FAILED: tmbed exited {out.returncode} on {species} shard {i}:\n"
                         f"{out.stderr.strip()[:2000]}")
            # Atomic: a truncated .pred must never be cached as a completed shard.
            dest.with_suffix(".tmp").rename(dest)
        dt = time.time() - s0
        rate = pos / max(time.time() - t0, 1e-9)
        say(f"    shard {pos}/{len(pending)}  {len(shards[i])} proteins in {dt:5.1f}s  "
            f"(~{(len(pending) - pos) / max(rate, 1e-9) / 60:.0f} min left)")
    key_path.write_text(key)

    feats: dict[str, dict] = {}
    labels: dict[str, str] = {}
    for i in range(len(shards)):
        f = shard_path(species, shard_size, i, pre)
        if f.exists():
            feats.update(LOC.parse_tmbed(f))
            labels.update(LOC.read_tmbed_labels(f))

    missing = [a for a in accs if a not in feats]
    if missing:
        sys.exit(f"FAILED: {len(missing)} {species} proteins have no TMbed result "
                 f"(first few: {', '.join(missing[:5])}).")
    # TMbed truncates nothing, so a length mismatch means the wrong sequence was scored.
    wrong = [a for a in accs if feats[a]["tmbed_length"] != len(seqs[a])]
    if wrong:
        sys.exit(f"FAILED: {len(wrong)} {species} label strings do not match their sequence length "
                 f"(first few: {', '.join(wrong[:5])}).")

    full = pd.DataFrame([{"uniprot_ac": a, **feats[a]} for a in accs])
    if not limit:
        full = M.reindex(full, species)   # canonical row order -- see src/matrices.py
    out_path = (SCRATCH_DIR / f"smoke_tmbed_{species}.tsv") if limit \
        else (EVIDENCE_DIR / f"tmbed_{species}.tsv")
    full[TMBED_OUT_COLUMNS].to_csv(out_path, sep="\t", index=False)
    full[TMBED_ACCESSORY_COLUMNS].to_csv(
        (SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}tmbed_topology_{species}.tsv", sep="\t", index=False)
    pd.DataFrame({"uniprot_ac": accs, "labels": [labels[a] for a in accs]}).to_csv(
        (SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}tmbed_labels_{species}.tsv", sep="\t", index=False)

    return {
        "species": species, "n": len(full), "n_expected": len(df),
        "mean_cytoplasmic_fraction": round(float(full["cytoplasmic_fraction"].mean()), 4),
        "with_signal_peptide": int(full["has_signal_peptide"].sum()),
        "with_tm_helix": int((full["n_tm_helix"] > 0).sum()),
        "beta_barrels": int((full["n_tm_strand"] >= LOC.BETA_BARREL_MIN_STRANDS).sum()),
        "shards_computed": len(pending), "shards_total": len(shards),
        "model": "TMbed 1.0 (ProtT5-XL-U50)", "out_format": TMBED_OUT_FORMAT,
        "shard_size": shard_size, "threads": threads,
        "seconds": round(time.time() - t0, 1),
        "path": str(out_path.relative_to(REPO_ROOT)),
        "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }


# ---------------- checks

def _read(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path, sep="\t") if path.exists() else None


def marker_checks(species_list, limit) -> tuple[list[str], float | None]:
    """The textbook-marker panel. Returns (failure lines, overall agreement)."""
    hits, total, fails = 0, 0, []
    for sp in species_list:
        path = (SCRATCH_DIR / f"smoke_deeplocpro_{sp}.tsv") if limit \
            else (EVIDENCE_DIR / f"deeplocpro_{sp}.tsv")
        dlp = _read(path)
        if dlp is None:
            continue
        m = dlp.merge(P.load(sp)[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left")
        for gene, expect in MARKERS.get(sp, {}).items():
            got = m.loc[m["gene_name"] == gene, "localization"].tolist()
            if not got:
                say(f"  {sp:<14} {gene:<6} not present under that gene name  (skipped)")
                continue
            ok = all(g == expect for g in got)
            total += 1
            hits += int(ok)
            say(f"  {sp:<14} {gene:<6} -> {','.join(got):<21} expected {expect:<21}"
                f"{'ok' if ok else '<- FAIL'}")
            if not ok:
                fails.append(f"{sp}/{gene}: got {got}, expected {expect}")
    if not total:
        return [], None
    say(f"\n  markers agreeing: {hits}/{total} = {100 * hits / total:.1f}%   "
        f"(floor {100 * MARKER_MIN_AGREEMENT:.0f}%)")
    if not fails:
        return [], hits / total
    say("  disagreements are listed above; the floor is what gates the run, not any single gene.")
    return fails, hits / total


def barrel_checks(species_list, limit) -> list[str]:
    """TMbed strand counts against the textbook values. Enforces the CALL, prints the count.

    Gram-negatives get the marker panel. A Gram-positive gets the opposite invariant -- **zero**
    beta-barrels, because it has no outer membrane to hold one -- which is a real check rather than
    a panel that cannot apply.
    """
    pre = "smoke_" if limit else ""
    failures = []
    for sp in species_list:
        topo = _read((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}tmbed_topology_{sp}.tsv")
        if topo is None:
            continue
        if LOC.GRAM_GROUP[sp] == "positive":
            n = int((topo["n_tm_strand"] >= LOC.BETA_BARREL_MIN_STRANDS).sum())
            say(f"  {sp:<14} {n:>4} beta-barrels   Gram-positive, so expected 0 "
                f"(no outer membrane){'':<19}{'ok' if n == 0 else '<- FAIL'}")
            if n:
                failures.append(f"{sp}: {n} proteins with >= {LOC.BETA_BARREL_MIN_STRANDS} strands, "
                                "but a Gram-positive has no outer membrane")
            continue
        m = topo.merge(P.load(sp)[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left")
        for gene, textbook in {**BARREL_MARKERS, **NOT_A_BARREL}.items():
            expect_barrel = gene in BARREL_MARKERS
            got = m.loc[m["gene_name"] == gene, "n_tm_strand"].tolist()
            if not got:
                continue
            is_barrel = [int(g) >= LOC.BETA_BARREL_MIN_STRANDS for g in got]
            ok = all(b == expect_barrel for b in is_barrel)
            note = "barrel" if expect_barrel else "NOT a barrel (trimer, 4 strands/monomer)"
            say(f"  {sp:<14} {gene:<6} {','.join(str(int(g)) for g in got):>4} strands  "
                f"textbook {textbook:>2}  {note:<44}{'ok' if ok else '<- FAIL'}")
            if not ok:
                failures.append(f"{sp}/{gene}: {got} strands, expected "
                                f"{'>=' if expect_barrel else '<'}{LOC.BETA_BARREL_MIN_STRANDS}")
    return failures


def cross_check(species_list, limit) -> None:
    """Two independent models against each other -- reported, never enforced.

    The separation is the evidence they agree: TM helices concentrate in `cytoplasmic_membrane`
    (5.2 against 0.0 in `cytoplasm`), strands only in `outer_membrane`, signal peptides in
    `periplasm` and `outer_membrane`. A collapse in that separation means one of them is broken.

    `cytoplasmic_fraction` is deliberately NOT expected to be monotone across the six compartments:
    a beta-barrel has periodic cytoplasmic turns between its strands, so `outer_membrane` sits above
    a fully-exported soluble periplasmic protein. It is monotone among the soluble compartments only.
    """
    pre = "smoke_" if limit else ""
    for sp in species_list:
        dlp = _read((SCRATCH_DIR / f"smoke_deeplocpro_{sp}.tsv") if limit
                    else (EVIDENCE_DIR / f"deeplocpro_{sp}.tsv"))
        tm = _read((SCRATCH_DIR / f"smoke_tmbed_{sp}.tsv") if limit
                   else (EVIDENCE_DIR / f"tmbed_{sp}.tsv"))
        topo = _read((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}tmbed_topology_{sp}.tsv")
        if dlp is None or tm is None or topo is None:
            continue
        m = dlp.merge(tm, on="uniprot_ac").merge(topo, on="uniprot_ac")
        say(f"  {sp}")
        say(f"    {'compartment':<22} {'n':>6} {'cyto_frac':>10} {'tm_helix':>9} "
            f"{'strands':>8} {'signal_pep':>11}")
        for cls in LOC.LOC_CLASSES:
            g = m[m["localization"] == cls]
            if g.empty:
                say(f"    {cls:<22} {0:>6}          -         -        -           -")
                continue
            say(f"    {cls:<22} {len(g):>6} {g['cytoplasmic_fraction'].mean():>10.3f} "
                f"{g['n_tm_helix'].mean():>9.1f} {g['n_tm_strand'].mean():>8.1f} "
                f"{100 * g['has_signal_peptide'].mean():>10.1f}%")
        say()


# ---------------- main

def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES),
                    choices=list(DEFAULT_SPECIES),
                    help="species to predict (human is not possible: DeepLocPro is prokaryote-only)")
    ap.add_argument("--only", choices=["deeplocpro", "tmbed"],
                    help="run one predictor; the other's outputs are left untouched")
    ap.add_argument("--shard-size", type=int, default=DEFAULT_SHARD_SIZE, help="TMbed shard size")
    ap.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"],
                    help="DeepLocPro device; TMbed is CPU-only here (it gates GPU on CUDA)")
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, help="only the first N proteins per species (smoke test)")
    ap.add_argument("--refresh", action="store_true", help="discard caches and re-predict")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    tracks = [args.only] if args.only else ["deeplocpro", "tmbed"]

    rule("=")
    say("STAGE 03 - localization from sequence")
    rule("=")
    say("  in       : data/processed/proteomes/proteome_<species>.tsv  (uniprot_ac + sequence)")
    say(f"  out      : {EVIDENCE_DIR.relative_to(REPO_ROOT)}/"
        "{deeplocpro,tmbed}_<species>.tsv  (+ scratch/)")
    say("             then: python scripts/localization/merge.py  ->  localization_<species>.tsv")
    say("  method   : DeepLocPro 1.0 -> compartment  ·  TMbed -> topology, sequence only")
    say(f"  classes  : {', '.join(LOC.LOC_CLASSES)}")
    say(f"  species  : {', '.join(args.species)}")
    say(f"  tracks   : {', '.join(tracks)}")
    say(f"  device   : {args.device} (DeepLocPro)   threads: {args.threads} (TMbed, CPU-only)")
    if args.limit:
        say(f"  limit    : {args.limit} proteins per species (SMOKE TEST)")
    if args.dry_run:
        say("\n  --dry-run: nothing predicted, nothing written.")
        for sp in args.species:
            for t in tracks:
                say(f"    {sp:<14} -> {t}_{sp}.tsv")
        return
    say()

    rule()
    say("PREDICTORS")
    rule()
    python, tmbed, versions = ensure_gradi_loc()
    say(f"  gradi-loc: {python.parent}")
    say(f"  versions : {versions}")
    say("  weights  : ESM-2 650M (~2.5 GB) + ProtT5-XL-U50 (~2.3 GB), cached from a previous run")
    say()

    frames = {sp: species_frame(sp, args.limit) for sp in args.species}
    dlp_rows, tmbed_rows = [], []

    # Sequentially, never concurrently: v1 measured TMbed at ~17 min/shard while DeepLocPro
    # competed for the CPU, against ~4 min/shard with the machine to itself.
    if "deeplocpro" in tracks:
        rule()
        say("PREDICT - DeepLocPro (compartment)")
        rule()
        dlp_rows = [run_deeplocpro(sp, frames[sp], python, args.device, args.refresh, args.limit)
                    for sp in args.species]
        pd.DataFrame(dlp_rows).to_csv(
            EVIDENCE_DIR / f"{'smoke_' if args.limit else ''}deeplocpro_manifest.tsv",
            sep="\t", index=False)
        say()

    if "tmbed" in tracks:
        rule()
        say("PREDICT - TMbed (topology)")
        rule()
        tmbed_rows = [run_tmbed(sp, frames[sp], tmbed, args.shard_size, args.threads,
                                args.refresh, args.limit) for sp in args.species]
        pd.DataFrame(tmbed_rows).to_csv(
            EVIDENCE_DIR / f"{'smoke_' if args.limit else ''}tmbed_manifest.tsv", sep="\t", index=False)
        say()

    if dlp_rows:
        rule()
        say("COVERAGE - compartments")
        rule()
        say(f"  {'species':<14} {'n':>6} " + " ".join(f"{LOC.LOC_CLASS_ABBREV[c]:>7}"
                                                      for c in LOC.LOC_CLASSES))
        for sp in args.species:
            counts = pd.read_csv(
                EVIDENCE_DIR / f"{'smoke_' if args.limit else ''}deeplocpro_counts_{sp}.tsv", sep="\t")
            by = dict(zip(counts["localization"], counts["n"]))
            say(f"  {sp:<14} {int(counts['n'].sum()):>6} "
                + " ".join(f"{by.get(c, 0):>7,}" for c in LOC.LOC_CLASSES))
        say("\n  Coverage is 100% by construction: DeepLocPro always returns a call, so there is")
        say("  no `unknown` class. Zeros for periplasm/outer_membrane on a Gram-positive are")
        say("  structural, not missing data -- their probability was folded into extracellular.")
        say()

    rule()
    say("SPOT CHECKS")
    rule()
    marker_fails, marker_rate = ([], None)
    barrel_fails: list[str] = []
    if dlp_rows:
        marker_fails, marker_rate = marker_checks(args.species, args.limit)
        say()
    if tmbed_rows:
        barrel_fails = barrel_checks(args.species, args.limit)
        say()
    for r in dlp_rows:
        if r["truncated"]:
            say(f"  {r['species']:<14} {r['truncated']} sequences truncated at "
                f"{r['max_length']} residues before ESM-2 (localization signal is N-terminal)")
    for r in dlp_rows:
        if r["gram_group"] == "positive":
            say(f"  {r['species']:<14} Gram-positive remap moved mean "
                f"{r['remapped_mass']:.3f} probability mass (periplasm + outer membrane) into")
            say(f"  {'':<14} `extracellular`. Do NOT compare its extracellular share against a "
                "Gram-negative.")
    say()

    if dlp_rows and tmbed_rows:
        rule()
        say("CROSS-CHECK - DeepLocPro compartment vs TMbed topology (reported, not enforced)")
        rule()
        cross_check(args.species, args.limit)

    rule()
    say("OUTPUTS")
    rule()
    for r in dlp_rows + tmbed_rows:
        say(f"  {r['path']:<56} {(REPO_ROOT / r['path']).stat().st_size / 1e3:>8.1f} kB")
    say("  evidence/ + scratch/    deeplocpro_probabilities_<species>.tsv, deeplocpro_counts_<species>.tsv,")
    say("                tmbed_topology_<species>.tsv, tmbed_labels_<species>.tsv,")
    say("                {deeplocpro,tmbed}_manifest.tsv, and the two resumable caches")
    total_min = sum(r["seconds"] for r in dlp_rows + tmbed_rows) / 60
    say(f"\n  total {total_min:.1f} min")
    say()

    bad = [r for r in dlp_rows + tmbed_rows if r["n"] != r["n_expected"]]
    if bad:
        sys.exit("FAILED: " + ", ".join(
            f"{r['species']} has {r['n']} rows, expected {r['n_expected']}" for r in bad))
    if barrel_fails:
        sys.exit("FAILED beta-barrel checks:\n  " + "\n  ".join(barrel_fails)
                 + "\n  TMbed's strand counts are the axis's best-validated signal; a drift here"
                   "\n  means it is misconfigured. Do not trust this run.")
    if marker_rate is not None and marker_rate < MARKER_MIN_AGREEMENT:
        sys.exit(f"FAILED marker panel: {100 * marker_rate:.1f}% agreement is below the "
                 f"{100 * MARKER_MIN_AGREEMENT:.0f}% floor.\n  "
                 + "\n  ".join(marker_fails))
    rule("=")
    say("stage 03 complete.")
    rule("=")


if __name__ == "__main__":
    main()
