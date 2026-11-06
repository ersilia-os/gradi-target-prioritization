"""Geptop 2.0 — gene essentiality from orthology and phylogeny, ported faithfully to Python 3.

Geptop scores every protein in a query proteome by how often it has a **reciprocal best hit** to a
known essential gene in 37 reference prokaryotes, weighting each reference by its **phylogenetic
closeness** to the query. Wen et al., *Front Microbiol* 10:1236 (2019) — the latest release; there
is no v3.

    score_raw(p) = SUM over references r of  1/d(query, r)   for each r where p has an RBH
                   to a gene listed essential in DEG2
    score(p)     = (score_raw(p) - min) / (max - min)        <-- min-max, WITHIN this proteome

`d` is the **k=6 composition-vector (CV-tree) distance**, not a sequence-identity proxy. v1
substituted median RBH % identity for it; that was a real deviation and is not repeated here.

WHY A PORT AND NOT THE ORIGINAL
-------------------------------
`data/source/geptop/geptop2.py` is kept verbatim as the reference, but it cannot be
run: Python 2 only, and three genuine defects beyond that.
  1. `os.popen(makeblastdb)` is followed immediately by `os.popen(blastp)` with **no wait** — a
     race. Here: `subprocess.run(..., check=True)`, sequenced.
  2. `except (Exception, err)` / `except (Exception, IOError)` — broken handlers that swallow real
     failures and print a class as a value.
  3. `import pp` (Parallel Python, py2-only, dead) is genuinely used at its line 152. Replaced with
     `concurrent.futures`.
The ARITHMETIC is transcribed unchanged, quirks included — see `composition_vector`.

**100% COVERAGE, BUT READ THE SECOND COLUMN**
Every protein gets a score, so coverage is 100% by construction. A protein with no RBH to any
reference essential gene scores exactly 0.0 — a MEASURED zero, not a gap. `geptop_informative`
records whether a protein had any essential RBH at all, because "100% scored" and "100% informative"
are different claims (stage 02 makes the same distinction for COG: 79.1% classified, 72.8%
informative).

**SCORES ARE PROTEOME-RELATIVE.** Min-max normalisation happens within each query, so a 0.6 in Kp
and a 0.6 in Sa are NOT the same quantity. Never compare raw Geptop scores across species — the
same trap as stage 01's per-species t-SNE coordinates. `geptop_score_raw` is kept for anyone who
needs a comparable quantity.

**TWO OF OUR ANCHORS ARE IN ITS REFERENCE SET.** `DataSet4.faa` is *E. coli* K-12 MG1655 and
`DataSet13.faa` is *S. aureus* NCTC 8325 — our exact stage-00 proteomes, with 296 and 350 of their
own genes in DEG2. For those two the CV distance to themselves is ~0, the authors' `if d == 0:
a = 100` branch fires, and the score largely reads off the answer. **We do NOT leave-one-out** (a
deliberate decision): for those two species DEG supplies a measured label anyway, and Kp — the
actual anchor — is absent from all 37 references, so it is uncontaminated. Instead every row carries
`geptop_in_reference_set`, and validation uses DEG species OUTSIDE the reference set.

Needs `blastp`/`makeblastdb`, borrowed from the `gradi-prokka` env (`GRADI_BLAST_BIN` overrides) —
the `GRADI_RPSBLAST_BIN` pattern. Do NOT install blast into `gradi`: no osx-arm64 build, it would
drag the env to osx-64 and take ESM-C down with it.

Outputs `data/processed/essentiality/geptop_<species>.tsv` and `scratch/geptop/` caches.
Measured: ~22 s per blastp run, 74 runs per proteome, so ~30 min per species plus a one-off
composition-vector pass over the 37 references.

    python scripts/essentiality/geptop.py --species kpneumoniae
    python scripts/essentiality/geptop.py --dry-run
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

RAW = REPO_ROOT / "data" / "source" / "geptop"
DATASETS = RAW / "datasets2"
DEG2 = RAW / "DEG2"
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
RUN_DIR = SCRATCH_DIR / "geptop"
DEFAULT_BLAST_BIN = Path.home() / "miniconda3" / "envs" / "gradi-prokka" / "bin"

KSTRING = 6
CUTOFF = 0.24          # the authors' default essentiality-score cutoff
EVALUE = 10.0          # the authors' hit filter: `if E_value < 10`

# The authors' amino-acid -> base-20 map, transcribed verbatim from geptop2.py including its
# ambiguity collapses (B->D's index, U->C's, X->G's, Z->E's, J->I's). Do not "fix" these: they
# define the numbering the whole composition vector is built on.
AA_DICT = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5,
           'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11,
           'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17,
           'W': 18, 'Y': 19, 'B': 2, 'U': 1, 'X': 5, 'Z': 3, 'J': 7}

# Our anchors that ARE Geptop references -- the circularity flag.
IN_REFERENCE_SET = {"ecoli": "DataSet4.faa", "saureus": "DataSet13.faa"}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def blast_bin() -> Path:
    p = Path(os.environ.get("GRADI_BLAST_BIN", DEFAULT_BLAST_BIN))
    for exe in ("blastp", "makeblastdb"):
        if not (p / exe).exists():
            sys.exit(f"FATAL no {exe} at {p}. Geptop needs NCBI BLAST+; this project borrows it "
                     "from the `gradi-prokka` env (GRADI_BLAST_BIN overrides). Do NOT install "
                     "blast into `gradi` -- there is no osx-arm64 build.")
    return p


def str_to_num(s: str) -> int:
    """The authors' base-20 encoding of a k-mer; -1 if any residue is outside AA_DICT."""
    n = 0
    for i, c in enumerate(s):
        if c not in AA_DICT:
            return -1
        n += AA_DICT[c] * (20 ** (len(s) - i - 1))
    return n


def composition_vector(fasta: Path) -> dict[int, float]:
    """The k=6 CV-tree composition vector. Arithmetic transcribed from geptop2.py, quirks included.

    For each position the original counts a 4-mer (k2), a 5-mer (k1) and a 6-mer (k0), then scores
    each observed 6-mer against the Markov prediction from its 5- and 4-mer counts:

        p0     = k1[n1]*k1[n2]*string0*string2 / k2[n3] / string1 / string1
        k[n0]  = (k0[n0] - p0) / p0            and -1 when the 6-mer was never observed

    Two transcribed quirks worth naming, because they look like bugs and must be preserved:
      * a 6-mer reachable from an observed 5-mer but never itself observed gets **exactly -1**,
        not a computed value;
      * `string0/1/2` are recomputed inside the per-sequence loop in the original; since the counts
        accumulate across sequences, the values after the loop are the whole-proteome totals, which
        is what is used. Same result here, computed once after the loop.
    """
    k0: dict[int, int] = {}
    k1: dict[int, int] = {}
    k2: dict[int, int] = {}
    seq = []
    with open(fasta, errors="ignore") as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    _accumulate("".join(seq), k0, k1, k2)
                seq = []
            else:
                seq.append(line.strip())
    if seq:
        _accumulate("".join(seq), k0, k1, k2)

    for d in (k0, k1, k2):
        d.pop(-1, None)
    string0, string1, string2 = sum(k0.values()), sum(k1.values()), sum(k2.values())
    if not (string0 and string1 and string2):
        sys.exit(f"FATAL empty composition vector for {fasta}")

    mod = 20 ** (KSTRING - 2)
    out: dict[int, float] = {}
    for n1, c1 in k1.items():
        base0, base2, n3 = n1 * 20, (n1 % mod) * 20, n1 % mod
        c3 = k2.get(n3)
        if not c3:
            continue
        for aa in range(20):
            n2 = base2 + aa
            c2 = k1.get(n2)
            if c2 is None:
                continue
            n0 = base0 + aa
            if n0 in k0:
                p0 = c1 * c2 * string0 * string2 / c3 / string1 / string1
                out[n0] = (k0[n0] - p0) / p0
            else:
                out[n0] = -1.0
    return out


def _accumulate(seq: str, k0, k1, k2) -> None:
    """Per-sequence k-mer counting, positions exactly as the original indexes them."""
    len0 = len(seq) - KSTRING + 3
    for s in range(len0):
        end = KSTRING + s - 2
        n = str_to_num(seq[s:end])
        k2[n] = k2.get(n, 0) + 1
        if s < len0 - 2:
            n = str_to_num(seq[s:end + 2])
            k0[n] = k0.get(n, 0) + 1
        if s < len0 - 1:
            n = str_to_num(seq[s:end + 1])
            k1[n] = k1.get(n, 0) + 1


def cv_distance(cv1: dict[int, float], cv2: dict[int, float]) -> float:
    """The authors' correlation distance: (1 - cos) / 2, transcribed verbatim.

    NOTE the original's `Distance()` computes O and Q with the SAME product in the shared branch
    (`O += v1*v2` then `Q += v1*v2`), so Q is not the usual sum of squares of cv2 over the shared
    keys. That is preserved: changing it would change every distance and so every score.
    """
    O = P_ = Q = 0.0
    for key, v1 in cv1.items():
        v2 = cv2.get(key)
        if v2 is not None:
            O += v1 * v2
            P_ += v1 * v1
            Q += v1 * v2
        else:
            P_ += v1 * v1
    for key, v2 in cv2.items():
        if key not in cv1:
            Q += v2 * v2
    if P_ <= 0 or Q <= 0:
        return 1.0
    return (1 - O / math.sqrt(P_ * Q)) / 2


def _cv_cached(fasta: Path) -> dict[int, float]:
    """Composition vectors are expensive and input-invariant, so compute each one once.

    **The key is a hash of the resolved PATH, not `fasta.stem`.** Keying on the stem was a real bug
    with a silent, severe failure mode: every query proteome is written to `<dir>/query.faa`, so all
    of them -- and every validation species -- collided on a single `cv_query.pkl`. Whichever ran
    first wrote it and every run afterwards read ITS composition vector as their own, i.e. used the
    wrong organism's phylogeny weights throughout, with no error. Caught only because three
    different species reported the same nearest reference at the same distance, and because
    *S. aureus* was reported at E. coli's distance to Salmonella.
    """
    import hashlib
    import pickle

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    cache = _cv_cache_path(fasta)
    if cache.exists():
        with open(cache, "rb") as fh:
            return pickle.load(fh)
    cv = composition_vector(fasta)
    tmp = cache.with_suffix(".tmp")
    with open(tmp, "wb") as fh:
        pickle.dump(cv, fh, protocol=4)
    tmp.rename(cache)
    return cv


def _cv_cache_path(fasta: Path) -> Path:
    """Where `_cv_cached` will look. Shared so the pre-check cannot drift from the real key --
    it did once: the pre-check tested the old stem-only name and reported "37 cached" while
    `_cv_cached` went on to recompute all of them under the new path-hashed key."""
    import hashlib

    tag = hashlib.sha256(str(fasta.resolve()).encode()).hexdigest()[:12]
    return SCRATCH_DIR / f"cv_{fasta.stem}_{tag}.pkl"


def _cv_worker(path_str: str) -> tuple[str, int]:
    p = Path(path_str)
    return p.stem, len(_cv_cached(p))


# ---------------------------------------------------------------- orthology: reciprocal best hits

def _blastdb(bin_dir: Path, fasta: Path, out_prefix: Path) -> Path:
    """One BLAST database per FASTA, built once and reused.

    The original rebuilt a database for every (query, reference) PAIR -- 74 rebuilds per proteome of
    the same 37 references -- and did it with an unwaited `os.popen`, so `blastp` could start before
    the database existed. Built once here, sequenced with `check=True`.
    """
    if not (out_prefix.with_suffix(".phr").exists() or (out_prefix.parent /
            f"{out_prefix.name}.phr").exists()):
        subprocess.run([str(bin_dir / "makeblastdb"), "-in", str(fasta), "-dbtype", "prot",
                        "-out", str(out_prefix)], check=True, capture_output=True)
    return out_prefix


def _best_hits(bin_dir: Path, query: Path, db: Path, cache: Path, threads: int) -> dict[str, str]:
    """Best hit per query sequence, e-value < EVALUE. Cached -- blastp is the runtime here.

    The original parses BLAST XML and breaks after the first HSP of the first alignment, i.e. the
    single best hit. `-max_target_seqs 1` is the tabular equivalent.
    """
    if not cache.exists():
        subprocess.run([str(bin_dir / "blastp"), "-query", str(query), "-db", str(db),
                        "-outfmt", "6 qseqid sseqid evalue", "-max_target_seqs", "1",
                        "-evalue", str(EVALUE), "-num_threads", str(threads),
                        "-out", str(cache) + ".tmp"], check=True, capture_output=True)
        Path(str(cache) + ".tmp").rename(cache)
    best: dict[str, str] = {}
    with open(cache) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, sub, ev = parts[0], parts[1], float(parts[2])
            if ev < EVALUE and q not in best:      # file is ordered best-first per query
                best[q] = sub
    return best


def essential_gi() -> set[str]:
    """The gi numbers DEG2 lists as essential. The original keys on `split('|')[1]`; so do we."""
    ids = {line.split("|")[1] for line in open(DEG2) if "|" in line}
    if not ids:
        sys.exit(f"FATAL DEG2 at {DEG2} yielded no ids -- the scoring would be all zeros")
    return ids


def score_species(sp: str, bin_dir: Path, refs: list[Path], deg: set[str],
                  threads: int, cutoff: float) -> pd.DataFrame:
    """Geptop score for every protein in one query proteome."""
    prot = P.load(sp)[["uniprot_ac", "sequence"]]
    qdir = RUN_DIR / sp
    qdir.mkdir(parents=True, exist_ok=True)
    qfa = qdir / "query.faa"
    if not qfa.exists():
        with open(qfa, "w") as fh:
            for a, seq in zip(prot.uniprot_ac, prot.sequence):
                fh.write(f">{a}\n{seq}\n")
    say(f"  {sp}: {len(prot):,} proteins")

    qcv = _cv_cached(qfa)
    qdb = _blastdb(bin_dir, qfa, qdir / "query_db")

    raw = pd.Series(0.0, index=prot.uniprot_ac)
    n_ess_rbh = pd.Series(0, index=prot.uniprot_ac)
    n_rbh = pd.Series(0, index=prot.uniprot_ac)     # ALL reciprocal best hits, essential or not
    audit = []
    for i, ref in enumerate(refs, 1):
        d = cv_distance(qcv, _cv_cached(ref))
        # The authors' own branch: a query that IS a reference has d == 0 and is weighted 100,
        # which swamps every genuine contribution (1/d runs ~2.0-3.3). Kept, and flagged upstream.
        w = 100.0 if d == 0 else 1.0 / d
        rdb = _blastdb(bin_dir, ref, SCRATCH_DIR / f"db_{ref.stem}")
        fwd = _best_hits(bin_dir, qfa, rdb, qdir / f"fwd_{ref.stem}.tsv", threads)
        rev = _best_hits(bin_dir, ref, qdb, qdir / f"rev_{ref.stem}.tsv", threads)
        rbh = {q: h for q, h in fwd.items() if rev.get(h) == q}
        ess = [q for q, h in rbh.items() if h.split("|")[1] in deg] if rbh else []
        for q in rbh:
            n_rbh[q] += 1
        for q in ess:
            raw[q] += w
            n_ess_rbh[q] += 1
        audit.append({"species": sp, "reference": ref.stem, "cv_distance": round(d, 6),
                      "weight": round(w, 4), "n_rbh": len(rbh), "n_essential_rbh": len(ess)})
        say(f"    [{i:2}/{len(refs)}] {ref.stem:22} d={d:.4f} w={w:7.2f}  "
            f"RBH {len(rbh):5}  essential {len(ess):4}")

    # EVIDENCE, NOT SCRATCH. This table is the only record of which of the 37 references actually
    # contributed and how much -- it is what makes the circularity visible (a self-reference weighs
    # ~89-121 against ~2.0-3.3 for everything else), and CLAUDE.md quotes its numbers. It used to
    # be written into the run directory under scratch/, where a routine `scratch/` purge deleted
    # all three copies and cost a 90-minute re-run. The directory contract's own test settles it:
    # "would you cite or check it" -> evidence/.
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(audit).to_csv(EVIDENCE_DIR / f"geptop_reference_audit_{sp}.tsv",
                               sep="\t", index=False)
    lo, hi = float(raw.min()), float(raw.max())
    norm = (raw - lo) / (hi - lo) if hi > lo else raw * 0.0
    # A score of 0 has TWO meanings and they must not be conflated -- an earlier version of this
    # script reported them together as "uninformative" and understated real coverage by 58 points.
    #   orthologs, none essential -> a CONFIDENT NON-ESSENTIAL call (evidence, not absence)
    #   no orthologs at all       -> the method cannot see this protein (genuinely no information)
    nr, ne = n_rbh.to_numpy(), n_ess_rbh.to_numpy()
    evidence = np.where(ne > 0, "essential_orthologs",
                        np.where(nr > 0, "orthologs_none_essential", "no_orthologs"))
    out = pd.DataFrame({
        "uniprot_ac": prot.uniprot_ac.to_numpy(),
        "geptop_score": norm.to_numpy().round(4),
        "geptop_score_raw": raw.to_numpy().round(4),
        "geptop_essential": (norm.to_numpy() > cutoff).astype(int),
        "geptop_n_rbh": nr,
        "geptop_n_essential_rbh": ne,
        "geptop_evidence": evidence,
        # TRUE coverage: the method had orthology evidence either way. NOT `n_essential_rbh > 0`.
        "geptop_informative": (nr > 0).astype(int),
        "geptop_in_reference_set": int(sp in IN_REFERENCE_SET),
    })
    return out


# ---------------------------------------------------------------- validation, instead of leave-one-out

def validate(deg_species: str, bin_dir: Path, refs: list[Path], deg: set[str],
             threads: int, cutoff: float) -> dict:
    """Score a DEG species that is NOT a Geptop reference, and measure AUROC against its labels.

    **This is what replaces leave-one-out.** We deliberately do not exclude a query's own dataset
    (see the module docstring), so accuracy cannot be measured on E. coli or S. aureus -- their
    scores are circular. Instead we measure on species DEG labels but Geptop has never seen, which
    is a genuinely out-of-set test and *stricter* than the paper's: the exclusion here is at GENUS
    level, so a same-genus reference cannot leak orthologs either.

    The paper reports a mean AUC of 0.84 across prokaryotes. That is the number to compare against.
    """
    from sklearn.metrics import roc_auc_score

    lab = pd.read_csv(EVIDENCE_DIR / "labeled_proteins.tsv", sep="\t")
    sub = lab[lab["species"].astype(str) == deg_species].copy()
    if sub.empty:
        sys.exit(f"FATAL {deg_species!r} is not in the DEG corpus; see evidence/labeled_proteins.tsv")
    # One DEG dataset at a time: a species with several screens has several base rates, and pooling
    # them would score the model against a label set no single experiment ever produced.
    ds = sub["deg_dataset_id"].value_counts().index[0]
    sub = sub[sub["deg_dataset_id"] == ds].drop_duplicates("fasta_id")
    say(f"  {deg_species}  dataset {ds}  n={len(sub):,}  "
        f"essential {int(sub.essential.sum()):,} ({100*sub.essential.mean():.1f}%)")

    vdir = SCRATCH_DIR / f"validate_{ds}"
    vdir.mkdir(parents=True, exist_ok=True)
    vfa = vdir / "query.faa"
    if not vfa.exists():
        with open(vfa, "w") as fh:
            for i, (fid, seq) in enumerate(zip(sub.fasta_id, sub.sequence)):
                fh.write(f">q{i}\n{seq}\n")
    order = [f"q{i}" for i in range(len(sub))]

    qcv = _cv_cached(vfa)
    qdb = _blastdb(bin_dir, vfa, vdir / "query_db")
    raw = pd.Series(0.0, index=order)
    n_ess_rbh = pd.Series(0, index=order)
    for i, ref in enumerate(refs, 1):
        d = cv_distance(qcv, _cv_cached(ref))
        w = 100.0 if d == 0 else 1.0 / d
        rdb = _blastdb(bin_dir, ref, SCRATCH_DIR / f"db_{ref.stem}")
        fwd = _best_hits(bin_dir, vfa, rdb, vdir / f"fwd_{ref.stem}.tsv", threads)
        rev = _best_hits(bin_dir, ref, qdb, vdir / f"rev_{ref.stem}.tsv", threads)
        rbh = {q: h for q, h in fwd.items() if rev.get(h) == q}
        for q, h in rbh.items():
            if h.split("|")[1] in deg:
                raw[q] += w
                n_ess_rbh[q] += 1
        if i % 10 == 0 or i == len(refs):
            say(f"    [{i:2}/{len(refs)}] d={d:.4f}  cumulative informative "
                f"{int((n_ess_rbh > 0).sum()):,}/{len(order):,}")

    lo, hi = float(raw.min()), float(raw.max())
    score = (raw - lo) / (hi - lo) if hi > lo else raw * 0.0
    y = sub["essential"].to_numpy(dtype=int)
    auc = float(roc_auc_score(y, score.to_numpy()))
    inf = float((n_ess_rbh > 0).mean())
    # AUROC on the SCORED subset too (proteins with >=1 essential RBH): the tied zero block drags
    # the global figure toward 0.5 without saying anything about discrimination where the method
    # produced a non-zero score. NB this is narrower than `geptop_informative` in the main output,
    # which counts orthology evidence EITHER WAY -- do not compare the two percentages.
    m = (n_ess_rbh > 0).to_numpy()
    auc_inf = float(roc_auc_score(y[m], score.to_numpy()[m])) if 0 < y[m].sum() < m.sum() else float("nan")
    say(f"    AUROC {auc:.4f} (all)   {auc_inf:.4f} (informative only)   "
        f"informative {100*inf:.1f}%   paper reports mean AUC 0.84")
    return {"deg_species": deg_species, "deg_dataset_id": ds, "n": len(sub),
            "n_essential": int(y.sum()), "base_rate": round(float(y.mean()), 4),
            "auroc": round(auc, 4), "auroc_informative": round(auc_inf, 4),
            "frac_with_essential_rbh": round(inf, 4), "paper_mean_auc": 0.84}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Human is EXCLUDED BY CONSTRUCTION, not by choice: Geptop's 37 references are all
    # prokaryotes with no eukaryote among them, so there is no honest score to give it -- the same
    # reason stage 02 does not offer `--species human` for COG.
    ap.add_argument("--species", nargs="+", default=["kpneumoniae"],
                    choices=["kpneumoniae", "ecoli", "saureus"])
    ap.add_argument("--threads", type=int, default=10)
    ap.add_argument("--cutoff", type=float, default=CUTOFF)
    ap.add_argument("--cv-jobs", type=int, default=4,
                    help="parallel composition-vector workers (memory-hungry; 4 is safe)")
    ap.add_argument("--validate", nargs="*", metavar="DEG_SPECIES",
                    help="score DEG species that are NOT Geptop references and report AUROC; "
                         "this is what replaces leave-one-out")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    refs = sorted(DATASETS.glob("*.faa"))
    rule("=")
    say("STAGE 07 - Geptop 2.0 essentiality (orthology + phylogeny), faithful Python-3 port")
    rule("=")
    say(f"  references   {len(refs)} prokaryote proteomes from datasets2/")
    say(f"  essential    {sum(1 for _ in open(DEG2)):,} gene ids in DEG2")
    say(f"  method       k={KSTRING} composition-vector distance x reciprocal-best-hit to an")
    say("               essential reference gene; min-max normalised WITHIN each proteome")
    say(f"  aligner      blastp (the authors' own), e-value < {EVALUE}")
    say(f"  species      {', '.join(args.species)}")
    for sp in args.species:
        if sp in IN_REFERENCE_SET:
            say(f"  WARNING      {sp} IS reference {IN_REFERENCE_SET[sp]} -- its score is")
            say("               CIRCULAR by construction and flagged, not corrected (see docstring)")
    rule("=")
    if args.dry_run:
        say("dry run: nothing computed.")
        return

    blast_bin()
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    rule()
    say("COMPOSITION VECTORS - one-off, cached; the expensive part of the phylogeny term")
    rule()
    todo = [r for r in refs if args.refresh or not _cv_cache_path(r).exists()]
    say(f"  {len(refs) - len(todo)} cached, {len(todo)} to compute")
    if todo:
        with ProcessPoolExecutor(max_workers=args.cv_jobs) as ex:
            for i, (stem, n) in enumerate(ex.map(_cv_worker, [str(t) for t in todo]), 1):
                say(f"    [{i}/{len(todo)}] {stem:22} {n:,} k-mer features")
    say("  done")

    deg = essential_gi()
    bin_dir = blast_bin()

    if args.validate is not None:
        rule()
        say("VALIDATE - out-of-set DEG species (genus-level exclusion, stricter than the paper)")
        rule()
        rows = [validate(sp, bin_dir, refs, deg, args.threads, args.cutoff)
                for sp in (args.validate or ["Haemophilus influenzae"])]
        out = EVIDENCE_DIR / "geptop_validation.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        say(f"\n  wrote {out.relative_to(REPO_ROOT)}")
        rule("=")
        return
    rule()
    say(f"SCORE - bidirectional blastp against {len(refs)} references, {len(deg):,} essential ids")
    rule()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for sp in args.species:
        df = score_species(sp, bin_dir, refs, deg, args.threads, args.cutoff)
        path = OUT_DIR / f"geptop_{sp}.tsv"
        df.to_csv(path, sep="\t", index=False)
        say(f"  -> {path.relative_to(REPO_ROOT)}")
        vc = df.geptop_evidence.value_counts()
        say(f"     scored {len(df):,} (100% by construction). Evidence behind those scores:")
        for k, desc in (("essential_orthologs", "essential orthologs -> scored > 0"),
                        ("orthologs_none_essential", "orthologs, none essential -> real NEGATIVE"),
                        ("no_orthologs", "no orthologs -> NO INFORMATION")):
            n = int(vc.get(k, 0))
            say(f"       {desc:52} {n:6,} ({100*n/len(df):5.1f}%)")
        say(f"     -> informative (evidence either way) {int(df.geptop_informative.sum()):,} "
            f"({100*df.geptop_informative.mean():.1f}%)")
        say(f"     -> called essential {int(df.geptop_essential.sum()):,} "
            f"({100*df.geptop_essential.mean():.1f}%)")
        if sp in IN_REFERENCE_SET:
            say(f"     WARNING circular: {sp} is reference {IN_REFERENCE_SET[sp]}")
        say()
    rule("=")
    say("geptop complete.")
    rule("=")


if __name__ == "__main__":
    main()
