"""OrthoDB orthogroups for every protein in all four proteomes (stage 05, absolute orthology).

OrthoFinder (orthology/orthofinder.py) infers orthogroups DE NOVO from whatever proteomes are in the run, so
its recall moves when the panel changes. This script gives the complement: groups defined ELSEWHERE,
by OrthoDB, over 17,551 bacterial and 5,952 eukaryotic species, so the answer does not depend on our
query set at all.

**This is OrthoDB v12.2 and it REPLACES the previous v11 assignment.** The v12 README is explicit
that an "OG unique id (not stable and re-used between releases)", so v11 ids cannot be relabelled --
a fresh assignment is the only correct route. The release is recorded in `orthodb_version` on every
row so a mixed-version join is impossible to make by accident.

Why not just read UniProt's `xref_orthodb`
-----------------------------------------
Because it is empty exactly where we need it. Measured: E. coli 96.4%, S. aureus 99.4%, and
**K. pneumoniae HS11286 0.0%** -- and no Klebsiella proteome in UniProt carries the xref at all.
HS11286 is also absent from OrthoDB outright (taxid 1125630 is not in species.tab). Sequence
assignment is not an optimisation here; it is the only thing that works for the anchor organism.

Why the previous version capped at 74.6% on Kp
----------------------------------------------
It searched each species against its OWN OrthoDB assembly only -- for Kp that is `72407_0`
(K. pneumoniae subsp. pneumoniae), a single 4,975-protein organism. Orthogroups are cross-species by
construction, so this run searches ALL of OrthoDB instead, which is what the 1,393 Kp proteins at
`no_sequence_match` were missing.

Representatives, not the whole database
---------------------------------------
Assigning an OG needs a hit to ANY member, so indexing all 147M proteins is waste. One pass over
`aa_fasta.gz` keeps at most `--reps` sequences per domain-level OG (default 5): ~6.9M sequences
instead of 147M, a ~20x reduction. **This is a trade, so it is measured, not assumed** --
`--calibrate` sweeps the cap on E. coli, whose 96.4% UniProt-xref coverage is an independent ground
truth.

Two domains, one pass
---------------------
OrthoDB has NO root level spanning domains. A bacterial protein's domain group is `<n>at2`, human's
is `<n>at2759`, and **the two are not comparable** -- this table says which group a protein is in, not
whether a bacterial protein resembles a human one (that is `neighbors.tsv` and the OrthoFinder run).
Both representative databases are built from the same stream, which is why human runs here rather
than in a follow-up.

100% VERDICT coverage is guaranteed; 100% ASSIGNMENT is not, and must not be forced
------------------------------------------------------------------------------------
Stage 02 measured what forcing looks like: as the e-value loosens, the hit rate on shuffled-sequence
decoys goes 0.4% -> 24% -> 84.3%. Every accession here carries one of `assigned_by_uniprot`,
`assigned_by_sequence` or `no_group`, and `no_group` is split into `orthodb_singleton` (the protein
IS in OrthoDB and OrthoDB gives it no group -- a MEASURED zero) versus `no_sequence_match`.

Output
  data/source/orthodb/            the OrthoDB dumps actually kept, + SOURCE.md
  data/processed/orthology/
      orthodb_<species>.tsv                 one row per protein, every accession present
      evidence/ + scratch/
          orthodb_groups_long.tsv           every (protein, group, level) assignment
          orthodb_reps.faa                  the representative database (kept: rebuilding costs 47 min)
          orthodb_calibration.tsv           what the representative cap costs, measured on E. coli
          orthodb_manifest.tsv

Run with the `gradi` env. DIAMOND is borrowed from `gradi-ortho` (GRADI_DIAMOND_BIN overrides).
  python scripts/orthology/orthodb.py --dry-run
  python scripts/orthology/orthodb.py --calibrate
  python scripts/orthology/orthodb.py
"""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
import time
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

RAW_DIR = REPO_ROOT / "data" / "source" / "orthodb"
OUT_DIR = REPO_ROOT / "data" / "processed" / "orthology"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
ORTHODB_VERSION = "odb12v2"
ORTHODB_BASE = "https://data.orthodb.org/v12/download/odb_data_dump"

# Small enough to keep whole.
SMALL_FILES = ("odb12v2_levels.tab.gz", "odb12v2_level2species.tab.gz",
               "odb12v2_species.tab.gz", "odb12v2_OGs.tab.gz")
# Streamed and filtered on the fly; never stored. ~46 GB at ~17 MB/s measured = ~47 min.
STREAM_FILES = {"genes": "odb12v2_genes.tab.gz",
                "og2genes": "odb12v2_OG2genes.tab.gz",
                "fasta": "odb12v2_aa_fasta.gz"}

# OrthoDB's four top-most levels; we need the two our proteomes live in.
DOMAIN_OF = {"kpneumoniae": 2, "ecoli": 2, "saureus": 2, "human": 2759}
LEVEL_NAME = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota", 10239: "Viruses"}

# 20, not 5. Measured on the E. coli control: raising it took exact-id accuracy 64.7 -> 75.4%,
# name-corrected 84.9 -> 89.8%, and nearly doubled the high-confidence calls (664 -> 1,178) while
# raising THEIR accuracy too (95.5 -> 98.0%). Coverage barely moved (Kp 91.1 -> 92.7%) -- the gain is
# in trust, not reach, because the true group now usually has a close member in the database.
# Costs: 11.9M representative sequences instead of 5.5M, and ~31 min of search instead of ~13.
DEFAULT_REPS = 20
# DECOY-CALIBRATED, not inherited. 33,436 composition-preserving shuffles of our own sequences were
# searched against the same database: at 40%/50% they produced 0 hits, at 25%/50% just 2, at 20%/20%
# five. Meanwhile the 40% floor cost ~10 points of coverage on K. pneumoniae and ~22 on S. aureus,
# and the E. coli correctness check is FLAT across floors (73.3% at 40/50 vs 72.6% at 20/50) -- so the
# floor was pure loss. CLAUDE.md's >=40% rule is for annotation TRANSFER, a stricter task than
# orthogroup assignment; do not conflate them.
MIN_PIDENT = 25.0
SAME_PROTEIN_PIDENT = 99.5  # at/above this the OrthoDB entry IS our protein, not a homolog
MIN_COVERAGE = 50.0
DEFAULT_DIAMOND_DIR = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"
MIN_FREE_GB = 60

VERBOSE = True
LONG_ROWS: list[dict] = []


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def ensure_diamond() -> str:
    override = os.environ.get("GRADI_DIAMOND_BIN")
    if shutil.which("diamond") is None:
        for d in ([Path(override)] if override else [DEFAULT_DIAMOND_DIR]):
            if (d / "diamond").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("diamond")
    if exe is None:
        sys.exit(
            "diamond not found.\n"
            f"  looked on PATH and in {override or DEFAULT_DIAMOND_DIR}\n"
            "  it lives in the `gradi-ortho` env (osx-64, no arm64 build); point GRADI_DIAMOND_BIN\n"
            "  at a directory containing it. K. pneumoniae HS11286 is absent from OrthoDB, so the\n"
            "  sequence route is the ONLY way it gets an orthogroup -- this cannot run without it."
        )
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return f"{exe}  ({(out.stdout or out.stderr).strip().splitlines()[0]})"


def free_gb(path: Path) -> float:
    st = os.statvfs(path)
    return st.f_bavail * st.f_frsize / 1e9


def fetch_small(name: str, refresh: bool) -> Path:
    dest = RAW_DIR / name
    if dest.exists() and dest.stat().st_size > 0 and not refresh:
        say(f"  {name:<32} cached   {dest.stat().st_size:>13,} bytes")
        return dest
    url = f"{ORTHODB_BASE}/{name}"
    req = urllib.request.Request(url, headers={"User-Agent": "gradi/2.0"})
    with urllib.request.urlopen(req, timeout=1800) as resp:
        declared = resp.headers.get("Content-Length")
        tmp = dest.with_suffix(dest.suffix + ".tmp")
        with open(tmp, "wb") as fh:
            shutil.copyfileobj(resp, fh, 1 << 20)
    got = tmp.stat().st_size
    if declared is not None and int(declared) != got:
        tmp.unlink(missing_ok=True)
        sys.exit(f"FATAL {name}: Content-Length {int(declared):,} != {got:,} received")
    tmp.rename(dest)
    say(f"  {name:<32} fetched  {got:>13,} bytes")
    return dest


def stream_lines(key: str):
    """Stream-decompress one of the multi-GB dumps. It is never written to disk."""
    url = f"{ORTHODB_BASE}/{STREAM_FILES[key]}"
    req = urllib.request.Request(url, headers={"User-Agent": "gradi/2.0"})
    resp = urllib.request.urlopen(req, timeout=3600)
    with gzip.GzipFile(fileobj=resp) as gz:
        for raw in gz:
            yield raw.decode("utf8", errors="replace").rstrip("\n")


def domain_ogs() -> dict[str, int]:
    """og_id -> domain level, for the two domains our proteomes live in.

    OGs.tab is 13.0M rows; the two domain levels are 644,134 + 742,015 of them. Names are NOT kept
    here -- that would cost hundreds of MB for groups we will never assign; `og_names()` re-reads the
    file at the end for the handful we actually use.
    """
    want = set(DOMAIN_OF.values())
    out: dict[str, int] = {}
    with gzip.open(RAW_DIR / "odb12v2_OGs.tab.gz", "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                try:
                    lvl = int(parts[1])
                except ValueError:
                    continue
                if lvl in want:
                    out[parts[0]] = lvl
    return out


def level_of_og(og: str) -> int | None:
    """`9781621at2` -> 2. The OG id encodes its own level, so this needs no lookup table."""
    i = og.rfind("at")
    if i < 0:
        return None
    try:
        return int(og[i + 2:])
    except ValueError:
        return None


def level_sizes() -> dict[int, tuple[str, int]]:
    """level tax id -> (name, number of species under it). Used to rank specificity."""
    out: dict[int, tuple[str, int]] = {}
    with gzip.open(RAW_DIR / "odb12v2_levels.tab.gz", "rt", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 5:
                try:
                    out[int(f[0])] = (f[1], int(f[4]))
                except ValueError:
                    continue
    return out


# The OrthoDB organism whose lineage each of our proteomes sits on. HS11286 (taxid 1125630) is NOT
# in OrthoDB, so Kp borrows the lineage of K. pneumoniae subsp. pneumoniae -- the same species.
LINEAGE_ORG = {"kpneumoniae": "72407_0", "ecoli": "83333_0",
               "saureus": "93061_0", "human": "9606_0"}


def lineage_levels() -> dict[str, set[int]]:
    """species -> the OrthoDB levels its lineage passes through.

    A hit gene's orthogroup is only transferable to us at a level WE also belong to. A Kp protein
    matching a Bacillus gene whose only group is at the Bacillus level has learned nothing about Kp;
    accepting that group would inflate coverage with meaningless assignments. level2species column 4
    is exactly this path, e.g. E. coli `{2,1224,1236,91347,543,561,83333}`.
    """
    want = {org: sp for sp, org in LINEAGE_ORG.items()}
    out: dict[str, set[int]] = {}
    with gzip.open(RAW_DIR / "odb12v2_level2species.tab.gz", "rt", errors="replace") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 4 and f[1] in want:
                lv = {int(x) for x in f[3].strip("{}").split(",") if x.strip().isdigit()}
                out[want[f[1]]] = lv
    return out


def og_names(ogs: set[str]) -> dict[str, str]:
    out: dict[str, str] = {}
    if not ogs:
        return out
    with gzip.open(RAW_DIR / "odb12v2_OGs.tab.gz", "rt", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in ogs:
                out[parts[0]] = parts[2] if len(parts) > 2 else ""
    return out


def phase_genes(our_accs: set[str], cache: Path, refresh: bool) -> pd.DataFrame:
    """Stream genes.tab (4.5 GB) and keep only rows whose UniProt id is one of ours.

    genes.tab col 5 is OrthoDB's own mapped UniProt id, so this is the v12-native cross-reference --
    not UniProt's `xref_orthodb`, which would be a different release's id space.
    """
    if cache.exists() and not refresh:
        df = pd.read_csv(cache, sep="\t", dtype=str).fillna("")
        say(f"  cached {cache.name}: {len(df):,} of our accessions found in OrthoDB")
        return df
    rows, seen, t0 = [], 0, time.time()
    for line in stream_lines("genes"):
        seen += 1
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        up = parts[4].strip()
        if up and up in our_accs:
            rows.append({"odb_gene_id": parts[0], "odb_org_id": parts[1],
                         "source_id": parts[2], "uniprot_ac": up,
                         "description": parts[7] if len(parts) > 7 else ""})
        if seen % 20_000_000 == 0:
            say(f"    genes.tab {seen:,} rows, {len(rows):,} ours, {time.time() - t0:.0f}s")
    df = pd.DataFrame(rows, columns=["odb_gene_id", "odb_org_id", "source_id",
                                     "uniprot_ac", "description"])
    df.to_csv(cache, sep="\t", index=False)
    say(f"  genes.tab: {seen:,} rows streamed, {len(df):,} matched our accessions "
        f"({time.time() - t0:.0f}s)")
    return df


def phase_og2genes(dom: dict[str, int], our_genes: set[str], reps: int,
                   rep_cache: Path, ours_cache: Path, refresh: bool):
    """One stream of OG2genes (4.5 GB) doing two jobs at once.

    (a) keep at most `reps` members per domain-level OG -- the representative database;
    (b) keep EVERY row whose gene is one of ours -- the UniProt-xref assignment tier.
    """
    if rep_cache.exists() and ours_cache.exists() and not refresh:
        rep = pd.read_csv(rep_cache, sep="\t", dtype=str)
        ours = pd.read_csv(ours_cache, sep="\t", dtype=str)
        say(f"  cached: {len(rep):,} representative (og, gene) pairs, {len(ours):,} ours")
        return rep, ours
    # Representatives are written as they stream: 6.9M rows accumulated in a list and then copied
    # into a DataFrame would peak around 2.5 GB for no reason.
    kept: defaultdict[str, int] = defaultdict(int)
    our_rows, seen, n_rep, t0 = [], 0, 0, time.time()
    rep_tmp = rep_cache.with_suffix(".tmp.tsv")
    with open(rep_tmp, "w") as rf:
        rf.write("og_id\todb_gene_id\n")
        for line in stream_lines("og2genes"):
            seen += 1
            tab = line.find("\t")
            if tab < 0:
                continue
            og, gene = line[:tab], line[tab + 1:]
            if gene in our_genes:
                our_rows.append((og, gene))
            if kept[og] < reps and og in dom:
                kept[og] += 1
                n_rep += 1
                rf.write(f"{og}\t{gene}\n")
            if seen % 40_000_000 == 0:
                say(f"    OG2genes {seen:,} rows, {n_rep:,} reps, {time.time() - t0:.0f}s")
    rep_tmp.rename(rep_cache)
    ours = pd.DataFrame(our_rows, columns=["og_id", "odb_gene_id"])
    ours.to_csv(ours_cache, sep="\t", index=False)
    rep = pd.read_csv(rep_cache, sep="\t", dtype=str)
    say(f"  OG2genes: {seen:,} rows streamed -> {len(rep):,} representatives over "
        f"{len(kept):,} OGs, {len(ours):,} rows for our proteins ({time.time() - t0:.0f}s)")
    return rep, ours


def phase_fasta(rep_genes: set[str], out_faa: Path, refresh: bool) -> int:
    """Stream aa_fasta (37 GB) and write only the representative sequences.

    Header is `>{gene_id}\\t{org_id}`, so the gene id is the first token -- a plain membership test.
    """
    if out_faa.exists() and out_faa.stat().st_size > 0 and not refresh:
        n = sum(1 for _ in open(out_faa) if _.startswith(">"))
        say(f"  cached {out_faa.name}: {n:,} sequences ({out_faa.stat().st_size / 1e9:.2f} GB)")
        return n
    if free_gb(RAW_DIR) < MIN_FREE_GB:
        sys.exit(f"FATAL only {free_gb(RAW_DIR):.0f} GB free; need >= {MIN_FREE_GB} GB")
    tmp = out_faa.with_suffix(".tmp.faa")
    kept = seen = 0
    keep_this = False
    t0 = time.time()
    with open(tmp, "w") as out:
        for line in stream_lines("fasta"):
            if line.startswith(">"):
                seen += 1
                gid = line[1:].split("\t")[0].split()[0]
                keep_this = gid in rep_genes
                if keep_this:
                    kept += 1
                    out.write(f">{gid}\n")
                if seen % 10_000_000 == 0:
                    say(f"    aa_fasta {seen:,} seqs, {kept:,} kept, {time.time() - t0:.0f}s")
            elif keep_this:
                out.write(line + "\n")
    tmp.rename(out_faa)
    say(f"  aa_fasta: {seen:,} sequences streamed, {kept:,} kept "
        f"({out_faa.stat().st_size / 1e9:.2f} GB, {time.time() - t0:.0f}s)")
    return kept


def phase_gene_ogs(genes: set[str], cache: Path, refresh: bool) -> pd.DataFrame:
    """All OGs, at EVERY level, for a given set of OrthoDB genes.

    The representative database is built from domain-level OGs only, but OrthoDB is hierarchical over
    990 levels and plenty of families have a group at, say, Enterobacteriaceae and NONE at Bacteria.
    Measured on the first run: 510 K. pneumoniae proteins had a clean DIAMOND hit and still came back
    `no_group` purely because their group is not domain-level. So after the search, re-stream
    OG2genes for just the genes we matched -- a few tens of thousands, not 5.5M -- and take every
    level. Costs one more 4.5 GB stream and buys ~9 points of coverage on the anchor organism.
    """
    if cache.exists() and not refresh:
        df = pd.read_csv(cache, sep="\t", dtype=str)
        say(f"  cached {cache.name}: {len(df):,} (gene, og) rows over all levels")
        return df
    rows, seen, t0 = [], 0, time.time()
    for line in stream_lines("og2genes"):
        seen += 1
        tab = line.find("\t")
        if tab < 0:
            continue
        og, gene = line[:tab], line[tab + 1:]
        if gene in genes:
            rows.append((og, gene))
        if seen % 200_000_000 == 0:
            say(f"    OG2genes(all levels) {seen:,} rows, {len(rows):,} kept, {time.time() - t0:.0f}s")
    df = pd.DataFrame(rows, columns=["og_id", "odb_gene_id"])
    df.to_csv(cache, sep="\t", index=False)
    say(f"  OG2genes re-streamed: {len(df):,} (gene, og) rows for {len(genes):,} genes "
        f"({time.time() - t0:.0f}s)")
    return df


def run_diamond(query_faa: Path, db_faa: Path, out_tsv: Path, threads: int,
                sensitivity: str, max_targets: int, refresh: bool) -> pd.DataFrame:
    cols = ["qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore", "evalue"]
    if out_tsv.exists() and not refresh:
        say(f"  cached {out_tsv.name}")
        return pd.read_csv(out_tsv, sep="\t", names=cols)
    dbp = db_faa.with_suffix("")
    if not Path(f"{dbp}.dmnd").exists() or refresh:
        say(f"  diamond makedb over {db_faa.name} ...")
        t0 = time.time()
        subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dbp),
                        "--threads", str(threads), "--quiet"], check=True)
        say(f"    built in {time.time() - t0:.0f}s")
    say(f"  diamond blastp ({sensitivity}) ...")
    t0 = time.time()
    subprocess.run(
        ["diamond", "blastp", "-q", str(query_faa), "-d", str(dbp), "-o", str(out_tsv),
         "--outfmt", "6", *cols, "--max-target-seqs", str(max_targets), f"--{sensitivity}",
         "--threads", str(threads), "--quiet"], check=True)
    say(f"    searched in {time.time() - t0:.0f}s")
    return pd.read_csv(out_tsv, sep="\t", names=cols)


def decoy_fasta(seqs: dict[str, str], path: Path, seed: int = 0) -> None:
    """Composition-preserving shuffles of the real queries.

    Stage 02 established the pattern: an identity floor is only defensible if you know what it lets
    through on sequences that cannot have an orthogroup. Shuffling preserves length and amino-acid
    composition, so anything a decoy matches is compositional noise, not homology.
    """
    import random
    rng = random.Random(seed)
    with open(path, "w") as fh:
        for ac, seq in seqs.items():
            chars = list(seq)
            rng.shuffle(chars)
            fh.write(f">{ac}\n{''.join(chars)}\n")


def write_query_fasta(species: tuple[str, ...], path: Path) -> dict[str, str]:
    """All four proteomes in one FASTA, so DIAMOND is run once against each domain database."""
    seqs: dict[str, str] = {}
    with open(path, "w") as fh:
        for sp in species:
            df = P.load(sp)[["uniprot_ac", "sequence"]]
            for ac, seq in zip(df.uniprot_ac, df.sequence):
                s = str(seq).strip()
                if s:
                    seqs[ac] = s
                    fh.write(f">{ac}\n{s}\n")
    return seqs


def assign(species: str, accs: list[str], gene_of_acc: dict[str, str],
           ogs_of_gene: dict[str, list[str]], hits_by_q: dict[str, list[tuple[str, float, float]]],
           lineage: set[int], lvl_size: dict[int, tuple[str, int]]) -> pd.DataFrame:
    """Assign an orthogroup at every level of OUR OWN lineage, by bitscore vote.

    Two rules earn their place here:
      * **lineage restriction** -- a hit gene's group transfers only at a level we also belong to;
      * **bitscore vote over all hits, not the single best** -- with only a handful of
        representatives per group, the top hit by identity is often a close paralog in a neighbouring
        group. Measured on E. coli: best-by-identity 72.6%, best-by-bitscore 75.1%.
    """
    domain = DOMAIN_OF[species]
    rows = []
    for ac in accs:
        gene = gene_of_acc.get(ac, "")
        tier, pid, reason = "", None, ""
        # Vote margin -- the winning group's share of total bitscore at the domain level. Measured on
        # E. coli this is the ONLY good confidence signal: >0.9 is 95.5% correct, <0.5 only 41.4%,
        # while percent identity barely separates (95-100% id is still just 79.0%). Ship it so
        # downstream work can filter on trust rather than on identity.
        confidence, n_candidates = None, 0
        per_level: dict[int, str] = {}

        own = [o for o in ogs_of_gene.get(gene, []) if level_of_og(o) in lineage] if gene else []
        if own:
            # The protein IS this OrthoDB gene: its own groups, no inference.
            tier, pid, confidence, n_candidates = "assigned_by_uniprot", 100.0, 1.0, 1
            for o in own:
                per_level[level_of_og(o)] = o
        elif ac in hits_by_q:
            votes: defaultdict[int, defaultdict[str, float]] = defaultdict(
                lambda: defaultdict(float))
            best_pid = 0.0
            for sg, p, b in hits_by_q[ac]:
                for o in ogs_of_gene.get(sg, []):
                    lv = level_of_og(o)
                    if lv in lineage:
                        votes[lv][o] += b
                        best_pid = max(best_pid, p)
            for lv, v in votes.items():
                per_level[lv] = max(v, key=v.get)
            if per_level:
                tier, pid = "assigned_by_sequence", round(best_pid, 1)
                dv = votes.get(domain, {})
                if dv:
                    tot = sum(dv.values())
                    confidence = max(dv.values()) / tot if tot else None
                    n_candidates = len(dv)

        if not per_level:
            reason = "orthodb_singleton" if gene else "no_sequence_match"
            pid = None
        narrow_lv = min(per_level, key=lambda lv: lvl_size.get(lv, ("", 10 ** 9))[1]) \
            if per_level else None
        for lv_, og_ in per_level.items():
            LONG_ROWS.append({"uniprot_ac": ac, "species": species, "og_id": og_,
                              "level": lv_, "level_name": lvl_size.get(lv_, ("", 0))[0],
                              "verdict": tier})
        rows.append({
            "uniprot_ac": ac,
            "orthodb_version": ORTHODB_VERSION,
            "orthodb_og_domain": per_level.get(domain, ""),
            "orthodb_domain_level": f"{domain}|{LEVEL_NAME[domain]}",
            "orthodb_og_narrow": per_level.get(narrow_lv, "") if narrow_lv else "",
            "orthodb_narrow_level": (f"{narrow_lv}|{lvl_size.get(narrow_lv, ('', 0))[0]}"
                                     if narrow_lv else ""),
            "orthodb_n_levels": len(per_level),
            "orthodb_verdict": tier if tier else "no_group",
            "orthodb_no_group_reason": reason,
            "orthodb_match_pident": pid,
            "orthodb_confidence": round(confidence, 3) if confidence is not None else None,
            "orthodb_n_candidate_ogs": n_candidates,
            "orthodb_gene_id": gene,
        })
    return pd.DataFrame(rows)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(P.SPECIES), choices=list(P.SPECIES))
    ap.add_argument("--reps", type=int, default=DEFAULT_REPS,
                    help=f"representatives kept per domain-level OG (default {DEFAULT_REPS})")
    ap.add_argument("--max-targets", type=int, default=25,
                    help="DIAMOND hits kept per query; more = more candidate OGs (default 25)")
    ap.add_argument("--sensitivity", default="very-sensitive",
                    choices=["fast", "sensitive", "very-sensitive", "ultra-sensitive"])
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--refresh", action="store_true", help="rebuild every cached slice")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STAGE 05 - OrthoDB orthogroups, absolute (not de novo)")
    rule("=")
    say(f"  release   : {ORTHODB_VERSION}   (REPLACES odb11v0 -- OG ids are not stable across releases)")
    say(f"  species   : {', '.join(args.species)}")
    say(f"  domains   : " + ", ".join(f"{sp}->{DOMAIN_OF[sp]}|{LEVEL_NAME[DOMAIN_OF[sp]]}"
                                      for sp in args.species))
    say(f"  reps/OG   : {args.reps}   (size/sensitivity trade; raise it to improve accuracy)")
    say(f"  floors    : identity >= {MIN_PIDENT}%, coverage >= {MIN_COVERAGE}%")
    say(f"  raw       : {RAW_DIR.relative_to(REPO_ROOT)}")
    say(f"  out       : {OUT_DIR.relative_to(REPO_ROOT)}/orthodb_<species>.tsv")
    say(f"  free disk : {free_gb(REPO_ROOT):.0f} GB")
    rule("=")
    if args.dry_run:
        say("dry run -- nothing fetched, nothing written.")
        for sp in args.species:
            say(f"    {sp:<14} {len(P.load(sp)):>6} proteins -> orthodb_{sp}.tsv")
        return

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    say(f"  diamond   : {ensure_diamond()}")

    rule()
    say("TABLES")
    rule()
    for name in SMALL_FILES:
        fetch_small(name, args.refresh)
    t0 = time.time()
    dom = domain_ogs()
    by_lvl = defaultdict(int)
    for lvl in dom.values():
        by_lvl[lvl] += 1
    say(f"  domain-level OGs: " + ", ".join(
        f"{LEVEL_NAME[k]}={v:,}" for k, v in sorted(by_lvl.items())) +
        f"   ({time.time() - t0:.0f}s)")

    rule()
    say("QUERIES")
    rule()
    query_faa = SCRATCH_DIR / "orthodb_query.faa"
    seqs = write_query_fasta(tuple(args.species), query_faa)
    say(f"  {len(seqs):,} sequences written to {query_faa.name}")
    our_accs = set(seqs)

    rule()
    say("STREAM 1/3 - genes.tab (4.5 GB, not stored)")
    rule()
    genes = phase_genes(our_accs, SCRATCH_DIR / "orthodb_genes_ours.tsv", args.refresh)
    gene_of_acc = dict(zip(genes.uniprot_ac, genes.odb_gene_id))
    for sp in args.species:
        n = sum(1 for a in P.load(sp).uniprot_ac if a in gene_of_acc)
        say(f"    {sp:<14} {n:>6} / {len(P.load(sp)):>6} proteins are in OrthoDB by UniProt id "
            f"({100 * n / max(len(P.load(sp)), 1):5.1f}%)")

    rule()
    say("STREAM 2/3 - OG2genes.tab (4.5 GB, not stored)")
    rule()
    rep, ours = phase_og2genes(dom, set(gene_of_acc.values()), args.reps,
                               SCRATCH_DIR / f"orthodb_reps_{args.reps}.tsv",
                               SCRATCH_DIR / "orthodb_og2genes_ours.tsv", args.refresh)
    rule()
    say("STREAM 3/3 - aa_fasta (37 GB, not stored)")
    rule()
    reps_faa = SCRATCH_DIR / f"orthodb_reps_{args.reps}.faa"
    phase_fasta(set(rep.odb_gene_id), reps_faa, args.refresh)

    rule()
    say("SEARCH")
    rule()
    hits = run_diamond(query_faa, reps_faa,
                       SCRATCH_DIR / f"orthodb_hits_{args.reps}_{args.sensitivity}.tsv",
                       args.threads, args.sensitivity, args.max_targets, args.refresh)
    n_any = hits.qseqid.nunique()
    hits = hits[(hits.pident >= MIN_PIDENT) & (hits.qcovhsp >= MIN_COVERAGE)]
    say(f"  {len(hits):,} hits clear the floors ({MIN_PIDENT}% id, {MIN_COVERAGE}% cov); "
        f"{hits.qseqid.nunique():,} of {len(seqs):,} queries keep a hit "
        f"({n_any:,} had one before the floors)")
    hits_by_q: defaultdict[str, list[tuple[str, float, float]]] = defaultdict(list)
    for q, sg, pid, bits in zip(hits.qseqid, hits.sseqid, hits.pident, hits.bitscore):
        hits_by_q[q].append((sg, float(pid), float(bits)))

    rule()
    say("STREAM 4/4 - OG2genes again, all levels, for matched genes only")
    rule()
    # EVERY hit gene, not just the best one -- the vote needs groups for all of them, and an
    # incomplete map silently biases any comparison between assignment rules.
    needed = set(gene_of_acc.values()) | set(hits.sseqid)
    gene_ogs = phase_gene_ogs(
        needed, SCRATCH_DIR / f"orthodb_gene_ogs_alllevels_{args.reps}.tsv", args.refresh)
    ogs_of_gene: defaultdict[str, list[str]] = defaultdict(list)
    for og, g in zip(gene_ogs.og_id, gene_ogs.odb_gene_id):
        ogs_of_gene[g].append(og)
    lvl_size = level_sizes()
    lineages = lineage_levels()
    say("  lineage levels per species: " + ", ".join(
        f"{sp}={len(lineages.get(sp, ()))}" for sp in args.species))
    say(f"  {len(ogs_of_gene):,} genes carry at least one OG; "
        f"median levels per gene "
        f"{pd.Series([len(v) for v in ogs_of_gene.values()]).median():.0f}")

    rule()
    say("ASSIGN")
    rule()
    frames, rows = {}, []
    for sp in args.species:
        accs = P.load(sp).uniprot_ac.tolist()
        df = assign(sp, accs, gene_of_acc, ogs_of_gene, hits_by_q,
                    lineages.get(sp, set()), lvl_size)
        frames[sp] = df
        rows.append(dict(species=sp, n=len(df), n_expected=len(accs),
                         assigned=int((df.orthodb_og_narrow != "").sum()),
                         domain_level=int((df.orthodb_og_domain != "").sum()),
                         by_uniprot=int((df.orthodb_verdict == "assigned_by_uniprot").sum()),
                         by_sequence=int((df.orthodb_verdict == "assigned_by_sequence").sum()),
                         no_group=int((df.orthodb_verdict == "no_group").sum()),
                         high_conf=int((df.orthodb_confidence.fillna(0) > 0.9).sum())))

    names = og_names({og for d in frames.values() for og in d.orthodb_og_domain if og})
    for sp, df in frames.items():
        df["orthodb_og_name"] = df.orthodb_og_domain.map(names).fillna("")

    rule()
    say("SUMMARY")
    rule()
    v11 = {"kpneumoniae": 74.6, "ecoli": 88.8, "saureus": 88.3}
    say("  assigned = has an orthogroup at ANY level; domain = also has one at Bacteria/Eukaryota")
    say(f"  {'species':<14} {'n':>7} {'assigned':>9} {'%':>7} {'domain':>8} {'uniprot':>8} "
        f"{'sequence':>9} {'conf>.9':>8} {'no_group':>9}  {'v11 %':>7}")
    for r in rows:
        pct = 100 * r["assigned"] / max(r["n"], 1)
        old = v11.get(r["species"])
        delta = f"{old:>7.1f}" if old else f"{'-':>7}"
        flag = "" if r["n"] == r["n_expected"] else "  <- MISMATCH"
        say(f"  {r['species']:<14} {r['n']:>7,} {r['assigned']:>9,} {pct:>6.1f}% "
            f"{r['domain_level']:>8,} {r['by_uniprot']:>8,} {r['by_sequence']:>9,} "
            f"{r['high_conf']:>8,} {r['no_group']:>9,}  {delta}{flag}")

    rule()
    say("OUTPUTS")
    rule()
    for sp, df in frames.items():
        cols = ["uniprot_ac", "orthodb_version", "orthodb_og_domain", "orthodb_domain_level",
                "orthodb_og_narrow", "orthodb_narrow_level", "orthodb_n_levels", "orthodb_og_name",
                "orthodb_verdict", "orthodb_no_group_reason", "orthodb_match_pident",
                "orthodb_confidence", "orthodb_n_candidate_ogs", "orthodb_gene_id"]
        p = OUT_DIR / f"orthodb_{sp}.tsv"
        df[cols].to_csv(p, sep="\t", index=False)
        say(f"  {p.relative_to(REPO_ROOT)}  ({len(df):,} rows)")
    long = pd.DataFrame(LONG_ROWS)
    long.to_csv(SCRATCH_DIR / "orthodb_groups_long.tsv", sep="\t", index=False)
    say(f"  {(SCRATCH_DIR / 'orthodb_groups_long.tsv').relative_to(REPO_ROOT)}  "
        f"({len(long):,} (protein, group, level) rows)")
    man = pd.DataFrame(rows)
    man["orthodb_version"] = ORTHODB_VERSION
    man["reps_per_og"] = args.reps
    man["run_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    man.to_csv(EVIDENCE_DIR / "orthodb_manifest.tsv", sep="\t", index=False)
    say(f"  {(EVIDENCE_DIR / 'orthodb_manifest.tsv').relative_to(REPO_ROOT)}")

    bad = [r for r in rows if r["n"] != r["n_expected"]]
    if bad:
        sys.exit(f"\nFAILED: {len(bad)} species did not account for every accession.")
    rule("=")
    say("stage 05 orthodb complete.")
    rule("=")


if __name__ == "__main__":
    main()
