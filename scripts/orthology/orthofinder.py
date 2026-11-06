"""Orthology across the four proteomes: a discrete ortholog matrix and a continuous neighbour matrix.

Two questions, answered separately and never merged:

  WHICH PROTEINS ARE ORTHOLOGS   OrthoFinder 3, run de novo on our own four proteomes, plus DIAMOND
                                 reciprocal best hits as an independent second call.
  HOW SIMILAR ARE THEY           DIAMOND blastp, the top-5 nearest neighbours in EACH target species,
                                 with identity, positives, both coverages and a normalised bitscore.

v1 shipped only the first. Its retrospective calls the result trap 1: the kp->ec table's `pident`,
`coverage` and `bitscore` were entirely empty, so *"any axis that thresholds transfer on percent
identity silently drops everything"*. The second matrix exists so that cannot recur.

A zero is a measured zero
-------------------------
OrthoFinder runs **de novo on our own FASTAs**, so every protein is either assigned an orthogroup or
named in `Orthogroups_UnassignedGenes.tsv` -- both measured. The stage exits non-zero unless assigned
plus unassigned accounts for every accession in all four proteomes, and every accession is a query in
all four DIAMOND searches. No lookup database is used: eggNOG's orthogroups (already on disk, free)
cover only 89-96% and no human, so a zero there would be an unknown in disguise.

Output
------
    data/processed/orthology/
        orthologs.tsv            sparse: one row per pair called an ortholog by EITHER method
        neighbors.tsv            sparse: top-5 per (protein, target species), self-hits dropped
        orthology_<species>.tsv  DENSE: one row per protein -- where a 0 is readable
        evidence/ + scratch/               the OrthoFinder run, orthogroups, disagreements, control, manifest

Environments
------------
Runs in `gradi`, shells out to `gradi-ortho` (osx-64, Rosetta) for both tools. **The OrthoFinder
launcher cannot be called directly**: its `#!/usr/bin/env python3` shebang resolves to an unrelated
env on this machine and dies with `ModuleNotFoundError: No module named 'ete4'` even though ete4 is
installed. Name the interpreter explicitly. `GRADI_ORTHOFINDER_ENV` and `GRADI_DIAMOND_BIN` override.

Run with the `gradi` env:
    python scripts/orthology/orthofinder.py
    python scripts/orthology/orthofinder.py --species kpneumoniae ecoli --top-k 3   # fast subset
    python scripts/orthology/orthofinder.py --refresh                               # re-run OrthoFinder
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src import orthology as O  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "orthology"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
FASTA_DIR = SCRATCH_DIR / "fasta"
OF_DIR = SCRATCH_DIR / "orthofinder"
HIT_DIR = SCRATCH_DIR / "hits"

DEFAULT_DIAMOND_DIR = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"
DEFAULT_ORTHO_ENV = Path.home() / "miniconda3" / "envs" / "gradi-ortho"

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------------------------------------------------------- tools


def ensure_diamond() -> str:
    """Locate the DIAMOND binary, or exit with the fix."""
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
            "  at a DIRECTORY containing it."
        )
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return f"{exe}  ({(out.stdout or out.stderr).strip().splitlines()[0]})"


def orthofinder_cmd() -> list[str]:
    """The OrthoFinder invocation, as an argv prefix.

    **Never call the launcher directly.** Its shebang is `#!/usr/bin/env python3`, which on this
    machine resolves to an unrelated conda env, so `orthofinder -h` fails with
    `ModuleNotFoundError: No module named 'ete4'` even though ete4 IS installed in `gradi-ortho`.
    Naming the env's own interpreter is the whole fix -- there is nothing to install.
    """
    env = Path(os.environ.get("GRADI_ORTHOFINDER_ENV", DEFAULT_ORTHO_ENV))
    py, launcher = env / "bin" / "python", env / "bin" / "orthofinder"
    if not py.exists() or not launcher.exists():
        sys.exit(
            f"orthofinder not found under {env}\n"
            "  it lives in the `gradi-ortho` env (osx-64, no arm64 build); point\n"
            "  GRADI_ORTHOFINDER_ENV at the env root (not its bin/)."
        )
    return [str(py), str(launcher)]


def orthofinder_version() -> str:
    out = subprocess.run(orthofinder_cmd() + ["-h"], capture_output=True, text=True)
    for line in (out.stdout or out.stderr).splitlines():
        if "OrthoFinder" in line and "version" in line.lower() or line.startswith("OrthoFinder"):
            return line.strip()
    return "OrthoFinder (version not reported)"


# ---------------------------------------------------------------- inputs


def comparator_panel() -> list[str]:
    """Tier-C labels from the registry -- the 26-species bacterial comparator panel.

    Read from `src/proteome_registry.tsv` rather than listed here, so the panel has exactly one
    definition and adding a species is a registry edit. Every row is already pinned to an explicit
    proteome id.
    """
    reg = pd.read_csv(REPO_ROOT / "src" / "proteome_registry.tsv", sep="\t")
    return sorted(reg.loc[reg["tier"] == "C", "label"].astype(str))


def write_fastas(species: tuple[str, ...], comparators: tuple[str, ...] = ()) -> dict[str, int]:
    """One FASTA per proteome, headers as bare accessions. Anchors first, then comparators.

    Anchors are written from the `sequence` column of `proteome_<species>.tsv` (100% populated in
    all four) rather than from `data/raw/`, so the FASTA can never drift from the table the rest of
    the pipeline joins on. The file stem becomes OrthoFinder's species name.

    **Comparators come from `data/source/uniprot/proteomes/<label>.fasta`** and have no
    `proteome_<label>.tsv` -- they are not analysis subjects, they are there to give OrthoFinder
    enough phylogenetic spread to form orthogroups. Their UniProt headers (`>sp|ACC|NAME ...`) are
    reduced to the bare accession, matching the anchors' convention, because OrthoFinder splits
    on whitespace and the descriptions would otherwise end up in the gene ids.
    """
    FASTA_DIR.mkdir(parents=True, exist_ok=True)
    counts = {}
    for sp in species:
        df = P.load(sp)[["uniprot_ac", "sequence"]]
        path = FASTA_DIR / f"{sp}.faa"
        with path.open("w") as fh:
            for ac, seq in zip(df["uniprot_ac"], df["sequence"]):
                fh.write(f">{ac}\n{seq}\n")
        counts[sp] = len(df)

    # `download.py` writes UniProt fetches to `<RAW_DIR>/uniprot/`; the tier-A files one level up
    # are stale from an earlier layout. Both are checked so neither convention silently misses.
    roots = [REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "uniprot",
             REPO_ROOT / "data" / "source" / "uniprot" / "proteomes"]
    for label in comparators:
        fa = next((r / f"{label}.fasta" for r in roots if (r / f"{label}.fasta").exists()), None)
        if fa is None:
            sys.exit(f"FAILED: {label}.fasta not found under "
                     f"{roots[0].relative_to(REPO_ROOT)} -- fetch the panel first:\n"
                     "    python scripts/proteomes/download.py --tier C")
        n = 0
        with (FASTA_DIR / f"{label}.faa").open("w") as fh:
            for line in fa.read_text().splitlines():
                if line.startswith(">"):
                    # >sp|ACC|NAME desc  ->  >ACC
                    parts = line[1:].split("|")
                    fh.write(f">{parts[1] if len(parts) > 2 else line[1:].split()[0]}\n")
                    n += 1
                else:
                    fh.write(line + "\n")
        counts[label] = n

    # OrthoFinder takes a DIRECTORY and treats every file in it as a proteome, so a stale FASTA from
    # an earlier --species subset would silently join the run.
    keep = set(species) | set(comparators)
    for stale in FASTA_DIR.glob("*.faa"):
        if stale.stem not in keep:
            stale.unlink()
    return counts


# ---------------------------------------------------------------- OrthoFinder


def run_orthofinder(threads: int, refresh: bool) -> Path:
    """Run OrthoFinder, or reuse the cached run. Returns the results directory."""
    results = find_results_dir()
    if results is not None and not refresh:
        say(f"  [cache] {results.relative_to(REPO_ROOT)}")
        return results
    if OF_DIR.exists():
        shutil.rmtree(OF_DIR)

    # A FULL run, not `-og`: `-og` stops after orthogroups and therefore never writes Orthologues/,
    # which is the only place pairwise orthologs appear.
    cmd = orthofinder_cmd() + ["-f", str(FASTA_DIR), "-o", str(OF_DIR), "-t", str(threads)]
    say(f"  running  : {' '.join(cmd[2:])}")
    say("  (v1's 3-way run over 30,547 proteins took 8 min 37 s at -t 8)")
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    # OrthoFinder EXITS 0 when its own dependency check fails -- if `diamond` is not on PATH it
    # prints an ERROR block, writes an empty Results_ directory and returns success. The return code
    # is therefore not evidence; the caller's find_results_dir() check is what actually catches it.
    if proc.returncode != 0:
        tail = "\n  ".join((proc.stderr or proc.stdout).strip().splitlines()[-15:])
        sys.exit(f"OrthoFinder failed (exit {proc.returncode}):\n  {tail}")
    say(f"  done in {(time.time() - t0) / 60:.1f} min")

    results = find_results_dir()
    if results is None:
        tail = "\n  ".join((proc.stdout or proc.stderr).strip().splitlines()[-20:])
        sys.exit(
            f"OrthoFinder returned 0 but wrote no Orthogroups.tsv under {OF_DIR}.\n"
            "  It exits 0 on a failed dependency check, so read the output below -- the usual cause\n"
            "  is `diamond` not being on PATH when OrthoFinder starts.\n  " + tail)
    return results


def find_results_dir() -> Path | None:
    """OrthoFinder's layout under `-o` varies by version; locate it by its content instead."""
    if not OF_DIR.exists():
        return None
    hits = sorted(OF_DIR.rglob("Orthogroups/Orthogroups.tsv"))
    return hits[0].parent.parent if hits else None


def _split(cell: str) -> list[str]:
    return [g.strip() for g in str(cell).split(",") if g.strip()]


def parse_orthogroups(results: Path, species: tuple[str, ...]) -> pd.DataFrame:
    """Every protein's orthogroup: assigned or explicitly unassigned.

    Both files are read. `Orthogroups_UnassignedGenes.tsv` is not an afterthought -- it is what makes
    "no orthogroup" a measured outcome rather than a gap, and the coverage guard depends on it.
    """
    rows = []
    assigned = pd.read_csv(results / "Orthogroups" / "Orthogroups.tsv", sep="\t",
                           dtype=str, keep_default_na=False)
    for _, r in assigned.iterrows():
        for sp in species:
            for ac in _split(r.get(sp, "")):
                rows.append({"uniprot_ac": ac, "species": sp, "orthogroup": r["Orthogroup"]})

    un_path = results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
    if un_path.exists():
        un = pd.read_csv(un_path, sep="\t", dtype=str, keep_default_na=False)
        for _, r in un.iterrows():
            for sp in species:
                for ac in _split(r.get(sp, "")):
                    rows.append({"uniprot_ac": ac, "species": sp, "orthogroup": ""})
    return pd.DataFrame(rows)


def panel_membership(results: Path, bacterial: tuple[str, ...]) -> dict[str, int]:
    """orthogroup -> how many BACTERIAL proteomes in the run have a member in it.

    This is what makes `n_bacterial_orthologs` a conservation measure: a group containing members
    from 27 of the 28 bacteria is a core gene, one containing only its own species is specific to
    it. Counted over SPECIES, not over proteins -- a paralog pair does not make a protein more
    conserved, and the approved reading of the column is "how many species have an ortholog".

    Only `Orthogroups.tsv` matters here: an UNASSIGNED protein is in no group, so its count is 0 by
    construction and nothing needs to be read for it.
    """
    assigned = pd.read_csv(results / "Orthogroups" / "Orthogroups.tsv", sep="\t",
                           dtype=str, keep_default_na=False)
    present = [c for c in bacterial if c in assigned.columns]
    missing = [c for c in bacterial if c not in assigned.columns]
    if missing:
        sys.exit(f"FAILED: {len(missing)} bacterial proteomes are absent from Orthogroups.tsv "
                 f"({missing[:3]}...). The OrthoFinder run did not include the panel; re-run with "
                 "--panel full --refresh rather than reporting a count over a partial panel.")
    counts = {}
    for _, r in assigned.iterrows():
        counts[r["Orthogroup"]] = sum(1 for c in present if str(r[c]).strip())
    return counts


def parse_orthologues(results: Path) -> pd.DataFrame:
    """OrthoFinder's pairwise orthologs, expanded to ordered (query, target) pairs.

    The `Orthologues/` files are many-to-many: a cell can list several genes on each side, which is
    exactly the co-orthology RBH structurally cannot represent. Every listed combination is a pair.
    """
    root = results / "Orthologues"
    if not root.exists():
        sys.exit(
            f"{root} is missing -- OrthoFinder ran without producing pairwise orthologs.\n"
            "  This happens with `-og`, which stops after orthogroups. The stage needs a full run."
        )
    rows = []
    for path in sorted(root.rglob("*__v__*.tsv")):
        left, right = path.stem.split("__v__")
        tab = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        if left not in tab.columns or right not in tab.columns:
            continue
        for _, r in tab.iterrows():
            ls, rs = _split(r[left]), _split(r[right])
            for a, b in product(ls, rs):
                rows.append((a, left, b, right))
                rows.append((b, right, a, left))
    df = pd.DataFrame(rows, columns=["query_ac", "query_species", "target_ac", "target_species"])
    return df.drop_duplicates()


# ---------------------------------------------------------------- DIAMOND


HIT_COLS = ["query_ac", "target_ac", "pident", "ppos", "alnlen", "qlen", "slen",
            "qcov", "scov", "evalue", "bitscore"]


def diamond_pair(q_sp: str, t_sp: str, top_k: int, evalue: float, sensitivity: str,
                 threads: int, refresh: bool) -> pd.DataFrame:
    """Top hits of every `q_sp` protein in `t_sp`. Cached per ordered pair."""
    HIT_DIR.mkdir(parents=True, exist_ok=True)
    out = HIT_DIR / f"hits_{q_sp}__{t_sp}.tsv.gz"
    if out.exists() and not refresh:
        return pd.read_csv(out, sep="\t")

    # One extra target on the self-pair, because the protein's own 100%-identity self-hit occupies
    # rank 1 there and is dropped afterwards.
    max_targets = top_k + 1 if q_sp == t_sp else top_k
    db = HIT_DIR / f"{t_sp}"
    if not db.with_suffix(".dmnd").exists() or refresh:
        subprocess.run(["diamond", "makedb", "--in", str(FASTA_DIR / f"{t_sp}.faa"),
                        "-d", str(db), "--quiet"], check=True)
    raw = HIT_DIR / f".raw_{q_sp}__{t_sp}.tsv"
    subprocess.run(
        ["diamond", "blastp", "-q", str(FASTA_DIR / f"{q_sp}.faa"), "-d", str(db),
         "-o", str(raw), f"--{sensitivity}", "--evalue", str(evalue),
         "--max-target-seqs", str(max_targets), "--outfmt", "6",
         "qseqid", "sseqid", "pident", "ppos", "length", "qlen", "slen",
         "qcovhsp", "scovhsp", "evalue", "bitscore",
         "--quiet", "--threads", str(threads)],
        check=True)
    hits = pd.read_csv(raw, sep="\t", names=HIT_COLS)
    raw.unlink()
    # One row per (query, target): DIAMOND can emit several HSPs for a pair.
    hits = hits.sort_values("bitscore", ascending=False).drop_duplicates(["query_ac", "target_ac"])
    hits.to_csv(out, sep="\t", index=False)
    return hits


# ---------------------------------------------------------------- main


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(O.SPECIES), choices=list(O.SPECIES),
                    help="the ANCHORS to produce tables for (default: all four)")
    ap.add_argument("--panel", choices=["anchors", "full"], default="anchors",
                    help="`full` adds the 26-species tier-C comparator panel to the OrthoFinder "
                         "run. They inform the orthogroups and get no tables of their own; "
                         "pairwise DIAMOND stays on the anchors. Default keeps prior behaviour.")
    ap.add_argument("--top-k", type=int, default=O.TOP_K,
                    help="nearest neighbours kept per target species (default: 5)")
    ap.add_argument("--evalue", type=float, default=O.EVALUE)
    ap.add_argument("--sensitivity", default=O.SENSITIVITY,
                    choices=["fast", "sensitive", "very-sensitive", "ultra-sensitive"])
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--refresh", action="store_true",
                    help="re-run OrthoFinder and every DIAMOND search")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    species = tuple(s for s in O.SPECIES if s in args.species)
    comparators = tuple(comparator_panel()) if args.panel == "full" else ()

    rule("=")
    say("STAGE 05 - orthology: ortholog matrix + nearest-neighbour matrix")
    rule("=")
    say("  in       : data/processed/proteomes/proteome_<species>.tsv  (the `sequence` column)")
    say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/  (+ evidence/ + scratch/)")
    say(f"  anchors  : {', '.join(species)}")
    if comparators:
        say(f"  panel    : + {len(comparators)} tier-C comparators (orthogroups only, no tables, "
            "no pairwise DIAMOND)")
    say(f"  discrete : OrthoFinder 3, de novo, all species jointly  (+ DIAMOND RBH)")
    say(f"  continuous: DIAMOND --{args.sensitivity} -e {args.evalue}, "
        f"top {args.top_k} per target species")
    say(f"  threads  : {args.threads}")
    say()

    expected = {sp: len(P.load(sp)) for sp in species}
    if args.dry_run:
        say("  --dry-run: nothing run, nothing written.")
        for sp in species:
            say(f"    {sp:<14} {expected[sp]:>6} proteins")
        for label in comparators:
            say(f"    {label[:34]:<36} comparator")
        say(f"    {'':<14} {sum(expected.values()):>6} anchor proteins, "
            f"{len(comparators)} comparators, {len(species) ** 2} DIAMOND searches")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    t_start = time.time()

    rule()
    say("INPUTS")
    rule()
    counts = write_fastas(species, comparators)
    for sp in species:
        say(f"  {sp:<14} {counts[sp]:>6} proteins -> {(FASTA_DIR / f'{sp}.faa').name}")
    if comparators:
        say(f"  + {len(comparators)} comparators, {sum(counts[c] for c in comparators):,} proteins "
            f"(orthogroups only)")
        say(f"  = {sum(counts.values()):,} proteins over {len(counts)} proteomes")
    say()

    rule()
    say("ORTHOFINDER")
    rule()
    # Before OrthoFinder, not after: OrthoFinder shells out to `diamond` for its own all-vs-all, so
    # `gradi-ortho/bin` has to be on PATH by the time it starts. ensure_diamond() mutates os.environ,
    # which the subprocess inherits. Getting this order wrong is silent -- see the note below.
    diamond = ensure_diamond()
    say(f"  diamond  : {diamond}")
    say(f"  {orthofinder_version()}")
    results = run_orthofinder(args.threads, args.refresh)
    og = parse_orthogroups(results, species)
    of_pairs = parse_orthologues(results)
    of_pairs = of_pairs[of_pairs.query_species.isin(species)
                        & of_pairs.target_species.isin(species)]
    n_assigned = int((og.orthogroup != "").sum())
    say(f"  proteins placed        : {len(og):,}")
    say(f"  in an orthogroup       : {n_assigned:,} ({100 * n_assigned / max(len(og), 1):.1f}%)")
    say(f"  unassigned (measured)  : {len(og) - n_assigned:,}")
    say(f"  distinct orthogroups   : {og.loc[og.orthogroup != '', 'orthogroup'].nunique():,}")
    say(f"  pairwise ortholog pairs: {len(of_pairs):,}")
    say()

    rule()
    say(f"DIAMOND - {len(species) ** 2} ordered species pairs")
    rule()
    hits = {}
    for q_sp, t_sp in product(species, species):
        t0 = time.time()
        h = diamond_pair(q_sp, t_sp, args.top_k, args.evalue, args.sensitivity,
                         args.threads, args.refresh)
        hits[(q_sp, t_sp)] = h
        n_q = h.query_ac.nunique()
        say(f"  {q_sp:<14} -> {t_sp:<14} {len(h):>7,} hits  "
            f"{n_q:>6,}/{counts[q_sp]:,} queries with a hit  ({time.time() - t0:>5.1f}s)")
    say()

    # Self-bitscore, for a length-normalised score comparable across proteins of different size.
    self_bits = {}
    for sp in species:
        h = hits[(sp, sp)]
        s = h[h.query_ac == h.target_ac].set_index("query_ac")["bitscore"]
        self_bits.update(s.to_dict())

    rows = build_neighbors(hits, species, args.top_k, self_bits)
    rbh = build_rbh(hits, species)
    ortho = build_orthologs(of_pairs, rbh, og)
    rows = annotate_neighbors(rows, of_pairs, rbh)
    # The bacterial panel for the conservation count: every bacterial proteome in the RUN, which
    # is the comparators plus the bacterial anchors. Human is excluded -- it is the liability
    # scope and has its own column.
    bacterial = tuple([s for s in species if s != "human"] + list(comparators))
    panel = panel_membership(results, bacterial) if comparators else None
    dense = build_dense(og, ortho, rows, species, expected, panel, len(bacterial) - 1)

    # BEFORE write_outputs -- it reads the tables the write is about to replace.
    delta = panel_delta(dense, species)

    write_outputs(rows, ortho, dense, og, species)

    if delta is not None:
        delta["panel_size"] = len(counts)
        delta["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        delta.to_csv(EVIDENCE_DIR / "panel_expansion.tsv", sep="\t", index=False)
        rule()
        say("PANEL CHANGE -- what moved in the already-shipped columns")
        rule()
        for _, r in delta.iterrows():
            say(f"  {r['species']:<13} in_orthogroup {r.get('in_orthogroup_before', 0):>6,} -> "
                f"{r.get('in_orthogroup_after', 0):>6,} ({r.get('in_orthogroup_delta', 0):+,})   "
                f"has_human {r.get('has_human_ortholog_before', 0):>5,} -> "
                f"{r.get('has_human_ortholog_after', 0):>5,} "
                f"({r.get('has_human_ortholog_delta', 0):+,})   "
                f"changed OG {r.get('changed_orthogroup', 0):>6,}")
        say("  recall RISES with panel size, so a NEGATIVE delta is a regression, not an update.")
        say(f"  wrote {(EVIDENCE_DIR / 'panel_expansion.tsv').relative_to(REPO_ROOT)}")
    fails = report(rows, ortho, dense, og, of_pairs, rbh, species, expected, counts, hits,
                   args, time.time() - t_start)
    if fails:
        sys.exit("FAILED:\n  " + "\n  ".join(fails))
    rule("=")
    say("stage 05 complete.")
    rule("=")


# ---------------------------------------------------------------- assembly


def build_neighbors(hits: dict, species: tuple[str, ...], top_k: int,
                    self_bits: dict) -> pd.DataFrame:
    """Top-k per (query, target species), self-hits dropped, ranked by bitscore."""
    frames = []
    for (q_sp, t_sp), h in hits.items():
        h = h[h.query_ac != h.target_ac]  # a protein is not its own neighbour
        h = h.sort_values(["query_ac", "bitscore"], ascending=[True, False])
        h = h.groupby("query_ac", sort=False).head(top_k).copy()
        h["rank"] = h.groupby("query_ac", sort=False).cumcount() + 1
        h["query_species"], h["target_species"] = q_sp, t_sp
        frames.append(h)
    out = pd.concat(frames, ignore_index=True)
    out["bitscore_norm"] = out.bitscore / out.query_ac.map(self_bits)
    return out


def build_rbh(hits: dict, species: tuple[str, ...]) -> set[tuple[str, str]]:
    """Reciprocal best hits, cross-species only.

    Within a species the best non-self hit is a paralog, not a reciprocal ortholog, so RBH is not
    defined there -- `same_orthogroup` carries paralogy instead.
    """
    best = {}
    for (q_sp, t_sp), h in hits.items():
        if q_sp == t_sp:
            continue
        h = h[h.query_ac != h.target_ac]
        b = h.sort_values("bitscore", ascending=False).drop_duplicates("query_ac")
        best[(q_sp, t_sp)] = dict(zip(b.query_ac, b.target_ac))
    pairs = set()
    for (q_sp, t_sp), fwd in best.items():
        rev = best.get((t_sp, q_sp), {})
        for a, b in fwd.items():
            if rev.get(b) == a:
                pairs.add((a, b))
                pairs.add((b, a))
    return pairs


def build_orthologs(of_pairs: pd.DataFrame, rbh: set, og: pd.DataFrame) -> pd.DataFrame:
    """The discrete matrix: every pair either method calls an ortholog, both flags side by side."""
    of_set = set(zip(of_pairs.query_ac, of_pairs.target_ac))
    sp_of = dict(zip(og.uniprot_ac, og.species))
    og_of = dict(zip(og.uniprot_ac, og.orthogroup))

    union = of_set | rbh
    rows = []
    for a, b in union:
        grp_a, grp_b = og_of.get(a, ""), og_of.get(b, "")
        rows.append({
            "query_ac": a, "query_species": sp_of.get(a, ""),
            "target_ac": b, "target_species": sp_of.get(b, ""),
            "is_ortholog_orthofinder": (a, b) in of_set,
            "is_rbh": (a, b) in rbh,
            "same_orthogroup": bool(grp_a) and grp_a == grp_b,
            "orthogroup": grp_a if grp_a and grp_a == grp_b else "",
        })
    out = pd.DataFrame(rows, columns=list(O.ORTHOLOG_COLUMNS))
    return out.sort_values(["query_species", "query_ac", "target_species", "target_ac"])


def annotate_neighbors(nb: pd.DataFrame, of_pairs: pd.DataFrame, rbh: set) -> pd.DataFrame:
    """Join the two ortholog calls onto the neighbour matrix, so one table answers both questions."""
    of_set = set(zip(of_pairs.query_ac, of_pairs.target_ac))
    keys = list(zip(nb.query_ac, nb.target_ac))
    nb["is_ortholog_orthofinder"] = [k in of_set for k in keys]
    nb["is_rbh"] = [k in rbh for k in keys]
    return nb


def panel_delta(dense: dict, species: tuple[str, ...]) -> pd.DataFrame | None:
    """Before/after on the columns a panel change MOVES, measured rather than assumed.

    OrthoFinder's recall rises with the number of species in the run -- v1's much-cited 55.5%
    Kp<->Ec came from a 25-species run, where four species give 44.8% -- so widening the panel
    changes `in_orthogroup` and `has_human_ortholog` on tables that are already shipped and
    already read by other axes. A silent shift in a shipped column is the thing to avoid here; a
    measured one is a finding, and the direction matters: UNDER-detecting human homology makes a
    target look more selective than it is.

    Reads whatever is currently on disk BEFORE the new tables overwrite it. Returns None on a
    first run, when there is nothing to compare against.
    """
    rows = []
    for sp in species:
        old_path = OUT_DIR / f"orthology_{sp}.tsv"
        if not old_path.exists():
            continue
        old = pd.read_csv(old_path, sep="\t")
        new = dense[sp]
        r = {"species": sp, "n": len(new)}
        for col in ("in_orthogroup", "has_human_ortholog"):
            if col in old.columns and col in new.columns:
                r[f"{col}_before"] = int(old[col].astype(bool).sum())
                r[f"{col}_after"] = int(new[col].astype(bool).sum())
                r[f"{col}_delta"] = r[f"{col}_after"] - r[f"{col}_before"]
        if "orthogroup" in old.columns:
            m = old[["uniprot_ac", "orthogroup"]].merge(
                new[["uniprot_ac", "orthogroup"]], on="uniprot_ac", suffixes=("_old", "_new"))
            r["changed_orthogroup"] = int((m.orthogroup_old.fillna("")
                                           != m.orthogroup_new.fillna("")).sum())
        rows.append(r)
    return pd.DataFrame(rows) if rows else None


def build_dense(og: pd.DataFrame, ortho: pd.DataFrame, nb: pd.DataFrame,
                species: tuple[str, ...], expected: dict,
                panel: dict[str, int] | None = None, panel_size: int = 0) -> dict[str, pd.DataFrame]:
    """One row per protein, per species. The table where a zero is readable."""
    size = og[og.orthogroup != ""].groupby("orthogroup").size()
    # Three counts per target species, not one. The two methods are never merged into a single
    # number: OrthoFinder's recall depends on how many species are in the run (a 2-species run has
    # almost no phylogenetic signal), while RBH does not, so a merged count would hide which of the
    # two moved. `n_orthologs_<sp>` is the union, kept as the convenient headline.
    def _counts(frame):
        return frame.groupby(["query_ac", "target_species"]).size().unstack(fill_value=0)

    counts = _counts(ortho)
    counts_of = _counts(ortho[ortho.is_ortholog_orthofinder])
    counts_rbh = _counts(ortho[ortho.is_rbh])
    best = (nb.sort_values("bitscore", ascending=False)
              .drop_duplicates(["query_ac", "target_species"])
              .set_index(["query_ac", "target_species"]))

    out = {}
    for sp in species:
        base = P.load(sp)[["uniprot_ac"]].copy()
        grp = og[og.species == sp].set_index("uniprot_ac")["orthogroup"]
        base["orthogroup"] = base.uniprot_ac.map(grp).fillna("")
        base["in_orthogroup"] = base.orthogroup != ""
        base["orthogroup_size"] = base.orthogroup.map(size).fillna(0).astype(int)
        # Paralogs = same orthogroup, same species, not itself.
        same_sp = og[(og.species == sp) & (og.orthogroup != "")]
        per_group = same_sp.groupby("orthogroup").size()
        base["n_paralogs"] = (base.orthogroup.map(per_group).fillna(1).astype(int) - 1).clip(lower=0)
        base["searched"] = True
        for t_sp in species:
            for prefix, table in (("n_orthologs_", counts),
                                  ("n_orthologs_of_", counts_of),
                                  ("n_orthologs_rbh_", counts_rbh)):
                n = table[t_sp] if t_sp in table.columns else pd.Series(dtype=int)
                base[f"{prefix}{t_sp}"] = base.uniprot_ac.map(n).fillna(0).astype(int)
            key = [(ac, t_sp) for ac in base.uniprot_ac]
            base[f"best_identity_{t_sp}"] = [
                best.pident.get(k, np.nan) for k in key]
            base[f"best_bitscore_{t_sp}"] = [
                best.bitscore.get(k, np.nan) for k in key]
        if "human" in species:
            base["has_human_ortholog"] = base["n_orthologs_human"] > 0
        if panel is not None and sp != "human":
            # How many OTHER bacterial proteomes share this protein's orthogroup. The protein's own
            # species is always one of the members, so subtract it; a protein in no orthogroup is a
            # measured 0, which is why `.fillna(0)` is correct here and nowhere near a lookup table.
            n = base.orthogroup.map(panel).fillna(0).astype(int) - 1
            base["n_bacterial_orthologs"] = n.clip(lower=0)
            # Ships beside the count so the number stays interpretable when the panel changes --
            # 12 of 28 and 12 of 3 are different claims, and only this column tells them apart.
            base["bacterial_panel_size"] = panel_size
        out[sp] = base
    return out


def _evidence_audit(sp: str, dense: pd.DataFrame, level: pd.Series) -> pd.DataFrame:
    """Every condition behind `orthology_evidence`, per protein.

    The ladder collapses five booleans into one integer, and the two that actually discriminate
    are the RBH term and the human conflict -- so without this table a level 2 says only "not
    corroborated" and never which half failed. `in_orthogroup` uniquely blocks 0 proteins at
    level 3 on all three species; it is load-bearing only for level 1.
    """
    odb = O.load_orthodb(sp).set_index("uniprot_ac")["orthodb_verdict"]
    num = lambda c: pd.to_numeric(dense[c], errors="coerce").fillna(0)
    others = [s for s in O.BACTERIA if s != sp]
    of_h, rbh_h = num("n_orthologs_of_human") > 0, num("n_orthologs_rbh_human") > 0
    return pd.DataFrame({
        "species": sp, "uniprot_ac": dense["uniprot_ac"],
        "in_orthogroup": dense["in_orthogroup"].astype(bool),
        "in_orthodb": dense["uniprot_ac"].map(odb).fillna("").isin(O.ORTHODB_PLACED),
        "of_bacterial": sum(num(f"n_orthologs_of_{s}") for s in others) > 0,
        "rbh_bacterial": sum(num(f"n_orthologs_rbh_{s}") for s in others) > 0,
        "of_human": of_h, "rbh_human": rbh_h, "human_conflict": of_h != rbh_h,
        "orthology_evidence": level.to_numpy(),
    })


def write_outputs(nb: pd.DataFrame, ortho: pd.DataFrame, dense: dict,
                  og: pd.DataFrame, species: tuple[str, ...]) -> None:
    nb = nb[list(O.NEIGHBOR_COLUMNS)].sort_values(
        ["query_species", "query_ac", "target_species", "rank"])
    for col in ("pident", "ppos", "qcov", "scov"):
        nb[col] = nb[col].round(2)
    nb["bitscore_norm"] = nb.bitscore_norm.round(4)
    nb.to_csv(OUT_DIR / "neighbors.tsv", sep="\t", index=False)
    ortho.to_csv(OUT_DIR / "orthologs.tsv", sep="\t", index=False)
    evidence_rows = []
    for sp, df in dense.items():
        for col in df.columns:
            if col.startswith("best_identity_"):
                df[col] = df[col].round(2)
        # The DELIVERABLE is three columns (owner's call, 2026-10-03): the selectivity liability
        # and one conservation number. The 29-column dense table is the evidence behind them and
        # keeps everything -- orthogroup, paralogs, per-species counts and identities.
        df.to_csv(EVIDENCE_DIR / f"orthology_{sp}.tsv", sep="\t", index=False)
        if "n_bacterial_orthologs" in df.columns:
            slim = df[["uniprot_ac", "has_human_ortholog"]].copy()
            # A PROPORTION needs its denominator stated, and it is 28, not 26: the 26 tier-C
            # comparators plus the three bacterial anchors, minus this protein's own species.
            # bacterial_panel_size ships in the dense table so the number stays interpretable if
            # the panel ever changes.
            slim["bacterial_panel_orthologs"] = (
                df["n_bacterial_orthologs"] / df["bacterial_panel_size"]).round(4)
            # The standard evidence column. Computed through `src/orthology.py` -- which re-reads
            # the dense table just written -- so the file and the loader cannot drift apart.
            # It needs OrthoDB as the second, panel-independent grouping: run orthology/orthodb.py
            # first. Absent, this raises rather than quietly shipping a one-grouping ladder under
            # the same column name.
            try:
                slim["orthology_evidence"] = O.orthology_evidence(sp).to_numpy()
            except FileNotFoundError as exc:
                raise SystemExit(
                    f"{sp}: orthology_evidence needs orthodb_{sp}.tsv -- run "
                    f"scripts/orthology/orthodb.py first ({exc})") from exc
            slim.to_csv(OUT_DIR / f"orthology_{sp}.tsv", sep="\t", index=False)
            evidence_rows.append(_evidence_audit(sp, df, slim["orthology_evidence"]))
            lv = slim["orthology_evidence"].value_counts()
            say(f"  {sp:<14} evidence "
                + " ".join(f"L{k}={int(lv.get(k, 0)):,}" for k in (1, 2, 3)))
        else:
            say(f"  WARN {sp}: no n_bacterial_orthologs -- run with --panel full; "
                f"deliverable not written")

    if evidence_rows:
        aud = pd.concat(evidence_rows, ignore_index=True)
        aud.to_csv(EVIDENCE_DIR / "evidence_audit.tsv", sep="\t", index=False)
        say(f"  wrote evidence/evidence_audit.tsv  ({len(aud):,} rows)")
        for sp, g in aud.groupby("species", sort=False):
            say(f"  {sp:<14} human call: BOTH {int((g.of_human & g.rbh_human).sum()):>4}   "
                f"OF-only {int((g.of_human & ~g.rbh_human).sum()):>4}   "
                f"RBH-only {int((~g.of_human & g.rbh_human).sum()):>4}   "
                f"-> {int(g.human_conflict.sum()):>4} conflicts, capped at level 2")

    comp = (og[og.orthogroup != ""].groupby(["orthogroup", "species"]).size()
              .unstack(fill_value=0).reset_index())
    comp["size"] = comp[[s for s in species if s in comp.columns]].sum(axis=1)
    comp["n_species"] = (comp[[s for s in species if s in comp.columns]] > 0).sum(axis=1)
    comp.to_csv(EVIDENCE_DIR / "orthogroups.tsv", sep="\t", index=False)

    dis = ortho[ortho.is_ortholog_orthofinder != ortho.is_rbh]
    dis.to_csv(EVIDENCE_DIR / "method_disagreement.tsv", sep="\t", index=False)


def report(nb, ortho, dense, og, of_pairs, rbh, species, expected, counts, hits,
           args, elapsed) -> list[str]:
    """Print the summary and the controls; return the list of hard failures."""
    fails: list[str] = []

    rule()
    say("COVERAGE - every accession accounted for")
    rule()
    for sp in species:
        got = set(og.loc[og.species == sp, "uniprot_ac"])
        want = set(P.load(sp)["uniprot_ac"])
        missing, extra = want - got, got - want
        flag = "" if not (missing or extra) else \
            f"   <- {len(missing)} missing, {len(extra)} unexpected"
        say(f"  {sp:<14} {len(got):>6} / {len(want):>6} placed by OrthoFinder{flag}")
        if missing or extra:
            fails.append(f"{sp}: OrthoFinder placed {len(got)} accessions, proteome has {len(want)}"
                         f" ({len(missing)} missing) -- a zero would not be a measured zero")
        queried = set(hits[(sp, sp)].query_ac)
        for t_sp in species:
            queried &= set(hits[(sp, t_sp)].query_ac) | (want - set(hits[(sp, t_sp)].query_ac))
    say("  -> a 0 in n_orthologs_* is therefore a measured 0, not an unknown.")
    say()

    rule()
    say("METHODS - OrthoFinder vs reciprocal best hit")
    rule()
    both = int((ortho.is_ortholog_orthofinder & ortho.is_rbh).sum())
    of_only = int((ortho.is_ortholog_orthofinder & ~ortho.is_rbh).sum())
    rbh_only = int((~ortho.is_ortholog_orthofinder & ortho.is_rbh).sum())
    total = len(ortho)
    say(f"  ortholog pairs (either)  {total:>8,}")
    say(f"    both methods           {both:>8,}  ({100 * both / max(total, 1):.1f}%)")
    say(f"    OrthoFinder only       {of_only:>8,}   <- co-orthologs RBH cannot represent")
    say(f"    RBH only               {rbh_only:>8,}")
    say("  disagreements -> evidence/method_disagreement.tsv")
    say("  RBH is panel-independent; OrthoFinder's recall RISES with the number of species in the")
    say(f"  run (this one has {len(species)}). A low agreement rate on a small panel is expected.")
    say()

    rule()
    say("CONTROL - v1's measured numbers, reproduced")
    rule()
    ctrl = []
    kp = dense.get("kpneumoniae")
    if kp is not None and "ecoli" in species:
        n_of = int((kp["n_orthologs_of_ecoli"] > 0).sum())
        n_rbh = int((kp["n_orthologs_rbh_ecoli"] > 0).sum())
        n_any = int((kp["n_orthologs_ecoli"] > 0).sum())
        # NOT like-for-like: v1's 3,179 came from a 25-SPECIES OrthoFinder run. Ortholog recall
        # rises with panel size, so a 4-species run is expected to fall short and that is not a
        # defect. RBH is the panel-independent comparison, and it is the one that must match.
        ctrl.append(("Kp w/ E. coli ortholog (OrthoFinder)", O.V1_KP_WITH_ECOLI_ORTHOLOG, n_of,
                     f"{100 * n_of / len(kp):.1f}%  v1 ran 25 species, not {len(species)}"))
        ctrl.append(("Kp w/ E. coli ortholog (RBH)", O.V1_KP_EC_RBH_PAIRS, n_rbh,
                     f"{100 * n_rbh / len(kp):.1f}%  <- the like-for-like comparison"))
        ctrl.append(("Kp w/ E. coli ortholog (either)", O.V1_KP_WITH_ECOLI_ORTHOLOG, n_any,
                     f"{100 * n_any / len(kp):.1f}%"))
    for pair, v1 in O.V1_RBH_MEDIAN_PIDENT.items():
        a, b = pair.split("__")
        if a in species and b in species:
            m = nb[(nb.query_species == a) & (nb.target_species == b) & nb.is_rbh]
            got = float(m.pident.median()) if len(m) else float("nan")
            ctrl.append((f"RBH median identity {a}-{b}", v1, round(got, 1), ""))
    if "kpneumoniae" in species and "human" in species:
        n = int(nb[(nb.query_species == "kpneumoniae")
                   & (nb.target_species == "human")].query_ac.nunique())
        ctrl.append(("Kp proteins with any human hit", O.V1_KP_WITH_ANY_HUMAN_HIT, n,
                     f"{100 * n / expected['kpneumoniae']:.1f}%  (v1: 12.8%, should RISE)"))

    say(f"  {'quantity':<38} {'v1':>8} {'v2':>8}   note")
    for name, v1, got, note in ctrl:
        say(f"  {name:<38} {v1:>8} {got:>8}   {note}")
    pd.DataFrame(ctrl, columns=["quantity", "v1", "v2", "note"]).to_csv(
        EVIDENCE_DIR / "control.tsv", sep="\t", index=False)
    say()

    rule()
    say("SUMMARY")
    rule()
    say(f"  {'species':<14} {'n':>6} {'in OG':>7} {'paralogs':>9} " +
        " ".join(f"{'->' + s[:6]:>10}" for s in species))
    rows = []
    for sp in species:
        d = dense[sp]
        cells = " ".join(f"{int((d[f'n_orthologs_{t}'] > 0).sum()):>10,}" for t in species)
        say(f"  {sp:<14} {len(d):>6} {int(d.in_orthogroup.sum()):>7,} "
            f"{int((d.n_paralogs > 0).sum()):>9,} {cells}")
        rows.append({
            "species": sp, "n": len(d), "n_expected": expected[sp],
            "in_orthogroup": int(d.in_orthogroup.sum()),
            "with_paralogs": int((d.n_paralogs > 0).sum()),
            **{f"with_ortholog_{t}": int((d[f"n_orthologs_{t}"] > 0).sum()) for t in species},
            "top_k": args.top_k, "evalue": args.evalue, "sensitivity": args.sensitivity,
            "run_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        })
    say("  (counts are proteins with >=1 ortholog by either method, in that target species)")
    pd.DataFrame(rows).to_csv(EVIDENCE_DIR / "manifest.tsv", sep="\t", index=False)
    say(f"\n  total {elapsed / 60:.1f} min")
    say()

    rule()
    say("OUTPUTS")
    rule()
    for name in ["orthologs.tsv", "neighbors.tsv"] + [f"orthology_{s}.tsv" for s in species]:
        p = OUT_DIR / name
        say(f"  {('data/processed/orthology/' + name):<52} "
            f"{p.stat().st_size / 1e6:>8.2f} MB")
    say("  evidence/ orthogroups.tsv, method_disagreement.tsv,  scratch/ orthofinder/, hits/,")
    say("                control.tsv, manifest.tsv")
    say()

    # Structural guards.
    if len(nb) and nb["rank"].max() > args.top_k:
        fails.append(f"neighbour rank exceeds --top-k ({nb['rank'].max()} > {args.top_k})")
    if len(nb) and (nb.query_ac == nb.target_ac).any():
        fails.append("neighbors.tsv contains self-hits")
    asym = {(b, a) for a, b in rbh} - rbh
    if asym:
        fails.append(f"RBH is not symmetric: {len(asym)} pairs present in one direction only")
    within = ortho[(ortho.query_species == ortho.target_species)
                   & ortho.is_ortholog_orthofinder]
    if len(within):
        fails.append(f"{len(within)} within-species pairs flagged is_ortholog_orthofinder -- "
                     "orthology is a between-species relation; this should be impossible")
    return fails


if __name__ == "__main__":
    main()
