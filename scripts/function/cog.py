"""Functional annotation -- broad functional classes per protein, from sequence.

Stage 02. Reads the stage-00 tables and assigns every bacterial protein a single COG functional
category letter: the coarse, controlled vocabulary the rest of the pipeline can filter and plot on.

COG is the first classification this stage produces and currently the only one, which is why the
stage is `02_function` while the artifact is `cog_<species>.tsv`. Anything else added here --
another scheme, another vocabulary -- gets its own `<scheme>_<species>.tsv` beside it and its own
`<scheme>_*` files in `evidence/ + scratch/`. Nothing in this stage is named after COG except COG's own files.

Why COG
-------
The 26 COG functional categories, grouped into four super-groups, are the standard broad functional
vocabulary in bacterial genomics. One letter per protein. That is the whole point -- InterPro gives
thousands of domains and PANTHER thousands of families, neither of which is *broad*.

Why not lift it from UniProt
----------------------------
Because the anchor organism does not have it. Measured on the stage-00 annotation layer:

    eggnog xref coverage:   Kp HS11286 0.00% (n=0)  |  Ec K-12 92.98%  |  Sa NCTC 8325 79.23%

Zero on K. pneumoniae. HS11286 is a dark TrEMBL proteome, so every accession-based route to COG
fails on exactly the species that matters most. A sequence-based one does not. This is the same
finding that forces DIAMOND-by-sequence everywhere else in this project.

v1 never computed COG at all. What it had was `functional_class` in `08a_webapp_export.py` -- a
hand-written substring keyword map over InterPro/PANTHER names, 11 ordered first-match-wins rules,
dumping 24.6% of Kp into `other`, and materialised only inside the webapp JSON. eggNOG-mapper was
identified as the fix in `legacy/docs/01_task_agnostic.md:240` and never built. This stage is that
gap closed, with a lighter tool.

The tool
--------
COGclassifier 2.0.0 (MIT). Protein FASTA in, category letter out:

    RPS-BLAST vs the CDD COG profile DB -> best hit -> CDD ID -> COG ID -> letter -> group

on **COG2024** (5,049 COGs; `cog_definition.tsv` and `cog_func_category.tsv` ship inside the
package). Downloads are ~199 MB (`Cog_LE.tar.gz` 191 MB + `cddid.tbl.gz` 7.4 MB), cached under
`data/source/cdd/`. eggNOG-mapper would have been ~20 GB unpacked for the same letter.

We use the library API, never its CLI or its Altair charts -- this repo plots with stylia. And we
drive `RpsBlast` directly rather than `CogClassifier.run()`, because the latter writes its hit table
to a tempdir; keeping it is what makes the run resumable.

rpsblast
--------
Not in `gradi`. It is borrowed from `gradi-prokka` (`rpsblast 2.17.0+`, x86_64 under Rosetta) by
prepending that env's `bin` to PATH -- the same pattern v1 used for `GRADI_DIAMOND_BIN`. Override
with `GRADI_RPSBLAST_BIN=<dir>`. Do NOT try to conda-install blast into `gradi`: there is no
osx-arm64 build, so it would drag the whole environment to osx-64.

Scope
-----
The three bacteria only -- 13,020 proteins. Human is out and cannot be otherwise: COG2024 is 2,296
genomes of **bacteria and archaea, zero eukaryotes**, so human would hit only the conserved core and
read as misleadingly sparse. Host/off-target comparison is orthology, a later stage's job.

Output
------
    data/processed/function/cog_<species>.tsv        one row per protein, keyed on uniprot_ac
        uniprot_ac cog_id cog_category cog_category_all cog_group cog_name cog_evalue cog_identity

    scratch/cog_rpsblast_<species>.tsv   raw outfmt-6 hits (the resumable cache)
    evidence/cog_counts_<species>.tsv     per-letter counts
    evidence/cog_control_ecoli.tsv        the ground-truth comparison, row by row
    evidence/cog_manifest.tsv

`cog_category` is exactly one letter. When a COG carries several (COG4862 is `KTN`), COG orders them
by importance and we take the first -- but the full string is kept in `cog_category_all` and the
count is reported, because this stage does not resolve ambiguity silently.

The E. coli control
-------------------
E. coli K-12 MG1655 (GCF_000005845.2) is itself one of the COG2024 reference genomes, so NCBI
publishes curated COG assignments for our exact anchor proteome. `--control` streams those rows out
of the 637 MB `cog-24.cog.csv` without landing the file, joins them on RefSeq WP_ accession and
reports agreement. Same spirit as stage 00's E. coli naming control: the one species where we can
check the machinery against a curated answer.

Run with the `gradi` env:
    python scripts/function/cog.py                                  # the three bacteria + control
    python scripts/function/cog.py --species saureus --limit 50     # smoke test
    python scripts/function/cog.py --refresh                        # re-run RPS-BLAST from scratch
"""

from __future__ import annotations

import argparse
import hashlib
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
from src import proteomes as P  # noqa: E402

RAW_DIR = REPO_ROOT / "data" / "source" / "cdd"
OUT_DIR = REPO_ROOT / "data" / "processed" / "function"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
DEFAULT_SPECIES = ("kpneumoniae", "ecoli", "saureus")
DEFAULT_EVALUE = 1e-2

# The one species in COG2024 that is also one of our anchors, and the URL of NCBI's own
# per-protein assignments. 637 MB, so it is streamed and filtered, never downloaded.
CONTROL_SPECIES = "ecoli"
CONTROL_ASSEMBLY = "GCF_000005845.2"          # E. coli K-12 MG1655
COG24_COG_CSV = "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/cog-24.cog.csv"
CONTROL_MIN_LETTER_AGREEMENT = 0.90

# rpsblast is not in `gradi`; borrow it from gradi-prokka unless told otherwise.
DEFAULT_RPSBLAST_DIR = Path.home() / "miniconda3" / "envs" / "gradi-prokka" / "bin"

OUT_COLUMNS = [
    "uniprot_ac", "cog_id", "cog_category", "cog_category_all",
    "cog_group", "cog_name", "cog_evalue", "cog_identity",
]

# Proteins whose category is not a matter of opinion. Checked after every run: if a Clp protease
# subunit is not `O` the CDD -> COG -> letter chain is broken, whatever the coverage table says.
SPOT_CHECKS = {
    "kpneumoniae": {"clpP": "O", "clpX": "O", "clpA": "O", "clpB": "O", "rpsA": "J", "rplB": "J"},
    "ecoli":       {"clpP": "O", "clpX": "O", "clpA": "O", "clpB": "O", "rpsA": "J", "rplB": "J"},
    "saureus":     {"clpP": "O", "clpX": "O", "clpC": "O", "clpB": "O", "rpsA": "J", "rplB": "J"},
}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------- rpsblast

def ensure_rpsblast() -> str:
    """Put rpsblast on PATH and return its version, or exit with the fix."""
    override = os.environ.get("GRADI_RPSBLAST_BIN")
    candidates = [Path(override)] if override else [DEFAULT_RPSBLAST_DIR]
    if shutil.which("rpsblast") is None:
        for d in candidates:
            if (d / "rpsblast").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("rpsblast")
    if exe is None:
        sys.exit(
            "rpsblast not found.\n"
            f"  looked on PATH and in {candidates[0]}\n"
            "  it lives in the `gradi-prokka` env; point GRADI_RPSBLAST_BIN at a directory\n"
            "  containing it, or `conda install -c bioconda blast` into a SEPARATE env --\n"
            "  installing blast into `gradi` would drag the whole env to osx-64."
        )
    out = subprocess.run([exe, "-version"], capture_output=True, text=True)
    version = (out.stdout or out.stderr).splitlines()[0].replace("rpsblast:", "").strip()
    return f"{exe}  (v{version})"


# ---------------- COG resources

def fetch_resources(refresh: bool) -> tuple[Path, Path]:
    """Download and unpack the CDD COG profile DB + the CDD->COG id table. Returns (db, cddid)."""
    from cogclassifier import const, utils

    RAW_DIR.mkdir(parents=True, exist_ok=True)
    cddid = utils.ftp_download(const.CDDID_TBL_FTP, RAW_DIR, overwrite=refresh)
    say(f"  cddid.tbl.gz        {cddid.stat().st_size / 1e6:>7.1f} MB")

    targz = utils.ftp_download(const.COG_LE_FTP, RAW_DIR, overwrite=refresh)
    say(f"  Cog_LE.tar.gz       {targz.stat().st_size / 1e6:>7.1f} MB")

    cog_le_dir = RAW_DIR / "Cog_LE"
    if refresh and cog_le_dir.exists():
        shutil.rmtree(cog_le_dir)
    if not cog_le_dir.exists():
        say(f"  unpacking -> {cog_le_dir.relative_to(REPO_ROOT)}/")
        shutil.unpack_archive(targz, cog_le_dir)

    db = cog_le_dir / "Cog"
    # The CDD COG set ships volume-split (Cog.00.*, Cog.01.*) behind a `Cog.pal` alias; a small
    # single-volume build would instead land a `Cog.pin`. `-db .../Cog` resolves either.
    if not ((cog_le_dir / "Cog.pal").exists() or (cog_le_dir / "Cog.pin").exists()):
        sys.exit(f"unpacked {cog_le_dir} but found neither Cog.pal nor Cog.pin -- corrupt download?")
    write_source_md(cddid, targz)
    return db, cddid


def write_source_md(cddid: Path, targz: Path) -> None:
    """Provenance for the two downloads, in the stage-00 SOURCE.md style."""
    from cogclassifier import const

    (RAW_DIR / "SOURCE.md").write_text(
        "# COG / CDD resources for stage 02\n\n"
        f"Fetched {datetime.now(timezone.utc).isoformat(timespec='seconds')} by scripts/function/cog.py\n\n"
        "| file | bytes | url |\n|---|---|---|\n"
        f"| `{cddid.name}` | {cddid.stat().st_size:,} | {const.CDDID_TBL_FTP} |\n"
        f"| `{targz.name}` | {targz.stat().st_size:,} | {const.COG_LE_FTP} |\n"
        f"| `Cog_LE/` | (unpacked) | RPS-BLAST profile DB, `-db Cog_LE/Cog` |\n\n"
        "The COG2024 definition and functional-category tables are NOT downloaded -- they ship\n"
        f"inside cogclassifier {getattr(__import__('cogclassifier'), '__version__', '?')}:\n\n"
        f"- `{const.COG_DEFINITION_FILE.name}` (cog-24.def.tab)\n"
        f"- `{const.COG_FUNC_CATEGORY_FILE.name}` (cog-24.fun.tab)\n\n"
        "Both originate from <https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/>.\n",
        encoding="utf-8",
    )


def load_cog_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """COG2024 definitions and the 26 functional categories, from the bundled package resources.

    Parsed through the package's own readers, not `pd.read_csv`: `cog_definition.tsv` is ragged
    (4 to 10 tab-separated fields per line, 5,043 of 5,049 rows at 7), which the C parser rejects.
    """
    from cogclassifier import const
    from cogclassifier.cog import CogDefinitionRecord, CogFuncCategoryRecord

    defs = pd.DataFrame(
        [(d.id, d.letter, d.cog_name) for d in CogDefinitionRecord(const.COG_DEFINITION_FILE).get_all()],
        columns=["cog_id", "cog_category_all", "cog_name"],
    )
    cats = pd.DataFrame(
        [(c.letter, c.group, c.desc) for c in CogFuncCategoryRecord(const.COG_FUNC_CATEGORY_FILE).get_all()],
        columns=["cog_category", "cog_group", "cog_category_desc"],
    )
    return defs, cats


# ---------------- per species

def cache_key(query: Path, evalue: float) -> str:
    """Identity of a hit table: which accessions were asked, at which e-value."""
    accs = sorted(ln[1:].strip() for ln in query.read_text().splitlines() if ln.startswith(">"))
    return hashlib.md5(("\n".join(accs) + f"|evalue={evalue}").encode()).hexdigest()


def write_query_fasta(species: str, limit: int | None, path: Path) -> int:
    """Write the query FASTA from the stage-00 table, headers = bare accessions.

    Deliberately NOT data/source/uniprot/proteomes/*.fasta: those carry `sp|A5A616|MGTS_ECOLI` headers, so
    every hit would need parsing back to an accession. Building it here makes QUERY_ID *be* the
    join key, and guarantees the query set is exactly the rows in the deliverable.
    """
    df = P.load(species)[["uniprot_ac", "sequence"]]
    if limit:
        df = df.head(limit)
    with path.open("w", encoding="utf-8") as f:
        for ac, seq in zip(df["uniprot_ac"], df["sequence"]):
            f.write(f">{ac}\n{seq}\n")
    return len(df)


def run_species(species, db, cddid, defs, cats, evalue, threads, limit, refresh) -> dict:
    from cogclassifier.blast import BlastAlignmentRecord, RpsBlast
    from cogclassifier.cog import CogCddIdTable, CogClassifyStats, CogDefinitionRecord, \
        CogFuncCategoryRecord
    from cogclassifier import const

    t0 = time.time()
    # A smoke test must never clobber a full run: with --limit everything is written under
    # scratch/smoke_*, so `cog_<species>.tsv` and the real hit cache are left alone.
    pre = "smoke_" if limit else ""
    hits_path = SCRATCH_DIR / f"{pre}cog_rpsblast_{species}.tsv"
    key_path = SCRATCH_DIR / f".{pre}cog_rpsblast_{species}.key"

    with tempfile.TemporaryDirectory() as tmp:
        query = Path(tmp) / f"{species}.faa"
        n_query = write_query_fasta(species, limit, query)
        say(f"  {species:<14} {n_query:>6} proteins")

        # Cache the hit table -- RPS-BLAST minutes are exactly what you want to keep across runs.
        # Keyed on (accession set, evalue), the stage-00 pattern: the file itself cannot tell you
        # what it covers, because a query with no hit leaves no row. A 50-protein smoke-test cache
        # silently serving a 2,889-protein run is precisely the v1 short-payload failure mode.
        key = cache_key(query, evalue)
        cached = False
        if hits_path.exists() and key_path.exists() and not refresh \
                and key_path.read_text().strip() == key:
            cached = True
            n_hit = len({ln.split("\t", 1)[0] for ln in
                         hits_path.read_text(encoding="utf-8").splitlines() if ln.strip()})
            say(f"  {'':<14} [cache] {hits_path.name}  {n_hit:,} queries with a hit")
        elif hits_path.exists() and not refresh:
            say(f"  {'':<14} [rerun] cached hits were built for a different query set or evalue")

        if cached:
            blast_rec = BlastAlignmentRecord(hits_path)
        else:
            say(f"  {'':<14} RPS-BLAST  evalue={evalue}  threads={threads} ...")
            blast_rec = RpsBlast(query, db, outfile=hits_path, evalue=evalue,
                                 thread_num=threads).run()
            key_path.write_text(key)

        stats = CogClassifyStats(
            query, blast_rec,
            CogFuncCategoryRecord(const.COG_FUNC_CATEGORY_FILE),
            CogDefinitionRecord(const.COG_DEFINITION_FILE),
            CogCddIdTable(cddid),
        )
        hits = stats.query_classify_df

    # Left-join onto the full proteome so unclassified proteins stay as rows, not silent absences.
    base = P.load(species)[["uniprot_ac"]]
    if limit:
        base = base.head(limit)

    hits = (hits[["QUERY_ID", "COG_ID", "EVALUE", "IDENTITY"]]
            .rename(columns={"QUERY_ID": "uniprot_ac", "COG_ID": "cog_id",
                             "EVALUE": "cog_evalue", "IDENTITY": "cog_identity"}))
    out = (base.merge(hits, on="uniprot_ac", how="left")
               .merge(defs, on="cog_id", how="left"))
    # cog_category is the FIRST letter; cog_category_all keeps the whole string. COG orders the
    # letters by importance, so "first" is COG's own choice, not ours -- but it is still a choice.
    out["cog_category"] = out["cog_category_all"].fillna("").str[:1]
    out = out.merge(cats[["cog_category", "cog_group"]], on="cog_category", how="left")
    out = out.fillna("")[OUT_COLUMNS]

    out_path = (SCRATCH_DIR / f"smoke_cog_{species}.tsv") if limit else (EVIDENCE_DIR / f"cog_{species}.tsv")
    out.to_csv(out_path, sep="\t", index=False)

    counts = (cats.merge(out["cog_category"].value_counts().rename("n"),
                         left_on="cog_category", right_index=True, how="left")
                  .fillna({"n": 0}))
    counts["n"] = counts["n"].astype(int)
    counts.insert(0, "species", species)
    counts.to_csv((SCRATCH_DIR if pre else EVIDENCE_DIR) / f"{pre}cog_counts_{species}.tsv", sep="\t", index=False)

    n = len(out)
    n_class = int((out["cog_category"] != "").sum())
    n_multi = int((out["cog_category_all"].str.len() > 1).sum())
    n_poor = int(out["cog_category"].isin(["R", "S"]).sum())
    return {
        "species": species, "n": n, "classified": n_class,
        "pct_classified": round(100 * n_class / n, 2),
        "informative": n_class - n_poor,
        "pct_informative": round(100 * (n_class - n_poor) / n, 2),
        "poorly_characterized_RS": n_poor, "multi_letter": n_multi,
        "unclassified": n - n_class, "seconds": round(time.time() - t0, 1),
        "cached": cached, "path": str(out_path.relative_to(REPO_ROOT)),
    }


# ---------------- the E. coli control

def run_control(defs: pd.DataFrame) -> dict | None:
    """Compare our E. coli calls against NCBI's own curated COG2024 assignments.

    E. coli K-12 MG1655 is one of the 2,296 COG2024 reference genomes, so this is a real
    ground truth for our exact anchor proteome -- not a second prediction.
    """
    ours_path = EVIDENCE_DIR / f"cog_{CONTROL_SPECIES}.tsv"
    if not ours_path.exists():
        say("  skipped: E. coli was not in this run")
        return None

    say(f"  streaming {COG24_COG_CSV}")
    say(f"  keeping only rows for {CONTROL_ASSEMBLY} (E. coli K-12 MG1655) -- 637 MB not landed")
    proc = subprocess.run(
        f"curl -sSf {COG24_COG_CSV} | grep ',{CONTROL_ASSEMBLY},' || true",
        shell=True, capture_output=True, text=True,
    )
    lines = [ln for ln in proc.stdout.splitlines() if ln.strip()]
    if not lines:
        say("  WARNING: no rows returned -- network problem or NCBI layout change. Control SKIPPED.")
        return None

    # Column 0 is the locus tag -- for MG1655 that is the b-number, which is exactly what stage 00
    # stores in locus_tags_ecoli.tsv, at 100% coverage. Column 2 is a RefSeq accession, but for this
    # genome NCBI uses NP_ while UniProt cross-references mostly WP_, so joining on it loses rows.
    # "Join on locus_tag, not gene name" (CLAUDE.md) turns out to apply to RefSeq too.
    ref = pd.DataFrame([ln.split(",") for ln in lines]).iloc[:, [0, 6]]
    ref.columns = ["locus_tag", "ref_cog_id"]
    ref = ref.drop_duplicates("locus_tag")
    ref = ref.merge(defs.rename(columns={"cog_id": "ref_cog_id",
                                         "cog_category_all": "ref_cat_all"})[
                        ["ref_cog_id", "ref_cat_all"]], on="ref_cog_id", how="left")
    ref["ref_category"] = ref["ref_cat_all"].fillna("").str[:1]
    say(f"  NCBI curated assignments for this genome: {len(ref):,} proteins")

    ours = pd.read_csv(ours_path, sep="\t", dtype=str, keep_default_na=False)
    tags = P.load_locus_tags(CONTROL_SPECIES)[["uniprot_ac", "locus_tag"]]
    merged = (ours.merge(tags, on="uniprot_ac", how="left")
                  .merge(ref, on="locus_tag", how="inner"))
    both = merged[(merged["cog_id"] != "") & (merged["ref_cog_id"] != "")]
    if both.empty:
        say("  WARNING: zero overlapping proteins -- control SKIPPED.")
        return None

    cog_agree = float((both["cog_id"] == both["ref_cog_id"]).mean())
    let_agree = float((both["cog_category"] == both["ref_category"]).mean())

    both.assign(cog_id_agrees=both["cog_id"] == both["ref_cog_id"],
                cog_category_agrees=both["cog_category"] == both["ref_category"]) \
        .to_csv(EVIDENCE_DIR / "cog_control_ecoli.tsv", sep="\t", index=False)

    say(f"  comparable proteins (both sides assigned): {len(both):,}")
    say(f"  COG id agreement       : {100 * cog_agree:6.2f}%")
    say(f"  COG category agreement : {100 * let_agree:6.2f}%   "
        f"(floor {100 * CONTROL_MIN_LETTER_AGREEMENT:.0f}%)")
    say("  row-by-row -> evidence/cog_control_ecoli.tsv")
    return {"n": len(both), "cog_agreement": cog_agree, "letter_agreement": let_agree}


def run_spot_checks(species_list: list[str], limit: int | None = None) -> list[str]:
    """Proteins whose category is not a matter of opinion. Returns a list of failures."""
    failures = []
    for sp in species_list:
        path = (SCRATCH_DIR / f"smoke_cog_{sp}.tsv") if limit else (EVIDENCE_DIR / f"cog_{sp}.tsv")
        if not path.exists():
            continue
        cog = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        names = P.load(sp)[["uniprot_ac", "gene_name"]]
        m = cog.merge(names, on="uniprot_ac", how="left")
        for gene, expect in SPOT_CHECKS.get(sp, {}).items():
            got = m.loc[m["gene_name"] == gene, "cog_category"].tolist()
            if not got:
                say(f"  {sp:<14} {gene:<6} not present under that gene name  (skipped)")
                continue
            ok = all(g == expect for g in got)
            say(f"  {sp:<14} {gene:<6} -> {','.join(got):<6} expected {expect}   "
                f"{'ok' if ok else '<- FAIL'}")
            if not ok:
                failures.append(f"{sp}/{gene}: got {got}, expected {expect}")
    return failures


# ---------------- main

def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES),
                    choices=list(DEFAULT_SPECIES),
                    help="species to classify (human is not possible: COG2024 has no eukaryotes)")
    ap.add_argument("--evalue", type=float, default=DEFAULT_EVALUE)
    ap.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, help="only the first N proteins per species (smoke test)")
    ap.add_argument("--refresh", action="store_true", help="re-download and re-run RPS-BLAST")
    ap.add_argument("--no-control", action="store_true", help="skip the E. coli ground-truth check")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and stop")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    import cogclassifier

    rule("=")
    say("STAGE 02 - functional annotation :: COG functional categories")
    rule("=")
    say("  in       : data/processed/proteomes/proteome_<species>.tsv  (uniprot_ac + sequence)")
    say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/cog_<species>.tsv  (+ evidence/ + scratch/)")
    say(f"  method   : RPS-BLAST vs CDD COG profiles -> COG2024 category "
        f"(cogclassifier {cogclassifier.__version__})")
    say(f"  species  : {', '.join(args.species)}")
    say(f"  evalue   : {args.evalue}    threads: {args.threads}")
    if args.limit:
        say(f"  limit    : {args.limit} proteins per species (SMOKE TEST)")
    if args.dry_run:
        say("\n  --dry-run: nothing fetched, nothing written.")
        for sp in args.species:
            say(f"    {sp:<14} -> cog_{sp}.tsv")
        return
    say()

    rule()
    say("RESOURCES")
    rule()
    say(f"  rpsblast : {ensure_rpsblast()}")
    db, cddid = fetch_resources(args.refresh)
    defs, cats = load_cog_tables()
    say(f"  COG2024  : {len(defs):,} COGs, {len(cats)} functional categories, "
        f"{cats['cog_group'].nunique()} groups")
    say()

    rule()
    say("CLASSIFY")
    rule()
    rows = [run_species(sp, db, cddid, defs, cats, args.evalue, args.threads,
                        args.limit, args.refresh) for sp in args.species]
    pd.DataFrame(rows).to_csv(
        (SCRATCH_DIR if args.limit else EVIDENCE_DIR) / ("smoke_cog_manifest.tsv" if args.limit else "cog_manifest.tsv"),
        sep="\t", index=False)
    say()

    rule()
    say("COVERAGE")
    rule()
    say(f"  {'species':<14} {'n':>6} {'classified':>18} {'informative':>18} "
        f"{'R+S':>10} {'multi-letter':>13}")
    for r in rows:
        say(f"  {r['species']:<14} {r['n']:>6} "
            f"{r['classified']:>9,} ({r['pct_classified']:>5.1f}%) "
            f"{r['informative']:>9,} ({r['pct_informative']:>5.1f}%) "
            f"{r['poorly_characterized_RS']:>10,} {r['multi_letter']:>13,}")
    say("\n  'informative' excludes R (general function prediction only) and S (function unknown):")
    say("  those are classified, but they are not an answer.")
    say()

    rule()
    say("SPOT CHECKS")
    rule()
    failures = run_spot_checks(args.species, args.limit)
    say()

    control = None
    if args.limit and not args.no_control:
        say("CONTROL skipped: --limit truncates the proteome, so agreement would be meaningless.\n")
    elif not args.no_control:
        rule()
        say("CONTROL - E. coli K-12 MG1655 vs NCBI's own COG2024 assignments")
        rule()
        control = run_control(defs)
        say()

    rule()
    say("OUTPUTS")
    rule()
    for r in rows:
        say(f"  {r['path']:<52} {(REPO_ROOT / r['path']).stat().st_size / 1e3:>8.1f} kB")
    say(f"  evidence/ + scratch/    cog_rpsblast_<species>.tsv, cog_counts_<species>.tsv, "
        f"cog_manifest.tsv{', cog_control_ecoli.tsv' if control else ''}")
    say()

    if failures:
        sys.exit("FAILED spot checks:\n  " + "\n  ".join(failures))
    if control and control["letter_agreement"] < CONTROL_MIN_LETTER_AGREEMENT:
        sys.exit(
            f"FAILED control: E. coli category agreement {100 * control['letter_agreement']:.2f}% "
            f"is below the {100 * CONTROL_MIN_LETTER_AGREEMENT:.0f}% floor.\n"
            "The CDD -> COG -> letter chain is wrong; do not trust this run."
        )
    rule("=")
    say("stage 02 complete.")
    rule("=")


if __name__ == "__main__":
    main()
