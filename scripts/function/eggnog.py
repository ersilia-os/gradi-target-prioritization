"""eggNOG-mapper v2 -- orthology-based functional annotation from sequence.

Stage 02, the annotation engine behind the GO-slim scheme. Runs eggNOG-mapper over each bacterial
proteome and writes its per-protein annotation as a first-class table.

Why this tool
-------------
Stage 02's first scheme (COG) tops out at 79.1 / 84.5 / 73.1% -- its own ceiling, since NCBI's own
curators reach only 81.6% on E. coli. The proteins it misses are short, uncharacterised, and for
Kp/Sa roughly 60% carry no Pfam and no InterPro either.

**The domain-based route cannot help there, and that was measured, not assumed:**

    interpro2go as a GO tier               filled 5 proteins out of 13,020
    InterPro live API vs UniProt's xref    0 of 15 GO-less Kp proteins would gain a GO term

UniProt's electronic GO is already InterPro2GO-derived, so InterProScan, Pfam+pfam2go and the
InterPro API all re-derive what stage 00 already has. Only a *different kind of evidence* can add
anything, and that is **orthology** -- hence eggNOG-mapper, which transfers annotation from
orthologous groups rather than from domain matches.

eggNOG-mapper is also the tool v1 identified as the right answer at
`legacy/docs/01_task_agnostic.md:240` and never built (`legacy/HISTORY.md:64`).

Environment
-----------
`gradi-emapper` (osx-64, Rosetta) -- eggnog-mapper has no osx-arm64 build, and installing it into
`gradi` would drag that whole env to osx-64 and take ESM-C with it. Same reasoning as `blast`.

Two packaging traps, both hit on the first run:

1. bioconda puts `diamond`/`mmseqs` in the env's `bin/`, but emapper looks for them inside
   `site-packages/eggnogmapper/bin/`. Symlinks are created by `install.sh`.
2. **`download_eggnog_data.py` is broken**: it fetches from `eggnogdb.embl.de`, which no longer
   resolves, and then prints "Finished" anyway -- every wget failed and the exit status was 0.
   The live host is `eggnog5.embl.de`. We fetch the three files ourselves and verify each against
   its Content-Length, because an exit code of 0 is not evidence of data.

Database
--------
~21 GB unpacked under `data/source/eggnog/` (eggnog.db ~12 GB, eggnog_proteins.dmnd ~9 GB,
eggnog.taxa.db). It is a public, re-derivable resource: **do not push it to eosvc.**

Output
------
    data/processed/function/eggnog_<species>.tsv     one row per protein, keyed on uniprot_ac
        uniprot_ac seed_ortholog evalue score eggnog_ogs max_annot_lvl cog_category
        description preferred_name gos ec kegg_ko kegg_pathway kegg_module brite cazy pfams

    scratch/eggnog_raw_<species>.tsv   the .emapper.annotations file as produced
    scratch/.eggnog_<species>.key      md5(accession set | emapper version) -- the cache key
    evidence/eggnog_manifest.tsv

Every protein gets a row; one eggNOG-mapper found no orthologous group for keeps its row with empty
strings. `gos` is consumed by `function/goslim.py`; the rest is here because the same run
produces it and it costs one `to_csv`.

Run with the `gradi` env (it shells out to `gradi-emapper`):
    python scripts/function/eggnog.py
    python scripts/function/eggnog.py --species saureus --limit 200   # smoke test
    python scripts/function/eggnog.py --refresh
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

DB_DIR = REPO_ROOT / "data" / "source" / "eggnog"
OUT_DIR = REPO_ROOT / "data" / "processed" / "function"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
DEFAULT_SPECIES = ("kpneumoniae", "ecoli", "saureus")
EMAPPER_ENV = Path.home() / "miniconda3" / "envs" / "gradi-emapper" / "bin"

# emapper's own column order, lower-cased. `#query` loses its hash when we parse the header.
RAW_COLUMNS = [
    "query", "seed_ortholog", "evalue", "score", "eggnog_ogs", "max_annot_lvl", "cog_category",
    "description", "preferred_name", "gos", "ec", "kegg_ko", "kegg_pathway", "kegg_module",
    "kegg_reaction", "kegg_rclass", "brite", "kegg_tc", "cazy", "bigg_reaction", "pfams",
]
OUT_COLUMNS = [
    "uniprot_ac", "seed_ortholog", "evalue", "score", "eggnog_ogs", "max_annot_lvl",
    "cog_category", "description", "preferred_name", "gos", "ec", "kegg_ko", "kegg_pathway",
    "kegg_module", "brite", "cazy", "pfams",
]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def emapper_bin() -> Path:
    override = os.environ.get("GRADI_EMAPPER_BIN")
    d = Path(override) if override else EMAPPER_ENV
    exe = d / "emapper.py"
    if not exe.exists():
        sys.exit(
            f"emapper.py not found at {exe}\n"
            "  create the env:  CONDA_SUBDIR=osx-64 conda create -n gradi-emapper "
            "-c conda-forge -c bioconda eggnog-mapper=2.1.15\n"
            "  or point GRADI_EMAPPER_BIN at the directory holding emapper.py"
        )
    return exe


def check_db() -> str:
    """The database must be present and whole. A missing file here is hours of confusion later."""
    need = {"eggnog.db": 5e9, "eggnog_proteins.dmnd": 5e9, "eggnog.taxa.db": 1e6}
    missing = []
    for name, floor in need.items():
        p = DB_DIR / name
        if not p.exists():
            missing.append(f"{name} (absent)")
        elif p.stat().st_size < floor:
            missing.append(f"{name} ({p.stat().st_size:,} bytes -- truncated)")
    if missing:
        sys.exit(
            "eggNOG database incomplete in " f"{DB_DIR}:\n  " + "\n  ".join(missing) + "\n\n"
            "  NOTE download_eggnog_data.py is broken -- it fetches from eggnogdb.embl.de, which no\n"
            "  longer resolves, and reports success anyway. Fetch from eggnog5.embl.de instead:\n"
            "    http://eggnog5.embl.de/download/emapperdb-5.0.2/{eggnog.db.gz,"
            "eggnog_proteins.dmnd.gz,eggnog.taxa.tar.gz}\n"
            "  and verify each against its Content-Length before decompressing."
        )
    return f"{sum((DB_DIR / n).stat().st_size for n in need) / 1e9:.1f} GB"


def write_source_md(version: str) -> None:
    (DB_DIR / "SOURCE.md").write_text(
        "# eggNOG-mapper database\n\n"
        f"Recorded {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
        "by scripts/function/eggnog.py\n\n"
        f"emapper: `{version}`\n\n"
        "| file | bytes |\n|---|---|\n"
        + "".join(f"| `{p.name}` | {p.stat().st_size:,} |\n"
                 for p in sorted(DB_DIR.glob("*")) if p.is_file() and p.suffix != ".md")
        + "\nSource: <http://eggnog5.embl.de/download/emapperdb-5.0.2/>\n\n"
        "**`download_eggnog_data.py` does not work**: it points at `eggnogdb.embl.de`, which no\n"
        "longer resolves, and prints `Finished.` even though every `wget` failed. These files were\n"
        "fetched directly from `eggnog5.embl.de` and each was checked against its Content-Length.\n\n"
        "This is a public, re-derivable database. **Do not upload it to eosvc.**\n",
        encoding="utf-8",
    )


def cache_key(accessions: list[str], version: str) -> str:
    return hashlib.md5(("\n".join(sorted(accessions)) + f"|{version}").encode()).hexdigest()


def run_species(species: str, exe: Path, version: str, cpu: int,
                limit: int | None, refresh: bool) -> dict:
    t0 = time.time()
    pre = "smoke_" if limit else ""
    raw_path = SCRATCH_DIR / f"{pre}eggnog_raw_{species}.tsv"
    key_path = SCRATCH_DIR / f".{pre}eggnog_{species}.key"
    out_path = (SCRATCH_DIR / f"smoke_eggnog_{species}.tsv") if limit else (EVIDENCE_DIR / f"eggnog_{species}.tsv")

    df = P.load(species)[["uniprot_ac", "sequence"]]
    if limit:
        df = df.head(limit)
    key = cache_key(df["uniprot_ac"].tolist(), version)

    cached = raw_path.exists() and key_path.exists() and not refresh \
        and key_path.read_text().strip() == key
    if cached:
        say(f"  {species:<14} {len(df):>6} proteins   [cache] {raw_path.name}")
    else:
        if raw_path.exists() and not refresh:
            say(f"  {species:<14} [rerun] cached run was for a different query set or version")
        with tempfile.TemporaryDirectory() as tmp:
            faa = Path(tmp) / f"{species}.faa"
            with faa.open("w") as f:
                for ac, seq in zip(df["uniprot_ac"], df["sequence"]):
                    f.write(f">{ac}\n{seq}\n")
            say(f"  {species:<14} {len(df):>6} proteins   running emapper (diamond, {cpu} cpu) ...")
            cmd = [str(exe), "-i", str(faa), "-o", species, "--output_dir", tmp,
                   "--temp_dir", tmp, "--data_dir", str(DB_DIR), "-m", "diamond",
                   "--itype", "proteins", "--cpu", str(cpu), "--override"]
            res = subprocess.run(cmd, capture_output=True, text=True)
            produced = Path(tmp) / f"{species}.emapper.annotations"
            if res.returncode != 0 or not produced.exists():
                say(res.stdout[-3000:])
                say(res.stderr[-3000:])
                sys.exit(f"emapper failed for {species} (exit {res.returncode})")
            shutil.copy(produced, raw_path)
            key_path.write_text(key)

    ann = pd.read_csv(raw_path, sep="\t", comment="#", header=None,
                      names=RAW_COLUMNS, dtype=str, keep_default_na=False)
    ann = ann.rename(columns={"query": "uniprot_ac"}).replace("-", "")
    out = df[["uniprot_ac"]].merge(ann, on="uniprot_ac", how="left").fillna("")
    out = out[OUT_COLUMNS]
    out.to_csv(out_path, sep="\t", index=False)

    n = len(out)
    return {
        "species": species, "n": n,
        "annotated": int((out["eggnog_ogs"] != "").sum()),
        "with_go": int((out["gos"] != "").sum()),
        "with_ko": int((out["kegg_ko"] != "").sum()),
        "with_ec": int((out["ec"] != "").sum()),
        "with_cog": int((out["cog_category"] != "").sum()),
        "with_name": int((out["preferred_name"] != "").sum()),
        "seconds": round(time.time() - t0, 1), "cached": cached,
        "path": str(out_path.relative_to(REPO_ROOT)),
    }


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES),
                    choices=list(DEFAULT_SPECIES))
    ap.add_argument("--cpu", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--limit", type=int, help="first N proteins per species (smoke test)")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    exe = emapper_bin()
    version = subprocess.run([str(exe), "--version"], capture_output=True, text=True)
    version = next((ln for ln in (version.stdout + version.stderr).splitlines()
                    if ln.startswith("emapper-")), "unknown").strip()

    rule("=")
    say("STAGE 02 - functional annotation :: eggNOG-mapper")
    rule("=")
    say("  in       : data/processed/proteomes/proteome_<species>.tsv  (uniprot_ac + sequence)")
    say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/eggnog_<species>.tsv  (+ evidence/ + scratch/)")
    say(f"  tool     : {version}")
    say(f"  db       : {DB_DIR.relative_to(REPO_ROOT)}/  ({check_db()})")
    say(f"  species  : {', '.join(args.species)}    cpu: {args.cpu}")
    if args.limit:
        say(f"  limit    : {args.limit} proteins per species (SMOKE TEST -> scratch/smoke_*)")
    if args.dry_run:
        say("\n  --dry-run: nothing run, nothing written.")
        return
    write_source_md(version)
    say()

    rule()
    say("ANNOTATE")
    rule()
    rows = [run_species(sp, exe, version, args.cpu, args.limit, args.refresh)
            for sp in args.species]
    pd.DataFrame(rows).to_csv(
        (SCRATCH_DIR if args.limit else EVIDENCE_DIR) / ("smoke_eggnog_manifest.tsv" if args.limit else "eggnog_manifest.tsv"),
        sep="\t", index=False)
    say()

    rule()
    say("COVERAGE")
    rule()
    say(f"  {'species':<14} {'n':>6} {'orthologous group':>19} {'GO':>15} {'KEGG ko':>15} "
        f"{'EC':>13} {'name':>15}")
    for r in rows:
        f = lambda k: f"{r[k]:>6,} ({100 * r[k] / r['n']:>4.1f}%)"  # noqa: E731
        say(f"  {r['species']:<14} {r['n']:>6} {f('annotated'):>19} {f('with_go'):>15} "
            f"{f('with_ko'):>15} {f('with_ec'):>13} {f('with_name'):>15}")
    say(f"\n  minutes: " + ", ".join(f"{r['species']} {r['seconds'] / 60:.1f}" for r in rows))
    say()

    rule("=")
    say("stage 02 (eggNOG-mapper) complete.")
    rule("=")


if __name__ == "__main__":
    main()
