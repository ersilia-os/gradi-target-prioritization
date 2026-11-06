"""One-shot: v2's external inputs move to `data/source/<provider>/`.

`data/raw/` today holds two unrelated things: the external data v2 fetches, bucketed by TASK, and
the frozen v1 archive, bucketed by organism. Two problems follow.

1. **A source is not owned by one task.** The UniProt proteomes feed every axis; bucketing them
   under `proteomes/` says otherwise. Sources are identified by WHERE THEY CAME FROM, which is
   stable, so `data/source/` is bucketed by PROVIDER.
2. **`data/raw/` cannot simply be renamed.** `legacy/data` is a symlink to `data/`, and the frozen
   v1 scripts reference `data/raw/{ecoli,human,kpneumoniae,legacy,other}` directly. Renaming the
   tree would break the v1 pipeline, which is meant to stay runnable.

So only the v2 buckets move out; `data/raw/` stays as the v1 archive that `PROVENANCE.md` indexes.

  python tools/migrate_source.py --dry-run
  python tools/migrate_source.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# old path under data/ -> new path under data/. Provider names, not task names.
MOVES = {
    "raw/proteomes/uniprot":            "source/uniprot/proteomes",
    "raw/embeddings/prott5_reference":  "source/uniprot/prott5_precomputed",
    "raw/embeddings/ecoli":             "source/uniprot/proteome_fasta/ecoli",
    "raw/embeddings/kpneumoniae":       "source/uniprot/proteome_fasta/kpneumoniae",
    "raw/embeddings/saureus":           "source/uniprot/proteome_fasta/saureus",
    "raw/function/cog":                 "source/cdd",
    "raw/function/eggnog":              "source/eggnog",
    "raw/function/go":                  "source/go",
    "raw/function/sprofgo":             "source/sprofgo",
    "raw/orthology/orthodb":            "source/orthodb",
    "raw/essentiality/deg":             "source/deg",
    "raw/essentiality/geptop":          "source/geptop",
    "raw/essentiality/proteomes":       "source/ncbi/deg_proteomes",
}

# code references: ("raw" / "x" / "y") tuples and the prose form, longest first
REWRITES = [
    ('"raw" / "proteomes" / "uniprot"',           '"source" / "uniprot" / "proteomes"'),
    ('"raw" / "embeddings" / "prott5_reference"', '"source" / "uniprot" / "prott5_precomputed"'),
    ('"raw" / "function" / "cog"',                '"source" / "cdd"'),
    ('"raw" / "function" / "eggnog"',             '"source" / "eggnog"'),
    ('"raw" / "function" / "go"',                 '"source" / "go"'),
    ('"raw" / "function" / "sprofgo"',            '"source" / "sprofgo"'),
    ('"raw" / "orthology" / "orthodb"',           '"source" / "orthodb"'),
    ('"raw" / "essentiality" / "deg"',            '"source" / "deg"'),
    ('"raw" / "essentiality" / "geptop"',         '"source" / "geptop"'),
    ('"raw" / "essentiality" / "proteomes"',      '"source" / "ncbi" / "deg_proteomes"'),
    ('"raw" / "proteomes"',                       '"source" / "uniprot" / "proteomes"'),
    ('"raw" / "embeddings"',                      '"source" / "uniprot" / "proteome_fasta"'),
    ('"raw" / "function"',                        '"source"'),
    # prose, in docstrings and docs
    ("data/raw/proteomes/uniprot",   "data/source/uniprot/proteomes"),
    ("data/raw/embeddings/prott5_reference", "data/source/uniprot/prott5_precomputed"),
    ("data/raw/embeddings",          "data/source/uniprot/proteome_fasta"),
    ("data/raw/function/cog",        "data/source/cdd"),
    ("data/raw/function/eggnog",     "data/source/eggnog"),
    ("data/raw/function/go",         "data/source/go"),
    ("data/raw/function/sprofgo",    "data/source/sprofgo"),
    ("data/raw/orthology/orthodb",   "data/source/orthodb"),
    ("data/raw/essentiality/deg",    "data/source/deg"),
    ("data/raw/essentiality/geptop", "data/source/geptop"),
    ("data/raw/essentiality/proteomes", "data/source/ncbi/deg_proteomes"),
    ("data/raw/proteomes",           "data/source/uniprot/proteomes"),
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    print("== moving v2 source trees (data/raw/ keeps the v1 archive)")
    for old, new in MOVES.items():
        src, dst = REPO / "data" / old, REPO / "data" / new
        if not src.exists():
            print(f"  SKIP (absent) {old}")
            continue
        if dst.exists():
            sys.exit(f"REFUSING: {dst} already exists")
        print(f"  data/{old}  ->  data/{new}")
        if not a.dry_run:
            dst.parent.mkdir(parents=True, exist_ok=True)
            src.rename(dst)
    # tidy the now-empty task buckets
    if not a.dry_run:
        for d in ("raw/proteomes", "raw/embeddings", "raw/function", "raw/orthology",
                  "raw/essentiality"):
            p = REPO / "data" / d
            if p.is_dir() and not any(p.iterdir()):
                p.rmdir()
                print(f"  removed empty data/{d}")

    print("\n== rewriting references")
    targets = [p for p in list(REPO.glob("scripts/**/*.py")) + list(REPO.glob("src/*.py"))
               + list(REPO.glob("docs/*.md")) if "legacy" not in p.parts]
    targets += [REPO / f for f in ("CLAUDE.md", "README.md")]
    n = 0
    for p in targets:
        if not p.exists():
            continue
        t = orig = p.read_text()
        for old, new in REWRITES:
            t = t.replace(old, new)
        if t != orig:
            n += 1
            print(f"  {p.relative_to(REPO)}")
            if not a.dry_run:
                p.write_text(t)
    print(f"  {n} files changed")


if __name__ == "__main__":
    main()
