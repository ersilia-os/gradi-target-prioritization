"""One-shot: renumbered stages -> task folders.

Moves `scripts/NN_*.py` into `scripts/<task>/<name>.py`, renames the `data/` and `output/`
stage directories to match, and rewrites every reference. Run once, from the repo root.

Ordering is load-bearing: script FILENAMES are rewritten before DIRECTORY tokens, because
`07_proteomelm_embeddings.py` contains the directory token `07_proteomelm` as a prefix.

  python tools/migrate_layout.py --dry-run
  python tools/migrate_layout.py
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

# old script filename -> new path relative to scripts/
SCRIPTS = {
    "00_download_proteomes.py":            "proteomes/download.py",
    "00_proteomes_plots.py":               "plots/proteomes.py",

    "01_embeddings.py":                    "embeddings/esmc.py",
    "01_embeddings_prott5.py":             "embeddings/prott5.py",
    "07_proteomelm_embeddings.py":         "embeddings/proteomelm.py",
    "01_embeddings_projection.py":         "embeddings/projection.py",
    "01_embeddings_projection_plots.py":   "plots/projection.py",

    "02_function_cog.py":                  "function/cog.py",
    "02_function_eggnog.py":               "function/eggnog.py",
    "02_function_goslim.py":               "function/goslim.py",
    "02_function_matrix.py":               "function/matrix.py",
    "02_function_deepgo.py":               "function/deepgo.py",
    "02_function_plots.py":                "plots/function.py",

    "03_localization.py":                  "localization/predict.py",
    "03_localization_plots.py":            "plots/localization.py",

    "05_orthology.py":                     "orthology/orthofinder.py",
    "05_orthology_orthodb.py":             "orthology/orthodb.py",
    "05_orthology_plots.py":               "plots/orthology.py",

    "04_degradability.py":                 "degradability/predict.py",
    "04_degradability_enrichment.py":      "degradability/enrichment.py",
    "04_degradability_plots.py":           "plots/degradability.py",
    "04_degradability_prediction_plots.py":"plots/degradability_predictions.py",
    "04_degradability_top_plots.py":       "plots/degradability_top.py",
    "04_degradability_proteomelm.py":      "degradability/head_comparison.py",
    "04_degradability_regressor.py":       "degradability/regressor.py",

    "07_essentiality_labels.py":           "essentiality/labels.py",
    "07_essentiality_proteomes.py":        "essentiality/deg_proteomes.py",
    "07_essentiality_geptop.py":           "essentiality/geptop.py",
    "07_essentiality_merge.py":            "essentiality/merge.py",

    "06_ligandability_chembl.py":          "ligands/chembl.py",
    "06_ligandability_bindingdb.py":       "ligands/bindingdb.py",
    "06_ligandability_plots.py":           "plots/ligands.py",
}

# workers move to the task that owns them; tabpfn_cv stays top-level (transversal).
WORKERS = {
    "prott5.py":     "embeddings/workers/prott5.py",
    "deeplocpro.py": "localization/workers/deeplocpro.py",
    "lazyqsar_cv.py":"degradability/workers/lazyqsar_cv.py",
}

# stage directory token -> task token (data/raw, data/processed, output/results, output/plots)
DIRS = {
    "00_proteomes":     "proteomes",
    "01_embeddings":    "embeddings",
    "07_proteomelm":    "embeddings",   # MERGE: all three embedding types in one tree
    "02_function":      "function",
    "03_localization":  "localization",
    "04_degradability": "degradability",
    "05_orthology":     "orthology",
    "06_ligandability": "ligands",
    "07_essentiality":  "essentiality",
}

TREES = ["data/raw", "data/processed", "output/results", "output/plots"]

# Merging 07_proteomelm into embeddings collides on `accessory/manifest.tsv` and
# `accessory/shards`: both stages named them generically. Disambiguate by embedding type,
# matching the `shards_prott5_*` / `prott5_*` convention ProtT5 already uses. Applied BEFORE
# the merge; the code constants are patched in CONST_PATCHES below, in the same run.
PRE_MERGE_RENAMES = [
    ("data/processed/01_embeddings/accessory/shards",       "data/processed/01_embeddings/accessory/shards_esmc"),
    ("data/processed/01_embeddings/accessory/manifest.tsv", "data/processed/01_embeddings/accessory/esmc_manifest.tsv"),
    ("data/processed/07_proteomelm/accessory/shards",       "data/processed/07_proteomelm/accessory/shards_proteomelm"),
    ("data/processed/07_proteomelm/accessory/manifest.tsv", "data/processed/07_proteomelm/accessory/proteomelm_manifest.tsv"),
]

# (file, old, new) -- the code that builds the paths renamed above.
CONST_PATCHES = [
    ("scripts/embeddings/esmc.py",       'SHARD_DIR = ACC_DIR / "shards"',      'SHARD_DIR = ACC_DIR / "shards_esmc"'),
    ("scripts/embeddings/esmc.py",       'ACC_DIR / "manifest.tsv"',            'ACC_DIR / "esmc_manifest.tsv"'),
    ("scripts/embeddings/esmc.py",       'accessory/shards/',                   'accessory/shards_esmc/'),
    ("scripts/embeddings/esmc.py",       'accessory/manifest.tsv',              'accessory/esmc_manifest.tsv'),
    ("src/embeddings.py",                'ACCESSORY_DIR / "manifest.tsv"',      'ACCESSORY_DIR / "esmc_manifest.tsv"'),
    ("scripts/embeddings/proteomelm.py", 'SHARD_DIR = ACC_DIR / "shards"',      'SHARD_DIR = ACC_DIR / "shards_proteomelm"'),
    ("scripts/embeddings/proteomelm.py", 'f"{pre}manifest.tsv"',                'f"{pre}proteomelm_manifest.tsv"'),
    ("scripts/embeddings/proteomelm.py", 'accessory/shards/',                   'accessory/shards_proteomelm/'),
    ("scripts/embeddings/proteomelm.py", 'accessory/manifest.tsv',              'accessory/proteomelm_manifest.tsv'),
    ("src/proteomelm.py",                'ACCESSORY_DIR / "manifest.tsv"',      'ACCESSORY_DIR / "proteomelm_manifest.tsv"'),
]


def pre_merge(dry: bool) -> None:
    print("\n== disambiguating names that would collide in the merge")
    for old, new in PRE_MERGE_RENAMES:
        src, dst = REPO / old, REPO / new
        if not src.exists():
            print(f"  SKIP (absent) {old}")
            continue
        if dst.exists():
            sys.exit(f"REFUSING: {new} already exists")
        print(f"  {old} -> {new}")
        if not dry:
            src.rename(dst)


def patch_consts(dry: bool) -> None:
    print("\n== patching the path constants behind those renames")
    for rel, old, new in CONST_PATCHES:
        p = REPO / rel
        if not p.exists():
            if dry:
                print(f"  (dry run: {rel} not moved yet, cannot check)")
                continue
            sys.exit(f"REFUSING: {rel} missing -- run after the moves")
        text = p.read_text()
        if old not in text:
            if new in text:
                print(f"  already patched: {rel}: {old!r}")
                continue
            sys.exit(f"REFUSING: {rel} does not contain {old!r}")
        print(f"  {rel}: {old!r} -> {new!r}  (x{text.count(old)})")
        if not dry:
            p.write_text(text.replace(old, new))


def sh(cmd: list[str], dry: bool) -> None:
    print("  $", " ".join(cmd))
    if not dry:
        subprocess.run(cmd, cwd=REPO, check=True)


def move_dirs(dry: bool) -> None:
    print("\n== data/ and output/ stage directories")
    for tree in TREES:
        for old, new in DIRS.items():
            src, dst = REPO / tree / old, REPO / tree / new
            if not src.is_dir():
                continue
            if dst.exists():
                print(f"  merge {tree}/{old} -> {tree}/{new}")
                merge_tree(src, dst, dry)
            else:
                mv(f"{tree}/{old}", f"{tree}/{new}", dry)


def merge_tree(src: Path, dst: Path, dry: bool) -> None:
    """Move src's contents into dst, descending into directories that exist on both sides.

    Refuses on a FILE collision -- two different files with one name is a decision, not a move.
    Two directories with one name is not: it just means both stages wrote an `accessory/`.
    """
    for p in sorted(src.iterdir()):
        target = dst / p.name
        if not target.exists():
            print(f"    {p.relative_to(REPO)} -> {target.relative_to(REPO)}")
            if not dry:
                p.rename(target)
        elif p.is_dir() and target.is_dir():
            merge_tree(p, target, dry)
        else:
            sys.exit(f"REFUSING: {target} exists and is not a mergeable directory")
    if not dry:
        src.rmdir()


def _tracked(rel: str) -> bool:
    r = subprocess.run(["git", "ls-files", "--error-unmatch", rel], cwd=REPO, capture_output=True)
    return r.returncode == 0


def mv(old: str, new: str, dry: bool) -> None:
    """git mv when git knows the path (history follows), plain mv otherwise."""
    sh(["git", "mv", old, new] if _tracked(old) else ["mv", old, new], dry)


def move_scripts(dry: bool) -> None:
    print("\n== scripts/")
    for new in sorted({v.split("/")[0] for v in SCRIPTS.values()} |
                      {v.rsplit("/", 1)[0] for v in WORKERS.values()}):
        d = REPO / "scripts" / new
        print(f"  mkdir {d.relative_to(REPO)}")
        if not dry:
            d.mkdir(parents=True, exist_ok=True)
    for old, new in sorted(SCRIPTS.items()):
        if not (REPO / "scripts" / old).exists():
            print(f"  SKIP (absent) {old}")
            continue
        mv(f"scripts/{old}", f"scripts/{new}", dry)
    for old, new in sorted(WORKERS.items()):
        if not (REPO / "scripts" / "workers" / old).exists():
            print(f"  SKIP (absent) workers/{old}")
            continue
        mv(f"scripts/workers/{old}", f"scripts/{new}", dry)


def rewrite(dry: bool) -> None:
    """Rewrite references. Filenames FIRST, then directory tokens (prefix collision)."""
    print("\n== rewriting references")
    targets: list[Path] = []
    for pat in ("scripts/**/*.py", "src/**/*.py", "docs/**/*.md", "tools/**/*.py"):
        targets += [p for p in REPO.glob(pat)
                    if "legacy" not in p.parts and p.resolve() != Path(__file__).resolve()]
    targets += [REPO / f for f in ("CLAUDE.md", "README.md", "install.sh", "requirements.txt")]

    n_files = 0
    for p in targets:
        if not p.exists():
            continue
        text = original = p.read_text()

        # 1. script filenames (longest first, so a shorter name is never a prefix of a longer one)
        for old in sorted(SCRIPTS, key=len, reverse=True):
            text = text.replace(old, SCRIPTS[old])
        # workers: only the `workers/<name>` form, so a bare `prott5.py` elsewhere is untouched
        for old, new in WORKERS.items():
            text = text.replace(f"workers/{old}", new)

        # 2. directory tokens
        for old, new in DIRS.items():
            text = text.replace(f'"{old}"', f'"{new}"')      # Path(...) literals
            text = text.replace(f"/{old}/", f"/{new}/")       # prose paths in docstrings/md
            text = text.replace(f"data/raw/{old}", f"data/raw/{new}")
            text = text.replace(f"data/processed/{old}", f"data/processed/{new}")
            text = text.replace(f"output/results/{old}", f"output/results/{new}")
            text = text.replace(f"output/plots/{old}", f"output/plots/{new}")

        # 3. repo-root depth: scripts are one level deeper now (src/ and legacy/ unchanged)
        if p.suffix == ".py" and p.parts[len(REPO.parts)] == "scripts":
            text = text.replace("parents[1]", "parents[2]")

        if text != original:
            n_files += 1
            print(f"  rewrote {p.relative_to(REPO)}")
            if not dry:
                p.write_text(text)
    print(f"  {n_files} files changed")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.dry_run:
        print("DRY RUN -- nothing is written\n")
    move_scripts(a.dry_run)
    pre_merge(a.dry_run)
    move_dirs(a.dry_run)
    rewrite(a.dry_run)
    patch_consts(a.dry_run)
    print("\ndone.")


if __name__ == "__main__":
    main()
