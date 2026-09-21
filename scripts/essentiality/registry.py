"""Every essentiality dataset we ever found, and whether it became a training set.

THE QUESTION THIS TABLE ANSWERS: "we looked at a lot of screens -- which ones are we actually
training on, and why not the rest?" That was previously answerable only by reading prose
(`docs/essentiality_screens.md`) and Python dicts (`screens.py`'s `HELD_BACK`, `UNJOINABLE`). Most
of the datasets on disk are NOT training sets, and the reason each one is out is the expensive part
of the knowledge -- a positives-only list, a condition-dependent screen, an undeposited assembly, a
DESeq run that did not converge.

GENERATED, NEVER HAND-WRITTEN. Every fact here already exists somewhere authoritative:

    screens.py SCREENS / HELD_BACK / UNJOINABLE   the disposition and the reason
    evidence/screen_join_audit.tsv                counts, base rates, join rates
    evidence/deg_datasets.tsv                     DEG's 66 datasets and its own `retained` flag
    evidence/ogee_taxa.tsv                        OGEE's 87 taxa and their tier
    data/{raw,source}/**/SOURCE.md                the per-directory provenance
    training_sets/*.tsv                           what actually exists on disk

A second, hand-maintained copy of that is how a registry stops agreeing with the code it
describes. So this script derives it, and RECONCILES BOTH DIRECTIONS: every `training_set` row must
have a file in `training_sets/`, and every file there must have a row. Either mismatch exits
non-zero, which is what makes the table trustworthy rather than merely present.

Run with the `gradi` env:
    python scripts/essentiality/registry.py
    python scripts/essentiality/registry.py --dry-run
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from screens import HELD_BACK, SCREENS, UNJOINABLE  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
TRAINING_DIR = OUT_DIR / "training_sets"
REGISTRY = OUT_DIR / "dataset_registry.tsv"

# Where a dataset's own provenance lives. A dataset with no SOURCE.md is a gap worth seeing.
SOURCE_ROOTS = (
    REPO_ROOT / "data" / "raw" / "ecoli" / "essentiality",
    REPO_ROOT / "data" / "raw" / "kpneumoniae" / "essentiality",
    REPO_ROOT / "data" / "raw" / "other" / "essentiality",
)
# Only the providers this axis actually uses. Sweeping all of `data/source` pulled in `cdd`,
# `eggnog`, `orthodb`, `go` and `sprofgo`, which belong to the function and orthology axes -- they
# are not essentiality datasets and listing them as "held_back" was simply wrong.
SOURCE_PROVIDERS = ("deg", "ogee", "proteomelm", "geptop")

COLUMNS = ["dataset_id", "organism", "assay", "disposition", "reason", "produces_column",
           "n_labels", "n_positives", "base_rate", "join_rate", "features_from",
           "path_on_disk", "has_source_md"]

VERBOSE = True


def say(m: str = "") -> None:
    if VERBOSE:
        print(m, flush=True)


def rule(c: str = "-", w: int = 118) -> None:
    say(c * w)


def source_dirs() -> dict[str, Path]:
    """directory name -> path, for every staged dataset directory."""
    out = {}
    for root in SOURCE_ROOTS:
        if not root.exists():
            continue
        for d in sorted(root.iterdir()):
            if d.is_dir():
                out[d.name] = d
    for name in SOURCE_PROVIDERS:
        d = REPO_ROOT / "data" / "source" / name
        if d.is_dir():
            out[name] = d
    return out


def declared_path(rel: str) -> tuple[str, bool]:
    """(relative path, has SOURCE.md) for a DECLARED source directory.

    There is no substring fallback on purpose. The previous version guessed, and matched
    `essential_ecoli_bw25113_tradis_goodall` to `data/source/go` -- a wrong provenance link that
    reads as a correct one. An unknown path is left empty, which is visible.
    """
    if not rel:
        return "", False
    p = REPO_ROOT / rel
    return (rel if p.exists() else ""), (p / "SOURCE.md").exists()


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("registry.py -- every essentiality dataset, and whether it became a training set")
    rule("=")
    dirs = source_dirs()
    say(f"  {len(dirs)} staged dataset directories, "
        f"{sum(1 for p in dirs.values() if (p / 'SOURCE.md').exists())} with a SOURCE.md")

    audit = {}
    ja = EVIDENCE_DIR / "screen_join_audit.tsv"
    if ja.exists():
        audit = {r["column_name"]: r for _, r in pd.read_csv(ja, sep="\t").iterrows()}

    rows = []

    # ---- 1. the training sets, from screens.py's own registry
    for s in SCREENS:
        r = audit.get(s.column, {})
        rel, has_md = declared_path(s.source_dir)
        rows.append({
            "dataset_id": s.column, "organism": s.organism, "assay": s.assay,
            "disposition": "training_set", "reason": s.comment,
            "produces_column": s.column,
            "n_labels": r.get("num_total"), "n_positives": r.get("num_positives"),
            "base_rate": r.get("base_rate"), "join_rate": r.get("join_rate"),
            "features_from": s.features,
            "path_on_disk": rel,
            "has_source_md": has_md,
        })

    # ---- 2. downloaded, parseable, deliberately unused
    for k, why in HELD_BACK.items():
        rel, has_md = declared_path(f"data/raw/kpneumoniae/essentiality/{k}")
        if not rel:
            rel, has_md = declared_path(f"data/raw/ecoli/essentiality/{k}")
        rows.append({"dataset_id": k, "organism": "", "assay": "",
                     "disposition": "held_back", "reason": why, "produces_column": "",
                     "path_on_disk": rel, "has_source_md": has_md})

    # ---- 3. measured and refuted -- the most expensive knowledge here
    for k, why in UNJOINABLE.items():
        rows.append({"dataset_id": k, "organism": "", "assay": "",
                     "disposition": "refuted", "reason": why, "produces_column": "",
                     "path_on_disk": "", "has_source_md": False})

    # ---- 4. the bulk corpora: many datasets, used as a corpus rather than one per column
    deg = EVIDENCE_DIR / "deg_datasets.tsv"
    if deg.exists():
        d = pd.read_csv(deg, sep="\t")
        col = "retained" if "retained" in d.columns else None
        kept = int(d[col].sum()) if col else len(d)
        rows.append({"dataset_id": "deg_corpus", "organism": "38 species", "assay": "various",
                     "disposition": "corpus_only",
                     "reason": f"DEG: {len(d)} datasets indexed, {kept} retained after excluding "
                               "non-genome-wide methods and condition-specific screens. Feeds "
                               "`deg_<species>.tsv` for the exact-anchor strains and "
                               "`geptop.py`'s reference labels; not one training set per dataset.",
                     "produces_column": "deg_ess", "n_labels": len(d),
                     "path_on_disk": "data/source/deg",
                     "has_source_md": (REPO_ROOT / "data/source/deg/SOURCE.md").exists()})
    ot = EVIDENCE_DIR / "ogee_taxa.tsv"
    if ot.exists():
        t = pd.read_csv(ot, sep="\t")
        usable = t[(t.base_rate < 0.5) & (t.tier == "prokaryote_ncbi")]
        rows.append({"dataset_id": "ogee_corpus", "organism": f"{len(t)} taxa",
                     "assay": "various", "disposition": "training_set",
                     "reason": f"OGEE v3: {len(t)} taxa with decided E/NE labels, of which "
                               f"{len(usable)} are prokaryotic AND carry both classes -- 40 taxa "
                               "are positives-only (25 from the RB-TnSeq Fitness Browser, which "
                               "cannot see essential genes by construction). Trained as ONE "
                               "cross-species corpus with leave-species-out CV, not per taxon.",
                     "produces_column": "ogee_corpus",
                     "n_labels": int(usable.n_entries.sum()),
                     "n_positives": int(usable.n_essential.sum()),
                     "base_rate": round(usable.n_essential.sum() / usable.n_entries.sum(), 4),
                     "features_from": "per-taxon NCBI assemblies",
                     "path_on_disk": "data/source/ogee",
                     "has_source_md": (REPO_ROOT / "data/source/ogee/SOURCE.md").exists()})

    # ---- 5. staged directories no registry entry claims. NOT noise: these are datasets someone
    # fetched and nothing consumes, which is exactly what this table exists to surface.
    claimed = {r["path_on_disk"] for r in rows if r["path_on_disk"]}
    for name, p in dirs.items():
        rel = str(p.relative_to(REPO_ROOT))
        if rel in claimed or not (p / "SOURCE.md").exists():
            continue
        if "essentiality" not in rel and not rel.startswith("data/source"):
            continue
        txt = (p / "SOURCE.md").read_text()
        used = "**Used in v2?** — YES" in txt or "**Used in v2?** — yes" in txt
        first = next((ln for ln in txt.splitlines() if ln.startswith("**Used in v2?**")), "")
        rows.append({"dataset_id": name, "organism": "", "assay": "",
                     "disposition": "corpus_only" if used else "held_back",
                     "reason": first.replace("**Used in v2?** — ", "")[:300] or "see SOURCE.md",
                     "produces_column": "", "path_on_disk": rel, "has_source_md": True})

    t = pd.DataFrame(rows)
    for c in COLUMNS:
        if c not in t.columns:
            t[c] = pd.NA
    t = t[COLUMNS].sort_values(["disposition", "dataset_id"]).reset_index(drop=True)

    rule()
    say("DISPOSITION")
    rule()
    for k, g in t.groupby("disposition"):
        say(f"  {k:14s} {len(g):>3}   {', '.join(g.dataset_id.head(6))}"
            f"{' ...' if len(g) > 6 else ''}")
    say(f"  {'TOTAL':14s} {len(t):>3}   "
        f"{int(t.has_source_md.fillna(False).sum())} carry a SOURCE.md")

    # ---- reconcile, both directions
    rule()
    say("RECONCILIATION  -- the check that makes this table worth reading")
    rule()
    on_disk = {p.stem for p in TRAINING_DIR.glob("*.tsv")} if TRAINING_DIR.exists() else set()
    declared = set(t.loc[t.disposition == "training_set", "produces_column"].dropna())
    missing = declared - on_disk
    orphan = on_disk - declared
    say(f"  declared training sets : {len(declared)}")
    say(f"  files in {TRAINING_DIR.relative_to(REPO_ROOT)}/ : {len(on_disk)}")
    if missing:
        say(f"  DECLARED BUT ABSENT : {sorted(missing)}")
    if orphan:
        say(f"  ON DISK BUT UNDECLARED : {sorted(orphan)}")

    if a.dry_run:
        rule("=")
        say("dry run: nothing written.")
    else:
        t["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        t.to_csv(REGISTRY, sep="\t", index=False)
        rule()
        say(f"  wrote {REGISTRY.relative_to(REPO_ROOT)}  ({len(t)} rows x {t.shape[1]})")

    if missing or orphan:
        sys.exit("FAILED reconciliation: the registry and training_sets/ disagree. A registry "
                 "that does not match the files it describes is worse than none.")
    rule("=")
    say("registry reconciles with training_sets/ in both directions.")
    rule("=")


if __name__ == "__main__":
    main()
