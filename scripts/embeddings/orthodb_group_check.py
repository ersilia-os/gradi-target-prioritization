"""Do the authors' OrthoDB group vectors speak the same id space as OUR OrthoDB groups?

THIS IS A GATE, NOT A STAGE. ProteomeLM's functional encoding is the mean ESM-C embedding of a
protein's OrthoDB orthologous group, keyed on the OrthoDB group id. We hold group assignments at
`odb12v2`; the authors' pickles carry whatever release they trained on. OrthoDB's own README says
an OG id is "not stable and re-used between releases", and CLAUDE.md already forbids joining v11
ids to v12 ids.

WHY THIS MUST RUN FIRST. If the id spaces differ, every lookup misses and
`build_group_embeddings_for_proteome` falls back to each protein's own ESM-C vector -- which is
exactly our current `self` mode. The "fix" would then reproduce the numbers we already have, with
nothing in the output saying why. A silent no-op that looks like a completed experiment is the
worst available outcome, so it is checked for its own sake before any compute is spent.

The four `group_vectors_*.pkl` files are DISJOINT SIZE BANDS, not nested supersets -- the suffix is
a group-size threshold and the authors' loader MERGES every file whose threshold is >= the
requested `min_group_size`. So `_200` alone holds only the largest groups (1.4 GB) while `_0` adds
the long tail of small ones (18.1 GB). A low per-protein coverage from `_200` alone is therefore
expected and is NOT evidence of a release mismatch; the discriminator is whether the ids overlap at
all and whether the key FORMAT matches.

Run with the `gradi` env:
    python scripts/embeddings/orthodb_group_check.py
    python scripts/embeddings/orthodb_group_check.py --min-group-size 50
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import orthology as O  # noqa: E402

DB_DIR = REPO_ROOT / "data" / "source" / "proteomelm"
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "embeddings" / "evidence"
SPECIES = ("ecoli", "kpneumoniae", "saureus")
GROUP_COLS = ("orthodb_og_domain", "orthodb_og_narrow")

# `<int>at<level_taxid>` -- the odb12 format, e.g. 9802800at2 (Bacteria), 43881at91347.
OG_RE = re.compile(r"^(\d+)at(\d+)$")

VERBOSE = True


def say(m: str = "") -> None:
    if VERBOSE:
        print(m, flush=True)


def rule(c: str = "-", w: int = 100) -> None:
    say(c * w)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--min-group-size", type=int, default=200, choices=[0, 10, 50, 200])
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("orthodb_group_check.py -- is the authors' group-vector id space OUR id space?")
    rule("=")
    say(f"  vectors  {DB_DIR.relative_to(REPO_ROOT)}/group_vectors_*.pkl  "
        f"(min_group_size={a.min_group_size})")
    say(f"  ours     data/processed/orthology/orthodb_<species>.tsv  (odb12v2)")
    say("  gate     near-zero id overlap == a different OrthoDB release == STOP")
    rule()

    # The authors' own loader, deliberately: it owns the band-merging semantics and the
    # `(mean_embedding, group_size)` value shape. Re-implementing the unpickling would be
    # inventing a second opinion about their file format.
    #
    # BUT OUR OWN `scripts/embeddings/proteomelm.py` SHADOWS THE INSTALLED PACKAGE. Python puts a
    # script's own directory first on sys.path, so `import proteomelm` finds our single-file module
    # and `proteomelm.utils` dies with "not a package" -- which reads like a broken install rather
    # than a name collision. Drop this directory before importing.
    _here = str(Path(__file__).resolve().parent)
    sys.path[:] = [p for p in sys.path if p not in ("", ".", _here)]
    from proteomelm.utils.proteome import load_orthodb_group_vectors
    say("loading group vectors (this reads the whole pickle into RAM) ...")
    means = load_orthodb_group_vectors(str(DB_DIR), min_group_size=a.min_group_size)
    if not means:
        sys.exit(f"FATAL no group vectors loaded from {DB_DIR}. Download at least "
                 f"group_vectors_{a.min_group_size}.pkl first.")

    keys = list(means)
    sample = means[keys[0]]
    say(f"  {len(keys):,} groups loaded")
    say(f"  value shape {tuple(sample.shape)}  dtype {sample.dtype}")
    if tuple(sample.shape) != (1152,):
        sys.exit(f"FATAL expected a (1152,) mean embedding, got {tuple(sample.shape)} -- this is "
                 "not the ESM-C-dimension group table ProteomeLM expects.")

    # Key format, and the taxonomic LEVELS present. The level tells us which clades the authors'
    # groups were taken at, which is independent evidence about the release and the scope.
    parsed = [OG_RE.match(str(k)) for k in keys]
    n_ok = sum(1 for m in parsed if m)
    say(f"  key format `<int>at<taxid>`: {n_ok:,}/{len(keys):,} ({n_ok / len(keys):.1%})")
    say(f"  examples: {keys[:5]}")
    if n_ok < 0.9 * len(keys):
        say(f"  NOTE only {n_ok / len(keys):.1%} parse as OrthoDB ids -- inspect before trusting")
    levels = Counter(m.group(2) for m in parsed if m)
    say(f"  top levels: {dict(levels.most_common(8))}")

    keyset = set(map(str, keys))
    rule()
    say("OVERLAP WITH OUR GROUPS")
    rule()
    say(f"  {'species':12s} {'column':20s} {'our groups':>11} {'in table':>9} {'':>7}  "
        f"{'proteins':>9} {'covered':>8}")
    rows = []
    for sp in SPECIES:
        d = O.load_orthodb(sp)
        for col in GROUP_COLS:
            g = d[col].dropna().astype(str)
            g = g[g != ""]
            uniq = set(g)
            hit_ids = uniq & keyset
            # per-PROTEIN coverage is what actually matters for the encoding
            covered = int(g.isin(keyset).sum())
            say(f"  {sp:12s} {col:20s} {len(uniq):11,} {len(hit_ids):9,} "
                f"{len(hit_ids) / max(len(uniq), 1):7.1%}  {len(d):9,} "
                f"{covered / len(d):8.1%}")
            rows.append({"species": sp, "group_column": col,
                         "min_group_size": a.min_group_size,
                         "n_group_vectors": len(keys),
                         "our_distinct_groups": len(uniq),
                         "groups_in_table": len(hit_ids),
                         "group_overlap_frac": round(len(hit_ids) / max(len(uniq), 1), 4),
                         "n_proteins": len(d),
                         "proteins_covered": covered,
                         "protein_coverage_frac": round(covered / len(d), 4)})

    t = pd.DataFrame(rows)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    out = EVIDENCE_DIR / "orthodb_group_vector_overlap.tsv"
    t.to_csv(out, sep="\t", index=False)
    rule()
    say(f"  wrote {out.relative_to(REPO_ROOT)}")

    best = t["group_overlap_frac"].max()
    rule("=")
    if best < 0.02:
        say(f"GATE FAILED -- best group-id overlap is {best:.2%}.")
        say("The authors' vectors are keyed on a DIFFERENT OrthoDB release than our odb12v2 ids.")
        say("Escalating to a larger pickle will NOT help: the id space itself does not match.")
        say("Do not proceed to the encoding switch.")
        rule("=")
        sys.exit(1)
    say(f"GATE PASSED -- best group-id overlap {best:.2%}, "
        f"best per-protein coverage {t['protein_coverage_frac'].max():.2%}.")
    say("The id spaces are compatible. If per-protein coverage is low, escalate")
    say("--min-group-size 50 / 10 / 0 (each ADDS a band of smaller groups).")
    rule("=")


if __name__ == "__main__":
    main()
