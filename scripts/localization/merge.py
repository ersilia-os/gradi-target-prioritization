"""Stack the two localization predictors into the axis's single deliverable.

    data/processed/localization/localization_<species>.tsv
        uniprot_ac  localization  confidence  cytoplasmic_fraction  has_signal_peptide

**No `evidence` column.** Both predictors cover 100% of every proteome by construction, so it was
constant across all 13,020 proteins and carried no information. The completeness check it would
have encoded still runs -- this stage exits non-zero if either track is silent for any protein.

**One table, two predictors, side by side -- this script does not arbitrate between them.**
DeepLocPro answers *which compartment* and TMbed answers *how much of the chain faces the
cytoplasm*; neither derives from the other, and where they disagree that disagreement is the
information. Reducing them to one call would throw away the only cross-check this axis has, since
TMbed corroborates `extracellular` -- DeepLocPro's weakest class -- from outside that model.

**Which to believe, when you must choose: prefer `cytoplasmic_fraction` over `localization`.**
DeepLocPro always returns a call, so its 100% coverage is a property of the method, not evidence,
and it has no `unknown` class to fall back on.

**The Gram-positive trap, carried in the data.** On *S. aureus* DeepLocPro runs in `positive`
mode, which does not merely mask `periplasm` and `outer_membrane` -- it ADDS their probability
mass into `extracellular`. So Sa has four reachable classes, its extracellular count absorbs
whatever the model wanted to call periplasmic, and the absence of those two labels is structural
rather than missing data. `evidence/deeplocpro_probabilities_<species>.tsv` is the only place that
moved mass can be measured.

Recomputes nothing. Seconds, no model, no GPU -- `scripts/localization/predict.py` does the work
and writes the two per-predictor tables to `evidence/`. Run in `gradi`:

    python scripts/localization/merge.py
    python scripts/localization/merge.py --species ecoli -q
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import localization as LOC  # noqa: E402
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "localization"
EVIDENCE_DIR = OUT_DIR / "evidence"

COLUMNS = ["uniprot_ac", "localization", "confidence",
           "cytoplasmic_fraction", "has_signal_peptide"]

# **No `evidence` column, by decision.** Both predictors cover 100% of every proteome BY
# CONSTRUCTION -- DeepLocPro always returns a call and TMbed labels every residue -- so the column
# was constant for all 13,020 proteins and said nothing. The check it would have carried is not
# dropped, only un-materialised: `merge_species` still derives per-row coverage and the stage
# exits non-zero if either track is ever silent. Do not re-add it as a constant.


def say(msg: str = "") -> None:
    print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def merge_species(species: str, quiet: bool) -> pd.DataFrame:
    """One species. Both predictors are complete by construction, so this is an inner join over
    identical keys -- but it is written as an outer join and audited, because a silently dropped
    protein is exactly the failure the completeness rule exists to catch."""
    dlp = LOC.load_deeplocpro(species)
    tmb = LOC.load_tmbed(species)
    df = dlp.merge(tmb, on="uniprot_ac", how="outer")

    # Per-row coverage, derived rather than assumed. It is not written out -- see COLUMNS -- but
    # it is what the caller's guard reads, so a silent track fails the stage instead of shipping
    # a half-empty column that looks like a measurement.
    df.attrs["n_deeplocpro"] = int(df["localization"].notna().sum())
    df.attrs["n_tmbed"] = int(df["cytoplasmic_fraction"].notna().sum())

    df = M.reindex(df[COLUMNS], species)   # canonical row order -- see src/matrices.py

    if not quiet:
        n = len(df)
        say(f"  {species:14s} {n:>6,} proteins")
        say(f"    deeplocpro             {df.attrs['n_deeplocpro']:>6,}  "
            f"{100 * df.attrs['n_deeplocpro'] / n:>5.1f}%")
        say(f"    tmbed                  {df.attrs['n_tmbed']:>6,}  "
            f"{100 * df.attrs['n_tmbed'] / n:>5.1f}%")
        reach = sorted(df["localization"].dropna().unique())
        unreached = [c for c in LOC.LOC_CLASSES if c not in reach]
        say(f"    compartments reached   {len(reach)}/{len(LOC.LOC_CLASSES)}"
            + (f"   never called: {', '.join(unreached)}" if unreached else ""))
        say(f"    cytoplasmic_fraction   median {df['cytoplasmic_fraction'].median():.3f}"
            f"   signal peptide {100 * df['has_signal_peptide'].mean():.1f}%")
    return df


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Stack deeplocpro_<sp>.tsv and tmbed_<sp>.tsv into localization_<sp>.tsv.")
    ap.add_argument("--species", nargs="+", default=list(LOC.SPECIES), choices=list(LOC.SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()

    if not args.quiet:
        rule("=")
        say("localization/merge.py - one table, two predictors, side by side")
        rule("=")
        say(f"  in       : {EVIDENCE_DIR.relative_to(REPO_ROOT)}/"
            "{deeplocpro,tmbed}_<species>.tsv")
        say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/localization_<species>.tsv")
        say(f"  columns  : {', '.join(COLUMNS)}")
        say(f"  species  : {', '.join(args.species)}")
        rule()

    failures = []
    for species in args.species:
        df = merge_species(species, args.quiet)
        n_expected = len(P.load(species))
        if len(df) != n_expected:
            failures.append(f"{species}: wrote {len(df)} rows, proteome has {n_expected}")
        for track, column in (("deeplocpro", "localization"), ("tmbed", "cytoplasmic_fraction")):
            missing = int(df[column].isna().sum())
            if missing:
                failures.append(f"{species}: {track} is silent for {missing} proteins -- both "
                                f"predictors cover every protein by construction, so this is a "
                                f"broken run, not a sparse one")
        df.to_csv(OUT_DIR / f"localization_{species}.tsv", sep="\t", index=False)

    if not args.quiet:
        rule()
        say(f"  wrote {len(args.species)} table(s) to "
            f"{OUT_DIR.relative_to(REPO_ROOT)}/localization_<species>.tsv")
        say("  the two per-predictor tables stay in evidence/ -- this one recomputes nothing")

    if failures:
        say("")
        for f in failures:
            print(f"  FAIL {f}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
