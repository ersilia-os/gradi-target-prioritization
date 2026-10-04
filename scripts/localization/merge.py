"""Stack the two localization predictors into the axis's single deliverable.

    data/processed/localization/localization_<species>.tsv
        uniprot_ac  localization  cytoplasmic_fraction

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

COLUMNS = ["uniprot_ac", "localization", "cytoplasmic_fraction", "localization_evidence"]
# `confidence` is NOT shipped (owner's call, 2026-10-03). Byte-identical in
# `evidence/deeplocpro_<species>.tsv`. KNOW WHAT GOES WITH IT: DeepLocPro always returns a call,
# so `localization` reads equally authoritative for every protein, and `confidence` was the only
# thing saying otherwise -- 12-15% of calls sit below 0.7 and 2-4% below 0.5, i.e. the winning
# class holds less than half the probability mass. Join it back before trusting a single label.
# `has_signal_peptide` is NOT shipped (owner's call, 2026-10-03). It stays byte-identical in
# `evidence/tmbed_<species>.tsv` via `load_tmbed()`, and the console summary below still reports
# it. What it uniquely said, and `cytoplasmic_fraction` alone cannot: WHY a fraction is near zero
# -- exported rather than membrane-buried. That distinction is mechanistically live for
# degradability (a secreted protein transits the cytoplasm unfolded and IS degradable; a membrane
# protein never does and is protected), so `degradability/enrichment.py` joins it back explicitly.

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

    # CROSS-AXIS READ, and the only one in this stage: the evidence level needs curated GO
    # cellular-component terms, which are the one signal here that is NOT a sequence model.
    # Both predictors read the sequence; a GO CC annotation is a human assignment, so it is the
    # only genuinely external check this axis has. Absent -> fail loudly, never silently drop to a
    # two-signal ladder that would wear the same column name.
    fn_path = REPO_ROOT / "data" / "processed" / "function" / f"function_{species}.tsv"
    if not fn_path.exists():
        sys.exit(f"FATAL missing {fn_path.relative_to(REPO_ROOT)} -- run "
                 f"scripts/function/matrix.py first. localization_evidence needs its GO "
                 f"cellular-component terms; computing the level without them would ship a "
                 f"different quantity under the same name.")
    go = pd.read_csv(fn_path, sep="\t", keep_default_na=False).set_index("uniprot_ac")
    df["goslim_terms"] = df["uniprot_ac"].map(go["goslim_terms"])
    df["localization_evidence"] = LOC.localization_evidence(
        df["localization"], df["cytoplasmic_fraction"], df["has_signal_peptide"],
        df["confidence"], df["goslim_terms"])

    keep = COLUMNS + ["confidence", "has_signal_peptide"]
    df = M.reindex(df[keep], species)   # canonical row order -- see src/matrices.py

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
        lv = df["localization_evidence"].value_counts().sort_index()
        say(f"    evidence               " + "  ".join(f"L{k}={v:,}" for k, v in lv.items()))
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
        df[COLUMNS].to_csv(OUT_DIR / f"localization_{species}.tsv", sep="\t", index=False)

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
