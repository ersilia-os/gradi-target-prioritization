"""Two essentiality columns per protein -- `geptop_ess` (predicted) and `deg_ess` (measured).

    geptop_ess   continuous 0-1, ALWAYS present. Geptop 2.0's orthology+phylogeny score.
    deg_ess      0 / 0.5 / 1 / EMPTY. The fraction of DEG screens on this exact strain that called
                 the protein essential. Empty means unmeasured -- which is ALL of K. pneumoniae.

**`deg_ess` being three-valued is incidental, not structural**: each species happens to have exactly
2 screens today, so the fraction can only be 0, 0.5 or 1. A third screen would make it quarters.
Read it as "fraction of screens agreeing", never as a fixed three-level scale.

A convenience column `essentiality` merges them -- measurement where it exists, prediction where it
does not -- with `essentiality_source` saying which. Read the two columns directly if you want to
apply your own rule; the merged one bakes in a choice (see the caveats below).

This is the column a target ranking actually consumes. It follows stage 04's pattern exactly —
`<x>_measured` (the label, empty where unmeasured), `<x>_score` (a value for EVERY protein), and
`<x>_source` saying which you got — rather than silently overwriting one with the other.

WHAT IS MEASURED, AND WHAT IS NOT
---------------------------------
DEG covers **our exact anchor strains** for two of three species, and not the third:

    ecoli     DEG1018 (genetic footprinting, 614 essential) + DEG1019 (single-gene knockout, 296)
    saureus   DEG1017 (TMDH, 345) + DEG1061 (Tn-seq, 284)
    kpneumon. NOTHING -- Klebsiella is absent from DEG entirely, which is why Geptop exists

Joined **by sequence, not accession** (the house rule): 99.4% of E. coli and 96.9% of S. aureus DEG
rows match a proteome sequence exactly.

**TWO SCREENS OF THE SAME STRAIN DISAGREE BY 2x.** E. coli MG1655 is 614 essential by footprinting
and 296 by knockout; S. aureus NCTC 8325 is 345 by TMDH and 284 by Tn-seq. So "the measured label"
is not a single thing, and this script does NOT pretend otherwise: it emits `deg_n_datasets`,
`deg_n_essential`, `deg_essential_any` and `deg_essential_all`, and `--rule` chooses which drives
the merged column. `any` is the default because a gene called essential by one credible genome-wide
screen is a real finding; `all` is the conservative read. **Neither is "the truth".**

THREE CAVEATS ON THE MERGED COLUMN, STATED BECAUSE IT IS THE ONE PEOPLE WILL RANK ON
------------------------------------------------------------------------------------
1. **It mixes units.** A measured call is pinned to 1.0 / 0.0; a prediction is a continuous Geptop
   score. A measured essential therefore outranks every prediction by construction. That is
   defensible within a species -- measurement beats prediction -- but see (2).
2. **It is NOT comparable across species.** E. coli and S. aureus rankings are dominated by
   measured calls while K. pneumoniae is entirely predicted, so pooling the three and sorting would
   systematically favour the measured species. Rank within a species.
3. **For E. coli and S. aureus the PREDICTION half is circular anyway** -- both are Geptop
   references, contributing 53.2% and 58.3% of their own score. It barely matters here because
   those two are mostly measured, but never quote their Geptop score as a prediction.

Writes `data/processed/essentiality/essentiality_<species>.tsv`.
Run with the `gradi` env; seconds, no external tools.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import essentiality as Es  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
# DEG datasets measured on OUR EXACT anchor strains. Anything else would be a cross-strain
# transfer, which is a different claim and is not made here.
ANCHOR_DATASETS = {
    "ecoli": ["DEG1018", "DEG1019"],
    "saureus": ["DEG1017", "DEG1061"],
    "kpneumoniae": [],                 # Klebsiella is absent from DEG
}
SPECIES = ("kpneumoniae", "ecoli", "saureus")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def measured(sp: str, lab: pd.DataFrame) -> pd.DataFrame:
    """DEG calls for one species, mapped onto our proteome BY SEQUENCE and aggregated per protein."""
    prot = P.load(sp)[["uniprot_ac", "sequence"]].copy()
    ids = ANCHOR_DATASETS[sp]
    empty = pd.DataFrame({"uniprot_ac": prot.uniprot_ac, "deg_n_datasets": 0,
                          "deg_n_essential": 0})
    if not ids:
        return empty
    sub = lab[lab["deg_dataset_id"].astype(str).isin(ids)].copy()
    if sub.empty:
        return empty
    # sequence -> uniprot_ac. A duplicated sequence in the proteome would make this ambiguous, so
    # resolve to the first accession and report how many were affected rather than hiding it.
    seq2ac = prot.drop_duplicates("sequence").set_index("sequence")["uniprot_ac"]
    dup = int(prot.sequence.duplicated().sum())
    sub["uniprot_ac"] = sub["sequence"].astype(str).map(seq2ac)
    matched = sub["uniprot_ac"].notna()
    say(f"    {sp}: {len(sub):,} DEG rows over {len(ids)} dataset(s) -> "
        f"{int(matched.sum()):,} matched by exact sequence ({100*matched.mean():.1f}%)"
        + (f"; {dup} duplicated proteome sequences resolved to first accession" if dup else ""))
    sub = sub[matched]
    # Count DISTINCT DATASETS, not matching rows. `size` here was a real bug: several DEG proteins
    # can share a sequence and collapse onto one accession, which inflated `deg_n_datasets` to as
    # much as 20 when only 2 datasets exist. Within a dataset, a protein counts as essential if ANY
    # of the rows mapping to it is essential.
    per_ds = (sub.groupby(["uniprot_ac", "deg_dataset_id"])["essential"].max()
              .reset_index())
    agg = per_ds.groupby("uniprot_ac").agg(deg_n_datasets=("deg_dataset_id", "nunique"),
                                           deg_n_essential=("essential", "sum")).reset_index()
    dup_rows = len(sub) - len(per_ds)
    if dup_rows:
        say(f"      {dup_rows} DEG rows collapsed onto an accession already covered by the same "
            "dataset (shared sequences); counted once per dataset, not once per row")
    return empty.drop(columns=["deg_n_datasets", "deg_n_essential"]).merge(agg, on="uniprot_ac",
                                                                          how="left").fillna(
        {"deg_n_datasets": 0, "deg_n_essential": 0})


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--rule", default="any", choices=["any", "all"],
                    help="with several screens of one strain: 'any' screen calling it essential "
                         "(default) or 'all' of them (conservative)")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STAGE 07 - merged essentiality: measurement where it exists, prediction where it does not")
    rule("=")
    say(f"  rule         '{args.rule}' -- {'any screen' if args.rule=='any' else 'every screen'} "
        "calling a gene essential makes it essential")
    say("  measured     DEG, joined BY SEQUENCE, only on our exact anchor strains")
    say("  predicted    Geptop 2.0 score (accuracy 0.59-0.84 by organism; the paper quotes 0.84)")
    say("  NOTE         the merged column MIXES UNITS (measured pinned to 1/0, predictions")
    say("               continuous) and is NOT comparable across species -- rank within a species")
    rule("=")

    lab = Es.load_labels()
    rule()
    say("JOIN")
    rule()
    rows = []
    for sp in args.species:
        g = (Es.load_geptop(sp)[["uniprot_ac", "geptop_score", "geptop_evidence",
                                 "geptop_in_reference_set"]]
             .rename(columns={"geptop_score": "geptop_ess"}))
        m = measured(sp, lab)
        if not ANCHOR_DATASETS[sp]:
            say(f"    {sp}: no DEG dataset on this strain -- every row is a PREDICTION")
        df = g.merge(m, on="uniprot_ac", how="left").fillna({"deg_n_datasets": 0,
                                                             "deg_n_essential": 0})
        df["deg_n_datasets"] = df["deg_n_datasets"].astype(int)
        df["deg_n_essential"] = df["deg_n_essential"].astype(int)
        has = df["deg_n_datasets"] > 0
        any_ess = df["deg_n_essential"] > 0
        all_ess = has & (df["deg_n_essential"] == df["deg_n_datasets"])
        df["deg_essential_any"] = np.where(has, any_ess.astype(int), pd.NA)
        df["deg_essential_all"] = np.where(has, all_ess.astype(int), pd.NA)
        # The only GRADED experimental signal available: what fraction of the screens covering this
        # protein called it essential. With 2 screens that is 0, 0.5 or 1 -- crude, but it separates
        # "both screens agree" from "the screens disagree", which is a real distinction: on E. coli
        # 487 of 692 essential-by-any calls come from ONE screen only, and just 205 have both.
        df["deg_ess"] = np.where(
            has, (df["deg_n_essential"] / df["deg_n_datasets"].where(has, 1)).round(4), pd.NA)
        chosen = any_ess if args.rule == "any" else all_ess
        # The merged column: 1.0/0.0 where measured, the Geptop score where not.
        df["essentiality"] = np.where(has, chosen.astype(float), df["geptop_ess"]).round(4)
        df["essentiality_source"] = np.where(has, "measured", "predicted")
        df["essentiality_rule"] = np.where(has, args.rule, "geptop")
        # ---- DEG gets its own per-species file, symmetric with geptop_<sp>.tsv and ogee_<sp>.tsv.
        # One file per EVIDENCE SOURCE, each carrying its own provenance columns, and a headline
        # table that summarises them with one column each. That is the relationship
        # geptop_<sp>.tsv already had with the headline; DEG now matches it instead of having its
        # five provenance columns inlined into the summary.
        deg_cols = ["uniprot_ac", "deg_ess", "deg_n_datasets", "deg_n_essential",
                    "deg_essential_any", "deg_essential_all"]
        deg_out = OUT_DIR / f"deg_{sp}.tsv"
        df[deg_cols].to_csv(deg_out, sep="\t", index=False)
        say(f"    -> {deg_out.relative_to(REPO_ROOT)}  {len(df):,} rows, "
            f"{'all unmeasured' if not has.any() else f'{int(has.sum()):,} measured'}")

        # ---- the headline: ONE COLUMN PER SOURCE plus the merge. `ogee_ess` and `screens_mean`
        # are filled by their own scripts and are absent until those have run, which is why they
        # are not invented here -- an invented column is indistinguishable from a measured one.
        # ONE COLUMN PER EVIDENCE SOURCE. Each source's own file carries its provenance columns;
        # the headline carries the summary. `ogee_ess` and `screens_mean` are joined in if their
        # scripts have run -- absent rather than invented, because an invented column is
        # indistinguishable from a measured one.
        for src, col in (("ogee", "ogee_ess"), ("screens", None)):
            f = OUT_DIR / f"{src}_{sp}.tsv"
            if not f.exists():
                say(f"       {src}_{sp}.tsv absent -- headline omits it (run its script to add)")
                continue
            t = pd.read_csv(f, sep="\t")
            if col:
                df = df.merge(t[["uniprot_ac", col]], on="uniprot_ac", how="left")
            else:
                # The nine screen columns are one opinion each; the mean is the summary. They
                # correlate heavily (all read the same embedding), so this is NOT nine independent
                # votes -- rank on it, do not read it as a consensus count.
                cols = [c for c in t.columns if c != "uniprot_ac"]
                df = df.merge(t[["uniprot_ac"]].assign(
                    screens_mean=t[cols].mean(axis=1).round(4)), on="uniprot_ac", how="left")
        head = ["uniprot_ac", "geptop_ess", "deg_ess", "ogee_ess", "screens_mean", "essentiality",
                "essentiality_source", "essentiality_rule",
                "geptop_evidence", "geptop_in_reference_set"]
        df = df[[c for c in head if c in df.columns]]
        out = OUT_DIR / f"essentiality_{sp}.tsv"
        df.to_csv(out, sep="\t", index=False)
        n_meas = int(has.sum())
        say(f"    -> {out.relative_to(REPO_ROOT)}  {len(df):,} rows   "
            f"measured {n_meas:,} ({100*n_meas/len(df):.1f}%)   predicted {len(df)-n_meas:,}")
        if n_meas:
            say(f"       measured essential: any={int(any_ess.sum()):,}  all={int(all_ess.sum()):,}"
                f"   (the {args.rule} rule is in `essentiality`)")
        rows.append({"species": sp, "n": len(df), "n_measured": n_meas,
                     "n_predicted": len(df) - n_meas,
                     "n_measured_essential_any": int(any_ess.sum()),
                     "n_measured_essential_all": int(all_ess.sum()),
                     "deg_datasets": ",".join(ANCHOR_DATASETS[sp]) or "none",
                     "rule": args.rule})
    pd.DataFrame(rows).to_csv(OUT_DIR / "evidence" / "essentiality_merge_manifest.tsv",
                              sep="\t", index=False)
    rule("=")
    say("merged essentiality complete.")
    rule("=")


if __name__ == "__main__":
    main()
