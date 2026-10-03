"""One column per evidence source, plus the merged column a ranking consumes.

    essentiality_<species>.tsv
        uniprot_ac   geptop_ess   proteomelm_ess   screens_ess_mean
        essentiality   essentiality_source   geptop_evidence   geptop_in_reference_set

    deg_<species>.tsv
        uniprot_ac   deg_ess   deg_n_datasets   deg_n_essential
        deg_essential_any   deg_essential_all

**The MEASURED column lives in `deg_<species>.tsv`, not in the summary** (owner's call,
2026-10-03). `deg_ess` used to appear in both, byte-identical and in the same canonical row order;
the per-source file is the richer of the two, because it carries `deg_essential_any` and
`deg_essential_all` side by side and the summary could only carry whichever `--rule` picked. The
measurement is still *inside* `essentiality`, with `essentiality_source` reading `measured` --
nothing about the merged column changed.

**`ogee_ess` is not in the summary either** (owner's call, 2026-10-03), though `ogee_<species>.tsv`
and its scripts are untouched. The reason is NOT that it duplicates `proteomelm_ess` -- measured,
those two are rho 0.33-0.64 with only 296-384 of their top 500 shared, so they are genuinely
different opinions. It is that `ogee_ess` is the one most redundant with **`screens_ess_mean`** (rho
0.49 / 0.55 / 0.37, against ProteomeLM's 0.32 / 0.32 / 0.13), and `screens_ess_mean` is both what drives
the merged column and the better-validated of the two: AUROC 0.89-0.96 on the three measured Kp
screens, against OGEE's leave-species-out spread of 0.529-0.940 and a Kp top-decile cut of 0.471
where E. coli reads 0.861. **Re-adding it is one line in `head` below.**

**`essentiality_rule` is gone too, and it carried nothing**: it was a 1:1 function of
`essentiality_source` (`measured` -> `any`, `predicted_screens_ess_mean` -> `screens_ess_mean`). Which rule
a run actually used is a property of the RUN, not of a protein, and it is recorded per species in
`evidence/essentiality_merge_manifest.tsv` -- next to both the `any` and the `all` counts, so the
3.4x spread between them stays visible rather than being collapsed into one repeated string.

**`deg_ess` being three-valued is incidental, not structural**: each species happens to have exactly
2 screens today, so the fraction can only be 0, 0.5 or 1. A third screen would make it quarters.
Read it as "fraction of screens agreeing", never as a fixed three-level scale.

**THE TABLE NO LONGER SHIPS A VERDICT COLUMN** (owner's call, 2026-10-03). `essentiality` and
`essentiality_source` are still computed -- the console summary and the manifest use them -- but
they are not written, along with `geptop_evidence` and `geptop_in_reference_set`, which were
byte-identical duplicates of columns in `geptop_<species>.tsv`.

**So this axis now hands over three predictors side by side and no answer.** That is the real
consequence, and it is deliberate: `essentiality` MIXED UNITS -- a measured call pinned to 1.0/0.0
against a continuous prediction, so every measured essential outranked every prediction by
construction -- and on K. pneumoniae, the anchor, it was a verbatim copy of `screens_ess_mean` for all
5,728 rows. Whoever ranks now chooses which column to rank on, which was always the honest
instruction. Reconstruct the old column with
`np.where(deg_essential_any.notna(), deg_essential_any, screens_ess_mean)`.

The caveats below still apply to any merge you build yourself.

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
        # `essentiality` is computed AFTER the source columns are joined, because its fallback is
        # now `screens_ess_mean` -- which does not exist yet at this point in the function.
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

        # ---- the headline: ONE COLUMN PER SOURCE plus the merge. `ogee_ess` and `screens_ess_mean`
        # are filled by their own scripts and are absent until those have run, which is why they
        # are not invented here -- an invented column is indistinguishable from a measured one.
        # ONE COLUMN PER EVIDENCE SOURCE. Each source's own file carries its provenance columns;
        # the headline carries the summary. `ogee_ess` and `screens_ess_mean` are joined in if their
        # scripts have run -- absent rather than invented, because an invented column is
        # indistinguishable from a measured one.
        for src, col in (("ogee", "ogee_ess"), ("proteomelm_ess", "proteomelm_ess"),
                         ("screens", None)):
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
                    screens_ess_mean=t[cols].mean(axis=1).round(4)), on="uniprot_ac", how="left")
        # ---- THE MERGED COLUMN. Measured where we have it; otherwise the best PREDICTOR we have
        # measured, which is no longer Geptop.
        #
        # It was `geptop_ess` until 2026-10-03. Two measurements moved it, both on K. pneumoniae,
        # the only anchor where every predictor is honest (`geptop_in_reference_set == 0`):
        #
        #   1. On the three measured Kp screens the screens-trained transfer models reach
        #      AUROC 0.89-0.96 (Goodall->Kp, assay-matched), against Geptop's own validated
        #      0.59-0.81 and ProteomeLM-Ess's 0.82-0.95.
        #   2. `geptop_ess` leaves 3,799 of 5,728 Kp proteins (66.3%) tied at EXACTLY 0 -- so
        #      two-thirds of the anchor proteome was unranked in the headline column, and this
        #      axis is consumed by ranking. `screens_ess_mean` has 2,596 distinct values.
        #
        # Do NOT re-derive this from an E. coli comparison: `geptop_ess` is 53.2% self-derived on
        # E. coli and 58.3% on S. aureus (both are Geptop reference genomes), so it looks like the
        # best predictor there -- AUPR 0.977 against Keio -- and that number is circular.
        # `geptop_ess` keeps its own column and its own file; only the FALLBACK changed.
        fallback, rule_name = "geptop_ess", "geptop"
        if "screens_ess_mean" in df.columns and df["screens_ess_mean"].notna().any():
            fallback, rule_name = "screens_ess_mean", "screens_ess_mean"
        df["essentiality"] = np.where(has, chosen.astype(float), df[fallback]).round(4)
        df["essentiality_source"] = np.where(has, "measured", f"predicted_{rule_name}")
        df["essentiality_rule"] = np.where(has, args.rule, rule_name)

        # Everything above stays a WORKING column -- the console lines and the manifest below
        # are computed from them -- but the table ships FOUR: the key and the three predictors.
        # Each dropped column is derivable, which was checked per species before they went:
        #   deg_ess                  byte-identical to the column in `deg_<sp>.tsv`
        #   ogee_ess                 byte-identical to the column in `ogee_<sp>.tsv`
        #   geptop_evidence          byte-identical to the column in `geptop_<sp>.tsv`
        #   geptop_in_reference_set  byte-identical to the column in `geptop_<sp>.tsv`
        #   essentiality             where(measured, deg_essential_any, screens_ess_mean)
        #   essentiality_source      where(deg_essential_any.notna(), measured, predicted_...)
        #   essentiality_rule        a 1:1 function of essentiality_source
        head = ["uniprot_ac", "geptop_ess", "proteomelm_ess", "screens_ess_mean"]
        out = OUT_DIR / f"essentiality_{sp}.tsv"
        df[[c for c in head if c in df.columns]].to_csv(out, sep="\t", index=False)
        n_meas = int(has.sum())
        say(f"    -> {out.relative_to(REPO_ROOT)}  {len(df):,} rows   "
            f"measured {n_meas:,} ({100*n_meas/len(df):.1f}%)   predicted {len(df)-n_meas:,}")
        if n_meas:
            say(f"       measured essential: any={int(any_ess.sum()):,}  all={int(all_ess.sum()):,}"
                f"   (not shipped as a column -- see deg_{sp}.tsv)")
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
