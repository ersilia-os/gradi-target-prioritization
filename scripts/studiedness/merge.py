"""The studiedness deliverable: one complete, canonical table per species.

    data/processed/studiedness/studiedness_<species>.tsv

    uniprot_ac · n_papers_uniprot_own · n_papers_uniprot_prokaryotic · evidence

**It recomputes nothing.** `transfer.py` already scored every protein; this reshapes those tables
into the deliverable, enforces canonical row order through `src/matrices.py`, and runs the checks
that say whether the axis can be trusted. Seconds, no network, no DIAMOND.

TWO PAPER COUNTS, BECAUSE THEY ANSWER DIFFERENT QUESTIONS
-----------------------------------------------------------
`n_papers_uniprot_own` is the number of curated references naming THIS PROTEIN -- by the ligands axis's
`exact` rule, the union of PMIDs over accession, identical sequence and same species at >= 95%
identity, because a protein does not stop being itself between strains. `n_papers_uniprot_prokaryotic` is the
count on its best-studied prokaryotic SwissProt homolog.

**`own` CAN EXCEED `family`, and that is not a bug**: `own` unions every strain entry of this
protein while `family` reads ONE donor's count. On E. coli the medians are own 8 against family
6, because a well-studied organism's own multi-strain literature outweighs any single homolog. Keeping both is what makes **"dark in
K. pneumoniae, famous in E. coli"** readable off a single row -- and on this anchor that is the
normal case, not an edge case. `n_papers_uniprot_own` is near-flat on Kp and Sa *by design*: it is the
measurement of darkness.

**Both are plain integers on the same footing**, so the three species and the two columns compare
directly. There is no scaling and no blend: the composite 0-1 score that shipped first was
removed on 2026-09-22 for being uninterpretable -- see `src/studiedness.py` for the three
measurements that killed it. `src.studiedness.scaled()` derives a 0-1 version on the fly for
anyone combining this axis with the others.

A ZERO IS NOT A MISSING VALUE, AND THERE ARE TWO KINDS OF IT
--------------------------------------------------------------
`n_papers_uniprot_prokaryotic == 0` comes with one of two evidence tiers, and they are different claims:
`no_hit` means **nothing among 575,748 curated entries resembles this protein at all** -- the
strongest novelty signal the axis produces -- while `below_floor` means a distant relative exists
whose literature is simply too far away to carry. Both are answers, not gaps; never impute them.

WHAT THIS SCRIPT CHECKS
------------------------
1. **Completeness and order** -- every protein present, canonical order, no nulls in either score.
2. **Spot checks** -- named workhorses (rpoB, gyrB, ftsZ, clpP, dnaA, secA, groEL, rplB) must rank
   high on `n_papers_uniprot_prokaryotic`. A transfer joined to the wrong accessions would still produce a
   well-formed table, so the control has to be named genes. The bar is calibrated on what these
   genes actually reach, not on an intuition -- see `SPOT_PERCENTILE`.
3. **Is Unknome telling us anything new?** Spearman of `unknome_knownness` against
   `n_papers_uniprot_prokaryotic` on the proteins where both exist. Reported plainly either way -- if it is
   redundant it is recorded as such rather than quietly dropped.
4. **Is the family score just the own score?** If they agreed, the transfer would be doing no work.
   On Kp and Sa they must NOT agree, and the gap is the point of the axis.

Run with the `gradi` env, after `transfer.py` (and `unknome.py` for check 3).
  python scripts/studiedness/merge.py
  python scripts/studiedness/merge.py --species kpneumoniae -q
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import studiedness as S  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# Workhorses whose families are unambiguously well studied in every bacterium. A wrong join still
# produces plausible numbers, so the check is named genes rather than a row count.
SPOT_GENES = ("rpoB", "gyrB", "ftsZ", "clpP", "dnaA", "secA", "groEL", "rplB")

# CALIBRATED ON THE MEASUREMENT, NOT ON AN ASSUMPTION -- the lesson CLAUDE.md records for the
# essentiality controls, where a bar set at 0.50 separation would have failed the gold standard.
# Measured here: seven of the eight land at 97.8-98.5th percentile in all three species, but
# `rplB` sits at 76.6 (Ec) / 82.1 (Kp) because ribosomal protein L2 genuinely has ~40 papers
# against RNA polymerase's 350 -- it is well studied as part of the ribosome, not as itself.
# That is the axis being right, so the bar is set to catch a BROKEN JOIN rather than to enforce
# an intuition: a wrong join scatters these eight to random percentiles, and all eight clearing
# the 70th by chance is p ~ 7e-5.
SPOT_PERCENTILE = 70.0

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def build(species: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The deliverable plus the full transfer table it came from."""
    path = EVIDENCE_DIR / f"transfer_{species}.tsv"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/transfer.py first")
    tr = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    for c in ("n_papers_uniprot_own", "n_papers_uniprot_prokaryotic", "donor_pident", "donor_n_pubs",
              "n_papers_pubtator_prokaryotic"):
        if c in tr.columns:
            tr[c] = pd.to_numeric(tr[c], errors="coerce")

    # `evidence` is DELIBERATELY NOT IN THE DELIVERABLE (project owner, 2026-10-03). It still
    # ships in evidence/transfer_<species>.tsv, which is where to look when a 0 needs explaining:
    # `no_hit` (nothing in SwissProt resembles this protein -- the strongest novelty claim the
    # axis makes) and `below_floor` (a too-distant curated relative) both read 0 here.
    out = tr[["uniprot_ac", "n_papers_uniprot_own", "n_papers_uniprot_prokaryotic"]].copy()
    out = add_pubtator(out, tr, species)
    if out[["n_papers_uniprot_own", "n_papers_uniprot_prokaryotic"]].isna().any().any():
        n = int(out[["n_papers_uniprot_own", "n_papers_uniprot_prokaryotic"]].isna().any(axis=1).sum())
        sys.exit(f"FATAL {species}: {n:,} rows have a null count. Every protein must carry a "
                 "number -- an unmatched protein has 0 papers, an answer, not a null.")
    for c in ("n_papers_uniprot_own", "n_papers_uniprot_prokaryotic"):
        if (out[c] % 1 != 0).any():
            sys.exit(f"FATAL {species}: {c} is not integral -- it is a paper count, not a score.")
        out[c] = out[c].astype(int)

    # Canonical row order by construction, not by hope. reindex() refuses to invent a missing
    # protein, so an incomplete axis fails here rather than shipping a NaN row.
    order = [c for c in ("uniprot_ac", "n_papers_uniprot_own",
                         "n_papers_uniprot_prokaryotic", "n_papers_pubtator_prokaryotic")
             if c in out.columns]
    out = out[order]
    # NULLABLE INTEGER, not float. pandas needs float to hold the blanks, which writes `14.0`
    # for a paper count and makes the column read as noisy. Int64 keeps blanks blank.
    if "n_papers_pubtator_prokaryotic" in out.columns:
        out["n_papers_pubtator_prokaryotic"] = out["n_papers_pubtator_prokaryotic"].astype("Int64")
    out = M.reindex(out, species)
    M.assert_canonical(out["uniprot_ac"], species)
    return out, tr


def add_pubtator(out: pd.DataFrame, tr: pd.DataFrame, species: str) -> pd.DataFrame:
    """The second count: PubTator3 text-mined papers on the donor's NCBI GeneID.

    A SECOND COLUMN, NEVER MERGED INTO THE FIRST. Two literature definitions mean two columns;
    `max(curated, text-mined)` would switch definition per protein, which is what killed the 0-1
    composite on 2026-09-22.

    **KEYED ON NCBI GeneID, SO THE SPECIES IS ALREADY IN THE KEY** -- GeneID 947587 IS E. coli
    K-12 `ftsZ`. There is no species term to add and no gene-symbol ambiguity: the symbol route
    this replaced read **345,630** for Kp `crp` (human C-reactive protein), where the GeneID
    route reads **376**. No API calls either -- the counts come from `gene2pubtator3.gz`.

    **ITS DONOR IS CHOSEN BY PUBTATOR, NOT BY CURATED PAPERS, AND THAT IS HALF THE SIGNAL.**
    Held-out control: selecting on curated and reading PubTator scores **0.1810**; selecting on
    PubTator and reading PubTator scores **0.3552**, which beats the curated column's **0.3280**
    on identical folds. `pubtator_donor_*` in `evidence/transfer_<sp>.tsv` names which donor.

    **A BLANK IS NOT A ZERO.** Empty means no in-scope donor carried an NCBI GeneID, so the route
    could not look anything up -- unreachable, not unstudied. A 0 means PubTator annotated 36M
    abstracts and never linked that gene. Coverage of scored proteins: ~95%.
    """
    col = "n_papers_pubtator_prokaryotic"
    if col not in tr.columns:
        say(f"      NOTE {species}: transfer table predates the pubtator donor column -- "
            "omitted. Re-run scripts/studiedness/pubtator.py then transfer.py.")
        return out
    v = pd.to_numeric(tr[col], errors="coerce")
    out[col] = v
    n_val = int(v.notna().sum())
    say(f"      pubtator   {n_val:,} with a value ({100 * n_val / len(out):.1f}%), "
        f"{int((v == 0).sum()):,} measured zero, {int(v.isna().sum()):,} no GeneID donor")
    return out


def spot_check(species: str, out: pd.DataFrame) -> tuple[list[str], list[str]]:
    """Named workhorses must rank high on n_papers_uniprot_prokaryotic."""
    prot = P.load(species)[["uniprot_ac", "gene_name"]]
    d = out.merge(prot, on="uniprot_ac", how="left")
    lines, failed = [], []
    for gene in SPOT_GENES:
        hit = d[d["gene_name"].fillna("").str.lower() == gene.lower()]
        if hit.empty:
            lines.append(f"      {gene:<7} not named in this proteome")
            continue
        row = hit.iloc[0]
        pct = 100 * (d["n_papers_uniprot_prokaryotic"] < row["n_papers_uniprot_prokaryotic"]).mean()
        ok = pct >= SPOT_PERCENTILE
        if not ok:
            failed.append(f"{gene} ({pct:.1f}th)")
        lines.append(f"      {gene:<7} {int(row['n_papers_uniprot_prokaryotic']):>4} papers  "
                     f"{pct:5.1f}th pct"
                     + ("" if ok else f"   <-- BELOW THE {SPOT_PERCENTILE:.0f}th FLOOR"))
    return lines, failed


def unknome_agreement(species: str, out: pd.DataFrame) -> dict | None:
    """Does the GO-based second opinion add anything to the literature-based score?"""
    path = EVIDENCE_DIR / f"unknome_{species}.tsv"
    if not path.exists():
        return None
    u = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    u["unknome_knownness"] = pd.to_numeric(u["unknome_knownness"], errors="coerce")
    d = out.merge(u[["uniprot_ac", "unknome_knownness"]], on="uniprot_ac", how="left")
    both = d[d["unknome_knownness"].notna()]
    if len(both) < 100:
        return None
    rho = both["unknome_knownness"].corr(both["n_papers_uniprot_prokaryotic"], method="spearman")
    # Where Unknome is silent, is the family score also low? If so the two agree even on the
    # proteins Unknome cannot reach, and it adds nothing there either.
    silent = d[d["unknome_knownness"].isna()]
    return {"species": species, "n_both": len(both),
            "spearman_vs_family": round(float(rho), 4),
            "median_family_where_unknome_present": int(both.n_papers_uniprot_prokaryotic.median()),
            "median_family_where_unknome_absent": int(silent.n_papers_uniprot_prokaryotic.median())
            if len(silent) else -1,
            "n_unknome_absent": len(silent)}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STUDIEDNESS - the deliverable")
    rule("=")
    say(f"  in    {(EVIDENCE_DIR / 'transfer_<species>.tsv').relative_to(REPO_ROOT)}")
    say(f"  out   {(OUT_DIR / 'studiedness_<species>.tsv').relative_to(REPO_ROOT)}")
    say("        uniprot_ac · n_papers_uniprot_own · n_papers_uniprot_prokaryotic · n_papers_pubtator_prokaryotic")
    say("  recomputes nothing; reshapes, reindexes to canonical order, and checks")
    say("  n_papers_uniprot_own is near-flat on Kp and Sa BY DESIGN -- rank on n_papers_uniprot_prokaryotic")
    say("  n_papers_uniprot_prokaryotic == 0 is an ANSWER, not a gap -- `no_hit` (nothing resembles it) and")
    say("  `below_floor` (a too-distant relative) are different claims; never impute either")
    rule("=")
    if args.dry_run:
        say("dry run: nothing written.")
        return

    try:
        dfn = S.definition()
        say(f"  number  {dfn['quantity']}")
        say(f"          source {dfn['count_source']}   donor scope {dfn['donor_scope']}   "
            f"floor {dfn['identity_floor']:.0f}% id / {dfn['coverage_floor']:.0f}% cov")
    except FileNotFoundError:
        sys.exit("FATAL evidence/definition.tsv is missing -- "
                 "run scripts/studiedness/transfer.py first")

    rows, unk_rows, failures = [], [], []
    for sp in args.species:
        rule()
        say(f"{sp}")
        rule()
        out, tr = build(sp)
        path = OUT_DIR / f"studiedness_{sp}.tsv"
        out.to_csv(path, sep="\t", index=False)

        # from `tr`, not `out` -- `evidence` is no longer a deliverable column, but the
        # tier breakdown is still the most useful line in the run log.
        tiers = tr["evidence"].value_counts().to_dict()
        gap = float((out["n_papers_uniprot_prokaryotic"] - out["n_papers_uniprot_own"]).median())
        say(f"    {len(out):,} x {out.shape[1]}  ->  {path.relative_to(REPO_ROOT)}")
        say(f"      evidence   " + "   ".join(f"{k} {v:,}" for k, v in sorted(tiers.items())))
        say(f"      n_papers_uniprot_own    median {out.n_papers_uniprot_own.median():.0f}   "
            f"max {out.n_papers_uniprot_own.max():,}   distinct {out.n_papers_uniprot_own.nunique():,}   "
            f"at 0 {int((out.n_papers_uniprot_own == 0).sum()):,}")
        say(f"      n_papers_uniprot_prokaryotic median {out.n_papers_uniprot_prokaryotic.median():.0f}   "
            f"max {out.n_papers_uniprot_prokaryotic.max():,}   distinct {out.n_papers_uniprot_prokaryotic.nunique():,}   "
            f"at 0 {int((out.n_papers_uniprot_prokaryotic == 0).sum()):,}")
        say(f"      family - own, median {gap:+.0f} papers   "
            f"(the transfer's contribution; must be large where the anchor is dark)")

        lines, failed = spot_check(sp, out)
        say(f"      spot checks (must be >= {SPOT_PERCENTILE:.0f}th pct on family):")
        for line in lines:
            say(line)
        if failed:
            failures.append(f"{sp}: {', '.join(failed)}")

        unk = unknome_agreement(sp, out)
        if unk:
            unk_rows.append(unk)
            say(f"      unknome    spearman vs family {unk['spearman_vs_family']:+.3f} "
                f"over {unk['n_both']:,} proteins")
            say(f"                 median family where unknome is present "
                f"{unk['median_family_where_unknome_present']:.3f} / "
                f"absent {unk['median_family_where_unknome_absent']:.3f} "
                f"({unk['n_unknome_absent']:,} proteins)")
        else:
            say("      unknome    not available -- run scripts/studiedness/unknome.py")

        rows.append({"species": sp, "n": len(out),
                     **{f"tier_{k}": v for k, v in tiers.items()},
                     "median_own": int(out.n_papers_uniprot_own.median()),
                     "median_family": int(out.n_papers_uniprot_prokaryotic.median()),
                     "max_family": int(out.n_papers_uniprot_prokaryotic.max()),
                     "distinct_own": int(out.n_papers_uniprot_own.nunique()),
                     "distinct_family": int(out.n_papers_uniprot_prokaryotic.nunique()),
                     "median_family_minus_own": int(gap)})

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(EVIDENCE_DIR / "manifest.tsv", sep="\t", index=False)
    if unk_rows:
        pd.DataFrame(unk_rows).to_csv(EVIDENCE_DIR / "unknome_agreement.tsv",
                                      sep="\t", index=False)

    rule("=")
    say(f"  wrote {(EVIDENCE_DIR / 'manifest.tsv').relative_to(REPO_ROOT)}")
    if unk_rows:
        say(f"  wrote {(EVIDENCE_DIR / 'unknome_agreement.tsv').relative_to(REPO_ROOT)}")
    if failures:
        sys.exit("FATAL spot checks failed -- these families are unambiguously well studied and "
                 "must rank high, so a failure means the transfer joined to the wrong donors:\n  "
                 + "\n  ".join(failures))
    say("  all spot checks passed.")
    say("  studiedness complete.")
    rule("=")


if __name__ == "__main__":
    main()
