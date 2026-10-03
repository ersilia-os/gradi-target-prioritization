"""How studiedness relates to every other axis -- the measurement, not a decision.

    data/processed/studiedness/evidence/confounds.tsv

**It recomputes nothing and fits nothing.** Every axis is read through its own `src/` loader,
aligned on canonical row order, and correlated against `n_papers_uniprot_prokaryotic` and `n_papers_uniprot_own`.
Seconds, no network, no model.

WHY THIS EXISTS
----------------
Nothing consumes studiedness yet, and before anything does, the size of the confound has to be on
the table. `docs/ligands.md` already measured one corner of it: `n_papers_uniprot_prokaryotic` predicts "has a
measurable bacterial ligand" at **AUROC 0.83-0.86**, which was one of the reasons a ligandability
model was rejected -- it would have been a citation count wearing a druggability label. If the
same relationship holds on the other axes, a prioritized shortlist built by stacking them is
partly a ranking of how famous each protein already is, which is the opposite of what GraDi wants.

**This script decides nothing.** It writes the numbers down so the project owner can choose the
role: a filter, a novelty bonus, or a confound to regress out. That is the house pattern --
`ligands/bindingdb.py` measured an overlap and was not promoted; the lazy-qsar comparison was run
and rejected; both are on disk with their numbers.

TWO COLUMNS, NOT ONE, AND THEY MEAN DIFFERENT THINGS
------------------------------------------------------
`n_papers_uniprot_prokaryotic` is the quantity to worry about: it is the real score, and it carries signal for
almost every protein. `n_papers_uniprot_own` is near-flat on Kp (5,710 of 5,728 proteins read exactly 1)
and on Sa, so a correlation against it on those two species is close to meaningless and is
reported only to make that visible. **On E. coli `_own` is the live column** -- it is the only
anchor whose own literature is real -- so Ec is where the two can disagree informatively.

LENGTH IS REPORTED BESIDE EVERY BINARY ENDPOINT, BECAUSE IT IS THE GENERIC CONFOUND
-------------------------------------------------------------------------------------
Protein length already tracks degradability more strongly than the labels do (CLAUDE.md records
-0.61 against -0.33) and scores AUROC 0.66-0.71 on ligand precedent. A studiedness AUROC means
nothing without it alongside: if studiedness and length score the same, studiedness has not been
shown to carry anything length does not.

PR-AUC ACCOMPANIES EVERY AUROC
--------------------------------
House rule, and it binds here hardest of all -- these endpoints have base rates of 2-10%, exactly
where AUROC flatters a ranking and PR-AUC does not. The base rate ships in the same row, because
a PR-AUC is unreadable without it.

AN AUROC BELOW 0.5 IS A DIRECTION, NOT A FAILURE
--------------------------------------------------
Every AUROC here is studiedness predicting the endpoint as written, so `function/goslim_
unannotated` reads **0.153** -- i.e. studiedness predicts *being annotated* at 0.847. It is left
uninverted on purpose: flipping it would hide which way each relationship runs, and the endpoint
names say which direction is which. Read distance from 0.5, then the sign.

THE TIER STRATIFICATION IS THE CONTROL
----------------------------------------
A third of K. pneumoniae scores `n_papers_uniprot_prokaryotic == 0`, split between `no_hit` and `below_floor`.
A correlation driven entirely by that block is a correlation with "did DIAMOND find anything",
not with studiedness -- so every pair is recomputed over the SCORED proteins alone
(`swissprot_direct|close|homolog`). **Read the two rows together**: where `scored_rho` collapses
toward zero, the apparent relationship was the zero block.

Run with the `gradi` env, after merge.py. No process boundary.
  python scripts/studiedness/confounds.py
  python scripts/studiedness/confounds.py --species kpneumoniae -q
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402
from src import essentiality as E  # noqa: E402
from src import function as F  # noqa: E402
from src import ligandability as L  # noqa: E402
from src import localization as LOC  # noqa: E402
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import studiedness as S  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# The tiers where the axis actually produced a number. The other two (`below_floor`, `no_hit`)
# both score 0 -- answers, not gaps, but not a ranking either.
SCORED_TIERS = ("swissprot_direct", "swissprot_close", "swissprot_homolog")

STUDIEDNESS_COLS = ("n_papers_uniprot_prokaryotic", "n_papers_uniprot_own")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 100) -> None:
    say(char * width)


# ---------------------------------------------------------------------------- metrics


def spearman(x: pd.Series, y: pd.Series) -> float | None:
    """Spearman rho, or None where it is undefined (one side constant, or too few rows)."""
    ok = x.notna() & y.notna()
    if ok.sum() < 50:
        return None
    a, b = x[ok], y[ok]
    if a.nunique() < 2 or b.nunique() < 2:
        return None
    return float(a.corr(b, method="spearman"))


def auroc(y_true: np.ndarray, score: np.ndarray) -> float | None:
    """Rank-based AUROC, ties averaged. None where one class is empty."""
    ok = ~pd.isna(score)
    y, s = np.asarray(y_true)[ok], np.asarray(score, dtype=float)[ok]
    n_pos, n_neg = int(y.sum()), int((~y.astype(bool)).sum())
    if n_pos == 0 or n_neg == 0:
        return None
    ranks = pd.Series(s).rank(method="average").to_numpy()
    return float((ranks[y.astype(bool)].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def pr_auc(y_true: np.ndarray, score: np.ndarray) -> float | None:
    """Average precision. None where there are no positives."""
    ok = ~pd.isna(score)
    y, s = np.asarray(y_true)[ok].astype(bool), np.asarray(score, dtype=float)[ok]
    if y.sum() == 0:
        return None
    order = np.argsort(-s, kind="mergesort")
    y = y[order]
    tp = np.cumsum(y)
    precision = tp / np.arange(1, len(y) + 1)
    return float(precision[y].sum() / y.sum())


# ---------------------------------------------------------------------- the endpoints


def endpoints(species: str) -> list[dict]:
    """Every per-protein score the other axes ship, aligned to canonical order.

    Each entry is {axis, endpoint, kind, values}. `kind` is "continuous" (spearman only) or
    "binary" (spearman plus AUROC/PR-AUC against length). A missing axis is skipped with a
    warning rather than failing the run -- this is a measurement script, not a gate.
    """
    out: list[dict] = []
    order = M.canonical(species)

    def align(df: pd.DataFrame, what: str) -> pd.DataFrame | None:
        """Canonical order asserted, never assumed -- a misaligned join reads as a real result."""
        if df is None:
            return None
        df = df.set_index("uniprot_ac").reindex(order)
        M.assert_canonical(pd.Index(order), species, what)
        return df

    def add(axis: str, endpoint: str, kind: str, values) -> None:
        v = pd.Series(np.asarray(values, dtype=float), index=range(len(order)))
        if v.notna().sum() < 50 or v.nunique(dropna=True) < 2:
            say(f"   skip  {axis}/{endpoint}: constant or too few values")
            return
        out.append({"axis": axis, "endpoint": endpoint, "kind": kind, "values": v})

    # -- protein length: the generic confound, and the baseline every AUROC is read against
    try:
        prot = align(P.load(species)[["uniprot_ac", "sequence"]], "proteomes")
        add("proteomes", "length", "continuous", prot["sequence"].fillna("").str.len())
    except Exception as exc:                                        # noqa: BLE001
        say(f"   WARN proteomes: {exc}")

    # -- degradability: the two activator probabilities
    try:
        deg = align(D.load(species), "degradability")
        for act in ("adep4", "onc212"):
            add("degradability", f"{act}_prob", "continuous", deg[f"{act}_prob"])
    except Exception as exc:                                        # noqa: BLE001
        say(f"   WARN degradability: {exc}")

    # -- essentiality: predictions and the one measured column
    try:
        # `deg_ess` and `ogee_ess` are no longer in the summary table -- they ship in
        # `deg_<species>.tsv` and `ogee_<species>.tsv`. Joined back explicitly rather than left to
        # the `if col in` guard below, which would have dropped them from this analysis silently --
        # including the one MEASURED column, on the axis CLAUDE.md calls the unpriced confound.
        ess = align(E.load(species)
                    .merge(E.load_deg(species)[["uniprot_ac", "deg_ess"]], on="uniprot_ac", how="left")
                    .merge(E.load_ogee(species)[["uniprot_ac", "ogee_ess"]], on="uniprot_ac", how="left"),
                    "essentiality")
        # `essentiality` is gone from the table and is NOT reconstructed here: it was
        # where(measured, deg_essential_any, screens_ess_mean), i.e. a function of two columns this
        # loop already prices separately. Correlating it too would double-count them.
        for col in ("geptop_ess", "deg_ess", "ogee_ess", "screens_ess_mean"):
            if col in ess.columns:
                add("essentiality", col, "continuous", pd.to_numeric(ess[col], errors="coerce"))
    except Exception as exc:                                        # noqa: BLE001
        say(f"   WARN essentiality: {exc}")

    # -- ligands: the precedent counts, and the binary endpoint docs/ligands.md already measured
    try:
        prec = align(L.load(species), "ligands/ligands")
        for col in ("n_ligands", "n_ligands_bacterial", "n_assayed_bacterial"):
            add("ligands", col, "continuous", pd.to_numeric(prec[col], errors="coerce"))
        # The reproduction target: docs/ligands.md reports AUROC 0.83-0.86 for n_papers_uniprot_prokaryotic
        # against "has any measurable bacterial ligand". n_measured_* lives in the full table.
        full = align(L.load_full(species), "ligands/ligands_full")
        measured = pd.to_numeric(full["n_measured_bacterial"], errors="coerce")
        add("ligands", "has_measurable_bacterial_ligand", "binary", (measured.fillna(0) > 0))
        potent = pd.to_numeric(prec["n_ligands_bacterial"], errors="coerce")
        add("ligands", "has_potent_bacterial_ligand", "binary", (potent.fillna(0) > 0))
    except Exception as exc:                                        # noqa: BLE001
        say(f"   WARN ligands: {exc}")

    # -- function: annotation darkness. The sanity check -- studiedness MUST predict this, and a
    #    weak number here means the join is wrong, not that the axes are independent.
    try:
        go = align(F.load_goslim_matrix(species), "function/goslim")
        add("function", "goslim_unannotated", "binary", (go["evidence"] == "none"))
        cog = align(F.load_cog_matrix(species), "function/cog")
        add("function", "cog_unclassified", "binary", (cog["evidence"] == "none"))
    except Exception as exc:                                        # noqa: BLE001
        say(f"   WARN function: {exc}")

    # -- localization: one binary per compartment, plus TMbed's fraction
    try:
        loc = align(LOC.load(species), "localization")
        add("localization", "cytoplasmic_fraction", "continuous",
            pd.to_numeric(loc["cytoplasmic_fraction"], errors="coerce"))
        for compartment in sorted(loc["localization"].dropna().unique()):
            add("localization", f"is_{compartment.lower().replace(' ', '_')}", "binary",
                (loc["localization"] == compartment))
    except Exception as exc:                                        # noqa: BLE001
        say(f"   WARN localization: {exc}")

    return out


# ------------------------------------------------------------------------------ rows


def rows_for(species: str, stud: pd.DataFrame, items: list[dict]) -> list[dict]:
    scored_mask = stud["evidence"].isin(SCORED_TIERS).to_numpy()
    length = next((it["values"] for it in items
                   if it["axis"] == "proteomes" and it["endpoint"] == "length"), None)
    rows: list[dict] = []

    for it in items:
        for col in STUDIEDNESS_COLS:
            s = pd.Series(stud[col].to_numpy(dtype=float))
            v = it["values"]
            row = {
                "species": species,
                "axis": it["axis"],
                "endpoint": it["endpoint"],
                "kind": it["kind"],
                "studiedness_column": col,
                "n": int((s.notna() & v.notna()).sum()),
                "rho": spearman(s, v),
                "n_scored": int((scored_mask & s.notna().to_numpy()
                                 & v.notna().to_numpy()).sum()),
                "scored_rho": spearman(s[scored_mask], v[scored_mask]),
            }
            if it["kind"] == "binary":
                y = v.fillna(0).to_numpy().astype(bool)
                row["base_rate"] = round(float(y.mean()), 4)
                row["auroc"] = auroc(y, s.to_numpy())
                row["pr_auc"] = pr_auc(y, s.to_numpy())
                if length is not None:
                    row["auroc_length"] = auroc(y, length.to_numpy())
                    row["pr_auc_length"] = pr_auc(y, length.to_numpy())
            rows.append(row)
    return rows


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--dry-run", action="store_true", help="measure and print, write nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    out_path = EVIDENCE_DIR / "confounds.tsv"
    rule("=")
    say("studiedness/confounds.py -- how studiedness relates to every other axis")
    say(f"   in   data/processed/studiedness/studiedness_<sp>.tsv + every other axis's loader")
    say(f"   out  {out_path.relative_to(REPO_ROOT)}")
    say(f"   species: {', '.join(args.species)}")
    rule("=")

    all_rows: list[dict] = []
    for species in args.species:
        say(f"\n[{species}]")
        stud = S.load(species)
        M.assert_canonical(stud["uniprot_ac"], species, "studiedness")
        n_scored = int(stud["evidence"].isin(SCORED_TIERS).sum())
        say(f"   {len(stud):,} proteins, {n_scored:,} scored "
            f"({100 * n_scored / len(stud):.1f}%); the rest score 0 in two tiers")

        items = endpoints(species)
        say(f"   {len(items)} endpoints from "
            f"{len({it['axis'] for it in items})} axes")
        rows = rows_for(species, stud, items)
        all_rows.extend(rows)

        # the closing table: family only, the column that carries the signal
        fam = [r for r in rows if r["studiedness_column"] == "n_papers_uniprot_prokaryotic"]
        fam.sort(key=lambda r: -abs(r["rho"] or 0))
        say(f"   {'endpoint':<40} {'rho':>7} {'scored':>7} {'AUROC':>7} {'PR':>7} {'(len)':>7}")
        for r in fam:
            def f(x, nd=3):
                return "   --  " if x is None else f"{x:>7.{nd}f}"
            say(f"   {r['axis'] + '/' + r['endpoint']:<40} {f(r['rho'])} "
                f"{f(r['scored_rho'])} {f(r.get('auroc'))} {f(r.get('pr_auc'))} "
                f"{f(r.get('auroc_length'))}")

    df = pd.DataFrame(all_rows)
    for c in ("rho", "scored_rho", "auroc", "pr_auc", "auroc_length", "pr_auc_length"):
        if c in df.columns:
            df[c] = df[c].astype(float).round(4)

    rule()
    if args.dry_run:
        say(f"--dry-run: {len(df):,} rows measured, nothing written")
        return
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    say(f"wrote {out_path.relative_to(REPO_ROOT)}  ({len(df):,} rows)")
    rule("=")


if __name__ == "__main__":
    main()
