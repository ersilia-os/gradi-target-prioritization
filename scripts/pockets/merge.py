"""The structural ligandability deliverable: one complete, canonical table per species.

    data/processed/pockets/pockets_<species>.tsv

    uniprot_ac · p2rank_score · fpocket_score · p2rank_n_pockets · holo_identity ·
    pdb_n_structures · pdb_coverage · af_plddt

**It recomputes nothing.** `structures.py`, `predict.py`, `holo.py` and `pdb_coverage.py` did the
work; this reduces
their evidence tables to one row per protein, enforces canonical order through `src/matrices.py`,
and runs the checks that say whether the axis can be trusted. Seconds, no network.

THE COLUMNS
------------
* `p2rank_score` -- P2Rank's calibrated probability for the best ADMITTED pocket (0-1).
* `fpocket_score` -- fpocket's druggability score for the best ADMITTED pocket (0-1).
* `p2rank_n_pockets` -- admitted P2Rank pockets with probability >= 0.5. Counted with P2Rank, not
  fpocket: fpocket finds at least one pocket on essentially every protein (v1: 5,697 of 5,727
  Kp), so its count carries no signal. fpocket's raw counts stay in `evidence/pockets_<sp>.tsv`.
* `holo_identity` -- % identity (0-100) to the closest bacterial PDB chain with a drug-like ligand
  in the aligned binding site (`holo.py`). 0 = none clears the floors. Sequence-based, so it is
  defined even for a protein with no AlphaFold model.
* `pdb_n_structures` -- distinct PDB entries with a chain that IS this protein (>= 95% identity,
  >= 50% of the PDB chain aligned), ligand or not (`pdb_coverage.py`). For conserved
  enterobacterial proteins this includes other species' structures of a near-identical protein:
  Kp rpoB counts E. coli's 411 RNA polymerase entries. 0 = no structure.
* `pdb_coverage` -- fraction (0-1) of this protein's residues covered by those chains. SEQRES, so
  loops missing from the density still count as covered.
* `af_plddt` -- mean AlphaFold pLDDT: whether the pocket columns can be trusted at all.
There is no `evidence` column: read it off the two columns that carry it. `af_plddt` is NA exactly
when no model exists, and `holo_identity > 0` exactly when a drug-like bacterial co-crystal was
found -- which is the whole of what `pdb+af` / `af_only` / `pdb_only` / `none` used to say. What
those labels ALSO said, and these columns cannot, is whether the model came from AlphaFold DB or
from ESMFold (the 34 proteins AFDB does not cover -- see `esmfold.py`); that is `model_source` in
`evidence/alphafold_<species>.tsv`. `af_plddt` on an ESMFold row is ESMFold's pLDDT, same 0-100
scale, so the column stays comparable either way.

"Admitted" means the pocket's lining residues average pLDDT >= 70 -- the one place model
confidence enters (see `predict.py`).

ZERO IS NOT MISSING
--------------------
With a model, a protein with no admitted pocket gets **0** -- the tools looked and found nothing.
Without a model the three pocket columns and `af_plddt` are **NA**, and that NA is the statement:
"could not look" is a different claim from "looked, found nothing", and filling it with 0 is the v1
mistake (`legacy/HISTORY.md:199`). **So never `fillna(0)` these columns.**

WHAT THIS SCRIPT CHECKS
------------------------
1. **Completeness and order** -- every protein, canonical order, no NA in the pocket columns
   where a model exists, no NA in `holo_identity`, `pdb_n_structures` or `pdb_coverage` anywhere.
2. **Spot checks** -- antibacterial targets with a published drug co-crystal must show
   `holo_identity >= 95` on E. coli (folA/methotrexate, gyrB, rpoB/rifampicin, fabI, murA, lpxC,
   ampC, acrB, def), named genes because a broken join still yields a well-formed table.
3. **Does the predicted pocket agree with the measurement?** AUROC of `p2rank_score` and
   `fpocket_score` for separating proteins with their OWN drug-like co-crystal
   (`holo_identity >= 95`) from those with no holo evidence at all, and again at family level
   (`holo_identity > 0`, `*_fam`), which has more positives. The two columns are computed
   with no shared input, so agreement is recovered, not built in. Reported plainly either way.
4. **Reproduces v1?** Proteins with any P2Rank pocket (before admission) against v1's
   4,542 Kp / 3,589 Ec (`legacy/docs/ligandability_log.md`).

Run with the `gradi` env, after structures.py -> predict.py -> holo.py -> pdb_coverage.py.
  python scripts/pockets/merge.py
  python scripts/pockets/merge.py --species saureus -q
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

SPECIES = ("kpneumoniae", "ecoli", "saureus")
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
COLUMNS = ["uniprot_ac", "p2rank_score", "fpocket_score", "p2rank_n_pockets", "holo_identity",
           "pdb_n_structures", "pdb_coverage", "af_plddt"]
# `evidence` stays a WORKING column -- the completeness checks below read it -- but is not shipped
# (owner's call, 2026-10-03). It is a function of two columns that ARE shipped:
#     af_plddt.notna()      a model exists, so the pocket columns could be computed
#     holo_identity > 0     a drug-like bacterial co-crystal was found for the family
# Verified exactly reconstructible on all three species before removal. The ONE thing it adds is
# the model's provenance -- AlphaFold DB vs ESMFold for the 34 proteins AFDB does not cover -- and
# that belongs to `model_source` in `evidence/alphafold_<sp>.tsv`, which is where to read it.
P2RANK_POCKET = 0.5
EXACT_HOLO = 95.0

# E. coli targets with a drug co-crystal of their own in the PDB. Named genes, not a row count.
SPOT_GENES_EC = ("folA", "gyrB", "rpoB", "fabI", "murA", "lpxC", "ampC", "acrB", "def")
V1_P2RANK_ANY = {"kpneumoniae": 4542, "ecoli": 3589}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUROC with ties counted half -- no sklearn needed for one number."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank(method="average").to_numpy()
    rp = ranks[: len(pos)].sum()
    return float((rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def build(species: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    for name in (f"alphafold_{species}", f"pocket_list_{species}", f"holo_{species}",
                 f"pdb_{species}"):
        if not (EVIDENCE_DIR / f"{name}.tsv").exists():
            sys.exit(f"FATAL missing evidence/{name}.tsv -- run structures.py, predict.py, "
                     "holo.py, pdb_coverage.py first")
    af = pd.read_csv(EVIDENCE_DIR / f"alphafold_{species}.tsv", sep="\t")
    pk = pd.read_csv(EVIDENCE_DIR / f"pocket_list_{species}.tsv", sep="\t")
    ho = pd.read_csv(EVIDENCE_DIR / f"holo_{species}.tsv", sep="\t")
    pdb = pd.read_csv(EVIDENCE_DIR / f"pdb_{species}.tsv", sep="\t")

    adm = pk[pk["admitted"]]
    p2 = adm[adm["tool"] == "p2rank"].groupby("uniprot_ac")["score"]
    fp = adm[adm["tool"] == "fpocket"].groupby("uniprot_ac")["score"]
    n2 = adm[(adm["tool"] == "p2rank") & (adm["score"] >= P2RANK_POCKET)].groupby("uniprot_ac").size()

    d = af[["uniprot_ac", "af_status", "model_source", "af_plddt"]].copy()
    has_model = d["af_status"] == "model"
    d["p2rank_score"] = d["uniprot_ac"].map(p2.max())
    d["fpocket_score"] = d["uniprot_ac"].map(fp.max())
    d["p2rank_n_pockets"] = d["uniprot_ac"].map(n2)
    # A model with nothing admitted is a measured 0; no model stays NA.
    for c in ("p2rank_score", "fpocket_score", "p2rank_n_pockets"):
        d.loc[has_model, c] = d.loc[has_model, c].fillna(0)
        d.loc[~has_model, c] = np.nan
    d["p2rank_n_pockets"] = d["p2rank_n_pockets"].astype("Int64")
    d = d.merge(ho[["uniprot_ac", "holo_identity_bacterial"]], on="uniprot_ac", how="left")
    d = d.rename(columns={"holo_identity_bacterial": "holo_identity"})
    d = d.merge(pdb[["uniprot_ac", "pdb_n_structures", "pdb_coverage"]], on="uniprot_ac",
                how="left")
    holo = d["holo_identity"] > 0
    esm = d["model_source"].fillna("").str.startswith("esmfold")
    d["evidence"] = np.select(
        [has_model & ~esm & holo, has_model & ~esm & ~holo,
         has_model & esm & holo, has_model & esm & ~holo, ~has_model & holo],
        ["pdb+af", "af_only", "pdb+esmfold", "esmfold_only", "pdb_only"], default="none")

    out = M.reindex(d[COLUMNS + ["evidence"]], species)
    M.assert_canonical(out["uniprot_ac"], species)
    for c in ("holo_identity", "pdb_n_structures", "pdb_coverage"):
        if out[c].isna().any():
            sys.exit(f"FATAL {species}: {c} has NA -- every protein was searched")
    out["pdb_n_structures"] = out["pdb_n_structures"].astype(int)
    model_rows = out["evidence"].isin(["pdb+af", "af_only", "pdb+esmfold", "esmfold_only"])
    if out.loc[model_rows, ["p2rank_score", "fpocket_score", "af_plddt"]].isna().any().any():
        sys.exit(f"FATAL {species}: a protein with a model has an NA pocket column")
    return out, pk


def checks(species: str, out: pd.DataFrame, pk: pd.DataFrame) -> tuple[dict, list[str]]:
    failed = []
    prot = P.load(species)[["uniprot_ac", "gene_name"]]
    d = out.merge(prot, on="uniprot_ac")
    if species == "ecoli":
        for g in SPOT_GENES_EC:
            r = d[d["gene_name"] == g]
            v = float(r["holo_identity"].iloc[0]) if len(r) else float("nan")
            ok = v >= EXACT_HOLO
            say(f"      {g:<5} holo {v:5.1f}  p2rank {float(r['p2rank_score'].iloc[0]):.3f}  "
                f"fpocket {float(r['fpocket_score'].iloc[0]):.3f}" + ("" if ok else "   <-- FAIL"))
            if not ok:
                failed.append(f"{g} holo {v:.1f}")

    m = d[d["evidence"].isin(["pdb+af", "af_only", "pdb+esmfold", "esmfold_only"])]
    pos = m[m["holo_identity"] >= EXACT_HOLO]
    fam = m[m["holo_identity"] > 0]
    neg = m[m["holo_identity"] == 0]
    any_p2 = pk.loc[pk["tool"] == "p2rank", "uniprot_ac"].nunique()
    res = {
        "species": species, "n": len(out),
        "model": int(len(m)),
        "p2rank>=0.5": int((m["p2rank_score"] >= P2RANK_POCKET).sum()),
        "median_p2rank": float(m["p2rank_score"].median()),
        "median_fpocket": float(m["fpocket_score"].median()),
        "holo>0": int((out["holo_identity"] > 0).sum()),
        "holo>=95": int((out["holo_identity"] >= EXACT_HOLO).sum()),
        "pdb>0": int((out["pdb_n_structures"] > 0).sum()),
        "pdb_cov>=0.9": int((out["pdb_coverage"] >= 0.9).sum()),
        "auroc_p2rank": auroc(pos["p2rank_score"].to_numpy(), neg["p2rank_score"].to_numpy()),
        "auroc_fpocket": auroc(pos["fpocket_score"].to_numpy(), neg["fpocket_score"].to_numpy()),
        "auroc_p2rank_fam": auroc(fam["p2rank_score"].to_numpy(), neg["p2rank_score"].to_numpy()),
        "auroc_fpocket_fam": auroc(fam["fpocket_score"].to_numpy(),
                                   neg["fpocket_score"].to_numpy()),
        "p2rank_any": any_p2, "v1_p2rank_any": V1_P2RANK_ANY.get(species),
    }
    return res, failed


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say("pockets/merge.py -- pockets_<species>.tsv, the structural ligandability deliverable")
    say("  in : data/processed/pockets/evidence/{alphafold,pockets,holo,pdb}_<species>.tsv")
    say(f"  out: data/processed/pockets/pockets_<species>.tsv  ({len(COLUMNS)} columns, "
        "complete, canonical)")
    say("=" * 92)

    rows, failed = [], []
    for sp in args.species:
        say(f"\n[{sp}]")
        out, pk = build(sp)
        path = TASK_DIR / f"pockets_{sp}.tsv"
        out[COLUMNS].to_csv(path, sep="\t", index=False)   # `evidence` is working-only
        say(f"  wrote {path.relative_to(REPO_ROOT)}  ({len(out):,} rows, canonical)")
        for k, v in out["evidence"].value_counts().items():
            say(f"    {k:<9} {v:>6,}  {100 * v / len(out):5.1f}%")
        res, f = checks(sp, out, pk)
        rows.append(res)
        failed += [f"{sp}: {x}" for x in f]

    say("\n" + "-" * 92)
    say("SUMMARY")
    say(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    say("\n  auroc_*: predicted pocket score, own drug-like co-crystal (holo >= 95) vs no holo "
        "evidence; *_fam: any holo evidence (holo > 0) vs none. 0.5 = no agreement.")
    if failed:
        sys.exit("FATAL spot checks failed: " + "; ".join(failed))


if __name__ == "__main__":
    main()
