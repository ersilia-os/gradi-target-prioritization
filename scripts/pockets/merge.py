"""The structural ligandability deliverable: one complete, canonical table per species.

    data/processed/pockets/pockets_<species>.tsv

    uniprot_ac · p2rank_score · fpocket_score · n_ligands_pdb · n_ligands_alphafill ·
    n_pdb_structures · af_plddt

**It recomputes nothing.** `structures.py`, `predict.py`, `holo.py` and `pdb_coverage.py` did the
work; this reduces
their evidence tables to one row per protein, enforces canonical order through `src/matrices.py`,
and runs the checks that say whether the axis can be trusted. Seconds, no network.

THE COLUMNS
------------
* `p2rank_score` -- P2Rank's calibrated probability for the best ADMITTED pocket (0-1).
* `fpocket_score` -- fpocket's druggability score for the best ADMITTED pocket (0-1).
* **No pocket COUNT -- removed 2026-10-03 on the owner's instruction, after measuring it.**
  `p2rank_n_pockets` (pockets at probability >= 0.5) needed a cutoff the P2Rank authors never
  recommend (PrankWeb lists every pocket; the papers evaluate by rank). Without one the count is
  protein size: Spearman 0.84 with length, and AUROC against holo evidence WITHIN length deciles
  0.43-0.49 -- no signal at all. A sweep found signal only from 0.3 to 0.5 and never above
  `p2rank_score`'s, so the count added nothing. Every pocket stays in
  `evidence/pocket_list_<sp>.tsv` for anyone who wants to count them.
* `n_ligands_pdb` -- **MEASURED**: non-redundant drug-like ligands in this protein's OWN PDB
  structures (`holo.py`; chains at >= 95% identity, joined to BioLiP with no alignment of our own).
* `n_ligands_alphafill` -- **MODELLED**: non-redundant drug-like ligands AlphaFill transplanted
  onto its AlphaFold model, from donors at ~30% median identity.
  **The two are never summed**: they differ by an order of magnitude in reach and by a great deal
  in strength of evidence. Non-redundant = distinct Bemis-Murcko generic scaffolds; the raw code
  counts ship as `n_codes_*` in `evidence/ligand_counts_<sp>.tsv`.
  They replaced `holo_identity` on 2026-10-03: that column reported the % identity of the closest
  bacterial holo chain found by a DIAMOND search of our own -- a hand-rolled AlphaFill. See
  `alphafill.py` for the measurement that retired it.
* `n_pdb_structures` -- distinct PDB entries with a chain that IS this protein (>= 95% identity,
  >= 50% of the PDB chain aligned), ligand or not (`pdb_coverage.py`). For conserved
  enterobacterial proteins this includes other species' structures of a near-identical protein:
  Kp rpoB counts E. coli's 411 RNA polymerase entries. 0 = no structure.
  A partial structure IS counted: there is no query-coverage floor, so a one-domain construct
  counts (Sa gyrB: 2 entries covering 36% of the protein). Measured: 1-6% of entry matches cover
  < 50% of the protein. Chains aligned over < 30 aa are 15-62 matches per species, and only 3
  proteins per species are matched through them alone -- genuinely tiny proteins with real
  structures (Ec TnaC leader peptide in the ribosome, Sa phenol-soluble modulins).
  How MUCH of the protein they cover (`pdb_coverage`) was dropped from this table on the owner's
  instruction (2026-10-03); it stays per protein in `evidence/pdb_<sp>.tsv`.
* `af_plddt` -- mean AlphaFold pLDDT: whether the pocket columns can be trusted at all.
There is no `evidence` column: read it off the columns that carry it. `af_plddt` is NA exactly
when no model exists, and `n_ligands_pdb > 0` exactly when a drug-like ligand was seen on this
protein -- which is the whole of what `pdb+af` / `af_only` / `pdb_only` / `none` used to say. What
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
   where a model exists, no NA in the ligand counts or `n_pdb_structures` anywhere.
2. **Spot checks -- NAMED LIGANDS, not counts.** A broken join still yields a well-formed table
   of plausible counts, so the check is that specific molecules land on specific genes:
   methotrexate and trimethoprim on E. coli `folA`, novobiocin on S. aureus `gyrB`, and E. coli
   `ftsZ` carrying none -- a real protein with no drug co-crystal of its own.
3. **Does the predicted pocket agree with the measurement?** AUROC of `p2rank_score` and
   `fpocket_score` for separating proteins that HAVE a drug-like ligand in their own PDB
   structures (`n_ligands_pdb > 0`) from those that do not. The two sides share no input -- one is
   geometry on a predicted model, the other is what crystallographers actually found.
   **And `*_len`: the family AUROC WITHIN length deciles** -- added because the raw AUROC is
   confounded by size: big proteins are both crystallised more often and offer more surface for
   pockets, and length ALONE scores 0.66-0.69 against holo evidence. Measured: P2Rank 0.54-0.58
   within length (weak but real), fpocket 0.46-0.51 (nothing). Quote `*_len`, not the raw number.
4. **Reproduces v1?** Proteins with any P2Rank pocket (before admission; AlphaFold DB models only,
   since v1 had no ESMFold) against v1's
   4,542 Kp / 3,589 Ec (`legacy/docs/ligandability_log.md`).

Run with the `gradi` env, LAST: after structures.py -> predict.py -> pdb_coverage.py ->
alphafill.py -> holo.py.
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
from src import pockets as K  # noqa: E402
from src import proteomes as P  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
# What `build()` assembles from the evidence tables...
SOURCE_COLUMNS = ["uniprot_ac", "p2rank_score", "fpocket_score", "n_ligands_pdb",
                  "n_ligands_alphafill", "n_pdb_structures", "af_plddt"]
# ...and what ships: the same, plus the two standard columns DERIVED from them at the end of
# `build()`. The split matters -- slicing with COLUMNS before they are computed raises KeyError.
COLUMNS = SOURCE_COLUMNS + ["pockets_consensus", "pockets_evidence"]
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
# NAMED LIGANDS on named genes. A broken join still produces plausible counts, so the check has
# to be that specific molecules arrive: MTX/TOP = methotrexate and trimethoprim (DHFR), NOV =
# novobiocin (gyrase B), TCL = triclosan (FabI), KHS = a ClpP activator. E. coli ftsZ is the
# negative control -- a real protein with no drug co-crystal of its own.
SPOT_LIGANDS = {
    "ecoli": {"folA": {"MTX", "TOP"}, "ftsZ": set()},
    "saureus": {"gyrB": {"NOV"}, "fabI": {"TCL"}},
    "kpneumoniae": {"clpP": {"KHS"}},
}
V1_P2RANK_ANY = {"kpneumoniae": 4542, "ecoli": 3589}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def auroc_within_length(d: pd.DataFrame, col: str, label: str, bins: int = 10) -> float:
    """AUROC of `col` for `label > 0` inside length deciles, weighted by positives.

    Length alone separates ligand-bearing proteins from the rest at ~0.66-0.69, so a raw AUROC
    mostly measures size; this asks whether the column adds anything once proteins of similar
    length are compared. Quote this, not the raw number."""
    d = d.assign(_bin=pd.qcut(d["seq_length"], bins, labels=False, duplicates="drop"))
    vals, w = [], []
    for _, g in d.groupby("_bin"):
        pos, neg = g.loc[g[label] > 0, col], g.loc[g[label] == 0, col]
        if len(pos) >= 5 and len(neg) >= 5:
            vals.append(auroc(pos.to_numpy(), neg.to_numpy()))
            w.append(len(pos))
    return float(np.average(vals, weights=w)) if vals else float("nan")


def auroc(pos: np.ndarray, neg: np.ndarray) -> float:
    """Mann-Whitney AUROC with ties counted half -- no sklearn needed for one number."""
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    allv = np.concatenate([pos, neg])
    ranks = pd.Series(allv).rank(method="average").to_numpy()
    rp = ranks[: len(pos)].sum()
    return float((rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def ligands_of(species: str, gene: str, prot: pd.DataFrame) -> set[str]:
    """Ligand codes on this gene's own PDB structures, read from holo.py's long table."""
    accs = set(prot.loc[prot["gene_name"] == gene, "uniprot_ac"])
    long = pd.read_csv(EVIDENCE_DIR / f"ligands_pdb_{species}.tsv", sep="\t",
                       keep_default_na=False, na_values=[""])
    return set(long.loc[long["uniprot_ac"].isin(accs), "ligand"])


def build(species: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    for name in (f"alphafold_{species}", f"pocket_list_{species}", f"ligand_counts_{species}",
                 f"pdb_{species}"):
        if not (EVIDENCE_DIR / f"{name}.tsv").exists():
            sys.exit(f"FATAL missing evidence/{name}.tsv -- run structures.py, predict.py, "
                     "pdb_coverage.py, alphafill.py, holo.py first")
    af = pd.read_csv(EVIDENCE_DIR / f"alphafold_{species}.tsv", sep="\t")
    pk = pd.read_csv(EVIDENCE_DIR / f"pocket_list_{species}.tsv", sep="\t")
    lig = pd.read_csv(EVIDENCE_DIR / f"ligand_counts_{species}.tsv", sep="\t")
    pdb = pd.read_csv(EVIDENCE_DIR / f"pdb_{species}.tsv", sep="\t")

    adm = pk[pk["admitted"]]
    p2 = adm[adm["tool"] == "p2rank"].groupby("uniprot_ac")["score"]
    fp = adm[adm["tool"] == "fpocket"].groupby("uniprot_ac")["score"]

    d = af[["uniprot_ac", "af_status", "model_source", "af_plddt", "seq_length"]].copy()
    has_model = d["af_status"] == "model"
    d["p2rank_score"] = d["uniprot_ac"].map(p2.max())
    d["fpocket_score"] = d["uniprot_ac"].map(fp.max())
    # A model with nothing admitted is a measured 0; no model stays NA.
    for c in ("p2rank_score", "fpocket_score"):
        d.loc[has_model, c] = d.loc[has_model, c].fillna(0)
        d.loc[~has_model, c] = np.nan
    d = d.merge(lig[["uniprot_ac", "n_ligands_pdb", "n_ligands_alphafill"]], on="uniprot_ac",
                how="left")
    d = d.merge(pdb[["uniprot_ac", "n_pdb_structures"]], on="uniprot_ac", how="left")

    out = M.reindex(d[SOURCE_COLUMNS + ["seq_length"]], species)   # seq_length: working only
    M.assert_canonical(out["uniprot_ac"], species)
    for c in ("n_ligands_pdb", "n_ligands_alphafill", "n_pdb_structures"):
        if out[c].isna().any():
            sys.exit(f"FATAL {species}: {c} has NA -- every protein was searched")
        out[c] = out[c].astype(int)
    model_rows = out["af_plddt"].notna()
    if out.loc[model_rows, ["p2rank_score", "fpocket_score", "af_plddt"]].isna().any().any():
        sys.exit(f"FATAL {species}: a protein with a model has an NA pocket column")

    # The two standard columns, through `src/pockets.py` so the file and the loader cannot drift.
    out["pockets_consensus"] = K.consensus(out).round(6)
    out["pockets_evidence"] = K.evidence(out)
    return out, pk


def checks(species: str, out: pd.DataFrame, pk: pd.DataFrame) -> tuple[dict, list[str]]:
    failed = []
    prot = P.load(species)[["uniprot_ac", "gene_name"]]
    d = out.merge(prot, on="uniprot_ac")
    for gene, want in SPOT_LIGANDS.get(species, {}).items():
        got = ligands_of(species, gene, prot)
        ok = (want <= got) if want else not got
        say(f"      {gene:<5} own-PDB ligands {len(got):>3}  "
            + (f"expect {sorted(want)}" if want else "expect none")
            + ("" if ok else "   <-- FAIL"))
        if not ok:
            failed.append(f"{gene}: wanted {sorted(want) or 'none'}, got {len(got)}")

    m = d[d["af_plddt"].notna()]
    pos = m[m["n_ligands_pdb"] > 0]
    neg = m[m["n_ligands_pdb"] == 0]
    # v1 only ever had AlphaFold DB models, so the reproduction check counts those alone.
    af_only = pk["model"].eq("alphafold_db_v6") if "model" in pk else True
    any_p2 = pk.loc[(pk["tool"] == "p2rank") & af_only, "uniprot_ac"].nunique()
    res = {
        "species": species, "n": len(out),
        "model": int(len(m)),
        "p2rank>=0.5": int((m["p2rank_score"] >= P2RANK_POCKET).sum()),
        "median_p2rank": float(m["p2rank_score"].median()),
        "median_fpocket": float(m["fpocket_score"].median()),
        "lig_pdb>0": int((out["n_ligands_pdb"] > 0).sum()),
        "lig_af>0": int((out["n_ligands_alphafill"] > 0).sum()),
        "pdb>0": int((out["n_pdb_structures"] > 0).sum()),
        "auroc_p2rank": auroc(pos["p2rank_score"].to_numpy(), neg["p2rank_score"].to_numpy()),
        "auroc_fpocket": auroc(pos["fpocket_score"].to_numpy(), neg["fpocket_score"].to_numpy()),
        "auroc_length": auroc(pos["seq_length"].to_numpy(), neg["seq_length"].to_numpy()),
        "auroc_p2rank_len": auroc_within_length(m, "p2rank_score", "n_ligands_pdb"),
        "auroc_fpocket_len": auroc_within_length(m, "fpocket_score", "n_ligands_pdb"),
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
    say("  in : evidence/{alphafold,pocket_list,ligand_counts,pdb}_<species>.tsv")
    say(f"  out: data/processed/pockets/pockets_<species>.tsv  ({len(COLUMNS)} columns, "
        "complete, canonical)")
    say("=" * 92)

    rows, failed, audit_rows = [], [], []
    for sp in args.species:
        say(f"\n[{sp}]")
        out, pk = build(sp)
        path = TASK_DIR / f"pockets_{sp}.tsv"
        out[COLUMNS].to_csv(path, sep="\t", index=False)   # `evidence` is working-only
        say(f"  wrote {path.relative_to(REPO_ROOT)}  ({len(out):,} rows, canonical)")

        lv = out["pockets_evidence"].value_counts()
        rho = float(out["pockets_consensus"].corr(out["seq_length"], method="spearman"))
        say("  evidence " + " ".join(f"L{k}={int(lv.get(k, 0)):,}" for k in (1, 2, 3))
            + f"   (1 predicted / 2 modelled / 3 measured)")
        # Printed EVERY run, not buried in the docs: this column is substantially protein size.
        say(f"  rho(pockets_consensus, length) = {rho:+.3f}   "
            f"-- the confound; p2rank_score alone is +0.69..+0.74")
        if rho > K.MAX_LENGTH_RHO:
            failed.append(f"consensus-length rho {rho:.3f} > {K.MAX_LENGTH_RHO}")
        audit_rows.append(pd.DataFrame({
            "species": sp, "uniprot_ac": out["uniprot_ac"],
            "n_pdb_structures": out["n_pdb_structures"],
            "n_ligands_pdb": out["n_ligands_pdb"],
            "n_ligands_alphafill": out["n_ligands_alphafill"],
            "af_plddt_isna": out["af_plddt"].isna(),
            "seq_length": out["seq_length"],
            "pockets_consensus": out["pockets_consensus"],
            "pockets_evidence": out["pockets_evidence"],
            "consensus_length_rho": round(rho, 4),
        }))
        res, f = checks(sp, out, pk)
        rows.append(res)
        failed += [f"{sp}: {x}" for x in f]

    if audit_rows:
        aud = pd.concat(audit_rows, ignore_index=True)
        aud.to_csv(EVIDENCE_DIR / "consensus_audit.tsv", sep="\t", index=False)
        say(f"\n  wrote evidence/consensus_audit.tsv  ({len(aud):,} rows)")

    say("\n" + "-" * 92)
    say("SUMMARY")
    say(pd.DataFrame(rows).to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    say("\n  auroc_*: pocket score separating proteins WITH a drug-like ligand in their own PDB\n"
        "  structures from those without. auroc_length: protein length ALONE as the score -- the\n"
        "  confound. *_len: the same within length deciles, the number to quote. "
        "0.5 = no agreement.")
    if failed:
        sys.exit("FATAL spot checks failed: " + "; ".join(failed))


if __name__ == "__main__":
    main()
