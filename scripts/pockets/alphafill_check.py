"""Does AlphaFill add holo evidence the PDB route misses? -- a MEASUREMENT, not a deliverable.

    data/processed/pockets/evidence/alphafill_comparison.tsv

The approved plan had `holo_identity` come from PDB (BioLiP) UNION AlphaFill. AlphaFill transplants
ligands from PDB-REDO homologs onto AlphaFold models, so it draws on the same experimental
structures BioLiP does; the question is whether its structural (rather than sequence) alignment
reaches proteins the DIAMOND route cannot. This answers it on the two species whose AlphaFill
responses v1 already cached (`data/source/alphafill/{kpneumoniae,ecoli}/`, AlphaFill 2.1.0/2.1.1),
using the SAME drug-like ligand classes (`evidence/ligand_classes.tsv`) and the same 40% floor.

Measured 2026-10-03: AlphaFill reaches **37 Kp / 21 Ec proteins** that the BioLiP route misses
even with no taxonomy restriction (61 / 42 against the bacterial-only column), out of 668 / 611.
The ligands behind that gain are mostly small molecules and crystallisation additives -- guanidine
(GAI) leads in both species, then ethylmercurithiosalicylate (EMT, a heavy-atom soak) -- with a few
real paromomycin (PAR) contacts on ribosomal proteins. So
AlphaFill is NOT part of the deliverable, and S. aureus was not fetched (~2,900 API calls for an
expected ~5% gain, most of it noise). AlphaFill also gives no donor taxonomy without a second
lookup, so its hits could not be held to the bacterial rule.

Codes AlphaFill transplants that BioLiP never lists (~110 per species) are counted and treated as
not drug-like -- they have no classification, and an unclassified ligand is not evidence.

Run with the `gradi` env, after `holo.py`. No network.
  python scripts/pockets/alphafill_check.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli")      # the two with cached AlphaFill responses
CACHE = REPO_ROOT / "data" / "source" / "alphafill"
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "pockets" / "evidence"
MIN_PIDENT = L.REMOTE_PIDENT

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def best_transplant(path: Path, druglike: set[str], known: set[str]) -> tuple[float, str, str, set]:
    """(best identity, ligand, pdb, unclassified codes) over drug-like transplants >= the floor."""
    if not path.exists():
        return float("nan"), "", "", set()
    d = json.loads(path.read_text()) or {}
    best, lig, pdb, unknown = 0.0, "", "", set()
    for hit in d.get("hits") or []:
        idt = 100 * float(hit["alignment"]["identity"])
        for t in hit.get("transplants") or []:
            c = t["compound_id"]
            if c not in known:
                unknown.add(c)
            if c in druglike and idt >= MIN_PIDENT and idt > best:
                best, lig, pdb = idt, c, hit["pdb_id"]
    return best, lig, pdb, unknown


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    VERBOSE = not ap.parse_args().quiet

    say("=" * 92)
    say("pockets/alphafill_check.py -- AlphaFill drug-like transplants vs the BioLiP holo route")
    say("  in : data/source/alphafill/<species>/*.json, evidence/{ligand_classes,holo_<sp>}.tsv")
    say("  out: data/processed/pockets/evidence/alphafill_comparison.tsv")
    say("=" * 92)

    cl = pd.read_csv(EVIDENCE_DIR / "ligand_classes.tsv", sep="\t", low_memory=False)
    druglike, known = set(cl.loc[cl["druglike"], "ligand"]), set(cl["ligand"])
    rows = []
    for sp in SPECIES:
        h = pd.read_csv(EVIDENCE_DIR / f"holo_{sp}.tsv", sep="\t")
        unknown_all: set[str] = set()
        n_resp = 0
        for acc, b_any, b_bact in zip(h["uniprot_ac"], h["holo_identity_any"],
                                      h["holo_identity_bacterial"]):
            best, lig, pdb, unk = best_transplant(CACHE / sp / f"{acc}.json", druglike, known)
            unknown_all |= unk
            n_resp += best == best  # not NaN
            rows.append({"species": sp, "uniprot_ac": acc, "alphafill_identity": best,
                         "alphafill_ligand": lig, "alphafill_pdb": pdb,
                         "holo_identity_any": b_any, "holo_identity_bacterial": b_bact})
        d = pd.DataFrame([r for r in rows if r["species"] == sp])
        af = d["alphafill_identity"].fillna(0) > 0
        say(f"\n[{sp}] {n_resp:,}/{len(d):,} cached AlphaFill responses; "
            f"{len(unknown_all):,} transplanted codes absent from BioLiP (treated as not drug-like)")
        say(f"  AlphaFill drug-like >= {MIN_PIDENT:.0f}%     {int(af.sum()):>5,}")
        say(f"  BioLiP route, any organism       {int((d['holo_identity_any'] > 0).sum()):>5,}")
        say(f"  BioLiP route, bacterial          {int((d['holo_identity_bacterial'] > 0).sum()):>5,}")
        say(f"  AlphaFill-only vs any organism   {int((af & (d['holo_identity_any'] == 0)).sum()):>5,}")
        say(f"  AlphaFill-only vs bacterial      {int((af & (d['holo_identity_bacterial'] == 0)).sum()):>5,}")
        only = d[af & (d["holo_identity_any"] == 0)]
        top = only["alphafill_ligand"].value_counts().head(8)
        say("  ligands behind the AlphaFill-only gain: "
            + ", ".join(f"{k}({v})" for k, v in top.items()))

    out = EVIDENCE_DIR / "alphafill_comparison.tsv"
    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    say(f"\nwrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
