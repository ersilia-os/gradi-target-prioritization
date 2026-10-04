"""Ligand precedent for any protein sequence -- three counts, in about a second.

    n_ligands_own / n_ligands_bacterial / n_ligands_human   POTENT, pChEMBL >= 6 (1 uM)
    n_assayed_own / n_assayed_bacterial / n_assayed_human   compounds tried at all
    (c) n_ligands_human      ligands on HUMAN orthologs

A QUERY TOOL over arbitrary input -- which is what separates it from `ligands/chembl.py`. That
script answers "what does our proteome have" and wants the 30.5 GB ChEMBL dump; this answers "what
about this sequence", for any sequence, from 82 MB of cached extracts.

(b) COUNTS MOLECULES, NOT TARGET-COMPOUND PAIRS. A compound tested against three homologous targets
counts ONCE -- the union of distinct `parent_molregno` over every target that passes the floors.
Summing per-target counts would inflate it, worse the wider the identity band.

(c) IS A LIABILITY, NOT A PRECEDENT. Never add it to (b). Measured on E. coli `clpP`: 136 bacterial
against **210 human** at 56.3% identity -- the human mitochondrial CLPP ortholog. That is a reason
to be careful with a target, not a reason to like it.

TWO THINGS THE NUMBERS DO NOT SAY:

  * `pchembl` is 100% populated in the extract, so these are POTENCY-MEASURABLE ligands only --
    `=` relations on IC50/EC50/Ki/Kd/Potency in nM. MIC and %-inhibition are absent BY
    CONSTRUCTION, so this tool does not say "has an antibiotic", and a ribosomal protein looking
    empty is a fact about assay type, not biology.
  * `n_ligands_bacteria_complex` is reported SEPARATELY and is not in (b). DNA gyrase is a
    `PROTEIN COMPLEX` in ChEMBL -- measured here, E. coli gyrB carries 666 single-protein ligands
    and **1,412 complex** ones. v1 dropped the complex track and made GyrA/GyrB look unliganded.

Run with the `gradi` env (DIAMOND from `gradi-ortho`; `GRADI_DIAMOND_BIN` overrides):
    python scripts/ligands/ligands.py --accession P0ABQ4
    python scripts/ligands/ligands.py --sequence MKTAYIAKQR...
    python scripts/ligands/ligands.py --fasta my_proteins.faa
    python scripts/ligands/ligands.py --species kpneumoniae      # -> ligands_<sp>.tsv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import matrices as M  # noqa: E402
from src import precedents as PR  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "ligands"
SPECIES = ("kpneumoniae", "ecoli", "saureus")

# THE DELIVERABLE IS FOUR COLUMNS PLUS THE KEY -- the project owner's choice. The full 13-column
# version is not discarded, it goes to evidence/precedents_full_<sp>.tsv, because two of the
# dropped columns answer questions a bare count cannot:
#
#   best_pident_bacterial         100% on KPC-2 and 87.9% on nfo are different claims, and the
#                                count alone does not distinguish "this protein" from "something
#                                that resembles it"
#   n_ligands_bacteria_complex   E. coli gyrB carries 666 single-protein ligands and 1,412 complex
#                                ones; Kp A0A0H3H0Y6 carries 1,410 complex and ZERO single. Absent
#                                from the root table, gyrase reads as unliganded -- which is
#                                exactly the v1 error this column was added to prevent.
#
# So the root table is the clean answer and evidence/ keeps the caveats reachable.
# The four counts the project owner asked for, plus the EVIDENCE label. The label is not an
# extra: CLAUDE.md's standing rule is that every axis ships one, because a 0 here is read by
# ~96% of every proteome and until now it meant three different things -- no homolog at all,
# a homolog nobody ever screened, and a homolog somebody screened that yielded nothing.
# The two standard columns are appended here rather than in `src/precedents.py`: that module is a
# per-sequence query tool with no species, and a within-species percentile rank is meaningless for
# one protein. They exist only in the batch deliverable.
DELIVERABLE_COLUMNS = ("uniprot_ac", *PR.DELIVERABLE,
                       "ligands_consensus", "ligands_evidence")

VERBOSE = True


def say(m: str = "") -> None:
    if VERBOSE:
        print(m, flush=True)


def rule(c: str = "-", w: int = 100) -> None:
    say(c * w)


def read_fasta(p: Path) -> dict[str, str]:
    seqs, cur, buf = {}, None, []
    for line in p.read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seqs[cur] = "".join(buf)
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


def control_row(sp: str, out: pd.DataFrame) -> dict:
    """Agreement with `chembl_<sp>.tsv` AT MATCHED SEMANTICS -- the check that validates the tool.

    The two disagree by default, and both reasons are deliberate rather than defects:

      * (b) EXCLUDES the complex track, which `remote_n_compounds` includes. Measured on Kp
        `A0A0H3H0Y6`: 1,410 complex-track ligands and ZERO single-protein ones -- a gyrase subunit.
      * the default counts every potency-measurable ligand, while `chembl.py` reports the
        pChEMBL >= 6 headline. 66 of 67 default-only proteins sit below 6.

    Put both back and they reconcile: measured 113 both-positive, 0 chembl-only, 99.9% agreement
    over all 5,728 Kp proteins. That is what this row records, so a future divergence is visible.
    """
    import numpy as np
    c = pd.read_csv(OUT_DIR / "evidence" / f"chembl_{sp}.tsv", sep="\t")
    o = PR.count(dict(zip(out.uniprot_ac, P.load(sp).set_index("uniprot_ac")
                          .loc[out.uniprot_ac, "sequence"].astype(str))),
                 {k: k for k in out.uniprot_ac},
                 organisms=PR.L.SPECIES_ORGANISM.get(sp), min_pchembl=6.0)
    o = o.rename(columns={"id": "uniprot_ac"})
    m = o.merge(c[["uniprot_ac", "remote_n_compounds"]], on="uniprot_ac")
    # SINGLE TRACK ONLY on both sides. `chembl.py:452` defines the remote pool as
    # `single & bact`, so `remote_n_compounds` excludes complexes -- adding ours to the comparison
    # and calling it "matched semantics" inflated the delta to +1,784 and hid the real residual.
    pb = m.n_ligands_bacterial > 0
    cb = m.remote_n_compounds.fillna(0) > 0
    # MAGNITUDE, NOT JUST PRESENCE. A boolean-only control reported 99.93% agreement while the
    # counts actually disagreed on 18 of 309 ligand-bearing Kp proteins by 1,019 ligands -- it was
    # structurally incapable of seeing the component_id->tid collapse that caused them. Comparing
    # totals would have caught it on the first run.
    mine = m.n_ligands_bacterial.fillna(0)
    theirs = m.remote_n_compounds.fillna(0)
    both = pb & cb
    delta = (mine[both] - theirs[both])
    return {"species": sp, "n": len(m),
            "both_positive": int(both.sum()),
            "precedents_only": int((pb & ~cb).sum()),
            "chembl_only": int((~pb & cb).sum()),
            "agreement": round(float((pb == cb).mean()), 4),
            "n_count_mismatch": int((delta != 0).sum()),
            "total_ligand_delta": int(delta.sum()),
            "max_abs_delta": int(delta.abs().max()) if len(delta) else 0,
            "matched_semantics": "pchembl>=6, SINGLE track both sides (chembl.py:452)"}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--sequence", help="one protein sequence")
    src.add_argument("--fasta", type=Path, help="a protein FASTA")
    src.add_argument("--accession", help="a UniProt accession from one of our four proteomes")
    ap.add_argument("--organism", default=None, metavar="NAME",
                    help="organism of the input (e.g. 'Klebsiella pneumoniae'). Enables the "
                         "species route of the exact count: same species at >=95%% identity is "
                         "THIS protein, not a homolog. --species sets it automatically.")
    src.add_argument("--species", nargs="+", choices=list(SPECIES),
                     help="score a whole anchor proteome -> ligands_<species>.tsv")
    ap.add_argument("--min-identity", type=float, default=PR.DEFAULT_MIN_IDENTITY,
                    help="identity floor for (b) and (c). Default 40, the house transfer floor; "
                         "95 is 'essentially this protein'.")
    ap.add_argument("--min-pchembl", type=float,
                    help="keep only ligands at or above this potency (6 is the house headline)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("-o", "--out", type=Path, help="write the table here instead of printing")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("precedents.py -- ligand precedent for a sequence: exact / bacterial / human")
    rule("=")
    say(f"  identity floor {a.min_identity:.0f}%   coverage floors q>={PR.L.MIN_QCOV:.0f} "
        f"s>={PR.L.MIN_SCOV:.0f} (imported, not restated)")
    say(f"  n_ligands_*     POTENT only: pChEMBL >= {PR.L.PCHEMBL_HEADLINE:g} (1 uM)")
    say("  n_assayed_*     compounds tried at all, whatever the outcome")
    say("  every count is DISTINCT MOLECULES over the union of the pool, never per-target sums")
    say("  *_bacterial INCLUDES the exact match -- nested, so never sum the two")
    say("  *_human is a LIABILITY -- never add it to the bacterial count")
    rule()

    if a.species:
        rows = []
        audit_rows = []
        for sp in a.species:
            d = P.load(sp)
            seqs = dict(zip(d["uniprot_ac"], d["sequence"].astype(str)))
            accs = {k: k for k in seqs}
            say(f"  {sp}: {len(seqs):,} proteins ...")
            out = PR.count(seqs, accs, organisms=PR.L.SPECIES_ORGANISM.get(sp),
                           min_identity=a.min_identity, min_pchembl=a.min_pchembl,
                           threads=a.threads)
            out = out.rename(columns={"id": "uniprot_ac"})
            out = M.reindex(out, sp)
            ev = OUT_DIR / "evidence"
            ev.mkdir(parents=True, exist_ok=True)
            out.to_csv(ev / f"precedents_full_{sp}.tsv", sep="\t", index=False)
            # Computed through `src/ligandability.py` so the file and the loader cannot drift.
            # The consensus ranks WITHIN this species, which is why it belongs here and not in
            # `precedents.count()`.
            out["ligands_consensus"] = L.consensus(out).round(6)
            out["ligands_evidence"] = L.evidence(out)
            p = OUT_DIR / f"ligands_{sp}.tsv"
            out[list(DELIVERABLE_COLUMNS)].to_csv(p, sep="\t", index=False)
            audit_rows.append(pd.DataFrame({
                "species": sp, "uniprot_ac": out["uniprot_ac"],
                "n_assayed_own": out["n_assayed_own"],
                "n_assayed_bacterial": out["n_assayed_bacterial"],
                "has_any_evidence": out["ligands_consensus"] > 0,
                "ligands_consensus": out["ligands_consensus"],
                "ligands_evidence": out["ligands_evidence"],
            }))
            lv = out["ligands_evidence"].value_counts()
            say(f"    evidence " + " ".join(f"L{k}={int(lv.get(k, 0)):,}" for k in (1, 2, 3))
                + f"   consensus: {int((out.ligands_consensus > 0).sum()):,} non-zero "
                f"({100 * (out.ligands_consensus > 0).mean():.1f}%)")
            l3 = out.loc[out.ligands_evidence == 3, "uniprot_ac"]
            if len(l3) <= 25:
                names = P.load(sp).set_index("uniprot_ac")["gene_name"]
                say(f"    assayed on themselves: "
                    + ", ".join(sorted(x for x in l3.map(names).fillna("?")) ) )
            say(f"    potent: exact {int((out.n_ligands_own > 0).sum()):5,}   "
                f"bacterial {int((out.n_ligands_bacterial > 0).sum()):5,}   "
                f"human {int((out.n_ligands_human > 0).sum()):5,}   |   "
                f"assayed: bacterial {int((out.n_assayed_bacterial.fillna(0) > 0).sum()):5,}")
            say(f"    -> {p.relative_to(REPO_ROOT)}  "
                f"({out.shape[0]} x {len(DELIVERABLE_COLUMNS)})   "
                f"full {out.shape[1]} cols -> evidence/precedents_full_{sp}.tsv")
            rows.append(control_row(sp, out))
        if audit_rows:
            ev = OUT_DIR / "evidence"
            ev.mkdir(parents=True, exist_ok=True)
            aud = pd.concat(audit_rows, ignore_index=True)
            aud.to_csv(ev / "consensus_audit.tsv", sep="\t", index=False)
            say(f"  wrote evidence/consensus_audit.tsv  ({len(aud):,} rows)")
        if rows:
            ev = OUT_DIR / "evidence"
            ev.mkdir(parents=True, exist_ok=True)
            cp = ev / "precedent_control.tsv"
            pd.DataFrame(rows).to_csv(cp, sep="\t", index=False)
            rule()
            say("CONTROL  -- agreement with the existing chembl axis at MATCHED semantics")
            rule()
            for r in rows:
                say(f"  {r['species']:12s} both {r['both_positive']:4d}  "
                    f"precedents-only {r['precedents_only']:3d}  chembl-only {r['chembl_only']:3d}"
                    f"   agreement {r['agreement']:.1%}")
            say(f"  wrote {cp.relative_to(REPO_ROOT)}")
        rule("=")
        return

    if a.accession:
        out = PR.for_accession(a.accession)
    else:
        seqs = read_fasta(a.fasta) if a.fasta else {"query": a.sequence.strip()}
        if not seqs:
            sys.exit(f"FATAL no sequences read from {a.fasta}")
        out = PR.count(seqs, None, organisms=a.organism,
                       min_identity=a.min_identity, min_pchembl=a.min_pchembl,
                       threads=a.threads)

    if a.out:
        out.to_csv(a.out, sep="\t", index=False)
        say(f"  wrote {a.out}")
    else:
        rule()
        for _, r in out.iterrows():
            say(f"  {r['id']}")
            def fmt(v):
                return "n/a" if pd.isna(v) else f"{int(v):,}"
            say(f"    this protein   {fmt(r['n_ligands_own']):>9s} potent / "
                f"{fmt(r['n_assayed_own']):>9s} assayed   via {r['exact_route']}"
                + (f", target {r['exact_target']}" if pd.notna(r["exact_target"]) else ""))
            say(f"    bacterial      {fmt(r['n_ligands_bacterial']):>9s} potent / "
                f"{fmt(r['n_assayed_bacterial']):>9s} assayed   "
                f"{r['n_targets_bacterial']} targets, best id {r['best_pident_bacterial']}%, "
                f"best pActivity {r['best_pactivity_bacterial']}   (INCLUDES this protein)")
            say(f"    human          {fmt(r['n_ligands_human']):>9s} potent / "
                f"{fmt(r['n_assayed_human']):>9s} assayed   "
                f"{r['n_targets_human']} targets, best id {r['best_pident_human']}%"
                + ("   <- LIABILITY: more human than bacterial"
                   if (pd.notna(r["n_ligands_human"]) and pd.notna(r["n_ligands_bacterial"])
                       and r["n_ligands_human"] > r["n_ligands_bacterial"]) else ""))
            if r["n_ligands_bacterial_complex"]:
                say(f"        + {int(r['n_ligands_bacterial_complex']):,} on bacterial PROTEIN "
                    "COMPLEXES, reported apart")
    rule("=")


if __name__ == "__main__":
    main()
