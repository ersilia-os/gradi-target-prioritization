"""ESMFold models for the proteins AlphaFold DB does not cover.

    data/processed/pockets/scratch/esmfold/<species>/<acc>.pdb       one model per protein
    data/processed/pockets/evidence/esmfold_<species>.tsv           one row per uncovered protein

WHY 34 PROTEINS HAD NO MODEL (measured 2026-10-03, all AFDB 404s)
------------------------------------------------------------------
* **21 E. coli micro-peptides of 8-15 aa** -- AlphaFold DB does not model anything that short.
* **3 E. coli formate dehydrogenases** (fdhF, fdnG, fdoG) -- each carries one **selenocysteine (U)**.
* **7 E. coli pseudogenes** (mdtQ, efeU, ybfI, yhdW, ycgI, ybfG, ypjI) -- UniProt writes the
  in-frame stop as **X**. Plus xtpA, a 2025 entry newer than AFDB v6.
* **2 giants** -- Kp irp1 (3,163 aa, the yersiniabactin NRPS) and Sa ebh (9,535 aa), above AFDB's
  2,700-residue limit for non-human proteomes.

**No other accession rescues any of them.** UniParc groups identical sequences across accessions;
none of the 34 has an AlphaFold DB model under ANY accession sharing its exact sequence, even where
hundreds do (pheM 1,062, fdhF 126). The exclusions are AFDB's rules, not accession gaps -- and
E. coli has the most because Swiss-Prot curates its leader peptides and pseudogene fragments,
which the Kp and Sa proteomes mostly do not list at all.

WHAT THIS DOES
---------------
Folds each protein up to `MAX_LEN` with **ESMFold v1** (Lin et al., Science 2023): Meta's
`facebook/esmfold_v1` weights through the Hugging Face `transformers` port, already in `gradi`
(transformers 4.48, torch 2.12), so no new env. CPU; the weights are 8.4 GB in the HF cache.

* **U -> C for prediction.** Selenocysteine is cysteine's selenium analogue and no structure
  predictor models it; this is the standard substitution. **X is kept** -- ESMFold has an unknown
  token -- and the model then has **no atoms at that position** (numbering preserved), so the
  residue is a gap: no pLDDT, never part of a pocket. Both are recorded per protein
  (`substitutions`).
* **pLDDT is written to the B-factor column on AlphaFold's 0-100 scale**, so the pLDDT >= 70 pocket
  admission means the same thing on both predictors. (The HF port reports 0-1; it is rescaled and
  the range asserted.)
* **The two giants are NOT folded -- the project owner's decision (2026-10-03).** `MAX_LEN` is
  AFDB's own 2,700. Measured CPU cost: 85 s at 221 aa, 1,232 s at 715 aa (~L^2.3), so the giants
  would need ~18 h even in 1,400-aa fragments. They keep NA pocket columns, labelled in `evidence`.
* **MPS was tested and is not faster**: same structure as CPU (CA RMSD 0.001 A, identical pLDDT)
  in 103 s against 85 s. The default device is CPU.

**This is a different predictor and is labelled as one everywhere.** ESMFold is less accurate than
AlphaFold2 on average, especially without close homologs. `structures.py` records `model_source`
per protein and the deliverable's `evidence` says `esmfold` rather than `af` for these rows.

Run with the `gradi` env, then re-run structures.py -> predict.py -> merge.py:
  python scripts/pockets/esmfold.py
  python scripts/pockets/esmfold.py --species ecoli -q
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
OUT_DIR = TASK_DIR / "scratch" / "esmfold"
MODEL = "facebook/esmfold_v1"
MAX_LEN = 2700   # AlphaFold DB's own limit for non-human proteomes

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def fold(model, tok, seq: str) -> tuple[str, np.ndarray]:
    """(PDB text with 0-100 pLDDT in the B-factors, CA pLDDT)."""
    import torch

    with torch.no_grad():
        ids = tok([seq], return_tensors="pt", add_special_tokens=False)["input_ids"]
        out = model(ids)
    lines = []
    for ln in model.output_to_pdb(out)[0].splitlines():
        if ln.startswith(("ATOM", "HETATM")):
            b = float(ln[60:66])
            if not 0.0 <= b <= 1.0 + 1e-6:
                raise ValueError(f"unexpected pLDDT {b} in ESMFold output -- scale changed?")
            ln = f"{ln[:60]}{100 * b:6.2f}{ln[66:]}"
        lines.append(ln)
    return "\n".join(lines) + "\n", ca_plddt(lines)


def ca_plddt(lines: list[str]) -> np.ndarray:
    return np.array([float(ln[60:66]) for ln in lines
                     if ln.startswith("ATOM") and ln[12:16].strip() == "CA"])


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("--max-len", type=int, default=MAX_LEN,
                    help=f"skip proteins longer than this (default {MAX_LEN}, AFDB's limit)")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say(f"pockets/esmfold.py -- {MODEL} for proteins without an AlphaFold DB model")
    say("  in : data/processed/pockets/evidence/alphafold_<species>.tsv (af_status != model)")
    say("  out: data/processed/pockets/scratch/esmfold/<species>/  +  evidence/esmfold_<species>.tsv")
    say(f"  U->C, X kept; proteins > {args.max_len} aa are skipped (owner's decision)")
    say("=" * 92)

    model = tok = None
    for sp in args.species:
        af = pd.read_csv(EVIDENCE_DIR / f"alphafold_{sp}.tsv", sep="\t")
        if "model_source" in af:
            uncovered = af["model_source"].fillna("") != "alphafold_db_v6"
        else:
            uncovered = af["af_status"] != "model"
        seqs = dict(P.load(sp)[["uniprot_ac", "sequence"]].values)
        todo = af.loc[uncovered, "uniprot_ac"].tolist()
        out_dir = OUT_DIR / sp
        out_dir.mkdir(parents=True, exist_ok=True)
        say(f"\n[{sp}] {len(todo)} proteins without an AlphaFold DB model")
        rows = []
        for acc in todo:
            raw = seqs[acc]
            subs = ";".join(f"{c}{i + 1}{'C' if c == 'U' else c}"
                            for i, c in enumerate(raw) if c not in "ACDEFGHIKLMNPQRSTVWY")
            seq = raw.replace("U", "C")
            row = {"uniprot_ac": acc, "seq_length": len(raw), "substitutions": subs}
            if len(seq) > args.max_len:
                row["status"] = "skipped_too_long"
                rows.append(row)
                say(f"  {acc:<11} {len(raw):>5} aa  skipped (> {args.max_len} aa)")
                continue
            path = out_dir / f"{acc}.pdb"
            t = time.time()
            if path.exists():
                ca = ca_plddt(path.read_text().splitlines())
                cached = True
            else:
                if model is None:
                    import torch
                    from transformers import AutoTokenizer, EsmForProteinFolding

                    torch.set_num_threads(args.threads)
                    tok = AutoTokenizer.from_pretrained(MODEL)
                    model = EsmForProteinFolding.from_pretrained(
                        MODEL, low_cpu_mem_usage=True).eval()
                    model.trunk.set_chunk_size(64)  # bounds pair-representation memory
                    say(f"  model loaded ({time.time() - t:.0f}s)")
                    t = time.time()
                text, ca = fold(model, tok, seq)
                path.write_text(text)
                cached = False
            # ESMFold writes NO atoms for an X (unknown) residue but keeps the numbering, so a
            # pseudogene's model has one CA fewer per X -- measured on mdtQ (X56: 477 CA for 478 aa).
            expected = len(seq) - seq.count("X")
            if len(ca) != expected:
                sys.exit(f"FATAL {acc}: model has {len(ca)} CA atoms, expected {expected} "
                         f"({len(seq)} aa minus {seq.count('X')} X)")
            row.update(status="model", mean_plddt=round(float(ca.mean()), 2),
                       seconds=None if cached else round(time.time() - t, 1))
            rows.append(row)
            say(f"  {acc:<11} {len(raw):>5} aa  pLDDT {row['mean_plddt']:5.1f}  "
                + ("cached" if cached else f"{row['seconds']:.1f}s")
                + (f"  [{subs}]" if subs else ""))
        out = EVIDENCE_DIR / f"esmfold_{sp}.tsv"
        pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
        say(f"  wrote {out.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
