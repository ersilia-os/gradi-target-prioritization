"""Functional annotation as two COMPLETE binary matrices -- the stage-02 deliverable.

Every axis in this project must end in a complete matrix: **one row per protein, always**, with an
unannotated protein present as an all-zero row rather than absent. This script converts stage 02's
three tables into exactly two:

    goslim_matrix_<species>.tsv    uniprot_ac + 97 GO-slim term columns + evidence
    cog_matrix_<species>.tsv       uniprot_ac + 26 COG letter columns   + evidence

**It recomputes nothing.** It reads `cog_<species>.tsv` and `goslim_<species>.tsv` and reshapes
them; no predictor runs, no database is needed. Seconds.

WHY BINARY AND MULTI-LABEL
--------------------------
The source tables carry BOTH a single chosen term (`goslim_mf`, `cog_category`) and the full set
(`goslim_mf_all`, `cog_category_all`). **The matrices are built from the full set**, because the
chosen term discards a great deal: measured, **34-41% of MF-annotated and 32-52% of BP-annotated
proteins carry more than one slim term** (max 7), and **12.4-12.6% of COG-classified proteins carry
more than one letter** (max 4). A one-hot of the chosen term would silently throw that away.

WHY THE FULL VOCABULARY, INCLUDING COLUMNS THAT ARE ALWAYS ZERO
---------------------------------------------------------------
All 97 GO-slim terms and all 26 COG letters appear as columns in every species file, so the three
stack without alignment. Some are never used and that is **information, not padding**:
  * 8 of 97 GO-slim terms never occur -- thylakoid, photosynthesis, extracellular matrix, protein
    tag, histone binding, ECM organization, nutrient reservoir, cell adhesion mediator: all
    eukaryote/plant concepts.
  * COG **`Y` (nuclear structure)** never occurs in any bacterium; *S. aureus* also lacks `W`, `Z`.
Dropping them would make the Gram-positive and Gram-negative files different shapes.

**A ZERO MEANS "NOT ANNOTATED", NOT "ABSENT".** For 26.3% of K. pneumoniae nothing is known at all.
`evidence` is what separates the two cases -- read it before treating a zero as a measured negative.

The COG column set comes from the **vendored** `data/source/cdd/cog_func_category.tsv`,
not from COGclassifier's installed package: a schema defined by a pip install is not reproducible.

Run with the `gradi` env.
  python scripts/function/matrix.py
  python scripts/function/matrix.py --species kpneumoniae -q
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import function as F  # noqa: E402
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "function"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
COG_VOCAB = REPO_ROOT / "data" / "source" / "cdd" / "cog_func_category.tsv"
SPECIES = ("kpneumoniae", "ecoli", "saureus")
ASPECTS = ("mf", "bp", "cc")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def cog_vocabulary() -> pd.DataFrame:
    """The 26 COG letters in the authors' own order (grouped), not alphabetical.

    Grouped order keeps related categories adjacent -- information storage, cellular processes,
    metabolism, poorly characterised -- which makes the matrix readable by eye. The file has NO
    header and NO trailing newline.
    """
    if not COG_VOCAB.exists():
        sys.exit(f"FATAL missing {COG_VOCAB}. It is vendored from COGclassifier's package "
                 "resources; see its SOURCE.md.")
    v = pd.read_csv(COG_VOCAB, sep="\t", header=None,
                    names=["letter", "group", "colour", "description"])
    if len(v) != 26 or set(v.letter) != set("ABCDEFGHIJKLMNOPQRSTUVWXYZ"):
        sys.exit(f"FATAL the COG vocabulary is not 26 letters A-Z: got {len(v)} "
                 f"({''.join(sorted(v.letter))})")
    return v


def goslim_vocabulary() -> pd.DataFrame:
    """The 97 GO-slim terms, ordered by aspect then id so the column order is stable."""
    t = F.load_goslim_terms().copy()
    if len(t) != 97:
        say(f"  NOTE the slim vocabulary has {len(t)} terms, not the expected 97")
    order = {"molecular_function": 0, "biological_process": 1, "cellular_component": 2}
    t["_a"] = t["aspect"].map(order).fillna(9)
    return t.sort_values(["_a", "go_id"]).drop(columns="_a").reset_index(drop=True)


def _split(field: str) -> list[str]:
    """Split a `;`-joined term list. Empty and whitespace-only fields give []."""
    return [x.strip() for x in str(field).split(";") if x.strip()]


def build_goslim(species: str, vocab: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per protein, 97 binary term columns, one evidence column."""
    src = F.load_goslim(species)
    terms = list(vocab.go_id)
    idx = {g: i for i, g in enumerate(terms)}
    M = np.zeros((len(src), len(terms)), dtype=np.int8)
    for r, row in enumerate(src.itertuples(index=False)):
        for aspect in ASPECTS:
            for go in _split(getattr(row, f"goslim_{aspect}_all")):
                j = idx.get(go)
                if j is None:
                    sys.exit(f"FATAL {species}: {go} is outside the 97-term slim vocabulary")
                M[r, j] = 1
    out = pd.DataFrame(M, columns=terms)
    out.insert(0, "uniprot_ac", src["uniprot_ac"].to_numpy())
    # `goslim_source` is already curated|eggnog|"" -- normalise the empty case to an explicit label
    out["evidence"] = np.where(src["goslim_source"].to_numpy() == "", "none",
                               src["goslim_source"].to_numpy())
    return out, src


def build_cog(species: str, vocab: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per protein, 26 binary letter columns, one evidence column."""
    src = F.load_cog(species)
    letters = list(vocab.letter)
    idx = {c: i for i, c in enumerate(letters)}
    M = np.zeros((len(src), len(letters)), dtype=np.int8)
    for r, field in enumerate(src["cog_category_all"].to_numpy()):
        for ch in str(field).strip():
            j = idx.get(ch)
            if j is None:
                sys.exit(f"FATAL {species}: COG letter {ch!r} is outside A-Z")
            M[r, j] = 1
    out = pd.DataFrame(M, columns=letters)
    out.insert(0, "uniprot_ac", src["uniprot_ac"].to_numpy())
    out["evidence"] = np.where(src["cog_category"].to_numpy() == "", "none", "cogclassifier")
    return out, src


def verify(species: str, mat: pd.DataFrame, src: pd.DataFrame, kind: str, vocab_cols: list) -> None:
    """Shape, alignment, and a ROUND-TRIP. A shape check passes on a wrong matrix; this does not."""
    prot = P.load(species)["uniprot_ac"].to_numpy()
    if len(mat) != len(prot):
        sys.exit(f"FATAL {kind}/{species}: {len(mat)} rows for {len(prot)} proteins")
    if not np.array_equal(mat["uniprot_ac"].to_numpy(), prot):
        sys.exit(f"FATAL {kind}/{species}: accessions differ from the proteome, or are reordered")
    if mat[vocab_cols].to_numpy().max(initial=0) > 1:
        sys.exit(f"FATAL {kind}/{species}: non-binary values present")

    # ---- the real test: rebuild the source strings from the matrix
    bad = 0
    if kind == "goslim":
        for r in range(len(mat)):
            got = {c for c in vocab_cols if mat.iat[r, mat.columns.get_loc(c)] == 1}
            want = set()
            for aspect in ASPECTS:
                want |= set(_split(src.iloc[r][f"goslim_{aspect}_all"]))
            if got != want:
                bad += 1
    else:
        for r in range(len(mat)):
            got = {c for c in vocab_cols if mat.iat[r, mat.columns.get_loc(c)] == 1}
            want = set(str(src.iloc[r]["cog_category_all"]).strip())
            if got != want:
                bad += 1
    if bad:
        sys.exit(f"FATAL {kind}/{species}: round-trip failed for {bad:,} proteins")
    say(f"      round-trip OK for all {len(mat):,} proteins")


PACKED_COLUMNS = ["uniprot_ac", "cog_categories", "goslim_terms"]


def pack(matrix: pd.DataFrame, vocab_cols: list[str], sep: str) -> list[str]:
    """Collapse a binary matrix back to one `;`-joined term list per row, vocabulary order."""
    arr = matrix[vocab_cols].to_numpy(dtype=bool)
    cols = np.array(vocab_cols, dtype=object)
    return [sep.join(cols[row]) for row in arr]


def build_packed(species: str, mats: dict, vocabs: dict) -> pd.DataFrame:
    """The deliverable: one row per protein, one packed column per scheme.

    **Derived FROM the matrices, not from the source**, so `verify()`'s round-trip against the
    long-form table still governs this table too -- and the reverse direction is asserted below,
    which makes the two forms provably interchangeable rather than merely both present.

    **Neither scheme carries an evidence column here** (owner's call, 2026-10-03). COG never had
    one to carry: `cogclassifier` is strictly 1:1 with "has a category" and `none` with "has not",
    measured on all three species. GO-slim's was real but thin -- `curated` vs `eggnog` separates
    only **322 proteins out of 13,020** (Kp 231 · Ec 75 · Sa 16) -- and it survives byte-identically
    in TWO places that stay on disk and stay audited: `evidence/goslim_matrix_<sp>.tsv` as
    `evidence`, and `evidence/goslim_<sp>.tsv` as `goslim_source`.

    **So the shipped table no longer says whether a GO term was UniProt-curated or inferred from an
    eggNOG orthogroup.** Read one of those two files before treating a term as curated.
    """
    cog, goslim = mats["cog"], mats["goslim"]
    out = pd.DataFrame({
        "uniprot_ac": cog["uniprot_ac"].to_numpy(),
        # COG letters are single characters and the canonical source packs them unseparated;
        # `;` here so one split rule reads both columns.
        "cog_categories": pack(cog, vocabs["cog"], ";"),
        "goslim_terms": pack(goslim, vocabs["goslim"], ";"),
    })

    # The reverse round-trip: re-expanding the packed columns must reproduce the matrices exactly.
    # A shape check passes on a wrong table; this does not.
    for kind, mat, cols, packed in (("cog", cog, vocabs["cog"], out["cog_categories"]),
                                    ("goslim", goslim, vocabs["goslim"], out["goslim_terms"])):
        idx = {t: i for i, t in enumerate(cols)}
        back = np.zeros((len(out), len(cols)), dtype=np.int8)
        for r, field in enumerate(packed):
            for term in _split(field):
                back[r, idx[term]] = 1
        if not np.array_equal(back, mat[cols].to_numpy(dtype=np.int8)):
            sys.exit(f"FATAL {species}/{kind}: packed column does not re-expand to its matrix")
    return out[PACKED_COLUMNS]


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    cvocab, gvocab = cog_vocabulary(), goslim_vocabulary()
    rule("=")
    say("STAGE 02 - functional annotation: one packed table, two matrices beside it")
    rule("=")
    say(f"  goslim   {len(gvocab)} term columns  (mf {int((gvocab.aspect=='molecular_function').sum())} · "
        f"bp {int((gvocab.aspect=='biological_process').sum())} · "
        f"cc {int((gvocab.aspect=='cellular_component').sum())}) + evidence")
    say(f"  cog      {len(cvocab)} letter columns (vendored vocabulary) + evidence")
    say("  built from the *_all columns, so MULTI-LABEL is preserved")
    say("  a zero means NOT ANNOTATED, not absent -- read `evidence`")
    rule("=")
    if args.dry_run:
        say("dry run: nothing written."); return

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp in args.species:
        rule()
        say(f"{sp}")
        rule()
        mats, vocabs = {}, {}
        for kind, builder, vocab_cols in (("goslim", build_goslim, list(gvocab.go_id)),
                                          ("cog", build_cog, list(cvocab.letter))):
            mat, src = builder(sp, gvocab if kind == "goslim" else cvocab)
            verify(sp, mat, src, kind, vocab_cols)
            mats[kind], vocabs[kind] = mat, vocab_cols
            # The matrices are EVIDENCE, not the deliverable (owner's call, 2026-10-03). They are
            # what carries the structural zeros -- a term a species cannot reach is a kept column
            # here and simply absent in the packed table, which cannot tell "impossible" from
            # "unannotated".
            path = EVIDENCE_DIR / f"{kind}_matrix_{sp}.tsv"
            mat.to_csv(path, sep="\t", index=False)
            annotated = int((mat["evidence"] != "none").sum())
            terms_per = mat[vocab_cols].to_numpy().sum(axis=1)
            dead = [c for c in vocab_cols if mat[c].sum() == 0]
            say(f"    {kind:<7} {mat.shape[0]:,} x {mat.shape[1]} -> "
                f"{path.relative_to(REPO_ROOT)}")
            say(f"      annotated {annotated:,} ({100*annotated/len(mat):.1f}%)   "
                f"all-zero rows {int((terms_per==0).sum()):,}   "
                f"mean terms/annotated protein {terms_per[terms_per>0].mean():.2f}")
            say(f"      columns always zero in this species: {len(dead)}")
            rows.append({"species": sp, "matrix": kind, "n": len(mat),
                         "n_columns": len(vocab_cols), "annotated": annotated,
                         "all_zero": int((terms_per == 0).sum()),
                         "dead_columns": len(dead),
                         "evidence": ";".join(f"{k}={v}" for k, v in
                                              mat.evidence.value_counts().items())})

        packed = M.reindex(build_packed(sp, mats, vocabs), sp)
        ppath = OUT_DIR / f"function_{sp}.tsv"
        packed.to_csv(ppath, sep="\t", index=False)
        n_cog = int((packed["cog_categories"] != "").sum())
        n_go = int((packed["goslim_terms"] != "").sum())
        say(f"    packed  {packed.shape[0]:,} x {packed.shape[1]} -> "
            f"{ppath.relative_to(REPO_ROOT)}")
        say(f"      cog_categories {n_cog:,} ({100*n_cog/len(packed):.1f}%)   "
            f"goslim_terms {n_go:,} ({100*n_go/len(packed):.1f}%)   "
            f"both empty {int(((packed.cog_categories=='') & (packed.goslim_terms=='')).sum()):,}")
        say("      re-expands to both matrices exactly (asserted)")

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(EVIDENCE_DIR / "matrix_manifest.tsv", sep="\t", index=False)
    rule("=")
    say(f"  wrote {(EVIDENCE_DIR / 'matrix_manifest.tsv').relative_to(REPO_ROOT)}")
    say("matrices complete.")
    rule("=")


if __name__ == "__main__":
    main()
