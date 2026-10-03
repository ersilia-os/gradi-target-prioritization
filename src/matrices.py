"""The canonical row order, and the check that every matrix obeys it.

**Every matrix in this project has the same rows, in the same order.** The order is the one
`data/processed/proteomes/proteome_<species>.tsv` is written in -- stage `proteomes/download.py`
is the unit of analysis, so it defines the index and everything else follows it.

Why an order and not just a set: once it holds, any two axes can be placed side by side with
`np.hstack` / `pd.concat(axis=1)` and no join at all. A join on `uniprot_ac` would give the same
answer, but only if someone remembers to write one -- and the failure mode of forgetting is not an
error, it is a silently misaligned column, which is the worst kind of bug this pipeline can have.
An embedding matrix has no accession column to join on in the first place: its rows are positional,
so alignment there is *only* ever by order.

Use `reindex()` when writing a matrix and `check()` (or `assert_canonical()`) when reading one.

    from src import matrices as M

    df = M.reindex(df, "kpneumoniae")            # write side: order it, and prove it is complete
    M.assert_canonical(df["uniprot_ac"], "kpneumoniae")   # read side: refuse to proceed if not

Completeness is checked in the same breath, because the two rules are one rule: a matrix that is
missing a protein cannot be in canonical order either. A protein an axis could not measure gets a
row of nulls or measured zeros -- never a missing row (CLAUDE.md, *Every axis ends in COMPLETE
matrices*).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from src import proteomes as P

REPO_ROOT = Path(__file__).resolve().parents[1]
SPECIES = ("kpneumoniae", "ecoli", "saureus")
KEY = "uniprot_ac"

_CACHE: dict[str, list[str]] = {}


def canonical(species: str) -> list[str]:
    """The accessions of `species`, in the order every matrix must use."""
    if species not in _CACHE:
        accs = P.load(species)[KEY].tolist()
        if len(set(accs)) != len(accs):
            raise ValueError(f"proteome_{species}.tsv has duplicate {KEY} -- the index is not a key")
        _CACHE[species] = accs
    return _CACHE[species]


def check(accessions, species: str) -> tuple[bool, str]:
    """(ok, message). Reports the FIRST thing wrong, most fundamental first."""
    ref = canonical(species)
    got = [str(a) for a in accessions]
    if got == ref:
        return True, f"canonical ({len(ref)} rows)"

    missing, extra = set(ref) - set(got), set(got) - set(ref)
    if missing or extra:
        bits = []
        if missing:
            bits.append(f"{len(missing)} missing (e.g. {sorted(missing)[:3]})")
        if extra:
            bits.append(f"{len(extra)} not in the proteome (e.g. {sorted(extra)[:3]})")
        return False, f"INCOMPLETE: {'; '.join(bits)}"
    if len(got) != len(ref):
        return False, f"duplicated rows: {len(got)} rows for {len(ref)} proteins"

    pos = next(i for i, (a, b) in enumerate(zip(got, ref)) if a != b)
    return False, f"MISORDERED: same {len(ref)} proteins, first differs at row {pos} ({got[pos]} != {ref[pos]})"


def assert_canonical(accessions, species: str, what: str = "matrix") -> None:
    ok, msg = check(accessions, species)
    if not ok:
        raise ValueError(f"{what} for {species} is not canonical -- {msg}")


def reindex(df: pd.DataFrame, species: str, key: str = KEY) -> pd.DataFrame:
    """Put `df` in canonical order. Raises rather than inventing rows for missing proteins.

    Deliberately NOT a reindex-with-fill: a silently NaN-filled row is a protein the axis never
    measured, dressed as one it measured as unknown. If an axis cannot cover the proteome it must
    say so where the rows are built, not here.
    """
    accs = df[key].astype(str)
    missing = set(canonical(species)) - set(accs)
    if missing:
        raise ValueError(
            f"cannot reindex {species}: {len(missing)} proteins absent from the frame "
            f"(e.g. {sorted(missing)[:3]}). Fill them where the rows are built, not here.")
    out = df.set_index(key).loc[canonical(species)].reset_index()
    out[key] = out[key].astype(str)
    return out


def reindex_arrays(accessions, matrix: np.ndarray, species: str) -> tuple[np.ndarray, np.ndarray]:
    """The same, for an (accessions, matrix) pair -- the shape embeddings are stored in."""
    accs = [str(a) for a in accessions]
    ref = canonical(species)
    missing = set(ref) - set(accs)
    if missing:
        raise ValueError(f"cannot reindex {species}: {len(missing)} proteins absent "
                         f"(e.g. {sorted(missing)[:3]})")
    where = {a: i for i, a in enumerate(accs)}
    take = np.array([where[a] for a in ref])
    return np.array(ref, dtype=object), matrix[take]


def audit(species: tuple[str, ...] = SPECIES) -> pd.DataFrame:
    """Check every matrix this project ships. One row per (species, matrix)."""
    D = REPO_ROOT / "data" / "processed"
    rows = []
    for sp in species:
        tables = [
            ("proteomes/proteome",          D / "proteomes"     / f"proteome_{sp}.tsv"),
            ("function/function",           D / "function"      / f"function_{sp}.tsv"),
            # The matrices are evidence now, not the deliverable -- still audited, because they
            # are what carries the structural zeros the packed table cannot express.
            ("function/goslim_matrix",      D / "function" / "evidence" / f"goslim_matrix_{sp}.tsv"),
            ("function/cog_matrix",         D / "function" / "evidence" / f"cog_matrix_{sp}.tsv"),
            ("localization/localization",   D / "localization"  / f"localization_{sp}.tsv"),
            ("degradability/degradability", D / "degradability" / f"degradability_{sp}.tsv"),
            ("essentiality/essentiality",   D / "essentiality"  / f"essentiality_{sp}.tsv"),
            ("embeddings/projection",       D / "embeddings"    / f"projection_{sp}.tsv"),
            ("studiedness/studiedness",     D / "studiedness"   / f"studiedness_{sp}.tsv"),
            # One file per essentiality EVIDENCE SOURCE, each complete and canonical. Reported
            # ABSENT until its script has run, which is information rather than a failure.
            ("essentiality/deg",            D / "essentiality"  / f"deg_{sp}.tsv"),
            ("essentiality/ogee",           D / "essentiality"  / f"ogee_{sp}.tsv"),
            ("essentiality/screens",        D / "essentiality"  / f"screens_{sp}.tsv"),
            # The ligands axis was canonical but UNAUDITED -- a regression there would not have
            # been caught by "N/N canonical". Added with the precedents table it now sits beside.
            ("ligands/chembl",              D / "ligands" / "evidence" / f"chembl_{sp}.tsv"),
            ("ligands/ligands",             D / "ligands"       / f"ligands_{sp}.tsv"),
        ]
        npzs = [
            ("embeddings/esmc",       D / "embeddings" / f"embeddings_{sp}.npz"),
            ("embeddings/prott5",     D / "embeddings" / f"prott5_{sp}.npz"),
            ("embeddings/proteomelm", D / "embeddings" / f"proteomelm_{sp}.npz"),
            # The same model under the functional encoding it was TRAINED with. A second
            # representation of the same proteins with the same shape, so it must be audited too --
            # an unaudited matrix is precisely the silent misalignment this module exists to catch.
            # Reported as ABSENT until `proteomelm.py --group-embeds orthodb` has been run for a
            # species, which is information rather than a failure.
            ("embeddings/proteomelm_orthodb",
             D / "embeddings" / f"proteomelm_{sp}_orthodb.npz"),
        ]
        for name, path in tables:
            if not path.exists():
                rows.append({"species": sp, "matrix": name, "status": "ABSENT", "detail": str(path)})
                continue
            accs = pd.read_csv(path, sep="\t", usecols=[0], dtype=str).iloc[:, 0]
            ok, msg = check(accs, sp)
            rows.append({"species": sp, "matrix": name, "status": "ok" if ok else "FAIL", "detail": msg})
        for name, path in npzs:
            if not path.exists():
                rows.append({"species": sp, "matrix": name, "status": "ABSENT", "detail": str(path)})
                continue
            z = np.load(path, allow_pickle=True)
            ok, msg = check(z["accessions"].astype(str), sp)
            rows.append({"species": sp, "matrix": name, "status": "ok" if ok else "FAIL", "detail": msg})
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = audit()
    with pd.option_context("display.width", 140, "display.max_colwidth", 70):
        print(df.to_string(index=False))
    bad = df[df.status != "ok"]
    print(f"\n{len(df) - len(bad)}/{len(df)} canonical" + (f"; {len(bad)} to fix" if len(bad) else ""))
