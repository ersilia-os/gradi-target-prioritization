"""Export the v2 per-axis tables to the browser payloads of the NEW site (`webapp/`).

Writes   webapp/data/ec.json, webapp/data/kp.json   (columnar: {columns, rows}, one row per protein)
Reads    only the `src/` loaders -- nothing is read by path and nothing is recomputed.

**This does NOT touch the live site.** The deployed app (`app/` on main, `legacy/app/` here) and its
`pages.yml` are left alone; the new site lives in `webapp/` until the owner decides how to publish it.

One row per protein, in the canonical proteome order, keyed on `uniprot_ac`. Every axis is exported
as the axis ships it: its `<axis>_consensus` (within-species 0-1 score, never comparable across
species) beside its `<axis>_evidence` (1-3). **No composite is written** -- the site's weight sliders
combine the consensus columns in the browser only.

Nulls stay null. A missing value is "not annotated / could not look", never 0: K. pneumoniae has no
measured essentiality (DEG indexes no Klebsiella), so its `deg_*` columns are exported empty and the
site says so. Sequences are NOT exported (2 MB per species; the viewer fetches models by accession), but `length` (amino acids) is: degradability scores fall with it, so it is a column in its own right.

Run with the `gradi` env:  python scripts/webapp/export.py [-q] [--species ecoli kpneumoniae]
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src import (  # noqa: E402
    degradability,
    essentiality,
    function,
    ligandability,
    localization,
    orthology,
    pockets,
    projections,
    proteomes,
    studiedness,
)

OUT = REPO_ROOT / "webapp" / "data"
SPECIES = {
    "ecoli": {"file": "ec.json", "prefix": "ec", "organism": "Escherichia coli", "strain": "K-12 MG1655", "proteome_id": "UP000000625"},
    "kpneumoniae": {"file": "kp.json", "prefix": "kp", "organism": "Klebsiella pneumoniae", "strain": "HS11286", "proteome_id": "UP000007841"},
}
# axis -> loader. `sequence` is dropped below.
AXES = [
    ("proteomes", proteomes.load),
    ("function", function.load),
    ("localization", localization.load),
    ("degradability", degradability.load),
    ("essentiality", essentiality.load),
    ("ligands", ligandability.load),
    ("pockets", pockets.load),
    ("orthology", orthology.load),
    ("studiedness", studiedness.load),
    ("projection", projections.load),
    ("deg", essentiality.load_deg),
]
DROP = {"sequence"}
DECIMALS = 4


def build(species: str, say) -> pd.DataFrame:
    base = proteomes.load(species)
    n = len(base)
    order = base["uniprot_ac"].tolist()
    df = base.drop(columns=[c for c in DROP if c in base.columns])
    # protein length in amino acids, from the reference-proteome sequence the whole project keys on
    seq = base["sequence"].astype(str)
    if (seq.str.len() == 0).any() or seq.str.contains(r"[^A-Za-z*]", regex=True).any():
        raise AssertionError(f"{species}: empty or non-alphabetic sequence in the proteome table")
    df.insert(df.columns.get_loc("protein_name") + 1, "length", seq.str.len().astype(int).values)
    for name, loader in AXES[1:]:
        part = loader(species)
        part = part.drop(columns=[c for c in DROP if c in part.columns])
        dup = [c for c in part.columns if c in df.columns and c != "uniprot_ac"]
        if dup:
            part = part.drop(columns=dup)
        before = len(df)
        df = df.merge(part, on="uniprot_ac", how="left", validate="one_to_one")
        if len(df) != before:
            raise AssertionError(f"{species}/{name}: merge changed the row count {before} -> {len(df)}")
        say(f"  + {name:14s} {part.shape[1] - 1:3d} columns")
    if df["uniprot_ac"].tolist() != order or len(df) != n:
        raise AssertionError(f"{species}: rows are not in canonical order")
    # First COG letter = the function filter. A protein with no COG is "none" (not annotated), which is
    # NOT a function class. Multi-letter COGs keep every letter in `cog_categories`; taking the first is a
    # display choice, recorded here and stated on the site.
    first = df["cog_categories"].fillna("").astype(str).str.split(";").str[0].str.strip()
    df.insert(df.columns.get_loc("cog_categories"), "cog", first.where(first != "", "none"))
    # names the browser code reads
    df = df.rename(columns={"uniprot_ac": "uniprot_accession", "gene_name": "gene", "protein_name": "name"})
    return df


def to_payload(df: pd.DataFrame, meta: dict) -> dict:
    cols = list(df.columns)
    rows = []
    for rec in df.itertuples(index=False, name=None):
        row = []
        for v in rec:
            if v is None or (isinstance(v, float) and np.isnan(v)) or v is pd.NA:
                row.append(None)
            elif isinstance(v, (np.floating, float)):
                row.append(round(float(v), DECIMALS))
            elif isinstance(v, (np.integer, np.bool_)):
                row.append(int(v))
            else:
                row.append(v)
        rows.append(row)
    return {
        "organism": meta["organism"],
        "strain": meta["strain"],
        "prefix": meta["prefix"],
        "proteome_id": meta["proteome_id"],
        "generated_at": date.today().isoformat(),
        "n": len(rows),
        "format": "columnar",
        "columns": cols,
        "rows": rows,
    }


def verify(df: pd.DataFrame, path: Path, n_cells: int = 2000, seed: int = 0) -> None:
    """Read the written file back and compare it with the frame it came from: shape, column order,
    per-column non-null counts, and a random sample of cells. A short or reordered payload is the
    silent failure this guards against."""
    back = json.loads(path.read_text(encoding="utf-8"))
    if back["columns"] != list(df.columns) or back["n"] != len(df) or len(back["rows"]) != len(df):
        raise AssertionError(f"{path.name}: shape or column order changed on the way to disk")
    for j, c in enumerate(back["columns"]):
        got = sum(r[j] is not None for r in back["rows"])
        if got != int(df[c].notna().sum()):
            raise AssertionError(f"{path.name}: column {c} has {got} values on disk, {int(df[c].notna().sum())} in the frame")
    rng = np.random.default_rng(seed)
    for _ in range(n_cells):
        i, j = int(rng.integers(len(df))), int(rng.integers(len(df.columns)))
        want, have = df.iloc[i, j], back["rows"][i][j]
        if want is None or (isinstance(want, float) and np.isnan(want)) or want is pd.NA:
            ok = have is None
        elif isinstance(want, (float, np.floating)):
            ok = abs(float(want) - have) <= 0.5 * 10 ** -DECIMALS + 1e-12
        else:
            ok = have == (int(want) if isinstance(want, (np.integer, np.bool_)) else want)
        if not ok:
            raise AssertionError(f"{path.name}: cell ({df.columns[j]}, row {i}) is {have!r} on disk, {want!r} in the frame")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    say = (lambda *_: None) if a.quiet else print

    print("scripts/webapp/export.py: v2 tables -> webapp/data/*.json (the live site is not touched)")
    OUT.mkdir(parents=True, exist_ok=True)
    for sp in a.species:
        meta = SPECIES[sp]
        say(f"\n{sp}")
        df = build(sp, say)
        payload = to_payload(df, meta)
        path = OUT / meta["file"]
        path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
        verify(df, path)
        say(f"  verified: shape, column order, per-column non-null counts and 2,000 random cells match the frame")
        say(f"  wrote {path.relative_to(REPO_ROOT)}: {payload['n']:,} rows x {len(payload['columns'])} columns, {path.stat().st_size / 1e6:.2f} MB")
        say(f"  {'column':34s}{'non-null':>10s}")
        for c in df.columns:
            say(f"  {c:34s}{100 * df[c].notna().mean():9.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
