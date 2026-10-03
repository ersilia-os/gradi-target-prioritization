"""AlphaFold DB models for the three bacterial proteomes -- the structures every pocket sits on.

    data/source/alphafold/<species>/AF-<acc>-F1-model_v6.cif       (cache, by provider)
    data/processed/pockets/evidence/alphafold_<species>.tsv           (one row per protein)

The evidence table is COMPLETE: every protein of the reference proteome has a row, and a protein
with no model gets `af_status = no_model` rather than a missing row. Columns:

    uniprot_ac · af_status · model_source · af_version · af_length · seq_length · af_plddt ·
    af_frac_ordered · af_frac_disordered

`model_source` is `alphafold_db_v6`, or -- for the 34 proteins AlphaFold DB does not cover --
`esmfold_v1` when `esmfold.py` has folded them (the two giants above 2,700 aa stay `no_model`). The `af_` prefix is kept
for the column names because the deliverable's `af_plddt` was named before ESMFold joined; on those
rows it is ESMFold's pLDDT, on the same 0-100 scale. ESMFold models are validated the same way:
their sequence must equal the proteome's after the recorded U -> C substitution.

`af_plddt` is the mean per-residue pLDDT read from the model itself (the B-factor column of the
CA atoms), not the API's summary field, so the number and the structure fed to the pocket tools
are guaranteed to be the same object. `af_frac_ordered` is the fraction of residues at pLDDT >= 70
(AlphaFold's own "confident" band) and `af_frac_disordered` the fraction below 50 ("very low").

THE MODEL MUST BE THIS PROTEIN. A cif is accepted only if its sequence equals the proteome's
sequence exactly; a mismatch is recorded as `sequence_mismatch` and the model is not used. An
accession that has been re-sequenced since AlphaFold DB was built would otherwise feed a pocket
predictor the wrong protein, and nothing downstream could tell.

v1 cached models for K. pneumoniae and E. coli (AFDB v6) in `data/processed/legacy/v1/`; they were
seeded into the cache and are re-validated here like any download. AFDB is still v6 (checked
2026-10-03), so they are current.

Run with the `gradi` env:
  python scripts/pockets/structures.py
  python scripts/pockets/structures.py --species saureus --workers 8 -q
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
AF_VERSION = 6
SOURCE_DIR = REPO_ROOT / "data" / "source" / "alphafold"
EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "pockets" / "evidence"
API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
FILE = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{v}.cif"
ESMFOLD_DIR = REPO_ROOT / "data" / "processed" / "pockets" / "scratch" / "esmfold"

PLDDT_ORDERED = 70.0     # AlphaFold's "confident" band starts here
PLDDT_DISORDERED = 50.0  # "very low" below this

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V", "SEC": "U",
    "PYL": "O", "UNK": "X",
}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def cif_path(species: str, acc: str) -> Path:
    return SOURCE_DIR / species / f"AF-{acc}-F1-model_v{AF_VERSION}.cif"


def esmfold_path(species: str, acc: str) -> Path:
    return ESMFOLD_DIR / species / f"{acc}.pdb"


def read_pdb_ca(path: Path) -> tuple[str, np.ndarray]:
    """(sequence, per-residue pLDDT) from the CA atoms of a PDB file (ESMFold, 0-100 B-factors).

    Indexed BY RESIDUE NUMBER: ESMFold writes no atoms for an X residue but keeps the numbering,
    so a gap is filled with 'X' and pLDDT NaN -- which is what makes the sequence check against
    the proteome (which has the X) pass, and keeps pocket residue numbers aligned."""
    by_num: dict[int, tuple[str, float]] = {}
    for ln in path.read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            by_num[int(ln[22:26])] = (THREE_TO_ONE.get(ln[17:20].strip(), "X"), float(ln[60:66]))
    n = max(by_num) if by_num else 0
    seq = "".join(by_num.get(i, ("X", np.nan))[0] for i in range(1, n + 1))
    pl = np.array([by_num.get(i, ("X", np.nan))[1] for i in range(1, n + 1)], dtype=float)
    return seq, pl


def model_file(species: str, acc: str) -> Path | None:
    """The structure pocket tools should see: the AFDB cif, else the ESMFold PDB, else None."""
    for p in (cif_path(species, acc), esmfold_path(species, acc)):
        if p.exists():
            return p
    return None


def model_ca(species: str, acc: str) -> tuple[str, np.ndarray, str]:
    """(model sequence, per-residue pLDDT, source). FileNotFoundError when there is no model."""
    p = model_file(species, acc)
    if p is None:
        raise FileNotFoundError(acc)
    if p.suffix == ".cif":
        seq, pl = read_ca(p)
        return seq, pl, "alphafold_db_v6"
    seq, pl = read_pdb_ca(p)
    return seq, pl, "esmfold_v1"


def read_ca(cif: Path) -> tuple[str, np.ndarray]:
    """(sequence, per-residue pLDDT) from the CA atoms of an AlphaFold cif.

    A minimal reader of `_atom_site` rather than biotite, because this runs over 13,000 files and
    only needs two columns -- and it makes the pLDDT provenance explicit (B_iso_or_equiv of CA)."""
    cols: list[str] = []
    seq, plddt = [], []
    in_loop = False
    with open(cif) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols.append(line.strip().split(".", 1)[1])
                in_loop = True
                continue
            if in_loop and cols and line.startswith(("ATOM", "HETATM")):
                f = line.split()
                if f[cols.index("label_atom_id")] != "CA":
                    continue
                seq.append(THREE_TO_ONE.get(f[cols.index("label_comp_id")], "X"))
                plddt.append(float(f[cols.index("B_iso_or_equiv")]))
            elif in_loop and cols and seq and not line.startswith(("ATOM", "HETATM")):
                break
    return "".join(seq), np.asarray(plddt, dtype=float)


def fetch_one(species: str, acc: str, session: requests.Session) -> str:
    """Download one model if absent. Returns a status string; never raises."""
    out = cif_path(species, acc)
    if out.exists() and out.stat().st_size > 0:
        return "cached"
    for attempt in range(4):
        try:
            r = session.get(API.format(acc=acc), timeout=30)
            if r.status_code == 404:
                return "no_model"
            r.raise_for_status()
            entries = [e for e in r.json() if e.get("uniprotAccession") == acc]
            if not entries:
                return "no_model"
            v = int(entries[0].get("latestVersion", AF_VERSION))
            if v != AF_VERSION:
                return f"version_{v}"
            c = session.get(FILE.format(acc=acc, v=AF_VERSION), timeout=60)
            c.raise_for_status()
            # An HTTP 200 is not evidence of data: assert it is a cif with atoms in it.
            if b"_atom_site." not in c.content:
                return "bad_payload"
            out.write_bytes(c.content)
            return "downloaded"
        except requests.RequestException:
            time.sleep(2 ** attempt)
    return "fetch_failed"


def run(species: str, workers: int) -> pd.DataFrame:
    prot = P.load(species)
    accs = prot["uniprot_ac"].tolist()
    seqs = dict(zip(prot["uniprot_ac"], prot["sequence"]))
    (SOURCE_DIR / species).mkdir(parents=True, exist_ok=True)

    todo = [a for a in accs if not cif_path(species, a).exists()]
    say(f"  {len(accs):,} proteins, {len(accs) - len(todo):,} cached, {len(todo):,} to fetch")
    fetch_status: dict[str, str] = {}
    if todo:
        with requests.Session() as s, ThreadPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(fetch_one, species, a, s): a for a in todo}
            for i, fut in enumerate(as_completed(futs), 1):
                fetch_status[futs[fut]] = fut.result()
                if i % 250 == 0 or i == len(todo):
                    say(f"    fetched {i:,}/{len(todo):,}")
        counts = pd.Series(fetch_status).value_counts()
        say("    fetch outcome: " + ", ".join(f"{k} {v:,}" for k, v in counts.items()))

    rows = []
    for acc in accs:
        row = {"uniprot_ac": acc, "seq_length": len(seqs[acc])}
        try:
            seq, pl, source = model_ca(species, acc)
        except FileNotFoundError:
            row["af_status"] = fetch_status.get(acc, "no_model")
            rows.append(row)
            continue
        expected = seqs[acc] if source == "alphafold_db_v6" else seqs[acc].replace("U", "C")
        row.update(model_source=source, af_length=len(seq),
                   af_version=AF_VERSION if source == "alphafold_db_v6" else None)
        if seq != expected:
            row["af_status"] = "sequence_mismatch"
        else:
            row.update(af_status="model",
                       af_plddt=round(float(np.nanmean(pl)), 2),
                       af_frac_ordered=round(float((pl[~np.isnan(pl)] >= PLDDT_ORDERED).mean()), 4),
                       af_frac_disordered=round(float(
                           (pl[~np.isnan(pl)] < PLDDT_DISORDERED).mean()), 4))
        rows.append(row)

    cols = ["uniprot_ac", "af_status", "model_source", "af_version", "af_length", "seq_length",
            "af_plddt", "af_frac_ordered", "af_frac_disordered"]
    df = M.reindex(pd.DataFrame(rows)[cols], species)
    df["af_version"] = df["af_version"].astype("Int64")
    df["af_length"] = df["af_length"].astype("Int64")
    return df


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say("pockets/structures.py -- AlphaFold DB v6 models, validated against the proteome")
    say(f"  in : data/processed/proteomes/proteome_<species>.tsv, {API.format(acc='<acc>')}")
    say(f"  out: {SOURCE_DIR.relative_to(REPO_ROOT)}/<species>/  +  "
        f"{EVIDENCE_DIR.relative_to(REPO_ROOT)}/alphafold_<species>.tsv")
    say("=" * 92)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    summary = []
    for sp in args.species:
        say(f"\n[{sp}]")
        df = run(sp, args.workers)
        out = EVIDENCE_DIR / f"alphafold_{sp}.tsv"
        df.to_csv(out, sep="\t", index=False)
        n = len(df)
        st = df["af_status"].value_counts()
        m = df[df["af_status"] == "model"]
        summary.append({
            "species": sp, "n": n, "model": int(st.get("model", 0)),
            "pct_model": 100 * st.get("model", 0) / n,
            "no_model": int(st.get("no_model", 0)),
            "mismatch": int(st.get("sequence_mismatch", 0)),
            "other": int(n - st.get("model", 0) - st.get("no_model", 0)
                         - st.get("sequence_mismatch", 0)),
            "median_plddt": m["af_plddt"].median(),
            "pct_plddt_ge70": 100 * (m["af_plddt"] >= PLDDT_ORDERED).mean(),
        })
        say(f"  wrote {out.relative_to(REPO_ROOT)}  ({n:,} rows, canonical)")
        for k, v in st.items():
            say(f"    {k:<18} {v:>6,}  {100 * v / n:5.1f}%")
        for k, v in m["model_source"].value_counts().items():
            say(f"      source {k:<22} {v:>6,}")

    say("\n" + "-" * 92)
    say("SUMMARY")
    s = pd.DataFrame(summary)
    say(s.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    if (s["pct_model"] < 90).any():
        sys.exit("FATAL fewer than 90% of a proteome has a usable model -- check the fetch.")


if __name__ == "__main__":
    main()
