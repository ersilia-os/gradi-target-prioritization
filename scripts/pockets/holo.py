"""What ligands have been seen bound to this protein? -- the measured half of the axis.

    data/processed/pockets/evidence/ligand_classes.tsv          one row per ligand code
    data/processed/pockets/evidence/ligands_pdb_<sp>.tsv        LONG: ligand x this protein's PDB chains
    data/processed/pockets/evidence/ligands_alphafill_<sp>.tsv  LONG: ligand x AlphaFill transplant
    data/processed/pockets/evidence/ligand_counts_<sp>.tsv      per protein, COMPLETE + canonical

Two counts, from two sources that are **never summed**:

* **`n_ligands_pdb`** -- drug-like ligands in **this protein's OWN** experimental structures.
  `pdb_coverage.py` already matched our proteins to PDB *chains* at >= 95% identity; those chains
  ARE this protein, so BioLiP's ligand rows attach by `(pdb, chain)` and **there is no alignment,
  no transfer and no binding-site test to make**.
* **`n_ligands_alphafill`** -- drug-like ligands AlphaFill transplanted onto this protein's
  AlphaFold model (`alphafill.py`), from donors at ~30% median identity. A model, not a
  measurement.

They differ by two orders of magnitude in reach (tens of proteins vs ~1,500), and a co-crystal of
this protein is not the same evidence as a remote structural transplant. Keeping them apart is the
project's standing rule; `src/pockets.py` never adds them.

**An earlier `holo_identity` column was removed on 2026-10-03.** It reported the % identity of the
closest bacterial holo chain found by our own DIAMOND search. It was bimodal rather than
continuous, identity is already known to be a weak axis here
(`ligands/transfer_calibration.py` found the bands uncalibratable), and above all it was a
hand-rolled AlphaFill -- see `alphafill.py` for the measurement that settled it.

NON-REDUNDANT = DISTINCT BEMIS-MURCKO GENERIC SCAFFOLDS
---------------------------------------------------------
Counting chemical-component codes overcounts: a congeneric series is one piece of evidence, not
six. The house definition is the ChEMBL axis's, and **`_scaffold_chunk` is imported from
`scripts/ligands/chembl.py` rather than restated**, so the two axes cannot drift. Acyclic ligands
share the single `""` bucket -- never dropped, never zero. Measured collapse against raw codes:
**25-27%**. Both numbers ship (`n_ligands_*` and `n_codes_*`), so the collapse stays a measurement.

WHAT COUNTS AS "DRUG-LIKE" -- established sources, one question each
---------------------------------------------------------------------
v1 used a hand-built denylist of ~300 codes that grew every time an artefact surfaced
(`legacy/src/ligandability.py:197-297`). It never converges. This replaces it with:

1. **Not an oligomer** -- BioLiP's own `peptide` / `rna` / `dna` classes.
2. **Parseable, and organic** -- a CCD SMILES RDKit reads, containing carbon (drops metal ions,
   Fe-S clusters, phosphate, sulfate).
3. **Not a polymer building block** -- the CCD's own `_chem_comp.type` says `*LINKING*` (amino
   acids, nucleotide units) or `*SACCHARIDE*`.
4. **Not a cofactor** -- the **PDBe cofactor classification** (28 classes, 364 codes; the
   successor to the retired EBI CoFactor database): heme, NAD(P), FAD/FMN, CoA, SAM, PLP, TPP,
   biotin, folate, molybdopterin...
5. **Not a nucleotide** -- a nucleobase N-glycosidically linked to a pentose carrying a
   5'-phosphate. sc-PDB keeps nucleotides as their own category; PDBbind excludes "simple
   cofactors such as ATP". Without this every ATPase looks liganded.
6. **Not a crystallisation artefact** -- **PLINDER's artefact badlist** (265 codes, curated April
   2024 against BioLiP/AF2/RFAA). BioLiP's own curation leaks: C8E, LDA, GOL, PEG and TRS survive
   it and were the best ligand of dozens of proteins before this was added.
7. **Not an endogenous metabolite** -- the **E. coli Metabolome Database** (3,760 metabolites),
   matched by InChIKey skeleton or ECMDB's own het code. A pocket holding pyruvate, glycerol or
   2-oxoglutarate is a substrate site, not a drug precedent.

**QED IS RECORDED AND NOT USED.** PLINDER's default is QED >= 0.2; on the antibiotics that define
bacterial ligandability it is wrong -- clorobiocin 0.086, novobiocin 0.184, rifampicin 0.109,
paromomycin 0.114, kanamycin 0.167 all fail it, because natural products are large and polar.
`druglike_qed` ships beside the headline so the cost is a number, not an argument.

**`in_biolip` is a flag, not a requirement.** AlphaFill transplants 119 codes BioLiP never lists;
for those, criterion 1 and BioLiP's curation cannot apply and only the rest do. The flag makes
that visible in the data rather than hidden in a filter.

Run with the `gradi` env, after `pdb_coverage.py` and `alphafill.py`. No network, ~2 min.
  python scripts/pockets/holo.py
  python scripts/pockets/holo.py --species saureus -q
"""

from __future__ import annotations

import argparse
import gzip
import importlib.util
import io
import json
import sys
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
SRC = REPO_ROOT / "data" / "source"
BIOLIP = SRC / "biolip" / "BioLiP.txt.gz"
CCD = SRC / "wwpdb" / "components.cif.gz"
COFACTORS = SRC / "pdbe" / "cofactors.json"
COFACTORS_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors"
PLINDER_ARTIFACTS = SRC / "plinder" / "artifacts_badlist.csv"
ECMDB = SRC / "ecmdb" / "ecmdb.json.zip"

EVIDENCE_DIR = REPO_ROOT / "data" / "processed" / "pockets" / "evidence"

QED_PLINDER = 0.2
ACYCLIC = "(acyclic)"   # chembl.py's shared bucket for ring-free compounds, named so it survives TSV
OLIGO = {"peptide", "rna", "dna"}
# nucleobase-N -- C1' of a pentose ring -- C4'-CH2-O-P (a 5'-phosphorylated nucleoside).
# Ring atoms are [R], not [R1]: in cyclic-di-GMP (C2E) the ribose atoms are also in the
# macrocycle, and an [R1] pattern missed it -- it was the most frequent "drug-like" code on Kp.
NUCLEOTIDE_SMARTS = "[n,N;R]-[C;R]1-[O;R]-[C;R](-[CH2]-O-P)-[C;R]-[C;R]-1"

BIOLIP_COLS = ["pdb", "chain", "resolution", "site", "ligand", "ligand_chain", "ligand_serial",
               "site_pdb", "site_renum", "cat_pdb", "cat_renum", "ec", "go", "aff_manual",
               "aff_moad", "aff_pdbbind", "aff_bindingdb", "uniprot", "pubmed", "ligand_resnum",
               "sequence"]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def read_tsv(path: Path, **kw) -> pd.DataFrame:
    """**`NA` is sodium's PDB chemical-component code**, and pandas reads it as a missing value by
    default -- it silently nulled 12% of the transplant rows the first time this ran. Empty fields
    are still NaN; the literal string "NA" survives."""
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""], **kw)


# --------------------------------------------------------------------------- sources
def load_biolip() -> pd.DataFrame:
    d = pd.read_csv(BIOLIP, sep="\t", header=None, names=BIOLIP_COLS, dtype=str,
                    keep_default_na=False, quoting=3)
    if len(d.columns) != 21 or d["sequence"].eq("").mean() > 0.01:
        sys.exit("FATAL BioLiP.txt does not have the documented 21 columns -- format changed?")
    return d


def load_cofactors() -> tuple[set[str], dict[str, str]]:
    if not COFACTORS.exists():
        COFACTORS.parent.mkdir(parents=True, exist_ok=True)
        COFACTORS.write_bytes(urllib.request.urlopen(COFACTORS_URL, timeout=60).read())
    data = json.loads(COFACTORS.read_text())
    code2class = {c: cls for cls, entries in data.items() for e in entries
                  for c in e.get("cofactors", [])}
    if len(data) < 20 or len(code2class) < 300:
        sys.exit(f"FATAL PDBe cofactor list looks truncated ({len(data)} classes, "
                 f"{len(code2class)} codes) -- delete {COFACTORS} and re-fetch.")
    return set(code2class), code2class


def load_artifacts() -> set[str]:
    codes = {ln.strip() for ln in PLINDER_ARTIFACTS.read_text().splitlines()
             if ln.strip() and not ln.startswith("#")}
    if len(codes) < 200:
        sys.exit(f"FATAL PLINDER artefact list has {len(codes)} codes -- expected 265")
    return codes


def load_ecmdb() -> tuple[set[str], set[str]]:
    """(InChIKey first blocks, PDB het codes) of every ECMDB metabolite."""
    import zipfile

    with zipfile.ZipFile(ECMDB) as z:
        recs = json.loads(z.read(z.namelist()[0]))
    if len(recs) < 3000:
        sys.exit(f"FATAL ECMDB has {len(recs)} records -- expected 3,760")
    keys = {r["moldb_inchikey"].replace("InChIKey=", "").split("-")[0]
            for r in recs if r.get("moldb_inchikey")}
    het = {r["het_id"].strip().upper() for r in recs if r.get("het_id")}
    return keys, het


def ccd_records(codes: set[str]) -> dict[str, dict]:
    """code -> {ccd_type, name, smiles, inchikey} from the Chemical Component Dictionary."""
    import biotite.structure.io.pdbx as pdbx

    f = pdbx.CIFFile.read(io.StringIO(gzip.open(CCD, "rt").read()))
    out = {}
    for code in codes:
        try:
            b = f[code]
        except KeyError:
            continue
        cc = b["chem_comp"]
        smiles, inchikey = None, None
        try:
            d = b["pdbx_chem_comp_descriptor"]
            rows = list(zip(d["type"].as_array(), d["program"].as_array(),
                            d["descriptor"].as_array()))
            ik = [str(s) for t, _, s in rows if t == "InChIKey"]
            inchikey = ik[0] if ik else None
            for want_t, want_p in (("SMILES_CANONICAL", "CACTVS"),
                                   ("SMILES_CANONICAL", "OpenEye"), ("SMILES", "")):
                hit = [s for t, p, s in rows if t == want_t and p.startswith(want_p)]
                if hit:
                    smiles = str(hit[0])
                    break
        except KeyError:
            pass
        out[code] = {"ccd_type": str(cc["type"].as_item()), "name": str(cc["name"].as_item()),
                     "smiles": smiles, "inchikey": inchikey}
    return out


def scaffold_fn():
    """`_scaffold_chunk` from the ChEMBL axis -- imported, never restated, so the two axes'
    definition of "non-redundant" cannot drift. `scripts/` is not a package, hence the spec load."""
    spec = importlib.util.spec_from_file_location(
        "ligands_chembl", REPO_ROOT / "scripts" / "ligands" / "chembl.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "_scaffold_chunk"):
        sys.exit("FATAL scripts/ligands/chembl.py has no _scaffold_chunk -- it was refactored; "
                 "this axis's non-redundancy must be re-aligned rather than guessed at")
    return mod._scaffold_chunk


# --------------------------------------------------------------------------- 1. vocabulary
def classify(codes: set[str], biolip_codes: set[str], counts: pd.Series) -> pd.DataFrame:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED

    RDLogger.DisableLog("rdApp.*")
    nuc = Chem.MolFromSmarts(NUCLEOTIDE_SMARTS)
    cof, cof_class = load_cofactors()
    artifacts = load_artifacts()
    ecmdb_keys, ecmdb_het = load_ecmdb()
    ccd = ccd_records(codes - OLIGO)

    rows = []
    for code in sorted(codes):
        r = {"ligand": code, "in_biolip": code in biolip_codes,
             "n_biolip_sites": int(counts.get(code, 0))}
        if code in OLIGO:
            r.update(is_oligo=True)
            rows.append(r)
            continue
        rec = ccd.get(code, {})
        mol = Chem.MolFromSmiles(rec["smiles"]) if rec.get("smiles") else None
        t = rec.get("ccd_type", "").upper()
        ik = rec.get("inchikey") or (Chem.MolToInchiKey(mol) if mol is not None else None)
        block = ik.split("-")[0] if ik else None
        route = ("inchikey" if block in ecmdb_keys else "het_id" if code in ecmdb_het else None)
        r.update(name=rec.get("name"), ccd_type=rec.get("ccd_type"), smiles=rec.get("smiles"),
                 inchikey=ik, is_oligo=False, parsed=mol is not None,
                 is_inorganic=(not any(a.GetSymbol() == "C" for a in mol.GetAtoms()))
                 if mol is not None else None,
                 is_polymer_unit=("LINKING" in t) or ("SACCHARIDE" in t),
                 is_cofactor=code in cof, cofactor_class=cof_class.get(code),
                 is_nucleotide=mol.HasSubstructMatch(nuc) if mol is not None else None,
                 is_artifact=code in artifacts,
                 is_metabolite=route is not None, metabolite_route=route)
        try:
            r["qed"] = round(QED.qed(mol), 4) if mol is not None else None
        except Exception:  # noqa: BLE001 -- some organometallics have no QED
            r["qed"] = None
        rows.append(r)

    d = pd.DataFrame(rows)
    for c in ("is_oligo", "in_biolip", "is_inorganic", "is_polymer_unit", "is_cofactor",
              "is_nucleotide", "parsed", "is_artifact", "is_metabolite"):
        d[c] = d[c].astype("boolean")
    # A ligand whose structure cannot be parsed cannot be shown to be organic -> not drug-like.
    d["druglike"] = (~d["is_oligo"].fillna(False) & d["parsed"].fillna(False)
                     & ~d["is_inorganic"].fillna(True) & ~d["is_polymer_unit"].fillna(False)
                     & ~d["is_cofactor"].fillna(False) & ~d["is_nucleotide"].fillna(False)
                     & ~d["is_artifact"].fillna(False)
                     & ~d["is_metabolite"].fillna(False)).astype(bool)
    d["druglike_qed"] = d["druglike"] & (pd.to_numeric(d["qed"]) >= QED_PLINDER).fillna(False)

    dl = d[d["druglike"]]
    scaf = scaffold_fn()(dl["smiles"].fillna("").tolist())
    # chembl.py returns "" for an acyclic compound -- one shared bucket, never dropped. An empty
    # string round-trips through TSV as NaN, and `nunique()` drops NaN, so every acyclic ligand
    # would silently vanish from the count. Name the bucket instead.
    d["scaffold"] = pd.Series([x if x else ACYCLIC for x in scaf],
                              index=dl.index).reindex(d.index)
    return d


def report(d: pd.DataFrame) -> None:
    say("  what each criterion removes (codes / BioLiP binding sites):")
    steps = [("all codes (BioLiP + AlphaFill)", pd.Series(True, index=d.index)),
             ("- oligo (peptide/rna/dna)", ~d["is_oligo"].fillna(False)),
             ("- unparseable SMILES", d["parsed"].fillna(False)),
             ("- inorganic (no carbon)", ~d["is_inorganic"].fillna(True)),
             ("- polymer unit / saccharide", ~d["is_polymer_unit"].fillna(False)),
             ("- PDBe cofactor", ~d["is_cofactor"].fillna(False)),
             ("- nucleotide", ~d["is_nucleotide"].fillna(False)),
             ("- PLINDER artefact list", ~d["is_artifact"].fillna(False)),
             ("- ECMDB metabolite", ~d["is_metabolite"].fillna(False))]
    keep = pd.Series(True, index=d.index)
    for label, m in steps:
        keep &= m.astype(bool)
        say(f"    {label:<32} {int(keep.sum()):>7,} codes  "
            f"{int(d.loc[keep, 'n_biolip_sites'].sum()):>9,} sites")
    q = d["druglike_qed"]
    say(f"    {'(variant) + QED >= 0.2':<32} {int(q.sum()):>7,} codes  "
        f"{int(d.loc[q, 'n_biolip_sites'].sum()):>9,} sites")
    say(f"  codes AlphaFill uses that BioLiP never lists: {int((~d['in_biolip']).sum()):,} "
        f"({int((~d['in_biolip'] & d['druglike']).sum()):,} of them drug-like)")


# --------------------------------------------------------------------------- 2. the two counts
def count(long: pd.DataFrame, species: str, name: str) -> pd.DataFrame:
    """Per protein: distinct scaffolds (the headline) and distinct codes (the comparison)."""
    accs = P.load(species)[["uniprot_ac"]]
    if long.empty:
        out = accs.assign(**{f"n_ligands_{name}": 0, f"n_codes_{name}": 0})
    else:
        g = long.groupby("uniprot_ac")
        out = accs.merge(pd.DataFrame({
            f"n_ligands_{name}": g["scaffold"].nunique(),
            f"n_codes_{name}": g["ligand"].nunique()}).reset_index(), on="uniprot_ac", how="left")
    for c in (f"n_ligands_{name}", f"n_codes_{name}"):
        out[c] = out[c].fillna(0).astype(int)
    return out


def pdb_ligands(species: str, bl: pd.DataFrame, classes: pd.DataFrame) -> pd.DataFrame:
    """Drug-like ligands on the PDB chains that ARE this protein. No alignment: see the header."""
    path = EVIDENCE_DIR / f"pdb_chains_{species}.tsv"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- run pockets/pdb_coverage.py")
    chains = read_tsv(path)
    dl = classes.loc[classes["druglike"], ["ligand", "scaffold"]]
    lig = bl[["pdb", "chain", "ligand"]].drop_duplicates().merge(dl, on="ligand")
    out = chains.merge(lig, on=["pdb", "chain"])
    return out[["uniprot_ac", "ligand", "scaffold", "pdb", "chain", "identity"]].drop_duplicates()


def alphafill_ligands(species: str, classes: pd.DataFrame) -> pd.DataFrame:
    path = EVIDENCE_DIR / f"transplants_{species}.tsv"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- run pockets/alphafill.py")
    t = read_tsv(path)
    dl = classes.loc[classes["druglike"], ["ligand", "scaffold"]]
    return t.merge(dl, on="ligand")


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say("pockets/holo.py -- ligands bound to this protein (own PDB chains) and transplanted onto "
        "it (AlphaFill)")
    say("  in : BioLiP, wwPDB CCD, PDBe cofactors, PLINDER artefacts, ECMDB, "
        "evidence/{pdb_chains,transplants}_<sp>.tsv")
    say("  out: evidence/{ligand_classes, ligands_pdb_<sp>, ligands_alphafill_<sp>, "
        "ligand_counts_<sp>}.tsv")
    say("  non-redundant = distinct Bemis-Murcko generic scaffolds (chembl.py's own function)")
    say("=" * 92)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    say("\n[1/3] the ligand vocabulary")
    bl = load_biolip()
    say(f"  BioLiP: {len(bl):,} ligand-site rows over {bl['pdb'].nunique():,} PDB entries")
    biolip_codes = set(bl["ligand"])
    af_codes: set[str] = set()
    for sp in SPECIES:
        p = EVIDENCE_DIR / f"transplants_{sp}.tsv"
        if p.exists():
            af_codes |= set(read_tsv(p, usecols=["ligand"])["ligand"].dropna())
    say(f"  AlphaFill transplant codes: {len(af_codes):,}")
    if not af_codes:
        sys.exit("FATAL no AlphaFill codes -- run pockets/alphafill.py first")
    classes = classify(biolip_codes | af_codes, biolip_codes, bl.groupby("ligand").size())
    classes.to_csv(EVIDENCE_DIR / "ligand_classes.tsv", sep="\t", index=False)
    report(classes)

    say("\n[2/3] ligands per protein")
    summary = []
    for sp in args.species:
        pdb_long = pdb_ligands(sp, bl, classes)
        af_long = alphafill_ligands(sp, classes)
        pdb_long.to_csv(EVIDENCE_DIR / f"ligands_pdb_{sp}.tsv", sep="\t", index=False)
        af_long.to_csv(EVIDENCE_DIR / f"ligands_alphafill_{sp}.tsv", sep="\t", index=False)
        counts = count(pdb_long, sp, "pdb").merge(count(af_long, sp, "alphafill"), on="uniprot_ac")
        counts = M.reindex(counts, sp)
        counts.to_csv(EVIDENCE_DIR / f"ligand_counts_{sp}.tsv", sep="\t", index=False)
        summary.append({
            "species": sp, "n": len(counts),
            "pdb>0": int((counts["n_ligands_pdb"] > 0).sum()),
            "pdb_codes": int(counts["n_codes_pdb"].sum()),
            "pdb_scaffolds": int(counts["n_ligands_pdb"].sum()),
            "af>0": int((counts["n_ligands_alphafill"] > 0).sum()),
            "af_codes": int(counts["n_codes_alphafill"].sum()),
            "af_scaffolds": int(counts["n_ligands_alphafill"].sum()),
        })
        say(f"  [{sp}] own-PDB {len(pdb_long):,} ligand-chain rows, "
            f"AlphaFill {len(af_long):,} drug-like transplants")

    say("\n[3/3] SUMMARY -- scaffolds are the shipped count; codes show what redundancy costs")
    s = pd.DataFrame(summary)
    s["pdb_collapse_%"] = (100 * (1 - s["pdb_scaffolds"] / s["pdb_codes"].clip(lower=1))).round(1)
    s["af_collapse_%"] = (100 * (1 - s["af_scaffolds"] / s["af_codes"].clip(lower=1))).round(1)
    say(s.to_string(index=False))
    say("\n  pdb = this protein's OWN structures (>=95% identity chains) -- a measurement.")
    say("  af  = AlphaFill transplants from ~30%-identity donors -- a model. NEVER summed.")


if __name__ == "__main__":
    main()
