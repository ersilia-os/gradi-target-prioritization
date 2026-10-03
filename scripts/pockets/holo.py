"""Has a drug-like molecule been seen BOUND to this protein's family? -- the measured half of the axis.

    data/processed/pockets/evidence/ligand_classes.tsv     one row per BioLiP ligand code
    data/processed/pockets/evidence/holo_<species>.tsv     one row per protein, COMPLETE
    data/processed/pockets/scratch/holo/                   FASTA, DIAMOND db and raw hits

`holo_identity` (0-100) is the % identity to the closest BACTERIAL protein chain observed in the
PDB with a drug-like ligand bound, 100 being this protein's own co-crystal and 0 meaning no chain
clears the floors. It is a measurement of the family, not a prediction.

WHAT COUNTS AS "DRUG-LIKE" -- established sources, one question each
---------------------------------------------------------------------
v1 used a hand-built denylist of ~300 chemical-component codes that grew every time a new artefact
surfaced (`legacy/src/ligandability.py:197-297`). It never converges. This replaces it with:

1. **Biologically relevant, not a crystallisation artefact** -- membership of **BioLiP**
   (Yang, Roy & Zhang, NAR 2013), the semi-manually curated set of ligand-protein contacts with
   buffers, cryoprotectants and other artefacts removed. Every row this script reads is in BioLiP.
   **Plus PLINDER's artefact badlist** (Durairaj et al. 2024; 265 codes, curated April 2024
   "comparing to BIOLIP, AF2 and RFAA artifact list"), vendored at `data/source/plinder/`. BioLiP
   alone measurably leaks: C8E, LDA (detergents), GOL, PEG and TRS survive its curation and were
   the best "drug-like" ligand of dozens of proteins on the first run.
2. **Not an oligomer** -- BioLiP's own `peptide` / `rna` / `dna` classes.
3. **Not inorganic** -- no carbon atom (metal ions, Fe-S clusters, phosphate, sulfate).
4. **Not a polymer building block** -- the PDB Chemical Component Dictionary's own `_chem_comp.type`
   says `*LINKING*` (free amino acids, nucleotide monomers) or `*SACCHARIDE*`.
5. **Not a cofactor** -- the **PDBe cofactor classification** (`/pdbe/api/pdb/compound/cofactors`,
   28 classes, derived from the retired EBI CoFactor database): heme, NAD(P), FAD/FMN, CoA, SAM,
   PLP, TPP, biotin, folate, molybdopterin, ...
6. **Not a nucleotide** -- a nucleobase N-glycosidically linked to a pentose carrying a 5'-phosphate
   (ATP, ADP, GDP, AMP-PNP...). sc-PDB (Desaphy et al., NAR 2015) keeps nucleotides as their own
   ligand category, distinct from organic compounds; PDBbind excludes "simple cofactors such as
   ATP". The PDBe classes do not cover them, and they would make every ATPase look liganded.
7. **Not an endogenous metabolite** -- not in the **E. coli Metabolome Database** (ECMDB; Guo et al.
   NAR 2013, Sajed et al. NAR 2016; 3,760 metabolites), matched by InChIKey skeleton (first block,
   so protonation and stereo variants match) or by ECMDB's own PDB het code. A pocket seen holding
   pyruvate, glycerol, adenine or 2-oxoglutarate is a substrate site, not a drug precedent -- and on
   the first run those were 10-30% of all evidence. ECMDB contains no antibiotics (checked:
   fosfomycin, cycloserine, trimethoprim, rifampicin, novobiocin, ampicillin, tetracycline,
   chloramphenicol, kanamycin, ciprofloxacin, erythromycin, triclosan, methotrexate -- none), so it
   removes metabolites without touching the small xenobiotics that a size rule would.
   Chosen by the project owner on 2026-10-03 over three measured alternatives, all reported every
   run: keep metabolites, drop PLINDER Ro3 "fragments" (also drops fosfomycin and D-cycloserine),
   or QED >= 0.2 (below).

**QED IS RECORDED AND NOT USED -- measured, not assumed.** PLINDER's default filter is QED >= 0.2,
and it was the first choice here. On the antibiotics that define bacterial ligandability it is
wrong: clorobiocin 0.086, novobiocin 0.184 (GyrB), rifampicin 0.109 (RNA polymerase), paromomycin
0.114, kanamycin 0.167 (ribosome) all fail it -- natural products are large and polar, which QED
penalises by design. A QED cut would erase exactly the sites antibacterial drugs bind. The
QED-filtered variant is computed beside the headline (`holo_identity_qed`) so the effect is a
number, not an argument.

MATCHING -- by sequence, and by the BINDING SITE, not the whole chain
----------------------------------------------------------------------
DIAMOND against every BioLiP receptor chain that carries a drug-like ligand. A hit counts when
  * identity >= 40% (`src.ligandability.REMOTE_PIDENT`, the house transfer floor),
  * the PDB chain is >= 50% aligned (`src.ligandability.MIN_SCOV`), and
  * **every binding-site residue of the ligand lies inside the aligned span.**
There is deliberately NO query-coverage floor: co-crystals are very often domain constructs, and
the house `MIN_QCOV = 50` would reject them -- E. coli GyrB's clorobiocin complex (1kzn) is a
186-residue ATPase domain of an 804-residue protein, i.e. 23% query coverage. What has to be shown
is that this protein contains the SITE, and BioLiP gives the site residues, so that is tested
directly. **Bacterial** = every NCBI taxid SIFTS assigns the chain descends from Bacteria (2) --
the house rule for a bacterial bucket (true Bacteria, not "non-human").

Run with the `gradi` env; DIAMOND from `gradi-ortho` (`GRADI_DIAMOND_BIN` overrides).
  python scripts/pockets/holo.py
  python scripts/pockets/holo.py --species saureus -q
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402
from src.precedents import diamond_bin  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
SRC = REPO_ROOT / "data" / "source"
BIOLIP = SRC / "biolip" / "BioLiP.txt.gz"
CCD = SRC / "wwpdb" / "components.cif.gz"
SIFTS_TAX = SRC / "sifts" / "pdb_chain_taxonomy.tsv.gz"
TAXDUMP = SRC / "ncbi" / "taxonomy" / "taxdump.tar.gz"
COFACTORS = SRC / "pdbe" / "cofactors.json"
COFACTORS_URL = "https://www.ebi.ac.uk/pdbe/api/pdb/compound/cofactors"
PLINDER_ARTIFACTS = SRC / "plinder" / "artifacts_badlist.csv"
ECMDB = SRC / "ecmdb" / "ecmdb.json.zip"

TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
WORK = TASK_DIR / "scratch" / "holo"

MIN_PIDENT = L.REMOTE_PIDENT
MIN_SCOV = L.MIN_SCOV
EXACT_PIDENT = L.DIRECT_PIDENT
QED_PLINDER = 0.2
BACTERIA_TAXID = 2
OLIGO = {"peptide", "rna", "dna"}

# nucleobase-N -- C1' of a pentose ring -- C4'-CH2-O-P  (a 5'-phosphorylated nucleoside).
# Ring atoms are [R], not [R1]: in cyclic-di-GMP (C2E) the ribose atoms are also in the macrocycle,
# and an [R1] pattern missed it -- it was the single most frequent "drug-like" ligand on Kp.
NUCLEOTIDE_SMARTS = "[n,N;R]-[C;R]1-[O;R]-[C;R](-[CH2]-O-P)-[C;R]-[C;R]-1"

BIOLIP_COLS = ["pdb", "chain", "resolution", "site", "ligand", "ligand_chain", "ligand_serial",
               "site_pdb", "site_renum", "cat_pdb", "cat_renum", "ec", "go", "aff_manual",
               "aff_moad", "aff_pdbbind", "aff_bindingdb", "uniprot", "pubmed", "ligand_resnum",
               "sequence"]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


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
    lines = PLINDER_ARTIFACTS.read_text().splitlines()
    codes = {ln.strip() for ln in lines if ln.strip() and not ln.startswith("#")}
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
    """code -> {type, name, smiles, inchikey} from the Chemical Component Dictionary."""
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


def bacterial_taxa() -> set[int]:
    """Every NCBI taxid descending from Bacteria (2), from nodes.dmp."""
    with tarfile.open(TAXDUMP) as t:
        nodes = t.extractfile("nodes.dmp").read().decode()
    parent = {}
    for line in nodes.splitlines():
        f = line.split("\t|\t")
        parent[int(f[0])] = int(f[1])
    memo: dict[int, bool] = {1: False, BACTERIA_TAXID: True}

    def is_bact(t: int) -> bool:
        path = []
        while t not in memo:
            path.append(t)
            p = parent.get(t)
            if p is None or p == t:
                memo[t] = False
                break
            t = p
        v = memo[t]
        for x in path:
            memo[x] = v
        return v

    return {t for t in parent if is_bact(t)}


# --------------------------------------------------------------------------- 1. ligand classes
def classify_ligands(bl: pd.DataFrame) -> pd.DataFrame:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import QED

    RDLogger.DisableLog("rdApp.*")
    nuc = Chem.MolFromSmarts(NUCLEOTIDE_SMARTS)
    cof, cof_class = load_cofactors()
    artifacts = load_artifacts()
    ecmdb_keys, ecmdb_het = load_ecmdb()
    counts = bl.groupby("ligand").size()
    codes = set(counts.index) - OLIGO
    ccd = ccd_records(codes)

    rows = []
    for code in sorted(set(counts.index)):
        r = {"ligand": code, "n_biolip_sites": int(counts[code])}
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
        r.update(name=rec.get("name"), ccd_type=rec.get("ccd_type"), inchikey=ik,
                 is_oligo=False, parsed=mol is not None,
                 is_artifact=code in artifacts,
                 is_metabolite=route is not None, metabolite_route=route,
                 is_inorganic=(not any(a.GetSymbol() == "C" for a in mol.GetAtoms()))
                 if mol is not None else None,
                 is_polymer_unit=("LINKING" in t) or ("SACCHARIDE" in t),
                 is_cofactor=code in cof, cofactor_class=cof_class.get(code),
                 is_nucleotide=mol.HasSubstructMatch(nuc) if mol is not None else None)
        try:
            r["qed"] = round(QED.qed(mol), 4) if mol is not None else None
        except Exception:  # noqa: BLE001 -- some organometallics have no QED
            r["qed"] = None
        rows.append(r)

    d = pd.DataFrame(rows)
    for c in ("is_oligo", "is_inorganic", "is_polymer_unit", "is_cofactor", "is_nucleotide",
              "parsed", "is_artifact", "is_metabolite"):
        d[c] = d[c].astype("boolean")
    # A ligand whose structure cannot be parsed cannot be shown to be organic -> not drug-like.
    d["druglike"] = (~d["is_oligo"] & d["parsed"].fillna(False) & ~d["is_inorganic"].fillna(True)
                     & ~d["is_polymer_unit"].fillna(False) & ~d["is_cofactor"].fillna(False)
                     & ~d["is_nucleotide"].fillna(False) & ~d["is_artifact"].fillna(False)
                     & ~d["is_metabolite"].fillna(False)).astype(bool)
    d["druglike_qed"] = d["druglike"] & (pd.to_numeric(d["qed"]) >= QED_PLINDER).fillna(False)
    return d


def report_classes(d: pd.DataFrame) -> None:
    say("  ligand codes in BioLiP, and what each criterion removes (codes / binding sites):")
    steps = [("all BioLiP ligand codes", pd.Series(True, index=d.index)),
             ("- oligo (peptide/rna/dna)", ~d["is_oligo"]),
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
        say(f"    {label:<32} {int(keep.sum()):>7,} codes  {int(d.loc[keep, 'n_biolip_sites'].sum()):>9,} sites")
    q = d["druglike_qed"]
    say(f"    {'(variant) + QED >= 0.2':<32} {int(q.sum()):>7,} codes  {int(d.loc[q, 'n_biolip_sites'].sum()):>9,} sites")
    top = d[d["druglike"]].nlargest(12, "n_biolip_sites")
    say("    most frequent drug-like codes: "
        + ", ".join(f"{r.ligand}({r.n_biolip_sites})" for r in top.itertuples()))


# --------------------------------------------------------------------------- 2. holo chains
def holo_chains(bl: pd.DataFrame, classes: pd.DataFrame) -> pd.DataFrame:
    """One row per (pdb, chain, ligand site) with a drug-like ligand, plus taxonomy."""
    dl = set(classes.loc[classes["druglike"], "ligand"])
    dq = set(classes.loc[classes["druglike_qed"], "ligand"])
    h = bl[bl["ligand"].isin(dl)].copy()
    h["qed_pass"] = h["ligand"].isin(dq)

    tax = pd.read_csv(SIFTS_TAX, sep="\t", comment="#", dtype=str, keep_default_na=False)
    tax.columns = [c.strip() for c in tax.columns]
    tax["TAX_ID"] = pd.to_numeric(tax["TAX_ID"], errors="coerce")
    bact = bacterial_taxa()
    tax["bact"] = tax["TAX_ID"].isin(bact)
    per_chain = tax.groupby(["PDB", "CHAIN"])["bact"].agg(["all", "size"]).reset_index()
    per_chain.columns = ["pdb", "chain", "bacterial", "n_taxa"]
    h = h.merge(per_chain, on=["pdb", "chain"], how="left")
    h["taxonomy_known"] = h["n_taxa"].notna()
    h["bacterial"] = h["bacterial"].astype("boolean").fillna(False).astype(bool)
    return h


# --------------------------------------------------------------------------- 3. DIAMOND
def write_fasta(seqs: dict[str, str], path: Path) -> None:
    with open(path, "w") as fh:
        for k, s in seqs.items():
            fh.write(f">{k}\n{s}\n")


def run_diamond(species: str, subj_faa: Path, threads: int) -> pd.DataFrame:
    dmnd = WORK / "holo_chains.dmnd"
    if not dmnd.exists() or dmnd.stat().st_mtime < subj_faa.stat().st_mtime:
        subprocess.run([diamond_bin(), "makedb", "--in", str(subj_faa), "-d", str(dmnd),
                        "--quiet"], check=True)
    q = WORK / f"query_{species}.faa"
    prot = P.load(species)
    write_fasta(dict(zip(prot["uniprot_ac"], prot["sequence"])), q)
    out = WORK / f"hits_{species}.tsv"
    if not out.exists() or out.stat().st_mtime < dmnd.stat().st_mtime:
        fmt = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
               "evalue", "bitscore", "qlen", "slen"]
        subprocess.run([diamond_bin(), "blastp", "-q", str(q), "-d", str(dmnd), "-o", str(out),
                        "--sensitive", "-k", "2000", "-e", "1e-5", "--threads", str(threads),
                        "--quiet", "--outfmt", "6", *fmt], check=True)
    cols = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart", "send",
            "evalue", "bitscore", "qlen", "slen"]
    return pd.read_csv(out, sep="\t", header=None, names=cols)


def site_positions(s: str) -> list[int]:
    """'N32 E36 V57' -> [32, 36, 57]: BioLiP's renumbered site, 1-based on the receptor sequence."""
    out = []
    for tok in s.split():
        digits = tok[1:]
        if digits.lstrip("-").isdigit():
            out.append(int(digits))
    return out


def score_species(species: str, hits: pd.DataFrame, h: pd.DataFrame,
                  seq_id_of: dict[tuple[str, str], str]) -> pd.DataFrame:
    """Per-protein holo evidence from raw hits + the per-site BioLiP table."""
    hits = hits[hits["pident"] >= MIN_PIDENT].copy()
    hits["scov"] = 100 * (hits["send"] - hits["sstart"] + 1) / hits["slen"]
    hits = hits[hits["scov"] >= MIN_SCOV]

    sites = h.assign(sid=[seq_id_of[(p, c)] for p, c in zip(h["pdb"], h["chain"])])
    sites["pos"] = sites["site_renum"].map(site_positions)
    sites = sites[sites["pos"].map(len) > 0]
    sites["site_min"] = sites["pos"].map(min)
    sites["site_max"] = sites["pos"].map(max)
    keep = ["sid", "pdb", "chain", "ligand", "bacterial", "qed_pass", "site_min", "site_max"]
    # One identical sequence can sit in hundreds of PDB entries with the same ligand in the same
    # place (DHFR, lysozyme). Collapse to distinct sites before the join, or it explodes.
    sites = sites[keep].drop_duplicates(
        ["sid", "ligand", "bacterial", "qed_pass", "site_min", "site_max"])
    m = hits.merge(sites, left_on="sseqid", right_on="sid", how="inner")
    m = m[(m["site_min"] >= m["sstart"]) & (m["site_max"] <= m["send"])]

    rows = []
    prot = P.load(species)[["uniprot_ac"]]
    g = dict(tuple(m.groupby("qseqid")))
    for acc in prot["uniprot_ac"]:
        r = {"uniprot_ac": acc}
        d = g.get(acc)
        for scope, mask in (("bacterial", "bacterial"), ("any", None)):
            dd = d if d is None or mask is None else d[d[mask]]
            if dd is None or dd.empty:
                r[f"holo_identity_{scope}"] = 0.0
                continue
            best = dd.sort_values(["pident", "bitscore"], ascending=False).iloc[0]
            r[f"holo_identity_{scope}"] = float(best["pident"])
            r[f"best_pdb_{scope}"] = f"{best['pdb']}_{best['chain']}"
            r[f"best_ligand_{scope}"] = best["ligand"]
            r[f"n_ligands_{scope}"] = int(dd["ligand"].nunique())
            if scope == "bacterial":
                ex = dd[dd["pident"] >= EXACT_PIDENT]
                r["n_ligands_exact"] = int(ex["ligand"].nunique())
                r["ligands_exact"] = ";".join(sorted(ex["ligand"].unique()))
        dq = None if d is None else d[d["bacterial"] & d["qed_pass"]]
        r["holo_identity_qed"] = 0.0 if dq is None or dq.empty else float(dq["pident"].max())
        rows.append(r)

    out = pd.DataFrame(rows)
    for c in ("n_ligands_bacterial", "n_ligands_any", "n_ligands_exact"):
        out[c] = out[c].fillna(0).astype(int)
    out["ligands_exact"] = out["ligands_exact"].fillna("")
    cols = ["uniprot_ac", "holo_identity_bacterial", "best_pdb_bacterial",
            "best_ligand_bacterial", "n_ligands_bacterial", "n_ligands_exact", "ligands_exact",
            "holo_identity_any", "best_pdb_any", "best_ligand_any", "n_ligands_any",
            "holo_identity_qed"]
    return M.reindex(out[cols], species)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say("pockets/holo.py -- closest bacterial PDB chain with a drug-like ligand in the same site")
    say("  in : BioLiP, wwPDB CCD, PDBe cofactor classes, SIFTS chain taxonomy, NCBI taxdump")
    say("  out: data/processed/pockets/evidence/{ligand_classes,holo_<species>}.tsv")
    say(f"  floors: identity >= {MIN_PIDENT:.0f}%, PDB chain >= {MIN_SCOV:.0f}% aligned, every "
        "site residue inside the alignment; no query-coverage floor (domain constructs)")
    say("=" * 92)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    say("\n[1/3] classify every BioLiP ligand")
    bl = load_biolip()
    say(f"  BioLiP: {len(bl):,} ligand-site rows over {bl['pdb'].nunique():,} PDB entries")
    classes = classify_ligands(bl)
    classes.to_csv(EVIDENCE_DIR / "ligand_classes.tsv", sep="\t", index=False)
    report_classes(classes)

    say("\n[2/3] holo chains: BioLiP sites with a drug-like ligand, with taxonomy")
    h = holo_chains(bl, classes)
    seqs: dict[str, str] = {}
    seq_id_of: dict[tuple[str, str], str] = {}
    by_seq: dict[str, str] = {}
    for (p, c), s in h.groupby(["pdb", "chain"])["sequence"].first().items():
        sid = by_seq.setdefault(s, f"s{len(by_seq)}")
        seqs[sid] = s
        seq_id_of[(p, c)] = sid
    subj = WORK / "holo_chains.faa"
    write_fasta(seqs, subj)
    n_chain = len(seq_id_of)
    say(f"  {len(h):,} drug-like sites on {n_chain:,} chains ({len(seqs):,} distinct sequences)")
    say(f"  taxonomy known for {100 * h['taxonomy_known'].mean():.1f}% of sites; "
        f"bacterial {100 * h['bacterial'].mean():.1f}%")

    say("\n[3/3] DIAMOND each proteome against the holo chains")
    summary = []
    for sp in args.species:
        hits = run_diamond(sp, subj, args.threads)
        out = score_species(sp, hits, h, seq_id_of)
        path = EVIDENCE_DIR / f"holo_{sp}.tsv"
        out.to_csv(path, sep="\t", index=False)
        n = len(out)
        b = out["holo_identity_bacterial"]
        summary.append({
            "species": sp, "n": n,
            "bacterial>0": int((b > 0).sum()), "pct": 100 * (b > 0).mean(),
            "exact(>=95)": int((b >= EXACT_PIDENT).sum()),
            "any_org>0": int((out["holo_identity_any"] > 0).sum()),
            "qed_variant>0": int((out["holo_identity_qed"] > 0).sum()),
            "median_when>0": float(b[b > 0].median()) if (b > 0).any() else np.nan,
        })
        say(f"  [{sp}] {len(hits):,} raw hits -> wrote {path.relative_to(REPO_ROOT)} ({n:,} rows)")

    say("\n" + "-" * 92)
    say("SUMMARY -- proteins with a drug-like holo structure in the family")
    say(pd.DataFrame(summary).to_string(index=False, float_format=lambda x: f"{x:.1f}"))


if __name__ == "__main__":
    main()
