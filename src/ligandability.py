"""Shared helpers for the ligandability assessment stage (docs §2, scripts/06*).

Centralises what the 06a–06h scripts all need: organism config, proteome-accession
loading, cross-species ortholog expansion (reusing the 03a orthology table), drug-like
ligand curation (the 04c ignore-list), and AlphaFold model paths. Keying is always by
UniProt accession (project convention).

Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import csv
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]

# organism -> (proteome_id used in filenames, short prefix used in output filenames)
ORGANISMS: dict[str, tuple[str, str]] = {
    "kpneumoniae": ("UP000007841_HS11286", "kp"),
    "ecoli": ("UP000000625_EcoliK12", "ec"),
}

# Ortholog tiers in the 03a long table that are *bacterial* (evidence transferable as a
# proxy for the anchor protein). "human" is kept separate: it informs whether a family is
# druggable at all, but human activity is NOT selective for an antibacterial programme.
BACTERIAL_TIERS = {"klebsiella", "gram_negative", "bacteria"}
HUMAN_TIER = "human"


# --------------------------------------------------------------------------- paths
def proteome_tsv(organism: str) -> Path:
    pid, _ = ORGANISMS[organism]
    return REPO_ROOT / "data" / "raw" / organism / "proteome" / f"{pid}.tsv"


def proteome_fasta(organism: str) -> Path:
    pid, _ = ORGANISMS[organism]
    return REPO_ROOT / "data" / "raw" / organism / "proteome" / f"{pid}.fasta"


def acc_from_header(header: str) -> str:
    """UniProt accession from a FASTA header like 'tr|A0A0..|NAME' or 'sp|P12345|..' or 'ACC'."""
    h = header.lstrip(">").strip()
    parts = h.split("|")
    if len(parts) >= 2 and parts[0] in ("sp", "tr"):
        return parts[1]
    return h.split()[0]


# DIAMOND binary (osx-64 `gradi-ortho` env; runs directly under Rosetta). Override via env var.
import os as _os  # noqa: E402

DIAMOND_BIN = _os.environ.get(
    "GRADI_DIAMOND_BIN", str(Path.home() / "miniconda3/envs/gradi-ortho/bin/diamond")
)


def run_diamond_blastp(
    query_fasta: Path,
    target_fasta: Path,
    out_tsv: Path,
    threads: int = 8,
    min_id: float = 25.0,
    query_cover: float = 50.0,
    max_target_seqs: int = 50,
    evalue: float = 1e-5,
) -> Path:
    """blastp query vs a target FASTA via DIAMOND (in the gradi-ortho env).

    Output columns (outfmt 6): qseqid sseqid pident qcovhsp scovhsp bitscore. Builds the db in a
    temp dir. Returns out_tsv. Idempotent: skips if out_tsv already exists and is non-empty.
    """
    if out_tsv.exists() and out_tsv.stat().st_size > 0:
        return out_tsv
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        db = Path(td) / "targets"
        subprocess.run(
            [DIAMOND_BIN, "makedb", "--in", str(target_fasta), "-d", str(db), "--quiet"],
            check=True,
        )
        subprocess.run(
            [DIAMOND_BIN, "blastp",
             "-q", str(query_fasta), "-d", str(db), "-o", str(out_tsv),
             "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore",
             "--id", str(min_id), "--query-cover", str(query_cover),
             "--max-target-seqs", str(max_target_seqs), "--evalue", str(evalue),
             "-p", str(threads), "--quiet"],
            check=True,
        )
    return out_tsv


def load_diamond_hits(out_tsv: Path) -> pd.DataFrame:
    cols = ["qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore"]
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return pd.DataFrame(columns=cols)
    return pd.read_csv(out_tsv, sep="\t", names=cols)


def af_cif_path(organism: str, accession: str) -> Path:
    return (
        REPO_ROOT
        / "data"
        / "processed"
        / organism
        / "alphafold"
        / "cif"
        / f"AF-{accession}-F1-model_v6.cif"
    )


def orthologs_long_tsv(organism: str) -> Path:
    _, prefix = ORGANISMS[organism]
    return (
        REPO_ROOT
        / "data"
        / "processed"
        / "other"
        / "orthology"
        / f"{prefix}_orthologs_long.tsv"
    )


def results_dir(organism: str) -> Path:
    d = REPO_ROOT / "output" / "results" / organism
    d.mkdir(parents=True, exist_ok=True)
    return d


def processed_dir(organism: str, *parts: str) -> Path:
    d = REPO_ROOT / "data" / "processed" / organism
    for p in parts:
        d = d / p
    d.mkdir(parents=True, exist_ok=True)
    return d


# --------------------------------------------------------------------------- loaders
def load_accessions(organism: str) -> list[str]:
    """UniProt accessions of the reference proteome, in file order."""
    accs: list[str] = []
    with open(proteome_tsv(organism)) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc = (row.get("Entry") or "").strip()
            if acc:
                accs.append(acc)
    return accs


def load_genes(organism: str) -> dict[str, str]:
    """accession -> first gene symbol (for plot/spot-check labels)."""
    genes: dict[str, str] = {}
    with open(proteome_tsv(organism)) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            acc = (row.get("Entry") or "").strip()
            gn = (row.get("Gene Names") or "").strip()
            if acc:
                genes[acc] = gn.split()[0] if gn else ""
    return genes


def load_orthologs(organism: str) -> pd.DataFrame:
    """The 03a long ortholog table for this anchor.

    Columns: anchor_uniprot, tier, species, target_uniprot, method, pident, coverage,
    evalue, bitscore, is_human_homolog. Adds a `bucket` column in {bacterial, human}.
    """
    df = pd.read_csv(orthologs_long_tsv(organism), sep="\t", low_memory=False)
    df = df[df["target_uniprot"].notna() & (df["target_uniprot"] != "")].copy()
    df["bucket"] = df["tier"].map(
        lambda t: "human" if t == HUMAN_TIER else ("bacterial" if t in BACTERIAL_TIERS else "other")
    )
    return df


def ortholog_map(organism: str) -> dict[str, pd.DataFrame]:
    """anchor_uniprot -> its ortholog rows (DataFrame), for fast per-protein lookup."""
    df = load_orthologs(organism)
    return {acc: sub for acc, sub in df.groupby("anchor_uniprot")}


# --------------------------------------------------------------------------- ligand curation
# Two-tier ligand curation. There are two questions we answer separately:
#   (a) is this het-group a real bound *molecule* at all (not solvent/ion/buffer)?  -> `is_ligand`
#   (b) is that molecule a plausible *drug-like* ligand, i.e. NOT a ubiquitous cofactor or
#       nucleotide that binds thousands of proteins (ATP, NAD, FAD, heme, Fe-S clusters, …)?
#       -> `is_druglike_ligand` (the stricter tier, used as the headline signal everywhere).
#
# LIGAND_IGNORE = solvent / ions / buffer / cryo / crystallisation additives (the broad floor;
# mirrors scripts/04c_pdb_coverage.py:LIGAND_IGNORE).
LIGAND_IGNORE: set[str] = {
    "HOH", "DOD",
    # monatomic ions
    "NA", "K", "LI", "CL", "BR", "IOD", "F", "MG", "CA", "ZN", "MN", "FE", "FE2",
    "CU", "CU1", "NI", "CD", "CO", "HG", "BA", "SR", "CS", "RB", "AL", "GA",
    # buffer / cryo / crystallisation additives
    "SO4", "PO4", "PI", "NO3", "NH4", "ACT", "FMT", "GOL", "EDO", "PEG", "PGE",
    "PG4", "1PE", "P6G", "MPD", "DMS", "TRS", "EPE", "MES", "BME", "IMD", "BO3",
    "CAC", "ACY", "FLC", "CIT", "TLA", "MLI", "OXL", "SCN", "AZI", "PER", "DTT",
    "BTB", "MRD", "POL", "2PE", "12P", "15P",
}

# COFACTOR_IGNORE = real organic molecules that are nonetheless NOT useful drug-likeness evidence:
# promiscuous cofactors, nucleotides/analogues, flavins, CoA, SAM/SAH, hemes/porphyrins, PLP,
# thiamine, cobalamin, biotin, folates/biopterin, molybdopterin, Fe-S / metal clusters,
# glutathione. These bind a large fraction of the proteome and would inflate "ligandable" counts
# (the ATP/ADP/NAD/SF4 problem). Excluded from the drug-like tier, but still counted as broad
# ligands (see `is_ligand`).
COFACTOR_IGNORE: set[str] = {
    # adenine nucleotides & non-hydrolysable analogues
    "AMP", "ADP", "ATP", "ANP", "ACP", "AGS", "APC", "ADX", "A2P",
    # guanine nucleotides & analogues
    "GMP", "GDP", "GTP", "GNP", "GSP", "GCP", "5GP",
    # cytidine / uridine / thymidine nucleotides
    "CMP", "CDP", "CTP", "C5P", "UMP", "UDP", "UTP", "U5P", "TMP", "TYD", "TTP",
    # dinucleotide redox cofactors (NAD(P)/H)
    "NAD", "NAI", "NAP", "NDP", "NAH", "NAJ", "NBD",
    # flavins
    "FAD", "FMN", "FDA", "RBF",
    # coenzyme A & acyl-CoAs
    "COA", "ACO", "MLC", "SCA", "COZ",
    # S-adenosyl methionine / homocysteine
    "SAM", "SAH", "MTA",
    # hemes / porphyrins
    "HEM", "HEA", "HEB", "HEC", "HEV", "DHE", "HDD", "HNI", "HAS", "HCO", "HEO", "COH",
    # pyridoxal-5'-phosphate (vitamin B6)
    "PLP", "PMP",
    # thiamine diphosphate (vitamin B1)
    "TPP", "TDP",
    # cobalamin (vitamin B12)
    "B12", "COB", "CNC",
    # biotin
    "BTN", "BTI",
    # folate / tetrahydrobiopterin
    "FOL", "THF", "THG", "H4B",
    # molybdopterin
    "MGD", "2MD", "MTE",
    # iron-sulfur & other metal clusters
    "SF4", "FES", "F3S", "CFM", "CLF", "ICS", "NFU", "HC0", "HC1", "OEX", "CUA",
    # glutathione (reduced / oxidised)
    "GSH", "GDS",
}

# NONDRUG_IGNORE = other common het-groups that are real molecules but not drug-like evidence:
# free amino acids, simple sugars, polyamines, unknown/placeholder codes, and oxyanion / metal-
# fluoride transition-state mimics. Like COFACTOR_IGNORE, excluded from the drug-like tier only.
NONDRUG_IGNORE: set[str] = {
    # free amino acids (the 20 standard + ornithine / 2-aminobutyrate)
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "ORN", "ABU",
    # simple sugars / amino sugars (substrates & cryoprotectants)
    "GLC", "BGC", "GAL", "MAN", "BMA", "FUC", "FRU", "RIB", "XYP", "ARA",
    "SUC", "MAL", "NAG", "NDG", "NGA", "SIA",
    # polyamines (crystallisation additives)
    "SPD", "SPM", "PUT", "CAD",
    # unknown / placeholder atoms & ligands
    "UNX", "UNL", "UNK", "UNG",
    # oxyanions / metal-fluoride transition-state mimics
    "VO4", "WO4", "MOO", "BEF", "ALF", "AF3",
    # trivial small molecules / additives
    "URE", "ACE", "EOH", "MOH", "IPA",
    # fatty acids
    "PLM", "MYR", "OLA", "STE", "LAU", "DAO", "HXA", "DKA", "PAM", "ACD",
    # phospholipids / glycerolipids (membrane structures)
    "3PE", "PEE", "PEF", "PEH", "PEK", "PEV", "PEU", "PE9", "PGV", "PGW", "PGT",
    "PG0", "LHG", "PTY", "PC1", "PCW", "PSC", "CDL", "DPG", "PX4", "P6L", "17F",
    "G4P", "OLC", "MC3",
    # detergents (membrane-protein crystallisation)
    "LMT", "LDA", "BOG", "BNG", "HTG", "OGA", "C8E", "DDQ", "JEG", "F09", "UMQ",
    "D10", "UND", "9OD",
    # blocked / modified terminal residues seen as het-groups
    "FME", "MSE", "FOR", "NH2", "MLY",
    # lone atoms / diatomics / hydroxide
    "O", "O2", "OXY", "OH", "NO", "CMO", "CO2", "CO3",
    # free nucleobases / nucleosides
    "GUN", "ADE", "CYT", "URA", "THY", "XAN", "HYP", "ADN", "GAO",
    # pyrophosphate / condensed phosphates / oxyanions
    "POP", "PPV", "DPO", "2HP", "PPF", "SO3", "SO2", "SUL", "PO3",
    # organic solvents / cryoprotectants / polyols
    "DIO", "TBU", "BU1", "BU3", "PDO", "PGO", "BUD", "ETN", "IPH",
    # PEG / glycol-ether fragments
    "AE3", "AE4", "1PG", "PG5", "PG6", "XPE", "7PE", "211", "33O", "PE3", "PE4", "M2M",
    # additional buffers (Good's buffers, bis-tris, taurine, …)
    "B3P", "NHE", "TAU", "BES", "CXS", "T3A", "2NV", "MD1",
    # TCA-cycle / small carboxylic-acid metabolites & fragments
    "SIN", "MLA", "MLT", "AKG", "PYR", "BEZ", "FUM", "LAC", "ICT", "MAE", "OAA", "2OG",
    "TAR", "TLA",
    # more cryoprotectants / detergents / buffers surfaced by AlphaFill donor diversity
    "PEO", "HEZ", "LMR", "DXC", "TAM", "OCT", "PGR", "MPO", "BCN", "HED", "DTU",
    "DTV", "MXE", "EGL", "P4G", "PE8", "1BO", "DOX", "OES", "PG6", "EGC",
}

# Everything that disqualifies a real bound molecule from the strict drug-like tier.
_NOT_DRUGLIKE = LIGAND_IGNORE | COFACTOR_IGNORE | NONDRUG_IGNORE


def is_ligand(chem_id: str) -> bool:
    """Broad tier: True if a chem-comp id is a real bound molecule (not solvent/ion/buffer).

    Includes promiscuous cofactors/nucleotides — use `is_druglike_ligand` to exclude those."""
    if not chem_id:
        return False
    return chem_id.strip().upper() not in LIGAND_IGNORE


def is_druglike_ligand(chem_id: str) -> bool:
    """Strict tier: True if a chem-comp id is a plausible drug-like ligand — a real bound
    molecule that is NOT a ubiquitous cofactor/nucleotide (ATP, NAD, FAD, heme, Fe-S, …) nor a
    free amino acid / simple sugar / polyamine / unknown / oxyanion mimic."""
    if not chem_id:
        return False
    return chem_id.strip().upper() not in _NOT_DRUGLIKE


def _dedup_filter(chem_ids, keep) -> list[str]:
    seen, out = set(), []
    for c in chem_ids or []:
        c = (c or "").strip().upper()
        if c and c not in seen and keep(c):
            seen.add(c)
            out.append(c)
    return out


def ligands(chem_ids) -> list[str]:
    """Broad tier: de-duplicated real bound molecules (cofactors/nucleotides included)."""
    return _dedup_filter(chem_ids, is_ligand)


def druglike_ligands(chem_ids) -> list[str]:
    """Strict tier: de-duplicated drug-like ligands (cofactors/nucleotides excluded)."""
    return _dedup_filter(chem_ids, is_druglike_ligand)
