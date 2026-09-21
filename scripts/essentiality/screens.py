"""Published essentiality screens -> one training set per output column.

Nine measured screens plus one conservation label. Each stays SEPARATE: different strains, different
assays, different analyses. Merging them would produce a number with no referent -- the base rate
alone spans 0.060 to 0.194 across this set, and that spread is real information about method and
strain rather than noise to average away.

THREE FILTERS define the set, each on instruction from the project owner:

  both classes   a positives-only gene list is not a training set. Drops Ramage 2017 (KPNIH1, 424
                 genes, no measured negatives) and Paczosa 2020 (310 published hits).
  no conditions  drops Choe's M9 arm (minimal medium -- its extra 131 "essentials" are
                 auxotrophies), Rome 2026 (iron-depleted), Short 2020's serum screens, and
                 Bruchmann's 2hpi/6hpi in-host columns.
  no duplicates  DEG1019 is the Keio collection again under DEG's curation; PEC's copy is kept.

Everything filtered out is still downloaded and parseable and is listed in HELD_BACK with its
reason -- those are decisions, not gaps. `--list` prints them, along with the dead ends in
UNJOINABLE that were measured and refuted and must not be re-derived.

WHY THE IDENTIFIER WORK DOMINATES THIS FILE: v1 held much of this data and joined it by GENE SYMBOL,
losing most of it -- Ramage 212/424, ECL8 954/5,165, KPPR1 32 of 3,791. Its own retrospective calls
that the single most expensive recurring cost in v1. Every join here lands on a stable identifier or
an exact sequence, the rate is MEASURED every run into `evidence/screen_join_audit.tsv`, and a rate
below the screen's floor exits non-zero rather than shipping a quietly decimated training set.

    THE ONE PLACE A SYMBOL JOIN IS SAFE. Gene-name coverage is 100% on E. coli against 63.4% on
    K. pneumoniae, so the ~27% loss v1 documents is a KLEBSIELLA problem. Goodall 2018 ships symbols
    and nothing else; on E. coli that costs little -- which this script proves each run by reporting
    the rate rather than assuming it, and by trying the compendium's symbol->b-number index first.

Run with the `gradi` env:
    python scripts/essentiality/screens.py
    python scripts/essentiality/screens.py --list
    python scripts/essentiality/screens.py --screen essential_ecoli_k12_knockout --dry-run

Writes one file per training set to `data/processed/essentiality/training_sets/`, schema
`key, label, source_id, features_from` (+ screen-specific extras). `key` is NOT a UniProt
accession for the five non-anchor screens -- see the comment at the write site.
"""

from __future__ import annotations

import argparse
import re
import sys
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable
from urllib.parse import unquote

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

RAW = REPO_ROOT / "data" / "raw"
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
# THE CLEAN TRAINING SETS LIVE IN THEIR OWN FOLDER, not among the 40-odd audit tables in
# evidence/. They are the input to every model on this axis, so "which files are the training
# sets" should be answerable by listing one directory rather than by knowing a filename prefix.
TRAINING_DIR = OUT_DIR / "training_sets"
NCBI_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "ncbi"
KP_STRAINS = REPO_ROOT / "data" / "source" / "ncbi" / "kp_strains"

COMPENDIUM = RAW / "other" / "essentiality" / "enterobacteriaceae_tradis" / "giant-tab_final.tsv"
KEIO_B = "Locus: Escherichia coli BW25113 (Keio)"      # b-numbers
BW_LOCUS = "Locus: Escherichia coli BW25113"            # BW25113_#### tags

ECL8 = "kpneumoniae__ecl8__GCA_000315385.1"
KPPR1 = "kpneumoniae__kppr1__GCF_000742755.1"

# The clade whose conservation we ship, and the cut. Conservation is strongly BIMODAL, so the clade
# choice IS the binning choice. Measured over the 4,301 b-numbered rows at an 80/20 cut:
#     Enterobacteriaceae   core 264  mid 100  non 3937     <- only 100 in the middle
#     Gammaproteobacteria  core 205  mid 530  non 3566     <- chosen
#     Bacteria             core 133  mid 557  non 3611
# ONE COLUMN was asked for, so the middle folds into non-core: 530 genes essential in SOME genomes
# are called non-core, which is the price of a binary. The raw 0-100 percentage ships beside the
# call in the screen table, so the cut can move without a re-ingest.
CLADE_COLUMN = "Gammaproteobacteria %essential"
CORE_CUT = 80.0

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 104) -> None:
    say(char * width)


# ------------------------------------------------------------------- shared identifier machinery

def ecoli_bnumber_map() -> dict[str, str]:
    """b-number -> uniprot_ac. E. coli locus tags ARE b-numbers (4,402 of 4,403)."""
    lt = P.load_locus_tags("ecoli")
    m = lt[lt["locus_tag"].astype(str).str.match(r"^b\d{4}$", na=False)]
    return dict(zip(m["locus_tag"], m["uniprot_ac"]))


def ecoli_symbol_map() -> dict[str, str]:
    """normalised gene symbol -> uniprot_ac, from EVERY name E. coli carries, not just the primary.

    `gene_name` holds one preferred name, `gene_name_all` the alternatives. A screen may use either,
    so both are indexed and the preferred name wins a collision.
    """
    df = P.load("ecoli")
    out: dict[str, str] = {}
    for col in ("gene_name_all", "gene_name"):           # primary LAST so it overwrites
        if col not in df.columns:
            continue
        for ac, names in zip(df["uniprot_ac"], df[col].fillna("")):
            for n in re.split(r"[;,\s]+", str(names)):
                if n.strip():
                    out[n.strip().lower()] = ac
    return out


def ecoli_sequence_map() -> dict[str, str]:
    """exact protein sequence -> uniprot_ac. The house join when identifier namespaces disagree."""
    df = P.load("ecoli")
    return dict(zip(df["sequence"].astype(str), df["uniprot_ac"]))


def compendium(usecols: list[str] | None = None) -> pd.DataFrame:
    if not COMPENDIUM.exists():
        sys.exit(f"FATAL missing {COMPENDIUM}. It is the Enterobacteriaceae-TraDIS compendium "
                 "(PMID 39207104), staged by v1 under data/raw/other/essentiality/.")
    return pd.read_csv(COMPENDIUM, sep="\t", low_memory=False, usecols=usecols)


def ncbi_tag_to_protein(label: str) -> pd.DataFrame:
    """(locus_tag, old_locus_tag, protein_id) for one NCBI assembly.

    NCBI puts `old_locus_tag` on the GENE feature and `protein_id` on the CDS, linked by ID/Parent.
    Both must be walked or the map comes out empty -- which looks exactly like a screen that does
    not join.
    """
    gff = NCBI_DIR / f"{label}.gff"
    if not gff.exists():
        sys.exit(f"FATAL {gff} missing -- run scripts/proteomes/download.py --tier D --only {label}")
    genes, cds = {}, []
    for line in gff.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].rstrip().split(";") if "=" in kv)
        if f[2] in ("gene", "pseudogene") and a.get("ID"):
            genes[a["ID"]] = (unquote(a.get("locus_tag", "")), unquote(a.get("old_locus_tag", "")))
        elif f[2] == "CDS" and a.get("protein_id"):
            cds.append({"parent": a.get("Parent", ""), "protein_id": a["protein_id"]})
    d = pd.DataFrame(cds).drop_duplicates("parent")
    d["locus_tag"] = d["parent"].map(lambda p: genes.get(p, ("", ""))[0])
    d["old_locus_tag"] = d["parent"].map(lambda p: genes.get(p, ("", ""))[1])
    return d


def ncbi_tag_map(label: str) -> dict[str, str]:
    """Every tag form -> protein_id, so a screen may key on the old or the current namespace.

    GENBANK (GCA) vs REFSEQ (GCF) IS PER-PAPER, NOT A DEFAULT. Screens that predate the current
    RefSeq annotation need GenBank, because PGAP's re-annotation drops submitter tags outright:

        BN373 / ECL8    GenBank 4,930/5,048 (97.7%)   RefSeq 4,767 (94.4%)   RefSeq adds 0
        KPNIH1          GenBank   424/424  (100.0%)   RefSeq   418 (98.6%)

    263 of RefSeq's 281 ECL8 misses were tags ABSENT from its annotation, not proteinless genes.
    But Paczosa keys on `VK055_RS*`, which IS RefSeq, and KPPR1 below uses the GCF assembly for
    exactly that reason. Check which namespace the paper used; do not apply either rule by reflex.
    """
    d = ncbi_tag_to_protein(label)
    m: dict[str, str] = {}
    for col in ("locus_tag", "old_locus_tag"):
        for tags, pid in zip(d[col], d["protein_id"]):
            for t in str(tags).split(","):
                if t.strip():
                    m[t.strip()] = pid
    return m


def kp_strain_map(label: str) -> dict[str, str]:
    """locus_tag -> itself, for a strain staged by strain_proteomes.py.

    That table IS keyed on the strain's own tags and carries the sequence, so the screen, the
    embedding and the label all share one namespace and the map is the identity on what exists.
    Its real job is to REJECT a tag the staged proteome has no sequence for.
    """
    p = KP_STRAINS / f"{label}.tsv"
    if not p.exists():
        sys.exit(f"FATAL {p} missing -- run scripts/essentiality/strain_proteomes.py")
    d = pd.read_csv(p, sep="\t")
    return {t: t for t in d["locus_tag"].astype(str)}


def deg_rows(dataset_id: str) -> pd.DataFrame:
    """One DEG dataset out of the labelled corpus. Every row carries its protein SEQUENCE."""
    lp = pd.read_csv(EVIDENCE_DIR / "labeled_proteins.tsv", sep="\t")
    d = lp[lp["deg_dataset_id"] == dataset_id]
    if not len(d):
        sys.exit(f"FATAL no rows for {dataset_id} -- run scripts/essentiality/deg_proteomes.py")
    return d


# ---------------------------------------------------------------------------------- the parsers

def s_keio() -> tuple[pd.DataFrame, str]:
    """PEC's curation of the Keio arrayed knockout collection -- the only non-transposon assay here.

    Column 9 is `Class(1:essential 2:noessential 3:unknown)`; column 3 packs synonyms including the
    b-number (`b0001,ECK0001,JW4367`). Unknown is DROPPED, never made negative.
    """
    p = RAW / "ecoli" / "essentiality" / "pec" / "PECData.dat"
    d = pd.read_csv(p, sep="\t", header=None, low_memory=False, skiprows=1,
                    names=[f"c{i}" for i in range(13)])
    d = d[pd.to_numeric(d["c1"], errors="coerce") == 1]          # feature type 1 = protein-coding
    cls = pd.to_numeric(d["c9"], errors="coerce")
    d = d.assign(label=np.where(cls == 1, 1, np.where(cls == 2, 0, np.nan)))
    d["uniprot_ac"] = d["c3"].astype(str).str.extract(r"\b(b\d{4})\b")[0].map(ecoli_bnumber_map())
    return d[["uniprot_ac", "label", "c2"]].rename(columns={"c2": "source_id"}), "b-number"


def s_gerdes() -> tuple[pd.DataFrame, str]:
    """RETIRED -- Gerdes 2003 genetic footprinting (DEG1018). Parser kept, not shipped.

    DROPPED ON THE PROJECT OWNER'S CALL, 2026-09-21, and the two independent diagnostics agreed:

        ribosome recall   0.538   worst of the ten; Keio 0.774, Goodall 0.961
        grouped AUROC     0.733   worst by 0.14; the next lowest is 0.860
        grouped AUPR      0.341   worst by 0.15; the next lowest is 0.491

    The label polarity is CORRECT (8x enrichment over the dispensables), so this is not a parsing
    error -- it is a 2003 assay that is genuinely noisy in both directions. It over-calls overall
    (base rate 0.142, more than double Keio's on the same organism) while still missing 11 core
    ribosomal proteins that Keio calls essential (rplE, rplP, rplQ, rplV, rpmA, rpmC, rpmD, rpsC,
    rpsD, rpsG, rpsK).

    NOTE FOR ANYONE TOUCHING `merge.py`: DEG1018 is ALSO one of the two screens behind `deg_ess`
    on E. coli, and this retirement does NOT change that column. Since 490 of its 695 E. coli
    essential calls rest on a single screen, dropping DEG1018 there too would move `deg_ess`
    substantially -- a separate decision, not a consequence of this one.
    """
    d = deg_rows("DEG1018").copy()
    d["uniprot_ac"] = d["sequence"].astype(str).map(ecoli_sequence_map())
    return (d[["uniprot_ac", "essential", "locus_tag"]]
            .rename(columns={"essential": "label", "locus_tag": "source_id"}), "exact sequence")


def s_goodall() -> tuple[pd.DataFrame, str]:
    """Goodall 2018 TraDIS on BW25113 -- the standard modern E. coli reference.

    `Unclear` (155) is DROPPED: the authors saying the insertion profile does not separate. Making
    it negative would manufacture negatives out of ignorance.

    TWO ROUTES, the stable one first. The compendium's `Gene_Keio` is itself a symbol->b-number
    index built by its authors, so that route lands on a stable identifier; our proteome's symbol
    index is the fallback. Measured on the 4,158 labelled rows: compendium 96.4%, ours 93.2%,
    union 98.6%.
    """
    x = pd.ExcelFile(RAW / "ecoli" / "essentiality" / "goodall2018_Ec_BW25113" /
                     "mbo001183726st4.xlsx")
    d = x.parse(x.sheet_names[0], header=1)
    d = d.assign(label=np.where(d["Essential"] == True, 1,                         # noqa: E712
                                np.where(d["Non-essential"] == True, 0, np.nan)))  # noqa: E712
    sym = d["Gene"].astype(str).str.strip().str.lower()
    c = compendium(["Gene_Keio", KEIO_B])
    c = c[c[KEIO_B].astype(str).str.match(r"^b\d{4}$", na=False)]
    sym2b = {str(s).strip().lower(): b for s, b in zip(c["Gene_Keio"], c[KEIO_B])
             if str(s) not in ("nan", "-", "")}
    via_b = sym.map(sym2b).map(ecoli_bnumber_map())
    d["uniprot_ac"] = via_b.fillna(sym.map(ecoli_symbol_map()))
    d["join_route"] = np.where(via_b.notna(), "symbol->b-number", "symbol->proteome")
    return (d[["uniprot_ac", "label", "Gene", "Insertion Index Score", "join_route"]]
            .rename(columns={"Gene": "source_id", "Insertion Index Score": "insertion_index"}),
            "symbol->b-number")


def s_ghomi_bw() -> tuple[pd.DataFrame, str]:
    """RETIRED -- the compendium's DESeq run on BW25113 did not converge. Parser kept, not shipped.

    It looked like a free second opinion on Goodall's data under a uniform pipeline. The ribosome
    control says otherwise: it calls 6 of 53 ribosomal proteins essential, against 0.77 for Keio
    and 0.96 for Goodall on the same organism. The cause is in the file:

        padj is NaN for 2,310 of 4,256 genes (54.3%)     BN373/ECL8, same pipeline: 7.0%
        of the 189 genes with ZERO insertion sites,      BN373: 91 of 91 -> Reduced
            108 are called `Unchanged`

    A gene with no insertions is the most essential thing a transposon library can show. DESeq
    cannot compute a statistic for it, so those genes fall through to `Unchanged` -- the label is
    INVERTED for exactly the genes that matter most. ECL8 is unaffected, which is why that column
    ships and this one does not.

    Nothing is lost: this is Goodall's data, and `essential_ecoli_bw25113_tradis_goodall` carries
    the authors' own calls on it.
    """
    d = pd.read_csv(COMPENDIUM.parent / "BW25113.out.DESeq.tsv", sep="\t")
    bridge = compendium([KEIO_B, BW_LOCUS])
    bridge = bridge[(bridge[KEIO_B].astype(str).str.match(r"^b\d{4}$", na=False))
                    & (bridge[BW_LOCUS].astype(str) != "-")]
    tag2b = dict(zip(bridge[BW_LOCUS].astype(str), bridge[KEIO_B].astype(str)))
    d["uniprot_ac"] = d["id"].astype(str).map(tag2b).map(ecoli_bnumber_map())
    d["label"] = np.where(d["Essentiality"].astype(str) == "Reduced", 1,
                          np.where(d["Essentiality"].astype(str) == "Unchanged", 0, np.nan))
    return (d[["uniprot_ac", "label", "id", "log2FoldChange", "padj"]]
            .rename(columns={"id": "source_id"}), "BW25113_ -> b-number")


def s_choe() -> tuple[pd.DataFrame, str]:
    """Choe 2023 Tn-seq on BW25113 -- LB ARM ONLY.

    The sheet also carries an M9-glucose arm. That is a growth CONDITION: its extra 131 "essential"
    genes are auxotrophies, genes a cell needs only because the medium withholds their product.
    Condition-dependent data is out of scope for now, so the M9 column is read and discarded here
    rather than silently averaged in.
    """
    d = pd.read_excel(RAW / "ecoli" / "essentiality" / "choe2023_ecoli" /
                      "msystems.00896-22-s0002.xlsx", sheet_name="Table S1", header=[0, 1])
    lt = d[("Locus Tag", "Unnamed: 5_level_1")].astype(str).str.strip()
    ess = d[("LB medium", "Essentiality")].astype(str)
    out = pd.DataFrame({"source_id": lt,
                        "label": np.where(ess == "E", 1, np.where(ess == "NE", 0, np.nan))})
    out["uniprot_ac"] = lt.map(ecoli_bnumber_map())
    return out[["uniprot_ac", "label", "source_id"]], "b-number"


def _deg_strain(dataset_id: str) -> tuple[pd.DataFrame, str]:
    """A DEG dataset on a NON-anchor strain, keyed on its own `fasta_id`.

    `fasta_id` is verbatim the header of `data/source/ncbi/deg_proteomes/<id>.faa`
    (`lcl|HG941718.1_prot_CDN80371.1_1`), which is what the embedding script keys its rows on --
    so label and feature vector align with no join at all. Measured against the E. coli anchor by
    exact sequence, DEG1048 shares 13.7% and DEG1056 23.0%: these are genuinely different
    proteomes and must be embedded in their own namespace, not mapped onto K-12.
    """
    d = deg_rows(dataset_id).copy()
    return (d[["fasta_id", "essential", "locus_tag"]]
            .rename(columns={"fasta_id": "uniprot_ac", "essential": "label",
                             "locus_tag": "source_id"}), "DEG fasta_id")


def s_st131():
    """DEG1048 -- E. coli ST131 EC958, the pandemic multidrug-resistant ExPEC lineage."""
    return _deg_strain("DEG1048")


def s_o157h7():
    """DEG1056 -- E. coli O157:H7, EHEC. Highest base rate of the set (0.194).

    DEG ships 14.1% of this dataset's "sequences" as the literal string `Not available now.`;
    `labeled_proteins.tsv` already validates the amino-acid alphabet, so those never reach here.
    """
    return _deg_strain("DEG1056")


def s_ecl8() -> tuple[pd.DataFrame, str]:
    """The compendium's DESeq call on K. pneumoniae ECL8 (K2-ST375) -- the first Kp endpoint.

    Keyed on `BN373_#####`, the deposited GenBank assembly's `old_locus_tag`. The authors' own
    `ecl8_*` tags annotate an undeposited assembly and are unjoinable by all three routes -- see
    UNJOINABLE. Same organism, same TraDIS data, uniformly reprocessed; what differs is the
    analysis (DESeq vs the authors' bimodal fit) and the positive count (562 vs 373).
    """
    d = pd.read_csv(COMPENDIUM.parent / "BN373.out.DESeq.tsv", sep="\t")
    d["uniprot_ac"] = d["id"].astype(str).map(ncbi_tag_map(ECL8))
    d["label"] = np.where(d["Essentiality"].astype(str) == "Reduced", 1,
                          np.where(d["Essentiality"].astype(str) == "Unchanged", 0, np.nan))
    return (d[["uniprot_ac", "label", "id", "log2FoldChange", "padj"]]
            .rename(columns={"id": "source_id"}), "BN373_ old_locus_tag")


def _bruchmann(table: str, mapper) -> tuple[pd.DataFrame, str]:
    """Bruchmann 2021 TraDIS, INPUT POOL column only.

    `essential_In` is essentiality in the input library -- in vitro, before infection. The same
    sheets carry `2hpi` and `6hpi` in-host columns, which are condition-dependent and deliberately
    unused.

    'ambiguous' is DROPPED (as is the file's own misspelling 'ambigiuous', which is why the label
    is built by whitelisting the two decided states rather than by negating 'essential'). It is the
    authors saying they cannot tell, not a measured negative.
    """
    d = pd.read_csv(RAW / "kpneumoniae" / "essentiality" / "bruchmann2021" / table, header=1)
    e = d["essential_In"].astype(str).str.strip()
    d["label"] = np.where(e == "essential", 1, np.where(e == "nonessential", 0, np.nan))
    d["uniprot_ac"] = d["locus_tag"].astype(str).str.strip().map(mapper)
    return (d[["uniprot_ac", "label", "locus_tag"]].rename(columns={"locus_tag": "source_id"}),
            "locus_tag")


def s_rh201207():
    """Kp RH201207 -- a clinical ST258-adjacent isolate. Its sequences come from GCF_905477585.1,
    because the GenBank assembly carries gene features with no CDS and no proteins at all."""
    return _bruchmann("TableS5.csv", kp_strain_map("rh201207_bruchmann"))


def s_atcc43816():
    """Kp ATCC 43816 (KPPR1) -- the standard murine-pneumonia lab strain. RefSeq tags."""
    return _bruchmann("TableS6.csv", ncbi_tag_map(KPPR1))


def s_core() -> tuple[pd.DataFrame, str]:
    """Core vs non-core essentiality across 12 Gammaproteobacteria genomes.

    NOT a per-organism measurement: it is the fraction of the clade's genomes in which the gene is
    essential, a property of the gene FAMILY. That makes it the one label here a per-protein
    embedding is straightforwardly the right feature for -- the embedding encodes the family.

    The computation is the compendium's; only the binning is ours, and the raw percentage ships
    beside the call so the cut can be revisited without re-running anything.
    """
    d = compendium([KEIO_B, CLADE_COLUMN])
    d = d[d[KEIO_B].astype(str).str.match(r"^b\d{4}$", na=False)].copy()
    v = pd.to_numeric(d[CLADE_COLUMN], errors="coerce")
    d["pct_essential"] = v
    d["label"] = np.where(v >= CORE_CUT, 1, np.where(v.notna(), 0, np.nan))
    d["uniprot_ac"] = d[KEIO_B].astype(str).map(ecoli_bnumber_map())
    return (d[["uniprot_ac", "label", KEIO_B, "pct_essential"]]
            .rename(columns={KEIO_B: "source_id"}), "b-number")


# ------------------------------------------------------------------------------- the registry

@dataclass
class Screen:
    column: str          # the name this becomes in the output matrix
    comment: str         # one line, human-readable -- goes straight into the summary table
    organism: str
    assay: str
    features: str        # WHICH proteome's embeddings this set's rows live in
    parser: Callable
    # Where the raw files live, DECLARED rather than inferred. `registry.py` used to guess it by
    # substring, which matched `essential_ecoli_bw25113_tradis_goodall` to `data/source/go` --
    # a false positive is worse than a miss, because it looks like a resolved provenance link.
    source_dir: str = ""
    floor: float = 0.90  # below this the join is structurally broken; exit non-zero


SCREENS: list[Screen] = [
    Screen("essential_kpneumoniae_ecl8_tradis",
           "Kp ECL8 (K2-ST375), compendium DESeq reanalysis of TraDIS",
           "K. pneumoniae ECL8", "TraDIS/DESeq", ECL8, s_ecl8,
           source_dir="data/raw/other/essentiality/enterobacteriaceae_tradis"),
    Screen("essential_kpneumoniae_rh201207_tradis",
           "Kp RH201207 clinical isolate, Bruchmann TraDIS input pool",
           "K. pneumoniae RH201207", "TraDIS", "rh201207_bruchmann", s_rh201207,
           source_dir="data/raw/kpneumoniae/essentiality/bruchmann2021"),
    Screen("essential_kpneumoniae_atcc43816_tradis",
           "Kp ATCC 43816 (KPPR1), Bruchmann TraDIS input pool",
           "K. pneumoniae ATCC 43816", "TraDIS", KPPR1, s_atcc43816,
           source_dir="data/raw/kpneumoniae/essentiality/bruchmann2021"),
    Screen("essential_ecoli_k12_knockout",
           "Ec K-12, Keio arrayed single-gene knockouts -- the only non-transposon assay",
           "E. coli K-12", "arrayed knockout", "ecoli", s_keio, floor=0.95,
           source_dir="data/raw/ecoli/essentiality/pec"),
    Screen("essential_ecoli_bw25113_tradis_goodall",
           "Ec BW25113, Goodall 2018 TraDIS -- the standard modern reference",
           "E. coli BW25113", "TraDIS", "ecoli", s_goodall,
           source_dir="data/raw/ecoli/essentiality/goodall2018_Ec_BW25113"),
    Screen("essential_ecoli_bw25113_tnseq_choe",
           "Ec BW25113, Choe 2023 Tn-seq, LB arm only (the M9 arm is a condition)",
           "E. coli BW25113", "Tn-seq", "ecoli", s_choe,
           source_dir="data/raw/ecoli/essentiality/choe2023_ecoli"),
    Screen("essential_ecoli_st131_tradis",
           "Ec ST131 EC958, pandemic MDR ExPEC lineage -- a non-K-12 background",
           "E. coli EC958", "TraDIS", "DEG1048", s_st131,
           source_dir="data/source/deg"),
    Screen("essential_ecoli_o157h7_tnseq",
           "Ec O157:H7 EHEC pathotype -- highest base rate of the set",
           "E. coli O157:H7", "Tn-seq", "DEG1056", s_o157h7,
           source_dir="data/source/deg"),
    Screen("core_essential_gammaproteobacteria",
           "Core vs non-core across 12 clade genomes -- a gene-FAMILY property, not per-organism",
           "Gammaproteobacteria", "clade conservation", "ecoli", s_core,
           source_dir="data/raw/other/essentiality/enterobacteriaceae_tradis"),
]

# Downloaded and parseable, deliberately not used. Decisions, not gaps.
HELD_BACK = {
    "ramage2017_kpnih1":   "positives-only -- 424 listed genes, no measured negatives",
    "paczosa2020_kppr1":   "positives-only (310 hits) AND condition-dependent (in-vivo, neutropenic)",
    "short2020_serum":     "condition-dependent -- serum resistance, 4 Kp strains",
    "rome2026_njst258":    "condition-dependent -- iron-depleted medium",
    "choe2023_m9":         "condition-dependent -- minimal medium; its extra essentials are auxotrophies",
    "bruchmann_2hpi_6hpi": "condition-dependent -- in-host timepoints from the same TraDIS sheets",
    "deg1019_keio":        "duplicate -- the Keio collection again under DEG's curation",
    "gerdes2003_mg1655":   "retired on the owner's call -- ribosome recall 0.538, AUROC 0.733, "
                           "AUPR 0.341, worst of the ten on every diagnostic. Correctly "
                           "polarised, just a noisy 2003 assay. Still feeds `deg_ess` via "
                           "merge.py, which is untouched.",
}

# Measured and refuted. The evidence is in docs/essentiality_screens.md; do not re-derive these.
UNJOINABLE = {
    "eichelberger2024_ecl8": "annotates an UNDEPOSITED assembly: tag suffixes 475/5,165, "
                             "coordinates 0/5,074, symbols 3.7%. The ecl8 column carries the same "
                             "TraDIS measurement instead.",
    "cain2017_njst258":      "DOES NOT EXIST. The DOI in circulation resolves to an unrelated "
                             "neurovascular paper; no NJST258 TraDIS screen was ever published.",
    "jung2019_mh258":        "ENA project PRJEB31265 exists but has released zero records.",
    "ghomi_bw25113_deseq":   "the compendium's DESeq run on BW25113 did not converge -- padj NaN "
                             "for 54.3% of genes and 108 of 189 zero-insertion genes called "
                             "`Unchanged`, giving ribosome recall 0.113 against Goodall's 0.961 "
                             "on the same data. ECL8's run (7.0% NaN) is unaffected.",
    "fitness_browser_rbtnseq": "RB-TnSeq cannot see essential genes BY CONSTRUCTION -- no "
                               "insertions survive, so they are absent rather than extreme.",
    "bvbrc_essentiality":    "FBA PREDICTION, not measurement -- circular as an ML label.",
}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    names = [s.column for s in SCREENS]
    ap.add_argument("--screen", nargs="+", choices=names, default=names)
    ap.add_argument("--list", action="store_true", help="show every screen, held-back and dead end")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("scripts/essentiality/screens.py -- published screens -> one training set per output column")
    rule("=")
    say(f"  in     {COMPENDIUM.parent.relative_to(REPO_ROOT)}/ + data/raw/{{ecoli,kpneumoniae}}/essentiality/")
    say(f"  out    {TRAINING_DIR.relative_to(REPO_ROOT)}/<column>.tsv")
    say(f"  clade  {CLADE_COLUMN}  (core >= {CORE_CUT:.0f}%)")
    say("  join   stable identifiers or exact sequence only -- v1 used symbols and lost most of it")

    if a.list:
        rule()
        say(f"IN USE -- {len(SCREENS)} columns")
        rule()
        for s in SCREENS:
            say(f"  {s.column:42s} {s.comment}")
            say(f"  {'':42s}   {s.organism} | {s.assay} | features from {s.features}")
        rule()
        say("HELD BACK -- downloaded and parseable, deliberately unused")
        rule()
        for k, v in HELD_BACK.items():
            say(f"  {k:24s} {v}")
        rule()
        say("MEASURED AND REFUTED -- do not re-derive")
        rule()
        for k, v in UNJOINABLE.items():
            say(f"  {k:24s} {v}")
        rule("=")
        return

    rows, tables = [], {}
    rule()
    say(f"  {'column':42s} {'rows':>6} {'pos':>6} {'base':>7} {'join':>7}  key")
    rule()
    for s in SCREENS:
        if s.column not in a.screen:
            continue
        df, key = s.parser()
        labelled = df[df["label"].notna()]
        # one row per protein: a screen may list paralogs separately; keep the strongest call
        joined = (labelled[labelled["uniprot_ac"].notna()]
                  .sort_values("label", ascending=False).drop_duplicates("uniprot_ac"))
        rate = len(joined) / len(labelled) if len(labelled) else 0.0
        pos = int(joined["label"].sum())
        say(f"  {s.column:42s} {len(joined):6d} {pos:6d} {pos / len(joined):7.3f} "
            f"{rate:7.3f}  {key}")
        if rate < s.floor:
            sys.exit(f"FAILED {s.column}: join rate {rate:.3f} is below the {s.floor:.2f} floor. "
                     "Something structural changed -- do NOT ship a decimated training set.")
        if set(joined["label"].unique()) != {0.0, 1.0}:
            sys.exit(f"FAILED {s.column}: not both classes present "
                     f"({sorted(joined['label'].unique())})")
        rows.append({"dataset": s.column.removeprefix("essential_"),
                     "short_comment": s.comment,
                     "column_name": s.column,
                     "num_positives": pos,
                     "num_total": len(joined),
                     "base_rate": round(pos / len(joined), 4),
                     "join_rate": round(rate, 4),
                     "join_key": key,
                     "organism": s.organism,
                     "assay": s.assay,
                     "features_from": s.features})
        tables[s.column] = joined

    if a.dry_run:
        rule("=")
        say("dry run: nothing written.")
        return

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    TRAINING_DIR.mkdir(parents=True, exist_ok=True)
    rule()
    say("TRAINING SETS")
    rule()
    feats = {s.column: s.features for s in SCREENS}
    for name, df in tables.items():
        # UNIFORM SCHEMA: key, label, source_id, features_from, then any screen-specific extras.
        #
        # `key` was called `uniprot_ac`, and that name was ACTIVELY MISLEADING: for 5 of the 9
        # screens the values are not UniProt accessions at all but `CCN31837.1` (ECL8 EMBL),
        # `WP_038431262.1` (KPPR1 RefSeq), `KPNRH_00001` (a locus tag) or
        # `lcl|HG941718.1_prot_...` (a DEG FASTA header). Joining one of those to an anchor table
        # on `uniprot_ac` returns silence, not an error -- measured at 0 of 4,930 for ECL8.
        #
        # `features_from` names the proteome whose embedding matrix `key` indexes, so no consumer
        # has to infer it.
        df = df.rename(columns={"uniprot_ac": "key"})
        df["features_from"] = feats[name]
        lead = ["key", "label", "source_id", "features_from"]
        df = df[lead + [c for c in df.columns if c not in lead]]
        out = TRAINING_DIR / f"{name}.tsv"
        df.to_csv(out, sep="\t", index=False)
        say(f"  {out.relative_to(REPO_ROOT)}  ({len(df)} x {df.shape[1]})")
    audit = pd.DataFrame(rows)
    audit["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    audit.to_csv(EVIDENCE_DIR / "screen_join_audit.tsv", sep="\t", index=False)
    say(f"  {(EVIDENCE_DIR / 'screen_join_audit.tsv').relative_to(REPO_ROOT)}")
    rule("=")
    say(f"screens complete -- {len(tables)} training sets, every one with both classes.")
    rule("=")


if __name__ == "__main__":
    main()
