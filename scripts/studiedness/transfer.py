"""Studiedness by sequence: the best-characterised SwissProt homolog of every protein.

The anchors' own literature cannot carry this axis. Measured on UniProt 2026_03: **5,710 of 5,728
K. pneumoniae proteins carry exactly one PubMed id** (the genome paper, 7 distinct values over the
whole proteome) and 2,532 of 2,889 S. aureus proteins carry none. So the usable number is
transferred from homologs -- which is also the house rule, *"map by sequence, not by accession,
when reaching an external database"*.

    data/processed/studiedness/evidence/own_<species>.tsv        this accession's own signals
    data/processed/studiedness/evidence/transfer_<species>.tsv   the chosen donor + its counts
    data/processed/studiedness/evidence/definition.tsv           what the shipped number means
    data/processed/studiedness/evidence/route_comparison.tsv     SwissProt vs the free route
    data/processed/studiedness/evidence/donor_scope_comparison.tsv  bacteria-only vs any donor
    data/processed/studiedness/evidence/decoy_calibration.tsv    what each floor lets through
    data/processed/studiedness/evidence/floor_sensitivity.tsv    the identity floor, swept
    data/processed/studiedness/evidence/control_ecoli_heldout.tsv   the control that matters

WHY SWISSPROT AND NOT THE ORTHOLOG TABLE WE ALREADY HAVE
---------------------------------------------------------
Measured from `src/orthology.load`: the free four-species panel reaches an E. coli or human
ortholog for only **Kp 68.5% and Sa 46.6%**. *S. aureus* is Gram-positive, so less than half of it
has an E. coli ortholog at all -- the free route structurally abandons most of that proteome.
SwissProt is 575,748 reviewed entries over ~14,000 species and holds the Gram-positive literature
the panel cannot contain. The free route is still scored here, into `route_comparison.tsv`,
because it costs nothing and the *Identifier mapping* rule says to measure every route. **Nothing
is merged.**

THE DONOR IS THE MOST-CITED HIT ABOVE THE FLOOR, NOT THE CLOSEST
-----------------------------------------------------------------
The question is "how well studied is the best-characterised thing this protein is homologous to",
so among hits clearing the floor the winner is the one with the most literature (ties broken by
annotation score, then identity). v1 used the same rule. **Preferring the closest band first was
tried and is measurably worse** -- held-out control 0.4426 against 0.4666 at the 40% floor --
because the quantity being carried is literature volume, and the best-cited homolog predicts
"this family is studied" better than the nearest one does. The nearest hit ships alongside in
`nearest_*` columns so the choice is auditable and never silent.

DONOR SCOPE: NOT EUKARYOTIC, AND ALL THREE ARE ALWAYS MEASURED
----------------------------------------------------------------
`--donor-scope prokaryotic` (default) admits any donor whose taxonomic lineage does NOT contain
"Eukaryota (domain)" -- Bacteria, Archaea **and phages**. `bacteria` is strict Bacteria only;
`any` lets a human, plant or fungal homolog donate.

**Phages are why the default is not simply "Bacteria", and that was measured, not assumed.** Of
the 93 proteins a strict bacterial rule stranded with no donor at all, **70 (75%) lost theirs to
a virus** -- *Escherichia* phage lambda and P1. A K. pneumoniae prophage protein whose
best-characterised relative is a lambda protein is emphatically not novel; lambda is among the
most studied systems in biology, and scoring it 0 would be wrong in the one direction this axis
must not get wrong.

**The restriction exists because "most-cited homolog anywhere" has a predictable failure.** Under
`any`, five of K. pneumoniae's ten highest-scoring proteins took human donors -- HSPD1, HADHA,
CTPS1, LONP1, AFG3L2 -- so the claim became "this protein is well studied because its human
mitochondrial homolog is". That is defensible for *is this family understood?* and useless for
*is this a novel antibacterial target?*, which is the question this project asks.

It is also the mistake the ligands axis already measured and recorded: an unrestricted non-human
bucket gave **424** apparent potent K. pneumoniae proteins against a true **175**, with a rat
protein winning bacterial slots at ~30% identity. Hence lineage, not organism name and not
"is not human" -- strains are filed under their own taxids.

**Both scopes are computed on every run** and land in `evidence/donor_scope_comparison.tsv`, with
`n_papers_uniprot_bacteria` and `n_papers_uniprot_any` side by side in the transfer table and
`donor_scope` naming the shipped one. The held-out control is also run under both, which is the
arbiter: it asks which scope better recovers E. coli's own measured literature from donors it was
not allowed to see.

THE FLOOR: 40%, AND WHY NOT 25%
--------------------------------
Bands are nested the way the `ligands/` axis does it -- **>= 95% direct, >= 60% close, >= 40%
homolog** -- each also requiring >= 50% coverage on query and subject.

Two things were measured before settling on 40%:
  * **Decoys do not decide it.** Composition-preserving shuffles of our own 13,020 sequences match
    SwissProt at **0.0% at every floor down to 20%** (`decoy_calibration.tsv`), so the floor is
    not defending against spurious homology -- DIAMOND's e-value already does.
  * **It is a coverage-versus-fidelity trade-off, and it is monotonic** (`floor_sensitivity.tsv`):

        floor   Kp cov   Sa cov   held-out control rho
         25%     77.5%    68.9%          0.3265
         40%     66.0%    55.9%          0.4666      <-- shipped
         60%     56.5%    40.5%          0.5060

    25% buys 11 points of K. pneumoniae coverage for 0.14 of control. 40% is CLAUDE.md's standing
    annotation-transfer threshold and nothing measured here justifies departing from it. The whole
    sweep ships, so the floor can be re-chosen without re-running DIAMOND.

DO NOT MAXIMISE THE CONTROL -- IT HAS A CEILING WELL BELOW 1, AND CHASING IT BREAKS THE AXIS
----------------------------------------------------------------------------------------------
The control correlates the transferred family score against E. coli's OWN literature. Those are
deliberately different quantities: a protein with 3 papers of its own whose human homolog has 300
*should* score low on `own` and high on `family` -- that gap is the entire point of the axis. So
perfect agreement is not the target and rho cannot reach 1.

Worse, the control rises monotonically with the floor (0.59 at 95%) for a reason that is close to
circular: the exclusion removes *Escherichia* only, so at a high floor the surviving donors are
largely Salmonella and Shigella near-duplicates whose publication counts track E. coli's because
they are effectively the same proteins. **Tuning the floor to maximise rho would drive it to 95%,
where `family` collapses onto `own` and the axis stops doing anything.** The control is a
pass/fail sanity check that transfer carries real signal on held-out data -- not a score to
optimise.

`--very-sensitive`, because the error direction matters: **under-detecting homology makes a target
look more novel than it is**, which is the failure this axis must not make.

`--max-target-seqs 500`, AND THAT NUMBER IS LOAD-BEARING
---------------------------------------------------------
DIAMOND returns the top k hits **by bitscore**, i.e. the CLOSEST relatives -- but this stage wants
the **best-cited** one, which for a conserved protein is usually a distant model-organism entry.
At k=50 the candidate list for *S. aureus* GroEL was 50 Staphylococci and Bacilli, so the chosen
donor was *B. subtilis* GroEL with **9 papers** and E. coli GroEL never appeared at all. The
effect is systematic and points the wrong way: it understates studiedness precisely for the most
conserved -- and therefore most studied -- families. Every spot-check gene was at the cap
(`n_candidates == 50`), which is what exposed it. `n_candidates` ships per protein and the run
prints what fraction of queries are still capped; if that is large, raise k again.

THE CONTROL: E. COLI WITH EVERY E. COLI DONOR REMOVED
-------------------------------------------------------
E. coli is the only anchor whose own literature is real, so it is the only place the transfer
mechanism can be tested at all. The control recomputes `n_papers_uniprot_prokaryotic` for E. coli with
**every Escherichia donor struck out of SwissProt**, then correlates it against E. coli's own
measured `n_papers_uniprot_own`. That is a genuine held-out test of the whole mechanism -- the same
shape as `embeddings/prott5.py` validating against UniProt's published vectors. The run exits
non-zero if it falls below the floor.

**It reads LOWER than it used to, and that is the circularity being removed.** While the axis
shipped a 0-1 composite the control read 0.5411, but the blend put a 0.4-weighted annotation term
on both sides and own-vs-donor annotation score correlates at **spearman 0.950** -- so part of
that number was annotation agreeing with itself, not literature transferring. Counting papers
only gives 0.328, which is the honest figure.

Run with the `gradi` env, after `fetch.py` and `gene2pubmed.py`.
DIAMOND is borrowed from `gradi-ortho` (GRADI_DIAMOND_BIN overrides) -- it has no arm64 build and
must not be installed into `gradi`.
  python scripts/studiedness/transfer.py
  python scripts/studiedness/transfer.py --limit 200          # smoke, writes only scratch/
  python scripts/studiedness/transfer.py --donor-scope any     # unrestricted donors
  python scripts/studiedness/transfer.py --refresh --threads 8
"""

from __future__ import annotations

import argparse
import gzip
import os
import random
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import orthology as O  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import studiedness as S  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
LIT_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "literature"
SPROT_FASTA_GZ = LIT_DIR / "uniprot_sprot.fasta.gz"
SPROT_META = LIT_DIR / "swissprot_meta.tsv.gz"
G2P_COUNTS = OUT_DIR / "scratch" / "gene2pubmed_counts.tsv"
# PubTator3, both routes. Written by scripts/studiedness/pubtator.py; see its docstring for why
# GeneID and symbol are different datasets for bacteria. Both are MEASURED ALTERNATIVES: they are
# carried per donor and scored by the held-out control, and neither feeds donor selection or the
# deliverable. Absent files degrade to zeros with a printed note, exactly like gene2pubmed.
PT_GENEID_COUNTS = OUT_DIR / "scratch" / "pubtator_geneid_counts.tsv"
PT_OWN_PMIDS = OUT_DIR / "scratch" / "pubtator_own_pmids.tsv"
PT_SYMBOL_COUNTS = OUT_DIR / "scratch" / "pubtator_symbol_counts.tsv"
DEFAULT_DIAMOND_DIR = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

SPECIES = ("kpneumoniae", "ecoli", "saureus")
# Species-rank taxids for the three anchors, for the `exact` same-species route.
# Matched on the organism-name BINOMIAL, the same way src/precedents.py does it -- no taxonomy
# dump needed, and the two axes stay literally comparable.
ANCHOR_BINOMIAL = {"kpneumoniae": "Klebsiella pneumoniae",
                   "ecoli": "Escherichia coli",
                   "saureus": "Staphylococcus aureus"}
HIT_COLS = ["qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore", "evalue"]
FLOOR_SWEEP = (20.0, 25.0, 30.0, 40.0, 60.0, 95.0)
# Donor scopes, all computed every run. `prokaryotic` is "not eukaryotic" -- Bacteria, Archaea
# AND phages; see the docstring for why phages must be in.
SCOPES = ("prokaryotic", "bacteria", "any")

# The control floor, RE-DERIVED for the count-based definition on 2026-09-22 and then frozen.
#
# It was 0.45 while the axis shipped a 0-1 composite, and that number was INFLATED by a
# circularity: the blend put a 0.4-weighted annotation term on both sides of the control, and
# own-vs-donor annotation score correlates at spearman 0.950. Dropping the blend removed the
# circularity and the honest literature-only signal is 0.322-0.339 across donor scopes.
#
# 0.25 leaves headroom below that while still being far above chance -- the guard exists to catch
# a broken join (which would read ~0), not to certify a particular rho.
CONTROL_RHO_FLOOR = 0.25

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def ensure_diamond() -> str:
    override = os.environ.get("GRADI_DIAMOND_BIN")
    if shutil.which("diamond") is None:
        for d in ([Path(override)] if override else [DEFAULT_DIAMOND_DIR]):
            if (d / "diamond").exists():
                os.environ["PATH"] = f"{d}{os.pathsep}{os.environ['PATH']}"
                break
    exe = shutil.which("diamond")
    if exe is None:
        sys.exit(
            "diamond not found.\n"
            f"  looked on PATH and in {override or DEFAULT_DIAMOND_DIR}\n"
            "  it lives in the `gradi-ortho` env (osx-64, no arm64 build); point GRADI_DIAMOND_BIN\n"
            "  at a directory containing it. The anchors are bibliometrically dark, so the\n"
            "  sequence route is the ONLY way they get a literature number -- this cannot run\n"
            "  without it."
        )
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    return f"{exe}  ({(out.stdout or out.stderr).strip().splitlines()[0]})"


def split_ids(cell: object) -> list[str]:
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    return [x.strip() for x in str(cell).split(";") if x.strip()]


# ---------------------------------------------------------------- inputs

def gene2pubmed_counts() -> dict[str, int]:
    if not G2P_COUNTS.exists():
        say("  NOTE gene2pubmed counts absent -- donor literature will be UniProt-only, which "
            "is ~4x coarser on well-studied organisms.")
        say(f"       run scripts/studiedness/gene2pubmed.py to build "
            f"{G2P_COUNTS.relative_to(REPO_ROOT)}")
        return {}
    df = pd.read_csv(G2P_COUNTS, sep="\t", dtype={"geneid": str, "n_pubs": int})
    say(f"  gene2pubmed counts for {len(df):,} GeneIDs")
    return dict(zip(df["geneid"], df["n_pubs"]))


def pubtator_counts() -> tuple[dict[str, int], dict[str, int]]:
    """PubTator3 counts by GeneID and by gene symbol. Either may be absent."""
    out = []
    for path, key, label in ((PT_GENEID_COUNTS, "geneid", "pubtator geneid"),
                             (PT_SYMBOL_COUNTS, "symbol", "pubtator symbol")):
        if not path.exists():
            say(f"  NOTE {label} counts absent -- run scripts/studiedness/pubtator.py")
            out.append({})
            continue
        df = pd.read_csv(path, sep="\t", dtype={key: str, "n_pubs": int})
        say(f"  {label:<16} counts for {len(df):,} {key}s")
        out.append(dict(zip(df[key], df["n_pubs"])))
    return out[0], out[1]


def load_swissprot(counts: dict[str, int],
                   pt_geneid: dict[str, int] | None = None,
                   pt_symbol: dict[str, int] | None = None) -> pd.DataFrame:
    """SwissProt metadata with the union literature count per entry."""
    if not SPROT_META.exists():
        sys.exit(f"FATAL missing {SPROT_META.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    t0 = time.time()
    sp = pd.read_csv(SPROT_META, sep="\t", dtype=str, keep_default_na=False)
    sp = sp.rename(columns={"Entry": "donor_ac", "Annotation": "donor_annotation_score",
                            "Protein existence": "donor_protein_existence",
                            "PubMed ID": "lit_pubmed_id", "GeneID": "geneid",
                            "Organism (ID)": "donor_taxid", "Organism": "donor_organism",
                            "Gene Names (primary)": "donor_gene",
                            "Protein names": "donor_protein_name"})
    sp["donor_annotation_score"] = pd.to_numeric(sp["donor_annotation_score"],
                                                 errors="coerce").fillna(1.0)
    # THE SHIPPED COUNT IS SWISSPROT-CURATED ONLY: papers a UniProt curator read and used. One
    # consistent definition, which is the whole point -- the composite score it replaced could
    # not be read off the number. gene2pubmed is computed beside it as a measured alternative
    # (larger for 87-94% of donors, median 2.8x) but does NOT enter selection or the deliverable.
    sp["donor_n_pubs_uniprot"] = sp["lit_pubmed_id"].apply(lambda c: len(split_ids(c)))
    # For the `exact` own-literature rule: the PMID set and the organism binomial per donor.
    sp["donor_pmids"] = sp["lit_pubmed_id"].apply(lambda c: frozenset(split_ids(c)))
    sp["donor_binomial"] = sp["donor_organism"].apply(
        lambda o: " ".join(str(o).split()[:2]) if o else "")
    sp["donor_n_pubs_gene2pubmed"] = (
        sp["geneid"].apply(lambda c: sum(counts.get(g, 0) for g in split_ids(c)))
        if counts else pd.Series(0, index=sp.index))
    # PubTator3, both routes, same status as gene2pubmed: measured, carried, never selected on.
    # `_symbol_anyspecies` is named for its limitation -- it pools every organism with a gene of
    # that name, so it CANNOT be donor-scoped and must never be merged into a scoped count.
    # MAX over the donor's GeneIDs, not SUM: several GeneIDs on one entry are the same gene
    # filed under different loci, so summing double-counts the same papers.
    #
    # THE ZERO/MISSING SPLIT IS THE POINT. A donor WITH a GeneID that PubTator never linked in
    # 36M abstracts is a MEASURED 0 -- it is inside PubTator's universe and absent. A donor with
    # NO GeneID is NaN: unreachable, not unstudied. 86.0% of selected donors have a GeneID;
    # 14.9% of those read 0.
    def _pt_geneid(cell: str) -> float:
        ids = split_ids(cell)
        if not ids:
            return float("nan")
        return float(max(pt_geneid.get(g, 0) for g in ids))

    sp["donor_n_pubs_pubtator_geneid"] = (
        sp["geneid"].apply(_pt_geneid) if pt_geneid else pd.Series(np.nan, index=sp.index))
    sp["donor_n_pubs_pubtator_symbol_anyspecies"] = (
        sp["donor_gene"].apply(lambda g: pt_symbol.get(str(g).strip(), 0) if str(g).strip() else 0)
        if pt_symbol else pd.Series(0, index=sp.index))
    sp["donor_n_pubs"] = sp["donor_n_pubs_uniprot"]
    # TRUE Bacteria, by lineage -- not by organism name and not by "is not human". The ligands
    # axis measured what the loose rule costs: an unrestricted non-human bucket gave 424 apparent
    # potent Kp proteins against a true 175, with a rat protein winning bacterial slots at ~30%
    # identity. `lineage` carries "Bacteria (domain)" explicitly.
    if "Taxonomic lineage" not in sp.columns:
        sys.exit("FATAL the SwissProt metadata has no `Taxonomic lineage` column -- re-run "
                 "scripts/studiedness/fetch.py --only swissprot --refresh to add it.")
    lin = sp["Taxonomic lineage"]
    sp["donor_is_bacteria"] = lin.str.contains("Bacteria (domain)", regex=False, na=False)
    sp["donor_is_eukaryote"] = lin.str.contains("Eukaryota (domain)", regex=False, na=False)
    # PROKARYOTIC = not eukaryotic: Bacteria, Archaea and their PHAGES. Phages are the reason
    # this is not simply "Bacteria" -- measured, 70 of the 93 proteins a strict bacterial rule
    # stranded lost their donor to a virus, and a Kp prophage protein whose best-characterised
    # relative is a lambda or P1 protein is emphatically not novel.
    sp["donor_is_prokaryotic"] = ~sp["donor_is_eukaryote"]
    say(f"  donor pool     {len(sp):,} reviewed entries: "
        f"{int(sp.donor_is_bacteria.sum()):,} Bacteria · "
        f"{int(sp.donor_is_eukaryote.sum()):,} Eukaryota · "
        f"{int((sp.donor_is_prokaryotic & ~sp.donor_is_bacteria).sum()):,} archaea/viruses/other")
    say(f"  swissprot      {len(sp):,} reviewed entries in {time.time() - t0:.0f}s")
    say(f"  curated refs   median {sp.donor_n_pubs_uniprot.median():.0f}   "
        f"P99 {np.percentile(sp.donor_n_pubs_uniprot, 99):.0f}   "
        f"max {sp.donor_n_pubs_uniprot.max():,}  "
        f"({sp.loc[sp.donor_n_pubs_uniprot.idxmax(), 'donor_ac']})   <- NO cap; an earlier note "
        "claiming a ceiling at 58 was wrong")
    say(f"  gene2pubmed    median {sp.donor_n_pubs_gene2pubmed.median():.0f}   "
        f"max {sp.donor_n_pubs_gene2pubmed.max():,}   "
        "(measured alternative; NOT the shipped count)")
    return sp


# The three routes that name THIS protein, imported in spirit from `src/precedents.py` so the two
# axes mean the same thing by "this protein". See exact_own_pmids().
EXACT_PIDENT = 95.0


def pubtator_own_pmids() -> dict[str, set]:
    """PubTator PMID sets per GeneID, for `n_papers_pubtator_own`. Written by pubtator.py."""
    if not PT_OWN_PMIDS.exists():
        say("  NOTE pubtator own-PMID sets absent -- n_papers_pubtator_own will be empty. "
            "Run scripts/studiedness/pubtator.py --route own-pmids")
        return {}
    out: dict[str, set] = {}
    for line in PT_OWN_PMIDS.read_text().splitlines()[1:]:
        g, _, p = line.partition("\t")
        out[g] = set(p.split(";")) if p else set()
    say(f"  pubtator own   PMID sets for {len(out):,} GeneIDs")
    return out


def exact_own_geneids(species: str, hits: pd.DataFrame, sp_meta: pd.DataFrame,
                      prot: pd.DataFrame) -> dict[str, set]:
    """GeneIDs that name THIS protein -- its own, plus the same three `exact` routes.

    The PubTator twin of `exact_own_pmids`. Same rule, different key: literature is reached
    through NCBI GeneIDs rather than UniProt accessions, so the union is taken over the GeneIDs
    of every entry that names this protein.
    """
    gid_of = dict(zip(sp_meta["donor_ac"], sp_meta["geneid"]))
    bino = dict(zip(sp_meta["donor_ac"], sp_meta["donor_binomial"]))
    want = ANCHOR_BINOMIAL[species]
    out: dict[str, set] = {}
    for a, cell in zip(prot["uniprot_ac"], prot["geneid"]):
        got = set(split_ids(cell))
        if got:
            out[a] = got
    h = hits[hits["qseqid"].isin(set(prot["uniprot_ac"]))]
    h = h[(h["qcovhsp"] >= S.COVERAGE_FLOOR) & (h["scovhsp"] >= S.COVERAGE_FLOOR)]
    same = h[(h["pident"] >= EXACT_PIDENT) & h["donor_ac"].map(lambda a: bino.get(a) == want)]
    ident = h[(h["pident"] >= 100) & (h["qcovhsp"] >= 100) & (h["scovhsp"] >= 100)]
    for frame in (same, ident):
        for q, d in zip(frame["qseqid"], frame["donor_ac"]):
            out.setdefault(q, set()).update(split_ids(gid_of.get(d, "")))
    return out


def exact_own_pmids(species: str, hits: pd.DataFrame, sp_meta: pd.DataFrame,
                    self_pmids: dict[str, set]) -> dict[str, set]:
    """PMIDs naming THIS protein, by the ligands axis's `exact` rule.

    **A PROTEIN DOES NOT STOP BEING ITSELF BETWEEN STRAINS.** `src/precedents.py` already settled
    this for ligands: `exact` is the UNION of three routes -- accession, identical sequence, and
    same species at >= 95% identity -- because an accession match and a same-species 99% match
    are the same protein filed twice. Literature must use the same definition or the two axes
    disagree about what "this protein" is.

    **UNION OF PMIDs, NOT MAX.** The ligands axis pools DISTINCT molecules over the matched
    components; the analogue here is distinct PubMed ids, because two strain entries cite
    different papers. Measured: max gains S. aureus 776 proteins, the union gains **965**, and on
    E. coli the union moves the median from 5 to 8 where max does not move it at all.

    What it changes (2026-10-03): E. coli 2,388 proteins gain literature (median 5 -> 8, mean
    6.7 -> 10.7); **S. aureus proteins with ANY literature go 357 -> 1,049**, because NCTC 8325
    is the anchor while most S. aureus curation sits under Newman, USA300, Mu50 and N315. Kp
    gains only 266 -- K. pneumoniae is genuinely thin in SwissProt, which is the axis's premise.

    Identical sequence is detected as DIAMOND pident 100 with full coverage both ways rather than
    by string comparison, which avoids re-reading the 89 MB SwissProt fasta; at 100% identity and
    100% coverage on both sides the sequences are the same.
    """
    pm = dict(zip(sp_meta["donor_ac"], sp_meta["donor_pmids"]))
    want = ANCHOR_BINOMIAL[species]
    bino = dict(zip(sp_meta["donor_ac"], sp_meta["donor_binomial"]))
    h = hits[hits["qseqid"].isin(self_pmids)]
    h = h[(h["qcovhsp"] >= S.COVERAGE_FLOOR) & (h["scovhsp"] >= S.COVERAGE_FLOOR)]
    same = h[(h["pident"] >= EXACT_PIDENT)
             & h["donor_ac"].map(lambda a: bino.get(a) == want)]
    ident = h[(h["pident"] >= 100) & (h["qcovhsp"] >= 100) & (h["scovhsp"] >= 100)]
    out: dict[str, set] = {}
    for frame in (same, ident):
        for q, d in zip(frame["qseqid"], frame["donor_ac"]):
            got = pm.get(d)
            if got:
                out.setdefault(q, set()).update(got)
    return out


def own_exact(species: str, hits: pd.DataFrame, sp_meta: pd.DataFrame,
              accs: set[str]) -> dict[str, set]:
    """Convenience wrapper: build the `exact` PMID map for one species' accessions."""
    return exact_own_pmids(species, hits[hits["qseqid"].isin(accs)], sp_meta,
                           {a: set() for a in accs})


def own_signals(species: str, counts: dict[str, int],
                exact: dict[str, set] | None = None,
                pt_gids: dict[str, set] | None = None,
                pt_pmids: dict[str, set] | None = None) -> pd.DataFrame:
    """This accession's own curation -- a measurement of darkness on Kp and Sa."""
    path = LIT_DIR / f"anchor_{species}.tsv.gz"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    lit = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).rename(
        columns={"Entry": "uniprot_ac", "Reviewed": "reviewed", "Annotation": "annotation_score",
                 "Protein existence": "protein_existence", "PubMed ID": "lit_pubmed_id"})
    prot = P.load_full(species)[["uniprot_ac", "gene_name", "geneid"]]
    df = prot.merge(lit.drop(columns=["GeneID"], errors="ignore"), on="uniprot_ac", how="left")
    if df["annotation_score"].isna().any():
        sys.exit(f"FATAL {species}: {int(df['annotation_score'].isna().sum()):,} proteins missing "
                 "from the UniProt literature stream -- re-run fetch.py --refresh")
    df["annotation_score"] = pd.to_numeric(df["annotation_score"], errors="coerce").fillna(1.0)
    df["n_pubs_uniprot"] = df["lit_pubmed_id"].apply(lambda c: len(split_ids(c)))
    df["n_pubs_gene2pubmed"] = df["geneid"].apply(
        lambda c: sum(counts.get(g, 0) for g in split_ids(c))) if counts else 0
    # THE LIGANDS AXIS'S `exact` RULE, so both axes mean the same thing by "this protein".
    # Union of PMIDs over accession + identical sequence + same species >= 95%.
    df["_self"] = df["lit_pubmed_id"].apply(lambda c: set(split_ids(c)))
    if exact:
        df["n_papers_uniprot_own"] = [len(s0 | exact.get(a, set()))
                              for a, s0 in zip(df["uniprot_ac"], df["_self"])]
        n_up = int((df["n_papers_uniprot_own"] > df["n_pubs_uniprot"]).sum())
        say(f"  exact rule     {n_up:,} of {len(df):,} proteins gain literature from a "
            "same-species or identical-sequence entry")
    else:
        df["n_papers_uniprot_own"] = df["n_pubs_uniprot"]
    df = df.drop(columns=["_self"])
    # The PubTator twin of n_papers_uniprot_own: same `exact` rule, same union-of-PMIDs, but
    # reached through NCBI GeneIDs. **Expect it to be near-empty on Kp and Sa** -- PubTator's
    # gene vocabulary barely covers those proteomes' own GeneIDs (Kp 1.9%, Sa 14.8%, Ec 91.5%).
    # That emptiness is a second, independent measurement of darkness, not a defect, which is
    # the same reason n_papers_uniprot_own is kept despite being flat on Kp.
    if pt_gids is not None and pt_pmids is not None:
        df["n_papers_pubtator_own"] = [
            len(set().union(*(pt_pmids.get(g, set()) for g in pt_gids.get(a, ()))) )
            if pt_gids.get(a) else 0
            for a in df["uniprot_ac"]]
    else:
        df["n_papers_pubtator_own"] = 0
    return df[["uniprot_ac", "gene_name", "reviewed", "annotation_score", "protein_existence",
               "n_pubs_uniprot", "n_pubs_gene2pubmed", "n_papers_uniprot_own",
               "n_papers_pubtator_own"]]


# ---------------------------------------------------------------- DIAMOND

def write_query_fasta(species: tuple[str, ...], path: Path, limit: int | None) -> int:
    n = 0
    with open(path, "w") as fh:
        for sp in species:
            d = P.load(sp)
            if limit:
                d = d.head(limit)
            for ac, seq in zip(d["uniprot_ac"], d["sequence"]):
                fh.write(f">{ac}\n{seq}\n")
                n += 1
    return n


def swissprot_fasta(refresh: bool) -> Path:
    """DIAMOND needs an uncompressed fasta; keep it in scratch, it is re-derivable."""
    out = SCRATCH_DIR / "uniprot_sprot.fasta"
    if out.exists() and out.stat().st_size > 0 and not refresh:
        say(f"  cached {out.relative_to(REPO_ROOT)} ({out.stat().st_size / 1e6:.0f} MB)")
        return out
    if not SPROT_FASTA_GZ.exists():
        sys.exit(f"FATAL missing {SPROT_FASTA_GZ.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    say("  decompressing swissprot fasta ...")
    t0 = time.time()
    with gzip.open(SPROT_FASTA_GZ, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst, length=1 << 22)
    say(f"    {out.stat().st_size / 1e6:.0f} MB in {time.time() - t0:.0f}s")
    return out


def run_diamond(query: Path, db_faa: Path, out_tsv: Path, threads: int, sensitivity: str,
                max_targets: int, refresh: bool) -> pd.DataFrame:
    if out_tsv.exists() and out_tsv.stat().st_size > 0 and not refresh:
        say(f"  cached {out_tsv.relative_to(REPO_ROOT)}")
        return pd.read_csv(out_tsv, sep="\t", names=HIT_COLS)
    dbp = db_faa.with_suffix("")
    if not Path(f"{dbp}.dmnd").exists() or refresh:
        say(f"  diamond makedb over {db_faa.name} ...")
        t0 = time.time()
        subprocess.run(["diamond", "makedb", "--in", str(db_faa), "-d", str(dbp),
                        "--threads", str(threads), "--quiet"], check=True)
        say(f"    built in {time.time() - t0:.0f}s")
    say(f"  diamond blastp ({sensitivity}, max-target-seqs {max_targets}) ...")
    t0 = time.time()
    subprocess.run(
        ["diamond", "blastp", "-q", str(query), "-d", str(dbp), "-o", str(out_tsv),
         "--outfmt", "6", *HIT_COLS, "--max-target-seqs", str(max_targets),
         f"--{sensitivity}", "--threads", str(threads), "--quiet"], check=True)
    say(f"    searched in {time.time() - t0:.0f}s")
    if out_tsv.stat().st_size == 0:
        # Zero hits is a failure for the real queries and the DESIRED result for decoys, so the
        # caller decides; return an empty frame with the right columns either way.
        return pd.DataFrame(columns=HIT_COLS)
    return pd.read_csv(out_tsv, sep="\t", names=HIT_COLS)


def decoy_calibration(species: tuple[str, ...], db: Path, threads: int, sensitivity: str,
                      refresh: bool) -> pd.DataFrame:
    """What the identity floor lets through on sequences that CANNOT have a real homolog.

    Composition-preserving shuffles of our own proteins: same length, same amino-acid composition,
    no evolutionary signal. Anything a decoy matches at a given floor is compositional noise, and
    that is what makes a floor defensible rather than inherited. The pattern is
    `orthology/orthodb.py`'s, applied to this database.
    """
    faa = SCRATCH_DIR / "decoy_query.faa"
    if not faa.exists() or refresh:
        rng = random.Random(0)
        with open(faa, "w") as fh:
            for sp in species:
                for ac, seq in zip(*(P.load(sp)[c] for c in ("uniprot_ac", "sequence"))):
                    chars = list(str(seq))
                    rng.shuffle(chars)
                    fh.write(f">{ac}\n{''.join(chars)}\n")
    hits = run_diamond(faa, db, SCRATCH_DIR / "decoy_hits.tsv", threads, sensitivity, 5, refresh)
    n_q = sum(1 for _ in open(faa) if _.startswith(">"))
    rows = []
    for floor in FLOOR_SWEEP:
        keep = hits[(hits["pident"] >= floor)
                    & (hits["qcovhsp"] >= S.COVERAGE_FLOOR)
                    & (hits["scovhsp"] >= S.COVERAGE_FLOOR)]
        rows.append({"identity_floor": floor, "n_decoys": n_q,
                     "decoys_with_hit": int(keep["qseqid"].nunique()),
                     "false_positive_pct": round(100 * keep["qseqid"].nunique() / n_q, 4)})
    return pd.DataFrame(rows)


def deepen_capped(hits: pd.DataFrame, species: tuple[str, ...], db: Path, max_targets: int,
                  deep: int, threads: int, sensitivity: str, limit: int | None, pre: str,
                  refresh: bool) -> pd.DataFrame:
    """Re-search only the queries that hit the k cap, with a much larger k.

    A capped query could only choose its donor from its CLOSEST relatives, and for a huge family
    like GroEL those 500 are all Firmicutes and Actinobacteria -- the heavily-cited E. coli entry
    never enters the list. Raising k globally is not affordable (13,020 x 5,000 rows), but only
    ~12% of queries are capped, so they are re-searched on their own. Bounded, and it is the
    conserved families -- the ones this axis should be most confident about -- that it repairs.
    """
    per_q = hits.groupby("qseqid").size()
    capped = sorted(per_q[per_q >= max_targets].index)
    if not capped:
        say("  no query hit the cap; nothing to deepen")
        return hits
    say(f"  {len(capped):,} capped queries re-searched at --max-target-seqs {deep:,}")
    faa = SCRATCH_DIR / f"{pre}deep_query.faa"
    want = set(capped)
    n = 0
    with open(faa, "w") as fh:
        for sp in species:
            d = P.load(sp)
            if limit:
                d = d.head(limit)
            for ac, seq in zip(d["uniprot_ac"], d["sequence"]):
                if ac in want:
                    fh.write(f">{ac}\n{seq}\n")
                    n += 1
    deep_hits = run_diamond(faa, db, SCRATCH_DIR / f"{pre}deep_hits.tsv", threads, sensitivity,
                            deep, refresh)
    if deep_hits.empty:
        say("  NOTE the deep pass returned nothing; keeping the capped results")
        return hits
    still = int((deep_hits.groupby("qseqid").size() >= deep).sum())
    if still:
        say(f"  WARNING {still:,} of {n:,} are STILL capped at {deep:,} -- those families are "
            "larger than the deep k and remain biased toward close relatives; raise "
            "--deep-targets")
    else:
        say(f"  none of the {n:,} is still capped at {deep:,}: every query now has its complete "
            "homolog list, so the most-cited donor is the true one")
    return pd.concat([hits[~hits["qseqid"].isin(want)], deep_hits], ignore_index=True)


def parse_hits(hits: pd.DataFrame) -> pd.DataFrame:
    """SwissProt fasta headers are `sp|ACC|NAME`; keep the accession."""
    parts = hits["sseqid"].str.split("|", n=2, expand=True)
    if parts.shape[1] < 2:
        sys.exit(f"FATAL unexpected subject id format: {hits['sseqid'].iloc[0]!r}")
    hits = hits.copy()
    hits["donor_ac"] = parts[1]
    return hits


# ---------------------------------------------------------------- donor choice

def pubtator_donors(hits: pd.DataFrame, sp_meta: pd.DataFrame, identity_floor: float,
                    exclude_taxids: set[str] | None = None,
                    scope: str = "any") -> pd.DataFrame:
    """Best donor per query by PUBTATOR count, and that donor's PubTator count.

    **THE SELECTION RULE MUST MATCH THE VALUE READ, AND GETTING THIS WRONG COSTS HALF THE
    SIGNAL.** Reading PubTator counts off the donor `choose_donors` picked -- the one with most
    CURATED papers -- scores **0.1716** on the held-out control. Selecting on PubTator and
    reading PubTator scores **0.3552**, which also beats the shipped curated column's 0.3280.
    Same hits, same folds, same floor; only the selection rule differs. A mismatched rule looks
    like "this data source is weak" when it is actually "we asked the wrong donor".

    **The count is keyed on NCBI GeneID, so the species is already in the key** -- GeneID 947587
    IS E. coli K-12 `ftsZ`. There is no separate species term to add, and no API call: the counts
    come from `gene2pubtator3.gz` via `pubtator.py`.

    Donors with NO GeneID are not eligible (the count is NaN, not 0): they are unreachable, not
    unstudied. 86.0% of selected donors carry one.
    """
    h = hits.merge(sp_meta, on="donor_ac", how="inner")
    if scope == "bacteria":
        h = h[h["donor_is_bacteria"]]
    elif scope == "prokaryotic":
        h = h[h["donor_is_prokaryotic"]]
    if exclude_taxids:
        h = h[~h["donor_taxid"].isin(exclude_taxids)]
    h = h[(h["pident"] >= identity_floor)
          & (h["qcovhsp"] >= S.COVERAGE_FLOOR)
          & (h["scovhsp"] >= S.COVERAGE_FLOOR)
          & h["donor_n_pubs_pubtator_geneid"].notna()]
    if h.empty:
        return pd.DataFrame(columns=["uniprot_ac"])
    best = (h.sort_values(["qseqid", "donor_n_pubs_pubtator_geneid", "pident"],
                          ascending=[True, False, False], kind="mergesort")
             .drop_duplicates("qseqid", keep="first"))
    return (best[["qseqid", "donor_ac", "donor_gene", "donor_organism", "pident",
                  "donor_n_pubs_pubtator_geneid"]]
            .rename(columns={"qseqid": "uniprot_ac",
                             "donor_ac": "pubtator_donor_ac",
                             "donor_gene": "pubtator_donor_gene",
                             "donor_organism": "pubtator_donor_organism",
                             "pident": "pubtator_donor_pident",
                             "donor_n_pubs_pubtator_geneid": "n_papers_pubtator_prokaryotic"}))


def choose_donors(hits: pd.DataFrame, sp_meta: pd.DataFrame, identity_floor: float,
                  exclude_taxids: set[str] | None = None,
                  scope: str = "any") -> pd.DataFrame:
    """Best donor per query: most-cited above the floor, then annotation score, then identity.

    Also records the NEAREST hit (highest identity) so the choice is auditable -- a most-cited
    donor at 41% identity and a nearest donor at 98% are very different situations and the table
    must show both.
    """
    h = hits.merge(sp_meta, on="donor_ac", how="inner")
    if scope == "bacteria":
        h = h[h["donor_is_bacteria"]]
    elif scope == "prokaryotic":
        h = h[h["donor_is_prokaryotic"]]
    if exclude_taxids:
        h = h[~h["donor_taxid"].isin(exclude_taxids)]
    h = h[(h["pident"] >= identity_floor)
          & (h["qcovhsp"] >= S.COVERAGE_FLOOR)
          & (h["scovhsp"] >= S.COVERAGE_FLOOR)]
    if h.empty:
        return pd.DataFrame(columns=["uniprot_ac"])

    n_cand = h.groupby("qseqid").size().rename("n_candidates")
    # MOST-CITED above the floor, not the closest. Preferring the closest band first was tried
    # and is measurably WORSE on the held-out control (0.4426 against 0.4666 at the 40% floor):
    # the quantity being transferred is literature volume, and the best-cited homolog predicts
    # "this family is studied" better than the nearest one does.
    best = (h.sort_values(["qseqid", "donor_n_pubs", "donor_annotation_score", "pident"],
                          ascending=[True, False, False, False], kind="mergesort")
             .drop_duplicates("qseqid", keep="first"))
    near = (h.sort_values(["qseqid", "pident", "bitscore"],
                          ascending=[True, False, False], kind="mergesort")
             .drop_duplicates("qseqid", keep="first")
             [["qseqid", "donor_ac", "pident", "donor_organism", "donor_n_pubs"]]
             .rename(columns={"donor_ac": "nearest_ac", "pident": "nearest_pident",
                              "donor_organism": "nearest_organism",
                              "donor_n_pubs": "nearest_n_pubs"}))

    out = (best[["qseqid", "donor_ac", "donor_gene", "donor_organism", "donor_taxid",
                 "pident", "qcovhsp", "donor_n_pubs", "donor_n_pubs_uniprot",
                 "donor_n_pubs_gene2pubmed", "donor_n_pubs_pubtator_geneid",
                 "donor_n_pubs_pubtator_symbol_anyspecies",
                 "donor_annotation_score", "donor_protein_name"]]
           .rename(columns={"qseqid": "uniprot_ac", "pident": "donor_pident",
                            "qcovhsp": "donor_qcov"})
           .merge(near.rename(columns={"qseqid": "uniprot_ac"}), on="uniprot_ac", how="left")
           .merge(n_cand.rename_axis("uniprot_ac").reset_index(), on="uniprot_ac", how="left"))
    return out


def assemble(species: str, own: pd.DataFrame, donors: pd.DataFrame,
             any_hit: set[str] | None = None) -> pd.DataFrame:
    """One row per protein, every protein present, evidence naming the tier.

    `any_hit` splits the zero-score rows into two genuinely different situations: a protein with a
    SwissProt hit that fell below the floor still has a distant relative, while one with no hit at
    all is the strongest novelty claim the axis can make. Both score 0 -- the axis will not invent
    a number for them -- but they must be distinguishable, or a third of K. pneumoniae becomes one
    undifferentiated tie at the bottom of the ranking.
    """
    df = own.merge(donors, on="uniprot_ac", how="left")
    # Straight through: the donor's curated reference count IS the number. No scaling, no blend.
    df["n_papers_uniprot_prokaryotic"] = df["donor_n_pubs"].fillna(0).astype(int)
    below = (df["uniprot_ac"].isin(any_hit) if any_hit is not None
             else pd.Series(False, index=df.index))
    df["evidence"] = np.select(
        [df["donor_ac"].notna() & (df["donor_pident"] >= S.IDENTITY_DIRECT),
         df["donor_ac"].notna() & (df["donor_pident"] >= S.IDENTITY_CLOSE),
         df["donor_ac"].notna(),
         below],
        ["swissprot_direct", "swissprot_close", "swissprot_homolog", "below_floor"],
        default="no_hit")
    return df


# ---------------------------------------------------------------- the control

def ecoli_taxids(sp_meta: pd.DataFrame) -> set[str]:
    """Every Escherichia taxid in SwissProt -- strains are filed under their own ids.

    The same trap the ligands axis records: matching on one tax_id finds 65 E. coli targets while
    matching on the organism name finds 225, because K-12 lives at 83333 and not 562.
    """
    mask = sp_meta["donor_organism"].str.startswith("Escherichia", na=False)
    return set(sp_meta.loc[mask, "donor_taxid"])


def run_control(hits: pd.DataFrame, sp_meta: pd.DataFrame, own_ec: pd.DataFrame,
                scope: str = "any") -> tuple[pd.DataFrame, dict]:
    """E. coli family score with every Escherichia donor removed, vs its own measured score."""
    excl = ecoli_taxids(sp_meta)
    say(f"    excluding {len(excl):,} Escherichia taxids "
        f"({int(sp_meta['donor_taxid'].isin(excl).sum()):,} SwissProt entries) from the donor pool")
    ec_hits = hits[hits["qseqid"].isin(set(own_ec["uniprot_ac"]))]
    donors = choose_donors(ec_hits, sp_meta, S.IDENTITY_FLOOR, exclude_taxids=excl, scope=scope)
    in_scope = (ec_hits if scope == "any" else ec_hits[ec_hits["donor_ac"].isin(
                    set(sp_meta.loc[sp_meta[f"donor_is_{scope}"], "donor_ac"]))])
    df = assemble("ecoli", own_ec, donors, any_hit=set(in_scope["qseqid"]))
    # Only proteins the held-out transfer could actually score are informative about the
    # mechanism; a no_homolog row says nothing about how good the transfer is when it fires.
    scored = df[df["donor_ac"].notna()]
    rho = scored["n_papers_uniprot_prokaryotic"].corr(scored["n_papers_uniprot_own"], method="spearman")
    # Pearson on raw counts is meaningless here -- both are heavily skewed integers -- so it is
    # taken on log1p, which is the scale on which "twice as studied" is a constant step.
    pearson = np.log1p(scored["n_papers_uniprot_prokaryotic"]).corr(np.log1p(scored["n_papers_uniprot_own"]))
    stats = {"n": len(df), "n_scored": len(scored),
             "donor_scope": scope,
             "coverage_pct": round(100 * len(scored) / len(df), 1),
             "spearman": round(float(rho), 4), "pearson": round(float(pearson), 4),
             "floor": CONTROL_RHO_FLOOR, "excluded_taxids": len(excl)}

    # THE ALTERNATIVE COUNTS, SCORED ON THE IDENTICAL HELD-OUT FOLDS. Each is the same transfer
    # mechanism carrying a different literature definition, so the comparison is like-for-like and
    # the only honest basis for promoting one. `n` is reported per route because the GeneID-keyed
    # routes reach fewer donors -- only 42.1% of prokaryotic SwissProt entries carry a GeneID at
    # all -- and a higher rho over a smaller, better-curated subset is not a better axis.
    # The SELF-CONSISTENT PubTator transfer: donor selected on PubTator, PubTator read.
    # Reported separately from the mismatched read below, because the gap between them (0.3552
    # vs 0.1716) is the measurement that justifies having a second donor at all.
    ptd = pubtator_donors(ec_hits, sp_meta, S.IDENTITY_FLOOR, exclude_taxids=excl, scope=scope)
    if len(ptd):
        p = scored[["uniprot_ac", "n_papers_uniprot_own"]].merge(ptd, on="uniprot_ac", how="inner")
        if len(p) >= 50:
            stats["rho_pubtator_own_donor"] = round(float(
                p["n_papers_pubtator_prokaryotic"].corr(p["n_papers_uniprot_own"], method="spearman")), 4)
            stats["n_pubtator_own_donor"] = len(p)

    alt_cols = ["donor_n_pubs_gene2pubmed", "donor_n_pubs_pubtator_geneid",
                "donor_n_pubs_pubtator_symbol_anyspecies"]
    have = [c for c in alt_cols if c in sp_meta.columns]
    if have:
        # `scored` ALREADY carries donor_n_pubs_gene2pubmed (assemble keeps it), so a plain merge
        # suffixes both copies to _x/_y and the bare name vanishes -- a KeyError, which is the
        # good case; a silent suffix on only some columns would have scored the wrong series.
        # Take the donor-side values as authoritative and drop the left copies first.
        alt = (scored.drop(columns=[c for c in have if c in scored.columns])
                     .merge(sp_meta[["donor_ac", *have]], on="donor_ac", how="left"))
        for c in have:
            v = pd.to_numeric(alt[c], errors="coerce")
            # NOT `v > 0`. A measured zero is a data point -- for the GeneID route it means
            # PubTator never linked this gene across 36M abstracts, which is exactly the
            # "nobody writes about it" signal the axis wants. Dropping zeros scored only the
            # already-studied donors and understated the route by a wide margin.
            ok = v.notna()
            short = c.replace("donor_n_pubs_", "")
            stats[f"rho_{short}"] = (round(float(v[ok].corr(alt.loc[ok, "n_papers_uniprot_own"],
                                                            method="spearman")), 4)
                                     if ok.sum() >= 50 else None)
            stats[f"n_{short}"] = int(ok.sum())
    keep = ["uniprot_ac", "gene_name", "n_papers_uniprot_own", "n_papers_uniprot_prokaryotic", "evidence",
            "donor_ac", "donor_organism", "donor_pident", "donor_n_pubs"]
    return df[keep], stats


def free_route_coverage(species: str) -> dict:
    """What the on-disk four-species ortholog table reaches. Measured, never merged."""
    d = O.load_dense(species)
    cols = [c for c in d.columns
            if c.startswith("n_orthologs_") and ("ecoli" in c or "human" in c)
            and "_of_" not in c and "_rbh_" not in c]
    if not cols:
        return {"species": species, "route": "orthologs_4sp", "reachable_pct": float("nan")}
    tot = sum(pd.to_numeric(d[c], errors="coerce").fillna(0) for c in cols)
    return {"species": species, "route": "orthologs_4sp", "n": len(d),
            "reachable": int((tot > 0).sum()),
            "reachable_pct": round(100 * float((tot > 0).mean()), 1)}


# ---------------------------------------------------------------- main

def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--sensitivity", default="very-sensitive",
                    choices=["fast", "sensitive", "more-sensitive", "very-sensitive", "ultra-sensitive"])
    ap.add_argument("--donor-scope", default="prokaryotic", choices=list(SCOPES),
                    help="which SwissProt entries may donate literature. ALL THREE are computed "
                         "every run and compared into evidence/donor_scope_comparison.tsv; this "
                         "picks which one becomes the shipped n_papers_uniprot_prokaryotic")
    ap.add_argument("--deep-targets", type=int, default=5000,
                    help="k for the second pass over queries that hit --max-targets (0 disables)")
    ap.add_argument("--max-targets", type=int, default=500,
                    help="DIAMOND hits per query. MUST be large: DIAMOND ranks by bitscore, so a "
                         "small cap returns only close relatives and the well-cited model-organism "
                         "homolog never enters the candidate list (see the docstring)")
    ap.add_argument("--limit", type=int, default=None,
                    help="smoke test: first N proteins per species, writes only scratch/")
    ap.add_argument("--no-control", action="store_true")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    smoke = args.limit is not None
    pre = "smoke_" if smoke else ""
    dest = SCRATCH_DIR if smoke else EVIDENCE_DIR

    rule("=")
    say("STUDIEDNESS - transfer from the best-characterised SwissProt homolog")
    rule("=")
    say("  in    data/source/uniprot/literature/{swissprot_meta.tsv.gz,uniprot_sprot.fasta.gz}")
    say(f"        {G2P_COUNTS.relative_to(REPO_ROOT)}")
    say(f"  out   {dest.relative_to(REPO_ROOT)}/{pre}{{own,transfer}}_<species>.tsv")
    say("  the number is a PAPER COUNT -- curated SwissProt references, no scaling, no blend")
    say(f"  donor = MOST-CITED hit at >= {S.IDENTITY_FLOOR:.0f}% id / "
        f">= {S.COVERAGE_FLOOR:.0f}% cov  (>= {S.IDENTITY_DIRECT:.0f}% = direct)")
    say("  the anchors are dark: 5,710/5,728 Kp proteins have exactly 1 PubMed id")
    say(f"  control: E. coli with every Escherichia donor removed, rho floor "
        f"{CONTROL_RHO_FLOOR}")
    if smoke:
        say(f"  SMOKE TEST --limit {args.limit}: writing to scratch/ with a `{pre}` prefix")
    rule("=")
    if args.dry_run:
        say("dry run: nothing written.")
        return

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    say(f"  diamond   : {ensure_diamond()}")

    rule()
    say("INPUTS")
    rule()
    counts = gene2pubmed_counts()
    pt_geneid, pt_symbol = pubtator_counts()
    pt_pmids = pubtator_own_pmids()
    sp_meta = load_swissprot(counts, pt_geneid, pt_symbol)
    top = sp_meta.loc[sp_meta["donor_n_pubs_uniprot"].idxmax()]
    say(f"  THE NUMBER  n_papers_uniprot_prokaryotic = the donor's curated reference count. No scaling, no "
        "blend.")
    say(f"              most-curated entry in SwissProt: {top['donor_ac']} with "
        f"{int(top['donor_n_pubs_uniprot'])} refs -- there is no cap")

    rule()
    say("DIAMOND vs SwissProt")
    rule()
    query = SCRATCH_DIR / f"{pre}query.faa"
    n_q = write_query_fasta(tuple(args.species), query, args.limit)
    say(f"  {n_q:,} query sequences -> {query.relative_to(REPO_ROOT)}")
    db = swissprot_fasta(args.refresh)
    raw = run_diamond(query, db, SCRATCH_DIR / f"{pre}hits.tsv", args.threads,
                      args.sensitivity, args.max_targets, args.refresh)
    if raw.empty:
        sys.exit("FATAL diamond returned zero hits against 575,748 reviewed entries. That is not "
                 "a biological result -- check the database and the query fasta.")
    hits = parse_hits(raw)
    say(f"  {len(hits):,} hits over {hits.qseqid.nunique():,} of {n_q:,} queries "
        f"({100 * hits.qseqid.nunique() / n_q:.1f}% have any hit at all)")
    per_q = raw.groupby("qseqid").size()
    capped = int((per_q >= args.max_targets).sum())
    say(f"  {capped:,} queries returned the full {args.max_targets} hits "
        f"({100 * capped / max(len(per_q), 1):.1f}% capped) -- a capped query can only choose a "
        "donor from its CLOSEST relatives, which biases conserved families downward")
    if args.deep_targets and capped:
        raw = deepen_capped(raw, tuple(args.species), db, args.max_targets, args.deep_targets,
                            args.threads, args.sensitivity, args.limit, pre, args.refresh)
    hits = parse_hits(raw)
    say(f"  {len(hits):,} hits after deepening over {hits.qseqid.nunique():,} queries")

    rule()
    say("PER-SPECIES")
    rule()
    scope_acs = {sc: set(sp_meta.loc[sp_meta[f"donor_is_{sc}"], "donor_ac"])
                 for sc in ("bacteria", "prokaryotic")}
    rows, route_rows, scope_rows = [], [], []
    for sp in args.species:
        own = own_signals(sp, counts,
                          own_exact(sp, hits, sp_meta,
                                    set(P.load_full(sp)["uniprot_ac"])),
                          exact_own_geneids(sp, hits, sp_meta, P.load_full(sp)), pt_pmids)
        if args.limit:
            own = own.head(args.limit)
        sp_hits = hits[hits.qseqid.isin(set(own["uniprot_ac"]))]
        # BOTH scopes, always. The shipped one is a choice; the other is the measurement that
        # justifies it, and a reader who disagrees can swap columns without re-running DIAMOND.
        scoped = {}
        for scope in SCOPES:
            d_s = choose_donors(sp_hits, sp_meta, S.IDENTITY_FLOOR, scope=scope)
            # `any_hit` is SCOPE-SPECIFIC on purpose: under the bacterial scope, `below_floor`
            # must mean "has a bacterial hit below the floor" and `no_hit` "no bacterial hit at
            # all". Using the unrestricted hit set would file a protein whose only curated
            # relative is eukaryotic as though it had a distant bacterial one, which is a
            # different and much weaker novelty claim.
            in_scope = (sp_hits if scope == "any"
                        else sp_hits[sp_hits["donor_ac"].isin(scope_acs[scope])])
            scoped[scope] = assemble(sp, own, d_s, any_hit=set(in_scope["qseqid"]))
        df = scoped[args.donor_scope].copy()
        df["donor_scope"] = args.donor_scope
        # The PubTator column with its OWN donor, selected on PubTator. See pubtator_donors().
        ptd = pubtator_donors(sp_hits, sp_meta, S.IDENTITY_FLOOR, scope=args.donor_scope)
        df = df.merge(ptd, on="uniprot_ac", how="left")
        for scope in SCOPES:
            df[f"n_papers_uniprot_{scope}"] = scoped[scope]["n_papers_uniprot_prokaryotic"].to_numpy()
        other = "any"
        df[f"donor_ac_{other}"] = scoped[other]["donor_ac"].to_numpy()
        df[f"donor_organism_{other}"] = scoped[other]["donor_organism"].to_numpy()
        base = scoped["any"]
        for scope in SCOPES:
            cur = scoped[scope]
            scope_rows.append({
                "species": sp, "scope": scope, "n": len(df),
                "with_donor": int(cur["donor_ac"].notna().sum()),
                "lost_vs_any": int((base["donor_ac"].notna().to_numpy()
                                    & cur["donor_ac"].isna().to_numpy()).sum()),
                "different_donor_vs_any": int(
                    (cur["donor_ac"].fillna("") != base["donor_ac"].fillna("")).sum()),
                "median_family": int(cur["n_papers_uniprot_prokaryotic"].median()),
                "spearman_vs_any": round(float(cur["n_papers_uniprot_prokaryotic"].corr(
                    base["n_papers_uniprot_prokaryotic"], method="spearman")), 4),
                "shipped": scope == args.donor_scope,
            })

        own_path = dest / f"{pre}own_{sp}.tsv"
        tr_path = dest / f"{pre}transfer_{sp}.tsv"
        own.to_csv(own_path, sep="\t", index=False)
        tr_cols = ["uniprot_ac", "donor_ac", "donor_gene", "donor_organism", "donor_taxid",
                   "donor_pident", "donor_qcov", "donor_n_pubs", "donor_n_pubs_uniprot",
                   "donor_n_pubs_gene2pubmed", "donor_n_pubs_pubtator_geneid",
                   "donor_n_pubs_pubtator_symbol_anyspecies",
                   "n_papers_pubtator_prokaryotic", "pubtator_donor_ac", "pubtator_donor_gene",
                   "pubtator_donor_organism", "pubtator_donor_pident",
                   "donor_annotation_score", "donor_protein_name",
                   "nearest_ac", "nearest_pident", "nearest_organism", "nearest_n_pubs",
                   "n_candidates", "n_papers_uniprot_own", "n_papers_uniprot_prokaryotic",
                   "evidence", "donor_scope",
                   # NOT n_papers_uniprot_prokaryotic again -- the shipped column IS the
                   # prokaryotic scope (verified byte-identical), so listing it twice would
                   # emit a duplicate header. The other two scopes are the real alternatives.
                   "n_papers_uniprot_bacteria", "n_papers_uniprot_any",
                   f"donor_ac_{other}", f"donor_organism_{other}"]
        df[tr_cols].to_csv(tr_path, sep="\t", index=False)

        tiers = df["evidence"].value_counts().to_dict()
        covered = int(df["donor_ac"].notna().sum())
        # Unrankable ties are this axis's version of stage 04's saturated-probability failure.
        # A count has no ceiling to saturate against, but it can still collapse: if one value
        # held most of a proteome the ranking would be useless.
        # ZERO IS EXCLUDED from the tie check. A third of K. pneumoniae genuinely has no
        # in-scope homolog, and that group is already split into two documented tiers
        # (`no_hit` vs `below_floor`) -- it is an answer, not an unrankable accident. What the
        # guard must catch is the SCORED part of the ranking collapsing.
        scored_fam = df.loc[df["n_papers_uniprot_prokaryotic"] > 0, "n_papers_uniprot_prokaryotic"]
        tie_counts = (scored_fam.value_counts() if len(scored_fam)
                      else pd.Series({0: 0}, dtype=int))
        frac = float(tie_counts.iloc[0] / max(len(scored_fam), 1))
        if frac > S.MAX_TIE_FRACTION:
            sys.exit(
                f"FATAL {sp}: {tie_counts.iloc[0]:,} of {len(scored_fam):,} SCORED proteins "
                f"({100 * frac:.1f}%) all have {int(tie_counts.index[0])} papers, above the "
                f"{100 * S.MAX_TIE_FRACTION:.0f}% guard -- the ranking has collapsed into one "
                "tie and cannot be sorted.")
        say(f"  {sp}")
        say(f"    {len(df):,} proteins   with a donor {covered:,} "
            f"({100 * covered / len(df):.1f}%)   "
            + "   ".join(f"{k} {v:,}" for k, v in sorted(tiers.items())))
        fam = df["n_papers_uniprot_prokaryotic"]
        say(f"    n_papers_uniprot_own    median {df.n_papers_uniprot_own.median():.0f}   "
            f"max {df.n_papers_uniprot_own.max():,}   distinct {df.n_papers_uniprot_own.nunique():,}")
        say(f"    n_papers_uniprot_prokaryotic median {fam.median():.0f}   max {fam.max():,}   "
            f"distinct {fam.nunique():,}   zero for {int((fam == 0).sum()):,}")
        top3 = fam.value_counts().sort_index(ascending=False).head(3)
        say("      top of the ranking: "
            + "  ".join(f"{int(k)}p x{int(v)}" for k, v in top3.items())
            + f"   largest non-zero tie {int(tie_counts.iloc[0]):,} at "
              f"{int(tie_counts.index[0])}p ({100 * frac:.1f}% of scored)")
        top_org = df.loc[df.donor_organism.notna(), "donor_organism"].value_counts().head(3)
        for org, n in top_org.items():
            say(f"      donor organism  {n:5,}  {org}")
        say(f"    -> {tr_path.relative_to(REPO_ROOT)}")
        rows.append({"species": sp, "n": len(df), "with_donor": covered,
                     "with_donor_pct": round(100 * covered / len(df), 1),
                     **{f"tier_{k}": v for k, v in tiers.items()},
                     "median_own": int(df.n_papers_uniprot_own.median()),
                     "median_family": int(fam.median()),
                     "max_family": int(fam.max()),
                     "distinct_family": int(fam.nunique()),
                     "largest_tie_pct": round(100 * frac, 1)})
        route_rows.append({"species": sp, "route": "swissprot_diamond", "n": len(df),
                           "reachable": covered,
                           "reachable_pct": round(100 * covered / len(df), 1)})
        if not smoke:
            route_rows.append(free_route_coverage(sp))

    if smoke:
        rule("=")
        say(f"smoke test complete; {len(rows)} species written to "
            f"{SCRATCH_DIR.relative_to(REPO_ROOT)}. No deliverable touched.")
        rule("=")
        return

    pd.DataFrame([{"quantity": "curated PubMed references on the best-studied in-scope "
                               "SwissProt homolog",
                   "count_source": "uniprot_lit_pubmed_id",
                   "not_used": "ncbi_gene2pubmed (measured alternative, see gene2pubmed_gain.tsv)",
                   "donor_scope": args.donor_scope,
                   "identity_floor": S.IDENTITY_FLOOR, "identity_close": S.IDENTITY_CLOSE,
                   "identity_direct": S.IDENTITY_DIRECT, "coverage_floor": S.COVERAGE_FLOOR,
                   "swissprot_entries": len(sp_meta),
                   "max_curated_refs_in_swissprot": int(sp_meta.donor_n_pubs_uniprot.max()),
                   "scaled_reference_if_needed": S.SCALE_REFERENCE}]).to_csv(
        EVIDENCE_DIR / "definition.tsv", sep="\t", index=False)
    pd.DataFrame(route_rows).to_csv(EVIDENCE_DIR / "route_comparison.tsv", sep="\t", index=False)
    pd.DataFrame(scope_rows).to_csv(EVIDENCE_DIR / "donor_scope_comparison.tsv",
                                    sep="\t", index=False)
    say("")
    say(f"  DONOR SCOPE (shipping `{args.donor_scope}`; all three computed, all three in the "
        "transfer table)")
    for r in scope_rows:
        mark = "  <-- shipped" if r["shipped"] else ""
        say(f"    {r['species']:<12} {r['scope']:<12} donor {r['with_donor']:5,}   "
            f"lose-vs-any {r['lost_vs_any']:4,}   different {r['different_donor_vs_any']:4,}   "
            f"median {r['median_family']:.3f}   rho-vs-any {r['spearman_vs_any']:.3f}{mark}")
    say("")
    say("  ROUTE COMPARISON (nothing is merged; the free route is measured, not used)")
    for r in route_rows:
        say(f"    {r['species']:<12} {r['route']:<18} {r.get('reachable', 0):5,} "
            f"({r['reachable_pct']}%)")

    rule()
    say("DECOY CALIBRATION - what each floor lets through on shuffled sequences")
    rule()
    decoys = decoy_calibration(tuple(args.species), db, args.threads, args.sensitivity,
                               args.refresh)
    decoys.to_csv(EVIDENCE_DIR / "decoy_calibration.tsv", sep="\t", index=False)
    for r in decoys.itertuples(index=False):
        say(f"  floor {r.identity_floor:5.0f}%   {r.decoys_with_hit:5,} of {r.n_decoys:,} "
            f"composition-preserving shuffles match   ({r.false_positive_pct}%)")
    say(f"  -> {(EVIDENCE_DIR / 'decoy_calibration.tsv').relative_to(REPO_ROOT)}")

    rule()
    say("FLOOR SENSITIVITY - coverage AND the held-out control, at every floor")
    rule()
    say("  a looser floor is only worth taking if the CONTROL does not degrade; coverage alone")
    say("  would argue for no floor at all")
    excl = ecoli_taxids(sp_meta)
    own_by_sp = {sp: own_signals(sp, counts,
                                 own_exact(sp, hits, sp_meta,
                                           set(P.load_full(sp)["uniprot_ac"])),
                                 exact_own_geneids(sp, hits, sp_meta, P.load_full(sp)), pt_pmids)
                 for sp in args.species}
    sweep = []
    for floor in FLOOR_SWEEP:
        row = {"identity_floor": floor}
        for sp in args.species:
            own = own_by_sp[sp]
            d = choose_donors(hits[hits.qseqid.isin(set(own["uniprot_ac"]))], sp_meta, floor,
                              scope=args.donor_scope)
            cov = 0 if d.empty else int(d["uniprot_ac"].nunique())
            row[f"{sp}_pct"] = round(100 * cov / len(own), 1)
        if "ecoli" in args.species:
            ec = own_by_sp["ecoli"]
            dc = choose_donors(hits[hits.qseqid.isin(set(ec["uniprot_ac"]))], sp_meta, floor,
                               exclude_taxids=excl, scope=args.donor_scope)
            cdf = assemble("ecoli", ec, dc)  # bands unused here; rho only
            sc = cdf[cdf["donor_ac"].notna()]
            row["control_spearman"] = (round(float(
                sc["n_papers_uniprot_prokaryotic"].corr(sc["n_papers_uniprot_own"], method="spearman")), 4)
                if len(sc) > 100 else float("nan"))
            row["control_n"] = len(sc)
        sweep.append(row)
        say(f"  floor {floor:5.0f}%   "
            + "   ".join(f"{sp} {row[f'{sp}_pct']:5.1f}%" for sp in args.species)
            + f"   control rho {row.get('control_spearman', float('nan')):.4f} "
              f"(n {row.get('control_n', 0):,})")
    pd.DataFrame(sweep).to_csv(EVIDENCE_DIR / "floor_sensitivity.tsv", sep="\t", index=False)
    say(f"  -> {(EVIDENCE_DIR / 'floor_sensitivity.tsv').relative_to(REPO_ROOT)}")
    say(f"  the shipped floor is {S.IDENTITY_FLOOR:.0f}% (src/studiedness.IDENTITY_FLOOR)")

    if not args.no_control and "ecoli" in args.species:
        rule()
        say("CONTROL - E. coli, every Escherichia donor struck out of SwissProt")
        rule()
        ec_own = own_signals("ecoli", counts,
                             own_exact("ecoli", hits, sp_meta,
                                       set(P.load_full("ecoli")["uniprot_ac"])),
                             exact_own_geneids("ecoli", hits, sp_meta,
                                               P.load_full("ecoli")), pt_pmids)
        both = {}
        for scope in SCOPES:
            both[scope] = run_control(hits, sp_meta, ec_own, scope=scope)
        # THE ARBITER between donor scopes: which one better recovers E. coli's own measured
        # literature from held-out donors? Reported for both, whichever is shipped.
        for scope, (_, st) in both.items():
            mark = "  <-- shipped" if scope == args.donor_scope else ""
            say(f"    {scope:<8} {st['n_scored']:,} of {st['n']:,} scored "
                f"({st['coverage_pct']}%)   spearman {st['spearman']}   "
                f"pearson {st['pearson']}{mark}")
        ctrl, stats = both[args.donor_scope]
        ctrl.to_csv(EVIDENCE_DIR / "control_ecoli_heldout.tsv", sep="\t", index=False)
        pd.DataFrame([st for _, st in both.values()]).to_csv(
            EVIDENCE_DIR / "control_ecoli_summary.tsv", sep="\t", index=False)
        say(f"    floor {CONTROL_RHO_FLOOR}")
        if stats["spearman"] < CONTROL_RHO_FLOOR:
            sys.exit(f"FATAL held-out transfer control failed: spearman {stats['spearman']} "
                     f"< {CONTROL_RHO_FLOOR}. The transfer does not recover E. coli's own "
                     "measured literature, so the family score cannot be trusted on Kp or Sa.")
        say("    control PASSED -- transfer recovers measured literature on held-out E. coli")

    pd.DataFrame(rows).to_csv(EVIDENCE_DIR / "transfer_manifest.tsv", sep="\t", index=False)
    rule("=")
    say(f"  wrote {(EVIDENCE_DIR / 'transfer_manifest.tsv').relative_to(REPO_ROOT)}")
    say("  next: scripts/studiedness/merge.py builds the deliverable")
    rule("=")


if __name__ == "__main__":
    main()
