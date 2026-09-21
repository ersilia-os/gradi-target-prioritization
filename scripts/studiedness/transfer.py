"""Studiedness by sequence: the best-characterised SwissProt homolog of every protein.

The anchors' own literature cannot carry this axis. Measured on UniProt 2026_03: **5,710 of 5,728
K. pneumoniae proteins carry exactly one PubMed id** (the genome paper, 7 distinct values over the
whole proteome) and 2,532 of 2,889 S. aureus proteins carry none. So the usable number is
transferred from homologs -- which is also the house rule, *"map by sequence, not by accession,
when reaching an external database"*.

    data/processed/studiedness/evidence/own_<species>.tsv        this accession's own signals
    data/processed/studiedness/evidence/transfer_<species>.tsv   the chosen donor + its counts
    data/processed/studiedness/evidence/scale.tsv                the fixed global scale anchor
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
`studiedness_family_bacteria` and `studiedness_family_any` side by side in the transfer table and
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
mechanism can be tested at all. The control recomputes `studiedness_family` for E. coli with
**every Escherichia donor struck out of SwissProt**, then correlates it against E. coli's own
measured `studiedness_own`. That is a genuine held-out test of the whole mechanism -- the same
shape as `embeddings/prott5.py` validating against UniProt's published vectors. The run exits
non-zero if it falls below the floor.

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
DEFAULT_DIAMOND_DIR = Path.home() / "miniconda3" / "envs" / "gradi-ortho" / "bin"

SPECIES = ("kpneumoniae", "ecoli", "saureus")
HIT_COLS = ["qseqid", "sseqid", "pident", "qcovhsp", "scovhsp", "bitscore", "evalue"]
FLOOR_SWEEP = (20.0, 25.0, 30.0, 40.0, 60.0, 95.0)
# Donor scopes, all computed every run. `prokaryotic` is "not eukaryotic" -- Bacteria, Archaea
# AND phages; see the docstring for why phages must be in.
SCOPES = ("prokaryotic", "bacteria", "any")

# The control floor. Set from the first full run and then frozen -- a floor tuned to each run
# is not a control. See the run log in docs/studiedness.md.
CONTROL_RHO_FLOOR = 0.45

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


def load_swissprot(counts: dict[str, int]) -> pd.DataFrame:
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
    n_uniprot = sp["lit_pubmed_id"].apply(lambda c: len(split_ids(c)))
    if counts:
        n_ncbi = sp["geneid"].apply(lambda c: sum(counts.get(g, 0) for g in split_ids(c)))
    else:
        n_ncbi = pd.Series(0, index=sp.index)
    # The union is bounded below by the max and above by the sum; the two sources overlap heavily
    # (gene2pubmed re-lists UniProt's curated papers), so summing would double-count. Max is the
    # honest choice, and both components ship so it can be undone.
    sp["donor_n_pubs_uniprot"] = n_uniprot
    sp["donor_n_pubs_gene2pubmed"] = n_ncbi
    sp["donor_n_pubs"] = np.maximum(n_uniprot, n_ncbi)
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
    say(f"  swissprot      {len(sp):,} reviewed entries in {time.time() - t0:.0f}s   "
        f"median n_pubs {sp.donor_n_pubs.median():.0f}   "
        f"P95 {np.percentile(sp.donor_n_pubs, 95):.0f}   max {sp.donor_n_pubs.max():,}")
    return sp


def own_signals(species: str, counts: dict[str, int]) -> pd.DataFrame:
    """This accession's own curation -- a measurement of darkness on Kp and Sa."""
    path = LIT_DIR / f"anchor_{species}.tsv.gz"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    lit = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).rename(
        columns={"Entry": "uniprot_ac", "Reviewed": "reviewed", "Annotation": "annotation_score",
                 "Protein existence": "protein_existence", "PubMed ID": "lit_pubmed_id"})
    prot = P.load(species)[["uniprot_ac", "gene_name", "geneid"]]
    df = prot.merge(lit.drop(columns=["GeneID"], errors="ignore"), on="uniprot_ac", how="left")
    if df["annotation_score"].isna().any():
        sys.exit(f"FATAL {species}: {int(df['annotation_score'].isna().sum()):,} proteins missing "
                 "from the UniProt literature stream -- re-run fetch.py --refresh")
    df["annotation_score"] = pd.to_numeric(df["annotation_score"], errors="coerce").fillna(1.0)
    df["n_pubs_uniprot"] = df["lit_pubmed_id"].apply(lambda c: len(split_ids(c)))
    df["n_pubs_gene2pubmed"] = df["geneid"].apply(
        lambda c: sum(counts.get(g, 0) for g in split_ids(c))) if counts else 0
    df["n_pubs"] = df[["n_pubs_uniprot", "n_pubs_gene2pubmed"]].max(axis=1)
    return df[["uniprot_ac", "gene_name", "reviewed", "annotation_score", "protein_existence",
               "n_pubs_uniprot", "n_pubs_gene2pubmed", "n_pubs"]]


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
                 "donor_n_pubs_gene2pubmed", "donor_annotation_score", "donor_protein_name"]]
           .rename(columns={"qseqid": "uniprot_ac", "pident": "donor_pident",
                            "qcovhsp": "donor_qcov"})
           .merge(near.rename(columns={"qseqid": "uniprot_ac"}), on="uniprot_ac", how="left")
           .merge(n_cand.rename_axis("uniprot_ac").reset_index(), on="uniprot_ac", how="left"))
    return out


def assemble(species: str, own: pd.DataFrame, donors: pd.DataFrame, p_ref: float,
             any_hit: set[str] | None = None) -> pd.DataFrame:
    """One row per protein, every protein present, evidence naming the tier.

    `any_hit` splits the zero-score rows into two genuinely different situations: a protein with a
    SwissProt hit that fell below the floor still has a distant relative, while one with no hit at
    all is the strongest novelty claim the axis can make. Both score 0 -- the axis will not invent
    a number for them -- but they must be distinguishable, or a third of K. pneumoniae becomes one
    undifferentiated tie at the bottom of the ranking.
    """
    df = own.merge(donors, on="uniprot_ac", how="left")
    df["studiedness_own"] = S.score(df["n_pubs"], df["annotation_score"], p_ref)
    fam_pubs = df["donor_n_pubs"].fillna(0)
    fam_ann = df["donor_annotation_score"].fillna(1.0)
    df["studiedness_family"] = np.where(df["donor_ac"].notna(),
                                        S.score(fam_pubs, fam_ann, p_ref), 0.0)
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
                p_ref: float, scope: str = "any") -> tuple[pd.DataFrame, dict]:
    """E. coli family score with every Escherichia donor removed, vs its own measured score."""
    excl = ecoli_taxids(sp_meta)
    say(f"    excluding {len(excl):,} Escherichia taxids "
        f"({int(sp_meta['donor_taxid'].isin(excl).sum()):,} SwissProt entries) from the donor pool")
    ec_hits = hits[hits["qseqid"].isin(set(own_ec["uniprot_ac"]))]
    donors = choose_donors(ec_hits, sp_meta, S.IDENTITY_FLOOR, exclude_taxids=excl, scope=scope)
    in_scope = (ec_hits if scope == "any" else ec_hits[ec_hits["donor_ac"].isin(
                    set(sp_meta.loc[sp_meta[f"donor_is_{scope}"], "donor_ac"]))])
    df = assemble("ecoli", own_ec, donors, p_ref, any_hit=set(in_scope["qseqid"]))
    # Only proteins the held-out transfer could actually score are informative about the
    # mechanism; a no_homolog row says nothing about how good the transfer is when it fires.
    scored = df[df["donor_ac"].notna()]
    rho = scored["studiedness_family"].corr(scored["studiedness_own"], method="spearman")
    pearson = scored["studiedness_family"].corr(scored["studiedness_own"])
    stats = {"n": len(df), "n_scored": len(scored),
             "donor_scope": scope,
             "coverage_pct": round(100 * len(scored) / len(df), 1),
             "spearman": round(float(rho), 4), "pearson": round(float(pearson), 4),
             "floor": CONTROL_RHO_FLOOR, "excluded_taxids": len(excl)}
    keep = ["uniprot_ac", "gene_name", "studiedness_own", "studiedness_family", "evidence",
            "donor_ac", "donor_organism", "donor_pident", "donor_n_pubs", "n_pubs"]
    return df[keep], stats


def free_route_coverage(species: str) -> dict:
    """What the on-disk four-species ortholog table reaches. Measured, never merged."""
    d = O.load(species)
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
    ap.add_argument("--scale-quantile", type=float, default=S.SCALE_QUANTILE,
                    help="quantile of the SwissProt publication distribution anchoring the log "
                         "scale (default 99; 95 saturates E. coli -- see src/studiedness.py)")
    ap.add_argument("--donor-scope", default="prokaryotic", choices=list(SCOPES),
                    help="which SwissProt entries may donate literature. ALL THREE are computed "
                         "every run and compared into evidence/donor_scope_comparison.tsv; this "
                         "picks which one becomes the shipped studiedness_family")
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
    sp_meta = load_swissprot(counts)
    pubs = sp_meta["donor_n_pubs"].to_numpy()
    p_ref = float(np.percentile(pubs, args.scale_quantile))
    if p_ref <= 0:
        sys.exit(f"FATAL the SwissProt P{args.scale_quantile:g} publication count is 0; the scale "
                 "would be undefined. Pass a higher --scale-quantile.")
    say(f"  SCALE     P{args.scale_quantile:g} = {p_ref:.0f} papers -> "
        f"log1p = {np.log1p(p_ref):.3f}")
    say("            one FIXED constant for all species and both columns (not a per-proteome "
        "percentile -- that is the geptop trap)")
    say("            distribution of reviewed-entry publication counts, for context:")
    say("              " + "   ".join(
        f"P{q:g} {np.percentile(pubs, q):,.0f}" for q in (50, 75, 90, 95, 99, 99.9))
        + f"   max {pubs.max():,.0f}")
    sat = float((S.literature_score(pubs, p_ref) >= 1.0).mean())
    say(f"            {100 * sat:.1f}% of SwissProt saturates at this anchor "
        f"(guard: {100 * S.MAX_SATURATED:.0f}%)")

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
        own = own_signals(sp, counts)
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
            scoped[scope] = assemble(sp, own, d_s, p_ref, any_hit=set(in_scope["qseqid"]))
        df = scoped[args.donor_scope].copy()
        df["donor_scope"] = args.donor_scope
        for scope in SCOPES:
            df[f"studiedness_family_{scope}"] = scoped[scope]["studiedness_family"].to_numpy()
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
                "median_family": round(float(cur["studiedness_family"].median()), 4),
                "spearman_vs_any": round(float(cur["studiedness_family"].corr(
                    base["studiedness_family"], method="spearman")), 4),
                "shipped": scope == args.donor_scope,
            })

        own_path = dest / f"{pre}own_{sp}.tsv"
        tr_path = dest / f"{pre}transfer_{sp}.tsv"
        own.to_csv(own_path, sep="\t", index=False)
        tr_cols = ["uniprot_ac", "donor_ac", "donor_gene", "donor_organism", "donor_taxid",
                   "donor_pident", "donor_qcov", "donor_n_pubs", "donor_n_pubs_uniprot",
                   "donor_n_pubs_gene2pubmed", "donor_annotation_score", "donor_protein_name",
                   "nearest_ac", "nearest_pident", "nearest_organism", "nearest_n_pubs",
                   "n_candidates", "studiedness_own", "studiedness_family", "evidence",
                   "donor_scope", "studiedness_family_prokaryotic",
                   "studiedness_family_bacteria", "studiedness_family_any",
                   f"donor_ac_{other}", f"donor_organism_{other}"]
        df[tr_cols].to_csv(tr_path, sep="\t", index=False)

        tiers = df["evidence"].value_counts().to_dict()
        covered = int(df["donor_ac"].notna().sum())
        # Unrankable ties are the failure mode this axis shares with stage 04's under-regularised
        # logistic: a score pinned at 1.0 over half a proteome cannot be sorted.
        for col in ("studiedness_own", "studiedness_family"):
            frac = float((df[col] >= 0.999).mean())
            if frac > S.MAX_SATURATED:
                sys.exit(
                    f"FATAL {sp}: {100 * frac:.1f}% of {col} is pinned at 1.0, above the "
                    f"{100 * S.MAX_SATURATED:.0f}% guard -- that many tied values cannot be "
                    f"ranked. Raise --scale-quantile (currently {args.scale_quantile:g}).")
        say(f"  {sp}")
        say(f"    {len(df):,} proteins   with a donor {covered:,} "
            f"({100 * covered / len(df):.1f}%)   "
            + "   ".join(f"{k} {v:,}" for k, v in sorted(tiers.items())))
        say(f"    studiedness_own    median {df.studiedness_own.median():.3f}   "
            f"distinct {df.studiedness_own.nunique():,}")
        say(f"    studiedness_family median {df.studiedness_family.median():.3f}   "
            f"distinct {df.studiedness_family.nunique():,}")
        top_org = df.loc[df.donor_organism.notna(), "donor_organism"].value_counts().head(3)
        for org, n in top_org.items():
            say(f"      donor organism  {n:5,}  {org}")
        say(f"    -> {tr_path.relative_to(REPO_ROOT)}")
        rows.append({"species": sp, "n": len(df), "with_donor": covered,
                     "with_donor_pct": round(100 * covered / len(df), 1),
                     **{f"tier_{k}": v for k, v in tiers.items()},
                     "median_own": round(float(df.studiedness_own.median()), 4),
                     "median_family": round(float(df.studiedness_family.median()), 4)})
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

    pd.DataFrame([{"p_ref_n_pubs": p_ref, "scale_quantile": args.scale_quantile,
                   "saturated_fraction_swissprot": round(sat, 4),
                   "w_literature": S.W_LITERATURE,
                   "w_annotation": S.W_ANNOTATION, "identity_floor": S.IDENTITY_FLOOR,
                   "identity_direct": S.IDENTITY_DIRECT, "coverage_floor": S.COVERAGE_FLOOR,
                   "swissprot_entries": len(sp_meta)}]).to_csv(
        EVIDENCE_DIR / "scale.tsv", sep="\t", index=False)
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
    own_by_sp = {sp: own_signals(sp, counts) for sp in args.species}
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
            cdf = assemble("ecoli", ec, dc, p_ref)  # bands unused here; rho only
            sc = cdf[cdf["donor_ac"].notna()]
            row["control_spearman"] = (round(float(
                sc["studiedness_family"].corr(sc["studiedness_own"], method="spearman")), 4)
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
        ec_own = own_signals("ecoli", counts)
        both = {}
        for scope in SCOPES:
            both[scope] = run_control(hits, sp_meta, ec_own, p_ref, scope=scope)
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
