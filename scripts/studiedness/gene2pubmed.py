"""Reduce NCBI gene2pubmed to the GeneIDs this project uses, and measure what it buys.

gene2pubmed is **83,245,982 rows over 43,745,777 distinct GeneIDs** -- far too many to carry
around, and we need a few hundred thousand of them: our three anchor proteomes plus every
SwissProt entry that could act as a literature donor. This script streams the file once, keeps
only that universe, and writes a small count table the rest of the axis reads.

    data/processed/studiedness/scratch/gene2pubmed_counts.tsv   geneid -> n_pubs  (the cache)
    data/processed/studiedness/evidence/gene2pubmed_<species>.tsv    per anchor protein
    data/processed/studiedness/evidence/gene2pubmed_gain.tsv         the measured verdict

IT DOES NOT RESCUE THE ANCHORS, AND THAT IS THE POINT OF MEASURING IT
----------------------------------------------------------------------
Measured 2026-09-21 against our own GeneIDs, before any of this was written:

    species   >=1 paper        distinct values   max   UniProt lit_pubmed_id for comparison
    Ec        4,175 (94.8%)    193               513   48 distinct, max 58  -- ~4x coarser
    Kp        5,342 (93.3%)    12                11    7 distinct           -- both dead
    Sa        427   (14.8%)    15                17    357 proteins with >=1 -- UniProt is better

So NCBI is dark on K. pneumoniae and S. aureus too, through a different door. What it does have
is finer RESOLUTION on well-studied organisms -- 189 distinct values against UniProt's 47 -- and
that is the measurement worth keeping.

**IT IS NOT THE SHIPPED COUNT -- that changed on 2026-09-22.** `n_papers_family` counts UniProt
curated references only, because "papers a curator read and used" is ONE consistent definition
applied to every row, whereas `max(curated, gene-linked)` switched definition per protein: Kp
`rpoB` took NCBI's 350 while ~11% of donors took SwissProt's number. One consistent definition
beat the bigger number.

This script therefore documents a **measured alternative** (the `interpro2go` precedent) rather
than feeding the deliverable. Its outputs stay in `evidence/`, both per-source counts ship in
`evidence/transfer_<species>.tsv`, and the 287 MB download is optional --
`fetch.py --only anchors swissprot unknome` skips it.

WHY THE JOIN IS FREE
--------------------
`proteome_<species>.tsv` already carries a `geneid` column (Kp 100% / Ec 95.1% / Sa 94.7%), and
the SwissProt metadata stream carries `xref_geneid`. No id-mapping round-trip anywhere.

Run with the `gradi` env, after `scripts/studiedness/fetch.py`.
  python scripts/studiedness/gene2pubmed.py
  python scripts/studiedness/gene2pubmed.py --refresh -q
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
LIT_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "literature"
GENE2PUBMED = REPO_ROOT / "data" / "source" / "ncbi" / "gene" / "gene2pubmed.gz"
COUNTS = SCRATCH_DIR / "gene2pubmed_counts.tsv"

SPECIES = ("kpneumoniae", "ecoli", "saureus")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def split_ids(cell: object) -> list[str]:
    """Split a `;`-joined UniProt list field. Empty, NaN and whitespace-only all give []."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    return [x.strip() for x in str(cell).split(";") if x.strip()]


def load_anchor_literature(species: str) -> pd.DataFrame:
    """One anchor's UniProt literature table, as fetched by fetch.py."""
    path = LIT_DIR / f"anchor_{species}.tsv.gz"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return df.rename(columns={"Entry": "uniprot_ac", "Reviewed": "reviewed",
                              "Annotation": "annotation_score",
                              "Protein existence": "protein_existence",
                              "PubMed ID": "lit_pubmed_id", "GeneID": "geneid"})


def load_swissprot_geneids() -> set[str]:
    """Every GeneID reachable from a reviewed entry -- the donor half of the universe."""
    path = LIT_DIR / "swissprot_meta.tsv.gz"
    if not path.exists():
        sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    ids: set[str] = set()
    n = 0
    for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False,
                             usecols=["GeneID"], chunksize=100_000):
        n += len(chunk)
        for cell in chunk["GeneID"]:
            ids.update(split_ids(cell))
    say(f"  swissprot      {n:,} entries -> {len(ids):,} distinct GeneIDs")
    return ids


def build_universe() -> tuple[set[str], dict[str, pd.DataFrame]]:
    """The GeneIDs worth counting: our anchors plus every SwissProt donor."""
    anchors = {}
    universe: set[str] = set()
    for sp in SPECIES:
        lit = load_anchor_literature(sp)
        prot = P.load_full(sp)[["uniprot_ac", "gene_name", "geneid"]]
        # The anchor stream and the stage-00 table are the same proteome; join to keep gene_name
        # for the spot checks, and prefer stage-00's geneid since that is the documented key.
        merged = prot.merge(lit.drop(columns=["geneid"]), on="uniprot_ac", how="left")
        if merged["annotation_score"].isna().any():
            missing = int(merged["annotation_score"].isna().sum())
            sys.exit(f"FATAL {sp}: {missing:,} proteins have no row in the UniProt literature "
                     "stream. The two fetches disagree -- re-run fetch.py --refresh.")
        anchors[sp] = merged
        gids = merged["geneid"].apply(split_ids)
        flat = {g for L in gids for g in L}
        universe |= flat
        say(f"  {sp:<14} {len(merged):,} proteins, "
            f"{int(gids.apply(bool).sum()):,} with a GeneID "
            f"({100 * gids.apply(bool).mean():.1f}%) -> {len(flat):,} distinct")
    universe |= load_swissprot_geneids()
    return universe, anchors


def stream_counts(universe: set[str], refresh: bool) -> dict[str, int]:
    """One pass over gene2pubmed, counting PMIDs for GeneIDs in the universe."""
    if COUNTS.exists() and not refresh:
        df = pd.read_csv(COUNTS, sep="\t", dtype={"geneid": str, "n_pubs": int})
        say(f"  cached {COUNTS.relative_to(REPO_ROOT)}: {len(df):,} GeneIDs")
        return dict(zip(df["geneid"], df["n_pubs"]))
    if not GENE2PUBMED.exists():
        sys.exit(f"FATAL missing {GENE2PUBMED.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    counts: dict[str, int] = {}
    seen = 0
    t0 = time.time()
    with gzip.open(GENE2PUBMED, "rt", errors="replace") as fh:
        header = fh.readline()
        if not header.startswith("#tax_id"):
            sys.exit(f"FATAL gene2pubmed header changed: {header!r}")
        for line in fh:
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            gid = parts[1]
            if gid in universe:
                counts[gid] = counts.get(gid, 0) + 1
            seen += 1
            if seen % 20_000_000 == 0:
                say(f"    {seen:,} rows, {len(counts):,} of our GeneIDs seen, "
                    f"{time.time() - t0:.0f}s")
    say(f"  streamed {seen:,} rows in {time.time() - t0:.0f}s")
    say(f"  {len(counts):,} of {len(universe):,} GeneIDs in our universe carry at least one paper "
        f"({100 * len(counts) / max(len(universe), 1):.1f}%)")
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    (pd.DataFrame({"geneid": list(counts), "n_pubs": list(counts.values())})
       .sort_values("geneid").to_csv(COUNTS, sep="\t", index=False))
    say(f"  wrote {COUNTS.relative_to(REPO_ROOT)}")
    return counts


def score_species(sp: str, anchor: pd.DataFrame, counts: dict[str, int]) -> tuple[pd.DataFrame, dict]:
    """Per-protein gene2pubmed count beside UniProt's, plus the union."""
    out = anchor[["uniprot_ac", "gene_name", "geneid"]].copy()
    out["n_pubs_uniprot"] = anchor["lit_pubmed_id"].apply(lambda c: len(split_ids(c)))
    out["n_pubs_gene2pubmed"] = anchor["geneid"].apply(
        lambda c: sum(counts.get(g, 0) for g in split_ids(c)))
    # The union is over PMIDs, not a sum of two counts -- gene2pubmed and UniProt overlap heavily
    # and adding them would double-count the same papers. gene2pubmed does not ship PMIDs here
    # (only counts), so the union is bounded below by the max and above by the sum; max is the
    # honest choice and is what downstream reads.
    out["n_pubs"] = out[["n_pubs_uniprot", "n_pubs_gene2pubmed"]].max(axis=1)
    gained = int(((out.n_pubs_gene2pubmed > 0) & (out.n_pubs_uniprot == 0)).sum())
    better = int((out.n_pubs_gene2pubmed > out.n_pubs_uniprot).sum())
    stats = {
        "species": sp, "n": len(out),
        "uniprot_ge1": int((out.n_pubs_uniprot > 0).sum()),
        "uniprot_distinct": int(out.n_pubs_uniprot.nunique()),
        "uniprot_max": int(out.n_pubs_uniprot.max()),
        "gene2pubmed_ge1": int((out.n_pubs_gene2pubmed > 0).sum()),
        "gene2pubmed_distinct": int(out.n_pubs_gene2pubmed.nunique()),
        "gene2pubmed_max": int(out.n_pubs_gene2pubmed.max()),
        "union_ge1": int((out.n_pubs > 0).sum()),
        "union_distinct": int(out.n_pubs.nunique()),
        "first_literature_from_ncbi": gained,
        "ncbi_higher_than_uniprot": better,
    }
    return out, stats


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--refresh", action="store_true", help="re-stream gene2pubmed")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STUDIEDNESS - NCBI gene2pubmed")
    rule("=")
    say("  in    data/source/ncbi/gene/gene2pubmed.gz  (83.2M rows, 43.7M distinct GeneIDs)")
    say("        data/source/uniprot/literature/{anchor_*,swissprot_meta}.tsv.gz")
    say(f"  out   {COUNTS.relative_to(REPO_ROOT)}")
    say(f"        {(EVIDENCE_DIR / 'gene2pubmed_<species>.tsv').relative_to(REPO_ROOT)}")
    say("  joins on the geneid column already in proteome_<species>.tsv -- no id mapping")
    say("  it does NOT rescue the dark anchors; it sharpens the DONOR ranking")
    rule("=")
    if args.dry_run:
        say("dry run: nothing written.")
        return

    rule()
    say("UNIVERSE")
    rule()
    universe, anchors = build_universe()
    say(f"  total {len(universe):,} distinct GeneIDs to count")

    rule()
    say("STREAM gene2pubmed")
    rule()
    counts = stream_counts(universe, args.refresh)

    rule()
    say("PER-SPECIES - gene2pubmed against UniProt's own reference list")
    rule()
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp in args.species:
        out, stats = score_species(sp, anchors[sp], counts)
        path = EVIDENCE_DIR / f"gene2pubmed_{sp}.tsv"
        out.to_csv(path, sep="\t", index=False)
        rows.append(stats)
        say(f"  {sp}")
        say(f"    uniprot      >=1 paper {stats['uniprot_ge1']:5,} "
            f"({100 * stats['uniprot_ge1'] / stats['n']:5.1f}%)   "
            f"distinct {stats['uniprot_distinct']:4,}   max {stats['uniprot_max']:4,}")
        say(f"    gene2pubmed  >=1 paper {stats['gene2pubmed_ge1']:5,} "
            f"({100 * stats['gene2pubmed_ge1'] / stats['n']:5.1f}%)   "
            f"distinct {stats['gene2pubmed_distinct']:4,}   max {stats['gene2pubmed_max']:4,}")
        say(f"    union        >=1 paper {stats['union_ge1']:5,} "
            f"({100 * stats['union_ge1'] / stats['n']:5.1f}%)   "
            f"distinct {stats['union_distinct']:4,}")
        say(f"    NCBI gives first literature to {stats['first_literature_from_ncbi']:,} proteins; "
            f"a higher count to {stats['ncbi_higher_than_uniprot']:,}")
        say(f"    -> {path.relative_to(REPO_ROOT)}")

    gain = pd.DataFrame(rows)
    gain.to_csv(EVIDENCE_DIR / "gene2pubmed_gain.tsv", sep="\t", index=False)
    rule("=")
    say(f"  wrote {(EVIDENCE_DIR / 'gene2pubmed_gain.tsv').relative_to(REPO_ROOT)}")
    say("  VERDICT: read gene2pubmed_gain.tsv. The anchors stay dark either way -- the value is")
    say("  in the donors, which transfer.py reads from the counts cache above.")
    rule("=")


if __name__ == "__main__":
    main()
