"""OGEE v3 -> an ML-ready essentiality label corpus, and an honest account of what is still missing.

THE POINT OF THIS SCRIPT IS THE TAXON REPORT, not the label table. OGEE gives labels for 89 taxa
and no sequences; ProteomeLM-Ess needs an embedding per protein. So the work between here and a
trained classifier is entirely IDENTIFIER MAPPING, and this script measures exactly how much of it
there is -- per taxon, per namespace -- before anyone downloads 89 proteomes.

WHAT THE PAPER DID, so we can say where we diverge (Malbranke 2026, Materials and Methods):
  * labels from OGEE v3; sequences recovered "by matching the gene names (gene IDs) provided by
    OGEE" against UniProt, NCBI, SGD and the Fitness Browser -> 87 taxonomic IDs
  * "relying on curated complete proteomes whenever possible" to avoid isoform/duplicate ambiguity;
    remaining duplicates merged, keeping both IDs
  * 83 of 87 genomes for training, HOLDING OUT S. cerevisiae and FOUR E. coli strains
  * split by clustering proteins ACROSS ALL GENOMES with MMseqs2 at 40% identity, then assigning
    whole clusters to train/val/test -- "at the protein level and not at the genome level, to avoid
    the model relying on sequence similarity between, say, two orthologs in similar genomes, as a
    shortcut"
  * head: 2-layer FC, hidden 2048, ReLU, dropout 0.5, two logits + softmax, inputs z-scored
    genome-wide; early stopping on validation. Best config ProteomeLM-L layer 8, AUC 0.93.

THREE THINGS THE PAPER DOES NOT SETTLE, all measured here because they change the corpus:

  1. THE NAMESPACE IS A PROPERTY OF THE DATASET BLOCK, NOT THE TAXON. `locus` is 100% populated but
     heterogeneous, and `Ref_db` is null for 158,340 of 236,144 rows. Measured: each
     (taxon, dataset) block is a SINGLE namespace at median purity 1.000, and 87 of 89 taxa use one
     namespace throughout. The two exceptions are real and matter -- E. coli K-12 carries PEC gene
     symbols in datasets 108/109 and b-numbers in 173/174/175, sharing ZERO ids, so a naive
     (taxid, locus) dedup counts every E. coli gene twice (9,496 "genes" for a 4,403-gene
     organism). Dedup must happen AFTER resolving to a common key, never before.
  2. DATASETS WITHIN A TAXON DISAGREE. 4,786 of 26,346 multi-dataset genes (18.2%) are called
     essential by one screen and not by another -- the same union-vs-intersection choice `deg_ess`
     records, and on E. coli it is a 1.8x spread (1,049 vs 584 essential). `--rule` makes it
     explicit rather than silent.
  3. TWO TAXA HAVE NO POSITIVES AT ALL. Taxids 187410 and 216596 carry 7,463 E/NE rows with zero
     `E`. An organism with no essential genes is not a biological finding; they are EXCLUDED, and
     training on them would teach the model that whole genomes are dispensable.

Run with the `gradi` env:
    python scripts/essentiality/ogee.py                 # audit only; no downloads
    python scripts/essentiality/ogee.py --rule any
    python scripts/essentiality/ogee.py --no-names      # skip the NCBI taxonomy lookup
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

OGEE = REPO_ROOT / "data" / "source" / "ogee" / "gene_essentiality.txt.gz"
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
TAXON_CACHE = SCRATCH_DIR / "ogee_taxon_names.json"

# Only the decided two-state calls. OGEE's other 18 values are conditional, ambiguous or
# dataset-specific (`C` conditional 11,831, `ND` not determined 1,728, `F` fitness-defect 1,638,
# `U` unknown 1,610, `ES`/`GD`/`S`/`DE`/`D`/`GA`, plus compound strings like
# `NE,E-infection,E-co-infection`). A conditional call is not an essentiality call, and the
# compound ones are in-host screens -- both out of scope for the same reason Choe's M9 arm is.
KEEP = ("E", "NE")

# Zero-positive taxa: see docstring point 3.
EXCLUDE_TAXA = {187410, 216596}

# Databases whose ids are eukaryotic model-organism accessions. These need per-database resolvers
# (SGD, PomBase, FlyBase, WormBase, ...) rather than the NCBI locus-tag route, and their proteomes
# are large -- so they are FLAGGED, not dropped, and reported as their own tier.
EUKARYOTE_DB = {"WormBase", "flybase", "ToxoDB", "PlasmoDB", "MGI", "SGD", "PomBase",
                "CGD", "DictyBase", "TAIR", "Ensembl", "ZFIN", "RGD", "AspGD", "TritrypDB"}

NCBI_TAXON = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/taxonomy/taxon/{}"
              "/dataset_report")

VERBOSE = True


def say(m: str = "") -> None:
    if VERBOSE:
        print(m, flush=True)


def rule(c: str = "-", w: int = 118) -> None:
    say(c * w)


# ------------------------------------------------------------------------ namespace classification

def namespace(s: str) -> str:
    """Which id space is this string in? Classified per DATASET BLOCK, where it is ~pure.

    Deliberately coarse: the point is to say WHICH RESOLVER a block needs, not to parse the id.
    """
    s = str(s)
    if re.match(r"^(WBGene|FBgn|ENS|Y[A-P][LR]\d|SP[AB][CP]|PF3D7|TGME49|LmjF|Tb\d|"
                r"ZDB-GENE|S0000|At\dg|AT\dG|Afu\dg|An\d\dg|AO\d{9}|EBG000)", s):
        return "model_organism_db_id"
    if re.fullmatch(r"\d{4,9}", s):
        return "ncbi_geneid"          # bare NCBI GeneID -- resolvable via gene2accession
    if re.match(r"^[A-Za-z][A-Za-z0-9]*_[A-Za-z0-9_]+$", s):
        return "locus_tag"            # PREFIX_NNNNN, incl. plasmid tags like SL1344_P3_0014
    if re.match(r"^[A-Za-z]{1,4}\d{3,5}(\.\d+)?[A-Za-z]?$", s):
        return "locus_tag"            # Rv0098, b4233, t0001, HI0220.1, PA0195.1, ACIAD3542
    if re.match(r"^[a-z]{3}[A-Z]?\d*$", s):
        return "gene_symbol"
    return "unclassified"


def load_ogee() -> pd.DataFrame:
    if not OGEE.exists():
        sys.exit(f"FATAL {OGEE} missing -- see data/source/ogee/SOURCE.md for the Wayback and "
                 "GitHub mirrors (the `id_` suffix on the Wayback URL is mandatory)")
    d = pd.read_csv(OGEE, sep="\t", low_memory=False)
    need = {"dataset", "taxaID", "locus", "essentiality", "pmid", "Ref_db"}
    missing = need - set(d.columns)
    if missing:
        sys.exit(f"FATAL OGEE columns changed: missing {sorted(missing)}")
    return d


def taxon_names(taxids: list[int], enabled: bool) -> dict[int, str]:
    """taxid -> scientific name, from NCBI, cached. Names are for the REPORT only.

    OGEE ships no organism name -- only `taxaID` -- so a human cannot read the corpus without
    this. Failure is non-fatal: an unnamed taxon is still a countable taxon.
    """
    cache: dict[str, str] = {}
    if TAXON_CACHE.exists():
        cache = json.loads(TAXON_CACHE.read_text())
    if not enabled:
        return {t: cache.get(str(t), "") for t in taxids}
    todo = [t for t in taxids if str(t) not in cache]
    if todo:
        say(f"  resolving {len(todo)} taxon names from NCBI (cached afterwards) ...")
    for i, t in enumerate(todo, 1):
        try:
            req = urllib.request.Request(NCBI_TAXON.format(t),
                                         headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=30) as fh:
                js = json.load(fh)
            rep = (js.get("reports") or [{}])[0].get("taxonomy", {})
            cache[str(t)] = rep.get("current_scientific_name", {}).get("name", "") or ""
        except (urllib.error.URLError, TimeoutError, ValueError, KeyError):
            cache[str(t)] = ""
        if i % 20 == 0:
            say(f"    {i}/{len(todo)}")
        time.sleep(0.12)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    TAXON_CACHE.write_text(json.dumps(cache, indent=0, sort_keys=True))
    return {t: cache.get(str(t), "") for t in taxids}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--rule", choices=["any", "all", "majority"], default="any",
                    help="how to combine datasets that disagree within a taxon (18.2%% of "
                         "multi-dataset genes do). `any` matches the permissive reading; `all` is "
                         "the conservative one and on E. coli they differ 1,049 vs 584.")
    ap.add_argument("--no-names", action="store_true", help="skip the NCBI taxonomy lookup")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("scripts/essentiality/ogee.py -- OGEE v3 label corpus + what it would take to use it")
    rule("=")
    say(f"  in    {OGEE.relative_to(REPO_ROOT)}")
    say(f"  out   {EVIDENCE_DIR.relative_to(REPO_ROOT)}/ogee_labels.tsv + ogee_taxa.tsv")
    say(f"  rule  {a.rule}   (datasets disagree on 18.2% of multi-dataset genes)")
    rule()

    raw = load_ogee()
    say(f"  {len(raw):,} rows, {raw.taxaID.nunique()} taxa, {raw.dataset.nunique()} datasets")
    d = raw[raw.essentiality.isin(KEEP)].copy()
    say(f"  {len(d):,} rows with a decided E/NE call  "
        f"({len(raw) - len(d):,} conditional/ambiguous dropped)")
    d = d[~d.taxaID.isin(EXCLUDE_TAXA)]
    say(f"  {len(d):,} rows after excluding the zero-positive taxa {sorted(EXCLUDE_TAXA)}")
    d["y"] = (d.essentiality == "E").astype(int)

    # namespace per (taxon, dataset) block, plus the purity that justifies treating it as one
    d["id_kind_row"] = d.locus.map(namespace)
    blk = (d.groupby(["taxaID", "dataset"])["id_kind_row"]
             .agg(id_kind=lambda s: s.value_counts().idxmax(),
                  id_purity=lambda s: round(s.value_counts(normalize=True).max(), 4))
             .reset_index())
    d = d.merge(blk, on=["taxaID", "dataset"], how="left")

    rule()
    say("NAMESPACE PER DATASET BLOCK  -- which resolver each block needs")
    rule()
    say(f"  {'namespace':24s} {'blocks':>7} {'taxa':>6} {'rows':>9} {'median purity':>14}")
    for k, g in d.groupby("id_kind"):
        nb = g.groupby(["taxaID", "dataset"]).ngroups
        say(f"  {k:24s} {nb:>7} {g.taxaID.nunique():>6} {len(g):>9,} "
            f"{g.id_purity.median():>14.3f}")

    # tier: which resolution route a taxon needs
    d["tier"] = "prokaryote_ncbi"
    d.loc[d.Ref_db.isin(EUKARYOTE_DB), "tier"] = "eukaryote_model_db"
    d.loc[d.id_kind == "model_organism_db_id", "tier"] = "eukaryote_model_db"

    # ---- labels, deduped WITHIN a namespace only (never across -- see docstring point 1)
    grp = d.groupby(["taxaID", "id_kind", "locus"])
    lab = grp.agg(n_datasets=("y", "size"), n_essential=("y", "sum"),
                  tier=("tier", "first"), ref_db=("Ref_db", "first"),
                  pmids=("pmid", lambda s: ";".join(sorted({str(x) for x in s})))).reset_index()
    if a.rule == "any":
        lab["label"] = (lab.n_essential > 0).astype(int)
    elif a.rule == "all":
        lab["label"] = (lab.n_essential == lab.n_datasets).astype(int)
    else:
        lab["label"] = (lab.n_essential * 2 > lab.n_datasets).astype(int)
    lab["disagrees"] = ((lab.n_essential > 0) & (lab.n_essential < lab.n_datasets)).astype(int)

    rule()
    say("LABELS  -- one row per (taxon, namespace, locus). NOT yet one row per PROTEIN")
    rule()
    say(f"  {len(lab):,} labelled entries, {int(lab.label.sum()):,} essential "
        f"({lab.label.mean():.1%})")
    say(f"  {int(lab.disagrees.sum()):,} carry a cross-dataset disagreement the `{a.rule}` rule "
        "resolved")
    for t, g in lab.groupby("tier"):
        say(f"  {t:22s} {len(g):>8,} entries  {g.taxaID.nunique():>3} taxa  "
            f"base {g.label.mean():.3f}")
    say("  NOTE: E. coli K-12 appears under TWO namespaces and is counted twice here. Collapsing "
        "requires resolving both to sequences first -- that is the missing step, not a bug.")

    # ---- per-taxon report: the actual deliverable
    names = taxon_names(sorted(lab.taxaID.unique().tolist()), not a.no_names)
    tax = (lab.groupby("taxaID")
              .agg(n_entries=("label", "size"), n_essential=("label", "sum"),
                   base_rate=("label", "mean"), n_namespaces=("id_kind", "nunique"),
                   id_kind=("id_kind", lambda s: "+".join(sorted(set(s)))),
                   tier=("tier", "first"), n_disagree=("disagrees", "sum"))
              .reset_index())
    tax["organism"] = tax.taxaID.map(names)
    tax["base_rate"] = tax.base_rate.round(4)
    ndat = d.groupby("taxaID").dataset.nunique().rename("n_datasets")
    tax = tax.merge(ndat, on="taxaID", how="left")
    tax = tax.sort_values("n_entries", ascending=False)

    rule()
    say("TAXA  -- top 20 by size. Full table in ogee_taxa.tsv")
    rule()
    say(f"  {'taxid':>8} {'organism':38s} {'entries':>8} {'ess':>6} {'base':>6} {'ds':>3} "
        f"{'namespace':22s} tier")
    for _, r in tax.head(20).iterrows():
        say(f"  {r.taxaID:>8} {str(r.organism)[:38]:38s} {r.n_entries:>8,} "
            f"{int(r.n_essential):>6,} {r.base_rate:>6.3f} {r.n_datasets:>3} "
            f"{r.id_kind[:22]:22s} {r.tier}")

    rule()
    say("BASE RATE SPREAD  -- check this before training; it is not a constant")
    rule()
    q = tax.base_rate.quantile([0, .25, .5, .75, 1]).round(3)
    say(f"  min {q[0]}  q25 {q[.25]}  median {q[.5]}  q75 {q[.75]}  max {q[1]}")
    say(f"  lowest : " + ", ".join(f"{r.organism or r.taxaID} {r.base_rate:.3f}"
                                   for _, r in tax.nsmallest(3, 'base_rate').iterrows()))
    say(f"  highest: " + ", ".join(f"{r.organism or r.taxaID} {r.base_rate:.3f}"
                                   for _, r in tax.nlargest(3, 'base_rate').iterrows()))

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    lab.to_csv(EVIDENCE_DIR / "ogee_labels.tsv", sep="\t", index=False)
    tax.to_csv(EVIDENCE_DIR / "ogee_taxa.tsv", sep="\t", index=False)
    rule()
    say(f"  wrote {(EVIDENCE_DIR / 'ogee_labels.tsv').relative_to(REPO_ROOT)}  ({len(lab):,} rows)")
    say(f"  wrote {(EVIDENCE_DIR / 'ogee_taxa.tsv').relative_to(REPO_ROOT)}  ({len(tax)} rows)")

    rule("=")
    say("WHAT IS STILL MISSING  -- no sequences are fetched by this script")
    rule("=")
    pro = tax[tax.tier == "prokaryote_ncbi"]
    euk = tax[tax.tier == "eukaryote_model_db"]
    say(f"  prokaryote_ncbi     {len(pro):>3} taxa, {pro.n_entries.sum():>7,} entries  "
        "-> taxid to assembly to protein FASTA + GFF, join on locus_tag")
    say(f"  eukaryote_model_db  {len(euk):>3} taxa, {euk.n_entries.sum():>7,} entries  "
        "-> needs a per-database resolver (SGD/PomBase/FlyBase/WormBase/...)")
    say("")
    say("  Then, per the paper: ESM-C per proteome, one ProteomeLM forward per proteome, MMseqs2")
    say("  clustering at 40% ACROSS ALL GENOMES for the split, and the 2-layer head.")
    rule("=")


if __name__ == "__main__":
    main()
