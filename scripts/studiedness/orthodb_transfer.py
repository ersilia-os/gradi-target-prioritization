"""REJECTED ALTERNATIVE -- literature via OrthoDB orthogroups, for proteins DIAMOND cannot reach.

**MEASURED AND REJECTED 2026-10-04. DO NOT SHIP THIS COLUMN.** Kept so the comparison is
reproducible, the way `gradi-lazyqsar` is kept for degradability. See docs/studiedness.md 3d.

**The transferred count carries no information.** Held-out E. coli, the anchor's genus struck out
of the member pool: reached-only spearman **-0.1068** (narrow -0.0852, domain -0.0059), while a
bare reached-vs-not flag scores **0.5058**. The apparent 0.4140 over all rows is presence/absence,
not literature. Filling those cells would be worse than leaving them blank -- a blank says
"unknown", a 4 says "four papers".

**And read `validate()` before trusting any number this prints.** The first version of the control
read **0.8351** and was pure self-correlation: a protein's orthogroup contains its OWN SwissProt
entry, so 98.3% of the "donors" were Escherichia. That is invisible in a correlation.

Literature for the proteins DIAMOND cannot reach, via OrthoDB orthogroup membership.

    data/processed/studiedness/scratch/orthodb_og_members.tsv    og -> SwissProt members (cache)
    data/processed/studiedness/evidence/orthodb_transfer_<sp>.tsv   per protein
    data/processed/studiedness/evidence/orthodb_control.tsv         the held-out validation

THE GAP THIS EXISTS TO FILL
-----------------------------
`transfer.py` finds a curated donor by sequence identity at a 40% floor. For **1,961 of 5,728
K. pneumoniae proteins (34%)** and **1,276 of 2,889 S. aureus (44%)** it finds nothing at all --
`no_hit` or `below_floor` -- so every literature column reads 0 for them and no choice of count or
corpus changes that. It is the axis's hard ceiling, and the only candidate that attacks it is
curated ORTHOLOGY rather than raw identity: **1,558 of those Kp proteins (79%) do carry an OrthoDB
orthogroup**, and a group can hold a well-studied member below 40% identity.

THIS IS A DIFFERENT MECHANISM, SO IT IS A DIFFERENT COLUMN
-------------------------------------------------------------
Group membership is not sequence identity. Mixing the two into one count would be the same
definition-switching that killed the 0-1 composite on 2026-09-22, so this writes its own column
with its own provenance and the held-out control decides whether it earns a place.

NARROW FIRST, THEN DOMAIN, AND THE LEVEL IS RECORDED
-------------------------------------------------------
OrthoDB assigns each protein a group at several levels. The narrow one is the most specific and the
better orthology claim, but **457 of the 1,558 Kp narrow groups sit at `570|Klebsiella`** -- a genus
whose members are all as dark as the query, so the narrow group is often empty of literature. The
domain group (`2|Bacteria`) always exists and always has well-studied members, at the cost of a much
weaker orthology claim. Both are tried, narrow first, and `orthodb_level_used` says which fired.

WHY THE STREAM IS FILTERED IN THIS ORDER
-------------------------------------------
A bacterial domain-level group can hold one member per species across 17,551 species. Collecting
members first and intersecting with SwissProt afterwards would hold tens of millions of gene ids in
memory. So `genes.tab` is streamed FIRST to learn which OrthoDB gene ids are SwissProt-backed (at
most 575,748), and `OG2genes` is then filtered against that set. Two passes, ~9 GB, nothing written
to disk -- the dumps are public and re-derivable and **must never go to eosvc**.

Run with the `gradi` env, after `transfer.py`. ~15 min cold, seconds from the cache.
  python scripts/studiedness/orthodb_transfer.py
  python scripts/studiedness/orthodb_transfer.py --species kpneumoniae -q
"""

from __future__ import annotations

import argparse
import gzip
import sys
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import orthology as O  # noqa: E402
from src import studiedness as S  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
LIT_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "literature"
SPROT_META = LIT_DIR / "swissprot_meta.tsv.gz"
MEMBERS = SCRATCH_DIR / "orthodb_og_members.tsv"
# Pass 1's result, cached separately. It is 19 minutes of streaming that does NOT depend on which
# groups are asked for, so widening the target set must not pay for it again.
GENE2AC = SCRATCH_DIR / "orthodb_gene2ac.tsv"

ORTHODB_BASE = "https://data.orthodb.org/v12/download/odb_data_dump"
GENES = "odb12v2_genes.tab.gz"
OG2GENES = "odb12v2_OG2genes.tab.gz"

SPECIES = ("kpneumoniae", "ecoli", "saureus")
NO_DONOR = ("no_hit", "below_floor")
# The genus to strike out of the MEMBER pool when validating, mirroring what the sequence control
# does to the donor pool. Without it the group simply returns the query's own SwissProt entry.
ANCHOR_GENUS = {"kpneumoniae": "Klebsiella", "ecoli": "Escherichia", "saureus": "Staphylococcus"}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 96) -> None:
    say(char * width)


def stream_lines(name: str):
    """Stream-decompress one of the multi-GB dumps. Never written to disk."""
    req = urllib.request.Request(f"{ORTHODB_BASE}/{name}",
                                 headers={"User-Agent": "gradi/2.0"})
    resp = urllib.request.urlopen(req, timeout=3600)
    with gzip.GzipFile(fileobj=resp) as gz:
        for raw in gz:
            yield raw.decode("utf8", errors="replace").rstrip("\n")


# ------------------------------------------------------------------ swissprot side


def swissprot_literature() -> pd.DataFrame:
    """Reviewed entries with their curated paper count and whether they are prokaryotic."""
    if not SPROT_META.exists():
        sys.exit(f"FATAL missing {SPROT_META.relative_to(REPO_ROOT)} -- run fetch.py first")
    m = pd.read_csv(SPROT_META, sep="\t", dtype=str, keep_default_na=False)
    m["n_pubs"] = m["PubMed ID"].map(lambda s: len([x for x in s.strip(";").split(";") if x]))
    # Same donor scope as transfer.py: not eukaryotic, so Bacteria, Archaea AND phages.
    m["prokaryotic"] = ~m["Taxonomic lineage"].str.contains("Eukaryota (domain)", regex=False,
                                                            na=False)
    say(f"  swissprot      {len(m):,} reviewed entries, "
        f"{int(m.prokaryotic.sum()):,} prokaryotic, "
        f"{int((m.n_pubs > 0).sum()):,} with at least one curated paper")
    return m[["Entry", "Organism", "Gene Names (primary)", "n_pubs", "prokaryotic"]].rename(
        columns={"Entry": "donor_ac", "Organism": "donor_organism",
                 "Gene Names (primary)": "donor_gene"})


# ---------------------------------------------------------------------- the targets


def targets() -> tuple[pd.DataFrame, set[str]]:
    """Every protein with no sequence donor, plus the held-out E. coli validation set."""
    rows = []
    for sp in SPECIES:
        t = S.load_transfer(sp)
        od = O.load_orthodb(sp)[["uniprot_ac", "orthodb_og_narrow", "orthodb_narrow_level",
                                 "orthodb_og_domain", "orthodb_verdict"]]
        m = t[["uniprot_ac", "evidence"]].merge(od, on="uniprot_ac", how="left")
        nod = m[m["evidence"].isin(NO_DONOR)].copy()
        nod["species"] = sp
        nod["track"] = "no_donor"
        rows.append(nod)
        with_og = int(nod["orthodb_verdict"].ne("no_group").sum())
        say(f"  {sp:<13} {len(nod):>5,} with no sequence donor, "
            f"{with_og:>5,} of them carry an OrthoDB group "
            f"({100 * with_og / max(len(nod), 1):.1f}%)")

    # The validation set: E. coli proteins that LOSE their donor when Escherichia is held out.
    # They are the only place this route can be tested, because they have measured own literature.
    # ALL E. coli proteins, not only the donorless ones. The validation needs every protein whose
    # own literature is measured, and restricting it to donorless proteins left 18 usable rows
    # once the genus exclusion below was applied.
    od = O.load_orthodb("ecoli")[["uniprot_ac", "orthodb_og_narrow", "orthodb_narrow_level",
                                  "orthodb_og_domain", "orthodb_verdict"]]
    v = od.copy()
    v["species"] = "ecoli"
    v["track"] = "control"
    v["evidence"] = "heldout"
    rows.append(v)
    say(f"  {'control':<13} {len(v):>5,} E. coli proteins -- the validation set, scored with "
        "Escherichia struck out of the MEMBER pool")

    tgt = pd.concat(rows, ignore_index=True)
    ogs = set()
    for c in ("orthodb_og_narrow", "orthodb_og_domain"):
        v = tgt[c].astype(str).str.strip()
        ogs |= set(v[v.ne("") & v.ne("nan") & v.ne("<NA>")])
    say(f"  targets        {len(tgt):,} proteins -> {len(ogs):,} distinct orthogroups to resolve")
    return tgt, ogs


# ------------------------------------------------------------------- the two streams


def og_members(ogs: set[str], sprot: set[str], refresh: bool) -> pd.DataFrame:
    """og_id -> the SwissProt accessions in it. Two streamed passes, nothing stored."""
    if MEMBERS.exists() and not refresh:
        df = pd.read_csv(MEMBERS, sep="\t", dtype=str, keep_default_na=False)
        say(f"  cached {MEMBERS.relative_to(REPO_ROOT)}: {len(df):,} og-member pairs")
        return df

    # PASS 1 -- which OrthoDB gene ids are SwissProt-backed? Column 5 is the UniProt accession.
    # Done FIRST so the second pass can be filtered against a set of at most 575,748.
    gene2ac: dict[str, str] = {}
    if GENE2AC.exists() and not refresh:
        for line in GENE2AC.read_text().splitlines()[1:]:
            g, _, a = line.partition("\t")
            gene2ac[g] = a
        say(f"  cached {GENE2AC.relative_to(REPO_ROOT)}: {len(gene2ac):,} SwissProt-backed genes")
    if not gene2ac:
        say(f"  streaming {GENES} (~4.5 GB) for the SwissProt-backed gene ids ...")
        rows = 0
        t0 = time.time()
        for line in stream_lines(GENES):
            rows += 1
            f = line.split("\t")
            if len(f) < 5:
                continue
            ac = f[4].strip()
            if ac and ac in sprot:
                gene2ac[f[0]] = ac
            if VERBOSE and rows % 20_000_000 == 0:
                say(f"    {rows:,} rows  {len(gene2ac):,} swissprot-backed  "
                    f"{time.time() - t0:.0f}s")
        say(f"  {rows:,} gene rows in {time.time() - t0:.0f}s -> "
            f"{len(gene2ac):,} SwissProt-backed OrthoDB genes")
        if len(gene2ac) < 10_000:
            sys.exit(f"FATAL only {len(gene2ac):,} SwissProt-backed genes -- the UniProt column "
                     "moved or the stream truncated. Expected hundreds of thousands.")
        SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
        with GENE2AC.open("w") as fh:
            fh.write("gene_id\tuniprot_ac\n")
            for g, a in gene2ac.items():
                fh.write(f"{g}\t{a}\n")
        say(f"  wrote {GENE2AC.relative_to(REPO_ROOT)}")

    # PASS 2 -- members of the target groups, keeping only the SwissProt-backed ones.
    say(f"  streaming {OG2GENES} (~4.5 GB) for the {len(ogs):,} target groups ...")
    out: list[tuple[str, str]] = []
    rows = 0
    t0 = time.time()
    for line in stream_lines(OG2GENES):
        rows += 1
        og, _, gid = line.partition("\t")
        if og in ogs:
            ac = gene2ac.get(gid)
            if ac:
                out.append((og, ac))
        if VERBOSE and rows % 50_000_000 == 0:
            say(f"    {rows:,} rows  {len(out):,} kept  {time.time() - t0:.0f}s")
    say(f"  {rows:,} membership rows in {time.time() - t0:.0f}s -> {len(out):,} og-member pairs")

    df = pd.DataFrame(out, columns=["og", "donor_ac"]).drop_duplicates()
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(MEMBERS, sep="\t", index=False)
    say(f"  wrote {MEMBERS.relative_to(REPO_ROOT)}")
    return df


# ------------------------------------------------------------------------ the scoring


def best_by_og(members: pd.DataFrame, lit: pd.DataFrame, drop_genus: str | None = None):
    """Most-cited prokaryotic SwissProt member per group. Same donor rule as transfer.py.

    `drop_genus` strikes a genus out of the MEMBER pool. **The validation is meaningless without
    it**: an E. coli protein's orthogroup contains that very protein's own SwissProt entry, so a
    first version of this control reported spearman 0.8351 of which **98.3% of the "donors" were
    Escherichia itself** -- self-correlation, not transfer. The sequence control removes
    Escherichia from the donor pool for exactly this reason; removing it here is the same rule
    applied to a different pool.
    """
    mem = members.merge(lit, on="donor_ac", how="inner")
    mem = mem[mem["prokaryotic"] & (mem["n_pubs"] > 0)]
    if drop_genus:
        mem = mem[~mem["donor_organism"].str.startswith(drop_genus)]
    return (mem.sort_values(["og", "n_pubs"], ascending=[True, False], kind="mergesort")
               .drop_duplicates("og", keep="first")
               .set_index("og"))


def score(tgt: pd.DataFrame, members: pd.DataFrame, lit: pd.DataFrame) -> pd.DataFrame:
    """Best-studied prokaryotic SwissProt member of the group. Narrow first, then domain."""
    best = best_by_og(members, lit)
    best_ctl = best_by_og(members, lit, drop_genus=ANCHOR_GENUS["ecoli"])
    say(f"  {len(best):,} of the target groups have a prokaryotic SwissProt member with papers "
        f"({len(best_ctl):,} once Escherichia is struck out, for the control track)")

    def pick(row) -> tuple:
        table = best_ctl if row["track"] == "control" else best
        for col, level in (("orthodb_og_narrow", "narrow"), ("orthodb_og_domain", "domain")):
            og = str(row[col]).strip()
            if og and og not in ("nan", "<NA>") and og in table.index:
                b = table.loc[og]
                return (int(b["n_pubs"]), b["donor_ac"], b["donor_gene"], b["donor_organism"],
                        og, level)
        return (0, "", "", "", "", "none")

    got = [pick(r) for _, r in tgt.iterrows()]
    out = tgt.copy()
    (out["n_papers_uniprot_orthodb"], out["orthodb_donor_ac"], out["orthodb_donor_gene"],
     out["orthodb_donor_organism"], out["orthodb_og_used"], out["orthodb_level_used"]) = zip(*got)
    return out


def validate(scored: pd.DataFrame) -> pd.DataFrame | None:
    """Held-out E. coli: does group-transferred literature predict measured own literature?"""
    v = scored[scored["track"] == "control"].copy()
    if len(v) < 50:
        return None
    own = S.load("ecoli")[["uniprot_ac", "n_papers_uniprot_own"]]
    v = v.merge(own, on="uniprot_ac", how="left").dropna(subset=["n_papers_uniprot_own"])
    hit = v[v["n_papers_uniprot_orthodb"] > 0]
    leak = hit["orthodb_donor_organism"].str.startswith(ANCHOR_GENUS["ecoli"]).sum()
    if leak:
        sys.exit(f"FATAL the control still has {int(leak):,} Escherichia donors -- the genus "
                 "exclusion is not being applied, and the correlation would be self-correlation.")
    rows = []
    for lbl, d in (("all", v), ("narrow", v[v.orthodb_level_used == "narrow"]),
                   ("domain", v[v.orthodb_level_used == "domain"])):
        if len(d) < 30:
            continue
        rows.append({"subset": lbl, "n": len(d),
                     "n_with_literature": int((d["n_papers_uniprot_orthodb"] > 0).sum()),
                     "spearman": round(float(d["n_papers_uniprot_orthodb"].corr(
                         d["n_papers_uniprot_own"], method="spearman")), 4)})
    say("")
    rule()
    say("CONTROL -- held-out E. coli proteins that LOSE their sequence donor")
    say("  these are the only place this route can be tested: no donor, but measured own papers")
    rule()
    say(f"  reached with literature: {len(hit):,} of {len(v):,} "
        f"({100 * len(hit) / max(len(v), 1):.1f}%)")
    for r in rows:
        say(f"    {r['subset']:<8} n={r['n']:>5,}  with literature {r['n_with_literature']:>5,}  "
            f"spearman {r['spearman']:>7.4f}")
    return pd.DataFrame(rows)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--refresh", action="store_true", help="re-stream the OrthoDB dumps")
    ap.add_argument("--dry-run", action="store_true", help="measure and print, write nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("studiedness/orthodb_transfer.py -- literature for the proteins DIAMOND cannot reach")
    say(f"   in   {ORTHODB_BASE}/{{{GENES},{OG2GENES}}}  (~9 GB, streamed, never stored)")
    say(f"   out  evidence/orthodb_transfer_<species>.tsv, evidence/orthodb_control.tsv")
    rule("=")

    lit = swissprot_literature()
    tgt, ogs = targets()
    members = og_members(ogs, set(lit["donor_ac"]), args.refresh)
    scored = score(tgt, members, lit)

    rule()
    say("REACH -- proteins that had NO literature of any kind, now given a group donor")
    rule()
    cols = ["uniprot_ac", "evidence", "orthodb_og_used", "orthodb_level_used",
            "orthodb_donor_ac", "orthodb_donor_gene", "orthodb_donor_organism",
            "n_papers_uniprot_orthodb"]
    for sp in args.species:
        d = scored[(scored["species"] == sp) & (scored["track"] == "no_donor")]
        if not len(d):
            continue
        got = int((d["n_papers_uniprot_orthodb"] > 0).sum())
        lv = d[d["n_papers_uniprot_orthodb"] > 0]["orthodb_level_used"].value_counts().to_dict()
        say(f"  {sp:<13} {got:>5,} of {len(d):>5,} reached "
            f"({100 * got / max(len(d), 1):>5.1f}%)   median {d.loc[d.n_papers_uniprot_orthodb > 0, 'n_papers_uniprot_orthodb'].median() if got else 0:>4.0f} papers"
            f"   by level: {lv}")
        if not args.dry_run:
            EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
            p = EVIDENCE_DIR / f"orthodb_transfer_{sp}.tsv"
            d[cols].to_csv(p, sep="\t", index=False)
            say(f"                -> {p.relative_to(REPO_ROOT)}")

    ctl = validate(scored)
    if ctl is not None and not args.dry_run:
        ctl.to_csv(EVIDENCE_DIR / "orthodb_control.tsv", sep="\t", index=False)
        say(f"  -> {(EVIDENCE_DIR / 'orthodb_control.tsv').relative_to(REPO_ROOT)}")
    rule("=")


if __name__ == "__main__":
    main()
