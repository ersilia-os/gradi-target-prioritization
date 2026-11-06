"""Unknome family knownness, joined through the PANTHER family id.

Unknome (Rocha, Jayaram, Stevens et al., PLoS Biol 21(8):e3002222, 2023) scores each orthologous
cluster by **how much is annotated about its best-annotated member** -- the largest weighted count
of GO terms over the cluster. That is a *family* property, which is exactly the right shape for
this project: our anchors are bibliometrically dark as individual accessions but their families
are often very well known.

    data/processed/studiedness/evidence/unknome_<species>.tsv

    uniprot_ac · panther_families · unknome_cluster_id · unknome_cluster_name · unknome_knownness
    · unknome_num_species · unknome_best_known_protein · unknome_evidence

USE THE CLUSTER TABLE. NEVER THE PER-PROTEIN TABLE.
---------------------------------------------------
Unknome publishes both. The per-protein download (`/download/prot_tsv`, 276 MB) reads
**`knownness = 0.000` for ALL 1,882 S. aureus entries** -- including `clpP`, `rpoB`, `ftsZ`,
`gyrB`, `dnaA`, `rplB` -- while those same proteins' clusters carry real values. Sa `clpP` sits in
`UKP00027`, the *same cluster* as Ec `clpP`, which reads 10.6. Measured 2026-09-21. Joining the
protein table would have silently zeroed a whole proteome while looking like a completed run.
So the join is: our `panther` xref -> the cluster table's `panther_id` -> its knownness.

WHY PANTHER AND NOT THE ACCESSION
----------------------------------
Unknome covers 143 species. E. coli (83333) and S. aureus (93061) are both there and are exactly
our anchor taxids -- but **there is no Klebsiella**. Going through the PANTHER family id reaches
all three species by the same route, and `panther` is already on disk from stage 00
(`data/processed/proteomes/evidence/annotation_<species>.tsv`), so nothing is downloaded but an
8.4 MB table.

COVERAGE IS CAPPED BY THE PANTHER XREF, AND THAT IS WHY THIS IS A COLUMN AND NOT THE AXIS
------------------------------------------------------------------------------------------
Measured: **Kp 71.2% · Ec 78.2% · Sa 65.1%** of proteins carry a PANTHER family at all, and
essentially every one that does joins (4,079 of 4,081 on Kp). The remaining ~29% of Kp has no
family to look up. The deliverable needs a score for **every** protein, so SwissProt homology
(`transfer.py`) carries the axis and this is an independent second opinion.

TWO CAVEATS TO READ BEFORE QUOTING A NUMBER
--------------------------------------------
1. Knownness counts **GO terms, not papers**. It measures annotation depth. That is what makes it
   genuinely independent of the literature route -- and also means the two are not interchangeable.
2. The clusters span 143 mostly-eukaryote species, so a bacterial protein's knownness can be set
   by a plant homolog: FtsZ's `best_known_protein` is `FTSZ1_ARATH`. Read
   `unknome_best_known_protein` before treating the number as bacterial evidence.

A protein with no PANTHER family gets an **empty row with `unknome_evidence = none`**, never a
zero -- a zero here would read as "nothing is known about this family", which is a different and
much stronger claim than "this protein has no PANTHER family".

Run with the `gradi` env, after `scripts/studiedness/fetch.py`.
  python scripts/studiedness/unknome.py
  python scripts/studiedness/unknome.py --species kpneumoniae -q
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
CLUSTERS = REPO_ROOT / "data" / "source" / "unknome" / "unknome_clusters.tsv"

SPECIES = ("kpneumoniae", "ecoli", "saureus")
PTHR = re.compile(r"PTHR\d+")

# Families whose knownness must be high in every species -- a join that silently produced the
# wrong ids would still look plausible, so the control is named genes, not a row count.
SPOT_GENES = ("rpoB", "gyrB", "ftsZ", "clpP", "dnaA", "secA", "groEL", "rplB")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def load_clusters() -> pd.DataFrame:
    """The Unknome cluster table, collapsed to one row per PANTHER family.

    A PANTHER family can carry more than one Unknome cluster (the clusters are finer). The family's
    knownness is the **max** over its clusters: knownness already means "the best-annotated member",
    so taking the best cluster is the same rule applied one level up.
    """
    if not CLUSTERS.exists():
        sys.exit(f"FATAL missing {CLUSTERS.relative_to(REPO_ROOT)} -- "
                 "run scripts/studiedness/fetch.py first")
    c = pd.read_csv(CLUSTERS, sep="\t", dtype=str, keep_default_na=False)
    c["knownness"] = pd.to_numeric(c["knownness"], errors="coerce")
    if c["knownness"].isna().any():
        sys.exit(f"FATAL {int(c['knownness'].isna().sum())} clusters have a non-numeric knownness")
    c["num_species"] = pd.to_numeric(c["num_species"], errors="coerce").fillna(0).astype(int)
    c = c.sort_values("knownness", ascending=False).drop_duplicates("panther_id", keep="first")
    say(f"  {len(c):,} PANTHER families   knownness median {c.knownness.median():.1f}  "
        f"max {c.knownness.max():.1f}  zero for {int((c.knownness == 0).sum()):,}")
    return c.set_index("panther_id")


def families(cell: object) -> list[str]:
    """PANTHER family ids from a UniProt xref cell, subfamilies (`:SF12`) folded into the family.

    The cell looks like `PTHR11066;PTHR11066:SF34;` -- the same family twice. Unknome keys on the
    family, so the subfamily suffix is dropped and the list de-duplicated in place.
    """
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    seen, out = set(), []
    for fam in PTHR.findall(str(cell)):
        if fam not in seen:
            seen.add(fam)
            out.append(fam)
    return out


def build(species: str, clusters: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """One row per protein, in the proteome's own order."""
    prot = P.load(species)[["uniprot_ac", "gene_name"]]
    ann = P.load_annotation(species)[["uniprot_ac", "panther"]]
    df = prot.merge(ann, on="uniprot_ac", how="left")
    if len(df) != len(prot):
        sys.exit(f"FATAL {species}: the annotation join changed the row count "
                 f"({len(prot):,} -> {len(df):,})")

    recs = []
    for cell in df["panther"]:
        fams = families(cell)
        hits = [(f, clusters.loc[f]) for f in fams if f in clusters.index]
        if not hits:
            recs.append(("", "", "", np.nan, 0, "", "no_family" if not fams else "no_cluster"))
            continue
        fam, best = max(hits, key=lambda kv: kv[1]["knownness"])
        recs.append((";".join(fams), best["cluster_id"], best["cluster_name"],
                     float(best["knownness"]), int(best["num_species"]),
                     best["best_known_protein_id"], "panther_cluster"))

    out = pd.DataFrame(recs, columns=[
        "panther_families", "unknome_cluster_id", "unknome_cluster_name", "unknome_knownness",
        "unknome_num_species", "unknome_best_known_protein", "unknome_evidence"])
    out.insert(0, "uniprot_ac", df["uniprot_ac"].to_numpy())

    has_family = int((df["panther"].fillna("").astype(str).str.contains("PTHR")).sum())
    joined = int(out["unknome_evidence"].eq("panther_cluster").sum())
    stats = {
        "species": species, "n": len(out),
        "with_panther": has_family,
        "joined": joined,
        "joined_pct": round(100 * joined / len(out), 1),
        "of_those_with_panther_pct": round(100 * joined / max(has_family, 1), 1),
        "median_knownness": round(float(out["unknome_knownness"].median()), 2)
        if joined else float("nan"),
        "knownness_zero": int((out["unknome_knownness"] == 0).sum()),
        "no_family": int(out["unknome_evidence"].eq("no_family").sum()),
        "no_cluster": int(out["unknome_evidence"].eq("no_cluster").sum()),
    }
    return out.merge(prot[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left"), stats


def spot_check(species: str, out: pd.DataFrame) -> list[str]:
    """Named genes whose families must be well known. A wrong join still looks plausible."""
    lines, bad = [], []
    scored = out[out["unknome_knownness"].notna()]
    if scored.empty:
        return ["      (nothing joined; no spot check possible)"]
    for gene in SPOT_GENES:
        hit = out[out["gene_name"].fillna("").str.lower() == gene.lower()]
        if hit.empty:
            lines.append(f"      {gene:<7} not named in this proteome")
            continue
        row = hit.iloc[0]
        k = row["unknome_knownness"]
        if pd.isna(k):
            lines.append(f"      {gene:<7} no PANTHER family")
            continue
        pct = 100 * (scored["unknome_knownness"] < k).mean()
        flag = "" if pct >= 50 else "   <-- LOW"
        if pct < 50:
            bad.append(gene)
        lines.append(f"      {gene:<7} knownness {k:6.1f}  {pct:5.1f}th pct  "
                     f"{row['unknome_cluster_id']} {row['unknome_best_known_protein']}{flag}")
    if bad:
        lines.append(f"      NOTE below median: {', '.join(bad)} -- knownness is a GO-term count "
                     "over mostly-eukaryote families, so a bacterial workhorse with a small "
                     "eukaryote family can sit low. Reported, not fatal.")
    return lines


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STUDIEDNESS - Unknome family knownness")
    rule("=")
    say(f"  in    {CLUSTERS.relative_to(REPO_ROOT)}  (cluster table -- NEVER the protein table)")
    say("        data/processed/proteomes/evidence/annotation_<species>.tsv  (panther xref)")
    say(f"  out   {(EVIDENCE_DIR / 'unknome_<species>.tsv').relative_to(REPO_ROOT)}")
    say("  knownness counts GO TERMS, not papers -- an independent second opinion")
    say("  no PANTHER family gives an EMPTY row, never a zero")
    rule("=")
    if args.dry_run:
        say("dry run: nothing written.")
        return

    clusters = load_clusters()
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for sp in args.species:
        rule()
        say(f"{sp}")
        rule()
        out, stats = build(sp, clusters)
        path = EVIDENCE_DIR / f"unknome_{sp}.tsv"
        out.drop(columns=["gene_name"]).to_csv(path, sep="\t", index=False)
        rows.append(stats)
        say(f"    {len(out):,} proteins   with a PANTHER family {stats['with_panther']:,} "
            f"({100 * stats['with_panther'] / stats['n']:.1f}%)")
        say(f"    joined to Unknome {stats['joined']:,} ({stats['joined_pct']}%)   "
            f"= {stats['of_those_with_panther_pct']}% of those with a family")
        say(f"    median knownness {stats['median_knownness']}   "
            f"knownness exactly 0 for {stats['knownness_zero']:,}   "
            f"no family {stats['no_family']:,}   family but no cluster {stats['no_cluster']:,}")
        for line in spot_check(sp, out):
            say(line)
        say(f"    -> {path.relative_to(REPO_ROOT)}")

    man = pd.DataFrame(rows)
    man.to_csv(EVIDENCE_DIR / "unknome_manifest.tsv", sep="\t", index=False)
    rule("=")
    say(f"  wrote {(EVIDENCE_DIR / 'unknome_manifest.tsv').relative_to(REPO_ROOT)}")
    say("  coverage is capped by the PANTHER xref, so this cannot be the axis on its own --")
    say("  transfer.py carries it. merge.py reports whether knownness adds anything.")
    rule("=")


if __name__ == "__main__":
    main()
