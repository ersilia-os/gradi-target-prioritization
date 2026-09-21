"""OGEE labels -> (protein sequence, label), for the 30-taxon prokaryotic corpus.

OGEE ships labels keyed on a locus identifier and NO sequences. This script supplies the sequences,
which is the entire gap between `ogee.py`'s label table and a trainable corpus. Output feeds an
`ogee_ess` column that sits BESIDE `geptop_ess` and `deg_ess` -- a third, independent opinion on
essentiality, never a replacement for either.

THE ASSEMBLY IS CHOSEN BY MEASUREMENT, NOT BY RULE. OGEE records a taxid, not an assembly, and the
annotation its authors keyed on may not be the current reference. Guessing is the failure mode:
`E. coli` ST131 (taxid 941322) has 37 complete assemblies and NONE flagged reference, while
`V. cholerae` C6706 (1124478) has ZERO complete assemblies. So for each taxon this script fetches
up to `--candidates` assemblies, builds each one's tag map, MEASURES the join rate against that
taxon's actual OGEE identifiers, and keeps the winner. The rate is written to the audit table for
every candidate tried, so a poor taxon is visibly poor rather than silently empty.

FOUR IDENTIFIER ROUTES, tried per taxon in order, because OGEE's `locus` namespace is a property of
the dataset block (see `ogee.py`):

  locus_tag       `PA14_00010`, `Rv0098`, `b4233`     -> GFF locus_tag, then old_locus_tag
  versioned tag   `HI0220.1`, `PA0195.1`              -> the same, with the `.N` suffix stripped
  ncbi_geneid     bare `199201`                       -> GFF Dbxref=GeneID:
  gene_symbol     `aaeA` (PEC's E. coli blocks only)  -> GFF gene / gene_synonym

`old_locus_tag` sits on the GENE feature while `protein_id` sits on the CDS, linked by ID/Parent --
walk both or the map comes out empty, which looks exactly like a taxon that does not join.

E. COLI K-12 MUST NOT BE DEDUPED ON `locus`. It carries two DISJOINT namespaces (PEC gene symbols
in datasets 108/109, b-numbers in 173/174/175, zero shared ids), so a naive dedup counts every gene
twice -- 9,496 entries for a 4,403-gene organism. Both namespaces are resolved separately and then
collapsed ON THE PROTEIN, which is what makes the count come out right.

Run with the `gradi` env:
    python scripts/essentiality/ogee_proteomes.py --dry-run
    python scripts/essentiality/ogee_proteomes.py --taxid 83333 --candidates 3
    python scripts/essentiality/ogee_proteomes.py
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from io import BytesIO
from pathlib import Path
from urllib.parse import unquote

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
TRAINING_DIR = OUT_DIR / "training_sets"
CACHE_DIR = REPO_ROOT / "data" / "source" / "ncbi" / "ogee_proteomes"

REPORT = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/taxon/{}/dataset_report"
          "?filters.assembly_level=complete_genome&filters.assembly_level=chromosome&page_size=40")
DOWNLOAD = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{}/download"
            "?include_annotation_type=PROT_FASTA&include_annotation_type=GENOME_GFF")

FLOOR = 0.80          # a taxon joining below this is REPORTED and excluded, never silently kept
MIN_TAXON_ROWS = 100  # below this a taxon cannot carry its own weight in the corpus

VERBOSE = True


def say(m: str = "") -> None:
    if VERBOSE:
        print(m, flush=True)


def rule(c: str = "-", w: int = 112) -> None:
    say(c * w)


def get_json(url: str, tries: int = 4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as fh:
                return json.load(fh)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if i == tries - 1:
                say(f"      ! {type(e).__name__} after {tries} tries")
                return None
            time.sleep(2 * (i + 1))
    return None


def candidate_assemblies(taxid: int, limit: int) -> list[str]:
    """Assembly accessions for a taxon, best-guess first.

    GENBANK BEFORE REFSEQ, deliberately: PGAP's RefSeq re-annotation drops submitter locus tags,
    and OGEE's identifiers predate the current annotation more often than not -- the same rule
    CLAUDE.md records for ECL8 and S. aureus COL. The ORDER is only a prior, though; the join rate
    decides.
    """
    js = get_json(REPORT.format(taxid))
    if not js or not js.get("reports"):
        return []
    scored = []
    for r in js["reports"]:
        acc = r.get("accession", "")
        info = r.get("assembly_info", {})
        cat = (info.get("refseq_category") or "").lower()
        rank = (0 if "reference" in cat else 1 if "representative" in cat else 2,
                0 if acc.startswith("GCA_") else 1,
                info.get("assembly_level", "") != "Complete Genome")
        scored.append((rank, acc))
    scored.sort()
    # GCA AND GCF OF THE SAME ASSEMBLY ARE BOTH KEPT. They are the same contigs but NOT the same
    # annotation, and the locus tags differ in exactly the way that decides this join: measured on
    # Synechococcus elongatus, OGEE keys on RefSeq `SYNPCC7942_RS00005` while the GenBank record
    # carries `Synpcc7942_0001`. An earlier version deduped the pair and scored that taxon 0.226.
    return [acc for _, acc in scored][:limit]


def fetch_assembly(acc: str) -> tuple[dict, dict] | None:
    """(tag_maps, protein_id -> sequence) for one assembly, cached as the raw zip.

    `tag_maps` is keyed by route: 'tag', 'geneid', 'symbol' -> {key: protein_id}.
    """
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    blob_path = CACHE_DIR / f"{acc}.zip"
    if blob_path.exists():
        blob = blob_path.read_bytes()
    else:
        try:
            with urllib.request.urlopen(DOWNLOAD.format(acc), timeout=600) as fh:
                blob = fh.read()
        except (urllib.error.URLError, TimeoutError) as e:
            say(f"      ! download failed for {acc}: {type(e).__name__}")
            return None
        blob_path.write_bytes(blob)
        time.sleep(1)
    try:
        z = zipfile.ZipFile(BytesIO(blob))
        gff = next(n for n in z.namelist() if n.endswith("genomic.gff"))
        faa = next(n for n in z.namelist() if n.endswith("protein.faa"))
    except (zipfile.BadZipFile, StopIteration):
        say(f"      ! {acc} has no GFF+protein.faa (annotation absent)")
        return None

    seqs, pid, buf = {}, None, []
    for line in z.read(faa).decode(errors="replace").splitlines():
        if line.startswith(">"):
            if pid:
                seqs[pid] = "".join(buf)
            pid, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if pid:
        seqs[pid] = "".join(buf)

    # old_locus_tag / gene / Dbxref live on the GENE feature; protein_id on the CDS. Walk both.
    genes: dict[str, dict] = {}
    maps = {"tag": {}, "geneid": {}, "symbol": {}, "norm": {}}
    for line in z.read(gff).decode(errors="replace").splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].rstrip().split(";") if "=" in kv)
        if f[2] in ("gene", "pseudogene") and a.get("ID"):
            genes[a["ID"]] = a
        elif f[2] == "CDS" and a.get("protein_id"):
            g = genes.get(a.get("Parent", ""), {})
            pid = a["protein_id"]
            if pid not in seqs:
                continue
            for src in (g, a):
                for key in ("locus_tag", "old_locus_tag"):
                    for t in unquote(src.get(key, "")).split(","):
                        if t.strip():
                            maps["tag"].setdefault(t.strip(), pid)
                            maps["norm"].setdefault(_norm(t), pid)
                for key in ("gene", "gene_synonym"):
                    for t in unquote(src.get(key, "")).split(","):
                        if t.strip():
                            maps["symbol"].setdefault(t.strip().lower(), pid)
                for x in unquote(src.get("Dbxref", "")).split(","):
                    if x.startswith("GeneID:"):
                        maps["geneid"].setdefault(x.split(":", 1)[1].strip(), pid)
    n_coll = len(maps["tag"]) - len({_norm(t) for t in maps["tag"]})
    if n_coll:
        say(f"      note {n_coll} locus tags collide once punctuation is stripped; exact matches "
            "still win, so this only affects the fallback route")
    return maps, seqs


def _norm(tag: str) -> str:
    """Locus tag with punctuation and case removed -- `HI0001` and `HI_0001` become one key.

    LOCUS TAG PUNCTUATION IS NOT SEMANTIC, and the mismatch is real: measured, OGEE keys
    H. influenzae on `HI0001` while the assembly writes `HI_0001`, H. pylori `HP0001` vs `HP_0001`,
    M. pneumoniae `MPN001` vs `MPN_001`. Three taxa and 4,199 labels were scoring 0.000 on a
    missing underscore.

    Used only AFTER exact matching, so it can never override a real hit, and collisions within a
    genome are counted and reported rather than silently resolved.
    """
    return re.sub(r"[^A-Za-z0-9]", "", str(tag)).lower()


def resolve(loci: pd.Series, maps: dict) -> pd.Series:
    """OGEE locus -> protein_id, trying every route. Reports nothing; the caller measures."""
    s = loci.astype(str).str.strip()
    out = s.map(maps["tag"])
    # versioned tags: HI0220.1 -> HI0220
    m = out.isna()
    if m.any():
        out[m] = s[m].str.replace(r"\.\d+$", "", regex=True).map(maps["tag"])
    m = out.isna()
    if m.any():
        out[m] = s[m].map(maps["geneid"])
    m = out.isna()
    if m.any():
        out[m] = s[m].str.lower().map(maps["symbol"])
    m = out.isna()
    if m.any():                    # last resort: punctuation-insensitive tag match
        out[m] = s[m].map(lambda x: maps["norm"].get(_norm(x)))
    return out


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--taxid", nargs="+", type=int, help="only these taxa")
    ap.add_argument("--candidates", type=int, default=3,
                    help="assemblies to try per taxon before keeping the best (default 3)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    lab_path, tax_path = EVIDENCE_DIR / "ogee_labels.tsv", EVIDENCE_DIR / "ogee_taxa.tsv"
    if not lab_path.exists():
        sys.exit(f"FATAL {lab_path} missing -- run scripts/essentiality/ogee.py first")
    lab = pd.read_csv(lab_path, sep="\t")
    tax = pd.read_csv(tax_path, sep="\t")

    # THE CORPUS DEFINITION, and it is the project's own filters, not OGEE's scope:
    # prokaryotes only (the eukaryotic taxa need per-database resolvers) and BOTH CLASSES present
    # (40 of 87 OGEE taxa are positives-only, 25 of them from the RB-TnSeq Fitness Browser, which
    # cannot see essential genes by construction).
    keep = tax[(tax.tier == "prokaryote_ncbi") & (tax.base_rate < 0.5)
               & (tax.n_entries >= MIN_TAXON_ROWS)]
    if a.taxid:
        keep = keep[keep.taxaID.isin(a.taxid)]
    names = dict(zip(tax.taxaID, tax.organism))

    rule("=")
    say("ogee_proteomes.py -- OGEE labels + protein sequences, for the ogee_ess corpus")
    rule("=")
    say(f"  taxa    {len(keep)} prokaryotic, both classes present, >= {MIN_TAXON_ROWS} entries")
    say(f"  labels  {int(keep.n_entries.sum()):,} entries, "
        f"{int(keep.n_essential.sum()):,} essential "
        f"(base {keep.n_essential.sum() / keep.n_entries.sum():.3f})")
    say(f"  cache   {CACHE_DIR.relative_to(REPO_ROOT)}/<accession>.zip")
    say(f"  floor   a taxon must join >= {FLOOR:.0%} of its labels or it is EXCLUDED and named")
    rule()
    if a.dry_run:
        say(f"  {'taxid':>8} {'organism':44s} {'entries':>8} {'base':>6}  namespace")
        for _, r in keep.iterrows():
            say(f"  {r.taxaID:>8} {str(r.organism)[:44]:44s} {r.n_entries:>8,} "
                f"{r.base_rate:>6.3f}  {r.id_kind}")
        rule("=")
        say("dry run: nothing fetched.")
        return

    say(f"  {'taxid':>8} {'organism':30s} {'assembly':18s} {'join':>7} {'n':>7} {'ess':>6}")
    rows, audit = [], []
    for _, r in keep.iterrows():
        t = int(r.taxaID)
        sub = lab[lab.taxaID == t]
        cands = candidate_assemblies(t, a.candidates)
        if not cands:
            audit.append({"taxaID": t, "organism": names.get(t), "accession": "",
                          "join_rate": 0.0, "n_labels": len(sub), "verdict": "no_assembly"})
            say(f"  {t:>8} {str(names.get(t))[:30]:30s} {'-':18s} {'0.000':>7}   NO ASSEMBLY")
            continue
        best = None
        for acc in cands:
            got = fetch_assembly(acc)
            if got is None:
                audit.append({"taxaID": t, "organism": names.get(t), "accession": acc,
                              "join_rate": 0.0, "n_labels": len(sub), "verdict": "fetch_failed"})
                continue
            maps, seqs = got
            pid = resolve(sub.locus, maps)
            rate = float(pid.notna().mean())
            audit.append({"taxaID": t, "organism": names.get(t), "accession": acc,
                          "join_rate": round(rate, 4), "n_labels": len(sub),
                          "verdict": "candidate"})
            if best is None or rate > best[0]:
                best = (rate, acc, pid, seqs)
            if rate >= 0.95:            # good enough; do not spend more downloads
                break
        if best is None:
            say(f"  {t:>8} {str(names.get(t))[:30]:30s} {'-':18s} {'0.000':>7}   ALL FETCHES FAILED")
            continue
        rate, acc, pid, seqs = best
        d = sub.assign(protein_id=pid.values)
        d = d[d.protein_id.notna()].copy()
        d["sequence"] = d.protein_id.map(seqs)
        d = d[d.sequence.notna() & d.sequence.str.len().gt(0)]
        # Collapse to one row per PROTEIN. This is what fixes E. coli's double count: its two
        # namespaces resolve to the same protein_id and merge here, never before.
        d = (d.sort_values("label", ascending=False)
              .drop_duplicates(["taxaID", "protein_id"], keep="first"))
        say(f"  {t:>8} {str(names.get(t))[:30]:30s} {acc:18s} {rate:>7.3f} {len(d):>7,} "
            f"{int(d.label.sum()):>6,}"
            + ("" if rate >= FLOOR else f"   <- BELOW {FLOOR:.0%} FLOOR, EXCLUDED"))
        for e in audit:
            if e["taxaID"] == t and e["accession"] == acc:
                e["verdict"] = "chosen" if rate >= FLOOR else "below_floor"
        if rate < FLOOR:
            continue
        d["organism"] = names.get(t)
        d["assembly"] = acc
        # UNIFORM TRAINING-SET SCHEMA: key, label, source_id, features_from -- the same shape as
        # the nine screen training sets, so everything in training_sets/ is interchangeable.
        # `key` is the assembly's protein_id and `features_from` is the taxid, because this
        # corpus's embeddings are computed PER TAXON (one FASTA and one ESM-C matrix each).
        d = d.rename(columns={"protein_id": "key", "locus": "source_id"})
        d["features_from"] = d["taxaID"].astype(str)
        rows.append(d[["key", "label", "source_id", "features_from", "taxaID", "organism",
                       "assembly", "id_kind", "sequence", "n_datasets", "n_essential",
                       "disagrees"]])

    if not rows:
        sys.exit("FATAL no taxon cleared the floor")
    out = pd.concat(rows, ignore_index=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    TRAINING_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(TRAINING_DIR / "ogee_corpus.tsv", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(EVIDENCE_DIR / "ogee_proteome_join.tsv", sep="\t", index=False)

    rule()
    say("CORPUS")
    rule()
    say(f"  {len(out):,} labelled proteins over {out.taxaID.nunique()} taxa, "
        f"{int(out.label.sum()):,} essential (base {out.label.mean():.3f})")
    say(f"  {out.sequence.nunique():,} distinct sequences  "
        f"({len(out) - out.sequence.nunique():,} shared across taxa -- expected, they are orthologs)")
    # A taxon is EXCLUDED only if no candidate was chosen. Reporting every taxon that had ANY
    # failed candidate listed E. coli K-12 as excluded when it had in fact joined at 0.921 on its
    # first candidate -- a reporting bug that looked like a data loss.
    chosen = {e["taxaID"] for e in audit if e["verdict"] == "chosen"}
    excl = sorted({(e["taxaID"], e["organism"], e["verdict"]) for e in audit
                   if e["taxaID"] not in chosen
                   and e["verdict"] in ("below_floor", "no_assembly", "fetch_failed")})
    if excl:
        say(f"  EXCLUDED: {len({e[0] for e in excl})} taxa")
        for _, org, why in excl:
            say(f"    {str(org)[:52]:52s} {why}")
    say(f"  wrote {(TRAINING_DIR / 'ogee_corpus.tsv').relative_to(REPO_ROOT)}")
    say(f"  wrote {(EVIDENCE_DIR / 'ogee_proteome_join.tsv').relative_to(REPO_ROOT)}")
    rule("=")


if __name__ == "__main__":
    main()
