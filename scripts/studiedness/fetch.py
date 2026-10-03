"""Download everything the studiedness axis reads, and verify every payload.

This axis asks **how much is already known about this protein?** -- needed in both directions,
because an uncharacterised target is a risk but also the novelty the collaboration is looking for.
Nothing here computes a score; it fetches the five inputs and writes a `SOURCE.md` beside each so
the run is reproducible after the bulk files are deleted.

    data/source/uniprot/literature/anchor_<species>.tsv.gz   annotation_score, PubMed ids, PE
    data/source/uniprot/literature/swissprot_meta.tsv.gz     the same for all reviewed entries
    data/source/uniprot/literature/uniprot_sprot.fasta.gz    the DIAMOND donor database
    data/source/ncbi/gene/gene2pubmed.gz                     NCBI's GeneID -> PMID links
    data/source/unknome/unknome_clusters.tsv                 PANTHER family knownness

WHY THE ANCHORS' OWN LITERATURE CANNOT CARRY THIS AXIS
------------------------------------------------------
Measured on 2026-09-21, which is the whole reason the axis is built on homology transfer:

    Kp HS11286   0.1% reviewed    5,710 of 5,728 have EXACTLY 1 PubMed id (the genome paper)
    Sa NCTC 8325 28.2% reviewed   2,532 of 2,889 have NONE; median 0
    Ec K-12      100% reviewed    median 5, max 58, 48 distinct values

99.7% of K. pneumoniae carries the identical count, so the column has no discriminative power at
all. This confirms and sharpens v1's finding (`legacy/HISTORY.md:60`). The number must come from
the best-characterised homolog; see `transfer.py`.

WHY SWISSPROT IS THE DONOR DATABASE
-----------------------------------
The free four-species ortholog table on disk reaches literature (an ortholog in Ec or human) for
only **Kp 68.5% and Sa 46.6%** -- *S. aureus* is Gram-positive, so less than half of it has an
E. coli ortholog at all. SwissProt is 575,748 reviewed entries over ~14,000 species and holds the
Gram-positive literature our four-species panel structurally cannot contain. Curated literature is
also, by definition, where reviewed entries are: an unreviewed entry has no curated references.

EVERY PAYLOAD IS VERIFIED, BECAUSE HTTP 200 IS NOT EVIDENCE OF DATA
-------------------------------------------------------------------
UniProt streams are checked against `x-total-results` from the matching search query -- an
authoritative count from the same server, not a guess. Plain files are checked against
`Content-Length`. A short read exits non-zero rather than writing a truncated file that would
silently shrink every downstream join.

Run with the `gradi` env. Nothing here needs DIAMOND.
  python scripts/studiedness/fetch.py
  python scripts/studiedness/fetch.py --only swissprot --refresh
  python scripts/studiedness/fetch.py --dry-run
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
SPROT_FASTA_URL = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
                   "knowledgebase/complete/uniprot_sprot.fasta.gz")
GENE2PUBMED_URL = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz"
UNKNOME_CLUSTER_URL = "https://unknome.mrc-lmb.cam.ac.uk/download/clust_tsv"

LIT_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "literature"
NCBI_GENE_DIR = REPO_ROOT / "data" / "source" / "ncbi" / "gene"
UNKNOME_DIR = REPO_ROOT / "data" / "source" / "unknome"

# The three bacteria. Human is out of scope for a bacterial target-prioritization axis.
SPECIES = ("kpneumoniae", "ecoli", "saureus")
PROTEOME_ID = {"kpneumoniae": "UP000007841", "ecoli": "UP000000625", "saureus": "UP000008816"}

# `xref_geneid` is the bridge to gene2pubmed and rides along free in the same call.
ANCHOR_FIELDS = ("accession", "reviewed", "annotation_score", "protein_existence",
                 "lit_pubmed_id", "xref_geneid")
SWISSPROT_FIELDS = ("accession", "reviewed", "annotation_score", "protein_existence",
                    "lit_pubmed_id", "xref_geneid", "organism_id", "organism_name",
                    "gene_primary", "protein_name", "lineage")
SWISSPROT_QUERY = "reviewed:true"

USER_AGENT = "gradi/2.0 (studiedness axis)"
TARGETS = ("anchors", "swissprot", "swissprot_fasta", "gene2pubmed", "unknome")

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def human(n: float) -> str:
    """Bytes as a short human string, so a 287 MB file does not print as 287114207."""
    for unit in ("B", "KB", "MB", "GB"):
        if abs(n) < 1024 or unit == "GB":
            return f"{n:,.1f} {unit}" if unit != "B" else f"{int(n)} B"
        n /= 1024
    return f"{n:.1f} GB"


def md5_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(url: str, dest: Path, refresh: bool, label: str) -> tuple[Path, int]:
    """Stream a URL to disk, asserting the byte count against Content-Length.

    Writes to a `.tmp` and renames only on success, so an interrupted run never leaves a
    half-file that a later run would treat as cached.
    """
    if dest.exists() and dest.stat().st_size > 0 and not refresh:
        size = dest.stat().st_size
        say(f"  {label:<22} cached    {human(size):>12}   {dest.relative_to(REPO_ROOT)}")
        return dest, size

    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    t0 = time.time()
    try:
        with urllib.request.urlopen(req, timeout=1800) as resp:
            declared = resp.headers.get("Content-Length")
            declared = int(declared) if declared else None
            got = 0
            with open(tmp, "wb") as fh:
                while True:
                    block = resp.read(1 << 20)
                    if not block:
                        break
                    fh.write(block)
                    got += len(block)
                    if declared and got % (64 << 20) < (1 << 20):
                        say(f"    {label}: {human(got)} / {human(declared)} "
                            f"({100 * got / declared:.0f}%)")
    except urllib.error.URLError as exc:
        tmp.unlink(missing_ok=True)
        sys.exit(f"FATAL {label}: download failed -- {exc}\n  url: {url}")

    # An HTTP 200 is not evidence of data. v1 was repeatedly bitten by short and empty payloads.
    if got == 0:
        tmp.unlink(missing_ok=True)
        sys.exit(f"FATAL {label}: server returned 0 bytes with no error.\n  url: {url}")
    if declared is not None and got != declared:
        tmp.unlink(missing_ok=True)
        sys.exit(f"FATAL {label}: short read -- got {got:,} bytes, Content-Length said "
                 f"{declared:,}.\n  url: {url}")
    tmp.replace(dest)
    say(f"  {label:<22} fetched   {human(got):>12}   {time.time() - t0:.0f}s   "
        f"{dest.relative_to(REPO_ROOT)}")
    return dest, got


def uniprot_total(query: str) -> int:
    """The authoritative row count for a UniProt query, from the server's own header."""
    url = f"{UNIPROT_SEARCH}?query={quote(query, safe=':*')}&format=tsv&size=1&fields=accession"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=300) as resp:
        total = resp.headers.get("x-total-results")
        release = resp.headers.get("x-uniprot-release")
        date = resp.headers.get("x-uniprot-release-date")
    if total is None:
        sys.exit(f"FATAL UniProt returned no x-total-results for {query!r}; cannot verify the "
                 "stream, and an unverified payload is not usable.")
    return int(total), (release or "?"), (date or "?")


def count_rows(path: Path) -> int:
    """Data rows in a gzipped TSV (header excluded)."""
    n = 0
    with gzip.open(path, "rt", errors="replace") as fh:
        for _ in fh:
            n += 1
    return max(n - 1, 0)


def fetch_uniprot_stream(query: str, fields: tuple[str, ...], dest: Path, refresh: bool,
                         label: str) -> dict:
    """Stream a UniProt TSV and assert its row count against the search endpoint's total."""
    expected, release, release_date = uniprot_total(query)
    url = (f"{UNIPROT_STREAM}?compressed=true&format=tsv"
           f"&query={quote(query, safe=':*')}&fields={','.join(fields)}")
    path, size = download(url, dest, refresh, label)
    observed = count_rows(path)
    if observed != expected:
        sys.exit(f"FATAL {label}: stream returned {observed:,} rows but UniProt's own search "
                 f"reports {expected:,} for {query!r}. Refusing a partial table -- re-run with "
                 "--refresh.")
    say(f"    verified {observed:,} rows against x-total-results (UniProt {release}, {release_date})")
    return {"url": url, "query": query, "rows": observed, "bytes": size,
            "release": release, "release_date": release_date, "md5": md5_of(path)}


def write_source_md(path: Path, title: str, body: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
    text = f"# {title}\n\n- fetched: {stamp}\n" + "\n".join(body) + "\n"
    path.write_text(text)
    say(f"  wrote {path.relative_to(REPO_ROOT)}")


def do_anchors(refresh: bool) -> list[dict]:
    rule()
    say("ANCHOR LITERATURE - annotation score, protein existence, PubMed ids")
    rule()
    meta = []
    for sp in SPECIES:
        info = fetch_uniprot_stream(
            f"proteome:{PROTEOME_ID[sp]}", ANCHOR_FIELDS,
            LIT_DIR / f"anchor_{sp}.tsv.gz", refresh, sp)
        expected = len(P.load(sp))
        if info["rows"] != expected:
            sys.exit(f"FATAL {sp}: UniProt returned {info['rows']:,} rows but "
                     f"proteome_{sp}.tsv has {expected:,}. The proteome moved under us; "
                     "re-run scripts/proteomes/download.py before this stage.")
        info["species"] = sp
        meta.append(info)
    write_source_md(
        LIT_DIR / "anchors.SOURCE.md", "uniprot literature -- anchor proteomes",
        [f"- fields: {','.join(ANCHOR_FIELDS)}", ""] +
        [f"- {m['species']}: {m['rows']:,} rows, {m['bytes']:,} bytes, md5 {m['md5']}, "
         f"UniProt {m['release']} ({m['release_date']})" for m in meta] +
        ["", "## URLs", ""] + [f"- {m['url']}" for m in meta] +
        ["", "## Why", "",
         "`proteomes/download.py` fetches neither annotation_score, protein_existence nor",
         "lit_pubmed_id, so this axis fetches them itself. Measured here: Kp is 0.1% reviewed",
         "with 5,710 of 5,728 proteins carrying exactly one PubMed id (the genome paper), so the",
         "column is kept as `n_papers_own` -- a measurement of darkness -- and the usable",
         "number is transferred from homologs by `transfer.py`."])
    return meta


def do_swissprot(refresh: bool) -> dict:
    rule()
    say("SWISSPROT - the donor database (metadata + sequences)")
    rule()
    meta = fetch_uniprot_stream(SWISSPROT_QUERY, SWISSPROT_FIELDS,
                                LIT_DIR / "swissprot_meta.tsv.gz", refresh, "swissprot meta")
    faa, faa_bytes = download(SPROT_FASTA_URL, LIT_DIR / "uniprot_sprot.fasta.gz",
                              refresh, "swissprot fasta")
    n_seq = 0
    with gzip.open(faa, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                n_seq += 1
    if n_seq != meta["rows"]:
        say(f"  NOTE fasta has {n_seq:,} sequences against {meta['rows']:,} metadata rows "
            "(isoforms and release skew); the transfer joins on accession, so this is reported "
            "rather than fatal.")
    else:
        say(f"    verified {n_seq:,} sequences, matching the metadata row count")
    write_source_md(
        LIT_DIR / "swissprot.SOURCE.md", "uniprot swissprot -- studiedness donor database",
        [f"- query: {SWISSPROT_QUERY}",
         f"- fields: {','.join(SWISSPROT_FIELDS)}",
         f"- metadata: {meta['rows']:,} rows, {meta['bytes']:,} bytes, md5 {meta['md5']}",
         f"- fasta: {n_seq:,} sequences, {faa_bytes:,} bytes, md5 {md5_of(faa)}",
         f"- UniProt release: {meta['release']} ({meta['release_date']})",
         "", "## URLs", "",
         f"- {meta['url']}", f"- {SPROT_FASTA_URL}",
         "", "## Why `lineage`", "",
         "The donor scope is restricted to true Bacteria by testing for 'Bacteria (domain)' in",
         "the taxonomic lineage -- the same rule `ligands/chembl.py` applies via ChEMBL's",
         "organism_class.l1. An organism NAME or a single tax_id cannot do this: strains are",
         "filed under their own taxids, and 'not human' is not the same as 'bacterial' (the",
         "ligands axis measured that mistake costing 424 apparent hits against a true 175).",
         "", "## Why reviewed-only", "",
         "Curated literature is by definition where reviewed entries are -- an unreviewed entry",
         "has no curated reference list. SwissProt spans ~14,000 species, which is what lets it",
         "reach the Gram-positive literature the project's own four-species ortholog panel",
         "cannot (that panel reaches Ec-or-human orthologs for only 46.6% of S. aureus).",
         "", "Deletable after `transfer.py` has run: every derived artifact is kept, and only",
         "re-running the DIAMOND search needs the fasta back."])
    meta["n_seq"] = n_seq
    return meta


def do_gene2pubmed(refresh: bool) -> dict:
    rule()
    say("NCBI gene2pubmed - GeneID -> PMID")
    rule()
    path, size = download(GENE2PUBMED_URL, NCBI_GENE_DIR / "gene2pubmed.gz", refresh, "gene2pubmed")
    with gzip.open(path, "rt", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
    if not header.startswith("#tax_id"):
        sys.exit(f"FATAL gene2pubmed: unexpected header {header!r}; expected '#tax_id\\tGeneID"
                 "\\tPubMed_ID'. The format changed -- check before parsing.")
    say(f"    header verified: {header}")
    write_source_md(
        NCBI_GENE_DIR / "gene2pubmed.SOURCE.md", "ncbi gene2pubmed",
        [f"- bytes: {size:,}", f"- md5: {md5_of(path)}",
         f"- columns: {header}",
         "- 83,245,982 rows over 43,745,777 distinct GeneIDs (measured 2026-09-21)",
         "", "## URLs", "", f"- {GENE2PUBMED_URL}",
         "", "## Why", "",
         "Joins on the `geneid` column already in proteome_<species>.tsv (Kp 100% / Ec 95.1% /",
         "Sa 94.7%), so no id-mapping round-trip is needed. Measured: it does NOT rescue the",
         "anchors (Kp 12 distinct values, Sa 14.8% coverage -- both still dark) but it resolves",
         "a well-studied organism ~4x more finely than UniProt's curated list (Ec 193 distinct",
         "values, max 513, against UniProt's 48 and max 58). That sharpens the DONOR ranking,",
         "but it is NOT the shipped count: n_papers_family counts UniProt curated refs only,",
         "one consistent definition per row. gene2pubmed is a measured alternative.",
         "", "Public and re-derivable -- never upload it to eosvc. Deletable after",
         "`gene2pubmed.py` has reduced it to the GeneIDs this project uses."])
    return {"bytes": size, "path": path}


def do_unknome(refresh: bool) -> dict:
    rule()
    say("UNKNOME - PANTHER family knownness")
    rule()
    path, size = download(UNKNOME_CLUSTER_URL, UNKNOME_DIR / "unknome_clusters.tsv",
                          refresh, "unknome clusters")
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        n = sum(1 for _ in fh)
    need = {"cluster_id", "panther_id", "knownness", "num_proteins", "num_species"}
    if not need.issubset(header):
        sys.exit(f"FATAL unknome: columns {sorted(need - set(header))} missing from {header}")
    say(f"    verified {n:,} clusters, {len(header)} columns")
    write_source_md(
        UNKNOME_DIR / "unknome.SOURCE.md", "unknome v3 -- cluster knownness",
        [f"- bytes: {size:,}", f"- md5: {md5_of(path)}", f"- clusters: {n:,}",
         f"- columns: {', '.join(header)}",
         "- citation: Rocha, Jayaram, Stevens et al., PLoS Biol 21(8):e3002222 (2023)",
         "", "## URLs", "", f"- {UNKNOME_CLUSTER_URL}",
         "",
         "## THE CLUSTER TABLE, NEVER THE PER-PROTEIN TABLE",
         "",
         "The sibling download (`/download/prot_tsv`, 276 MB) reads `knownness = 0.000` for ALL",
         "1,882 S. aureus entries -- including clpP, rpoB, ftsZ, gyrB, dnaA -- while those same",
         "proteins' clusters carry real values (Sa clpP is in UKP00027, the same cluster as Ec",
         "clpP, knownness 10.6). Joining the protein table would silently zero a whole proteome.",
         "Measured 2026-09-21. This file is the cluster table and is the one to use.",
         "",
         "## Species scope",
         "",
         "143 species. E. coli (83333) and S. aureus (93061) are both present and are exactly",
         "our anchor taxids, but there is NO Klebsiella -- Kp is reachable only through the",
         "PANTHER family id, which is how every species is joined here for consistency."])
    return {"bytes": size, "clusters": n}


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", nargs="+", default=list(TARGETS), choices=list(TARGETS),
                    help="fetch only these targets (default: all)")
    ap.add_argument("--refresh", action="store_true", help="re-download even if cached")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STUDIEDNESS - fetch")
    rule("=")
    say("  inputs   UniProt anchor literature · SwissProt metadata + fasta · NCBI gene2pubmed ·")
    say("           Unknome cluster knownness")
    say(f"  outputs  {LIT_DIR.relative_to(REPO_ROOT)}")
    say(f"           {NCBI_GENE_DIR.relative_to(REPO_ROOT)}")
    say(f"           {UNKNOME_DIR.relative_to(REPO_ROOT)}")
    say("  every payload is verified against an authoritative count; HTTP 200 is not evidence")
    say(f"  targets  {', '.join(args.only)}")
    rule("=")
    if args.dry_run:
        say("dry run: nothing fetched.")
        say("  approximate transfer: swissprot meta ~60 MB · fasta 94 MB · gene2pubmed 287 MB ·")
        say("                        unknome 8 MB · anchors < 5 MB")
        return

    free = shutil.disk_usage(REPO_ROOT).free / 1e9
    if free < 5:
        sys.exit(f"FATAL only {free:.1f} GB free; this stage needs ~1 GB plus room for DIAMOND.")

    t0 = time.time()
    if "anchors" in args.only:
        do_anchors(args.refresh)
    if "swissprot" in args.only or "swissprot_fasta" in args.only:
        do_swissprot(args.refresh)
    if "gene2pubmed" in args.only:
        do_gene2pubmed(args.refresh)
    if "unknome" in args.only:
        do_unknome(args.refresh)

    rule("=")
    say(f"fetch complete in {time.time() - t0:.0f}s.")
    say("  next: scripts/studiedness/gene2pubmed.py, then unknome.py, then transfer.py")
    rule("=")


if __name__ == "__main__":
    main()
