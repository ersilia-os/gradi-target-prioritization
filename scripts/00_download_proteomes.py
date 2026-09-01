"""Fetch the project's reference proteomes and build one maximally-identified table per species.

Stage 00. Four species, **one reference proteome each** -- that proteome is the unit of analysis for
every downstream stage:

  HS11286    UP000007841  K. pneumoniae   5,728  the anchor (KPHS_*)
  K-12       UP000000625  E. coli         4,403  co-equal second organism; 100% reviewed (b-numbers)
  NCTC 8325  UP000008816  S. aureus       2,889  the ClpP-activator organism (SAOUHSC_*)
  Human      UP000005640  H. sapiens     20,416  reviewed canonical, for selectivity

The problem this stage exists to solve
--------------------------------------
The anchors are badly under-named. Measured against UniProt 2026_02:

                       gene name     synonyms    locus tag
  Kp HS11286            18.4%          0.6%        100%
  Ec K-12              100.0%         45.5%        100%
  Sa NCTC 8325          28.2%          2.9%         98.4%

Coverage tracks curation almost exactly (Sa: 815 reviewed, 815 named). The *mappings* are fine --
RefSeq/GeneID are ~100% even on Kp -- it is specifically names and synonyms that are missing, and
those are what literature and database lookups key on. A gene-symbol join on the raw anchor loses
~27% of known essentials (see legacy/HISTORY.md trap 59).

NCBI is not the answer: the HS11286 RefSeq GFF sets `Name=KPHS_00010` for all 5,867 genes -- zero
real symbols (Sa yields only 219).

What is the answer is the rest of the species. The Kp species tree carries 26,956 entries WITH a
gene name; Sa carries 17,562. So names are filled in three ordered, labelled tiers, all of which
stay inside the species -- no cross-species inference:

  anchor            UniProt's own name on the anchor entry
  species_exact     identical sequence, same species          Kp 18.4 -> 40.4%, Sa 28.2 -> 39.5%
  species_uniref90  same UniRef90 cluster, same species        (>=90% identity, no alignment needed)
  none

Cross-species naming (inheriting an E. coli ortholog's name) is deliberately NOT done here: it is an
inference, it belongs downstream of real orthology, and it would drag a DIAMOND dependency into
stage 00. This script needs nothing but `requests` and `pandas`.

Ambiguity is never resolved silently. Where donors disagree the choice is the donor consensus, ties
abstain, and every candidate is kept in `gene_name_candidates` plus a row in `name_audit.tsv`.

Outputs
-------
  data/raw/00_proteomes/uniprot/<label>.{fasta,tsv}   exactly as fetched, plus SOURCE.md
  data/raw/00_proteomes/ncbi/<label>.{faa,gff}        tier D bridge strains
  data/processed/00_proteomes/
      proteome_<species>.tsv    THE deliverable -- four tables, one per species
      accessory/
          locus_tags_<species>.tsv  locus_tag + locus_tag_all
          annotation_<species>.tsv  the wide xref layer, free in the same request
          name_audit.tsv            every name fill: donor, candidates, contested, rule
          registry.tsv              the registry as actually fetched
          manifest.tsv              label, url, n, sha256, release, fetched_at
          .uniref90_<species>.json  the clustering cache

Run with the `gradi` env:
    python scripts/00_download_proteomes.py                      # tiers A + B
    python scripts/00_download_proteomes.py --dry-run            # show the plan, fetch nothing
    python scripts/00_download_proteomes.py --tier C --tier D    # comparator panel + bridge strains
    python scripts/00_download_proteomes.py --only saureus__nctc8325__UP000008816 --refresh
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import io
import re
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

REPO_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = REPO_ROOT / "src" / "proteome_registry.tsv"
RAW_DIR = REPO_ROOT / "data" / "raw" / "00_proteomes"
OUT_DIR = REPO_ROOT / "data" / "processed" / "00_proteomes"
# Supporting detail, not staging artifacts on the way to something else -- hence
# "accessory" rather than "intermediate". Keeps the top level to the four tables.
ACC_DIR = OUT_DIR / "accessory"

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
IDMAP_RUN = "https://rest.uniprot.org/idmapping/run"
IDMAP_STATUS = "https://rest.uniprot.org/idmapping/status"
IDMAP_UNIREF = "https://rest.uniprot.org/idmapping/uniref/results"
NCBI_DL = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{acc}/download"

# --- the two output schemas -------------------------------------------------------------------
# Identity: what you need to recognise a protein and find it in the literature / a database.
IDENTITY_FIELDS = [
    "accession", "reviewed", "gene_primary", "gene_synonym", "gene_oln", "gene_orf",
    "protein_name", "length", "sequence", "organism_id", "xref_refseq", "xref_geneid", "id",
]
# Annotation: free in the same request. `kegg`/`string` live here by explicit request -- promoting
# either to the identity table later is a column move, not a refetch.
ANNOTATION_FIELDS = [
    "accession", "xref_kegg", "xref_string", "xref_embl", "xref_eggnog", "xref_biocyc",
    "xref_interpro", "xref_pfam", "xref_panther", "go_id", "ec", "protein_families",
    "xref_pdb", "xref_alphafolddb",
]
# The donor pool only ever supplies names, so it stays narrow (and therefore fast).
DONOR_FIELDS = ["accession", "reviewed", "gene_primary", "gene_synonym", "sequence"]

# UniProt TSV header label -> our column name. UniProt returns display labels, not field ids.
RENAME = {
    "Entry": "uniprot_ac", "Entry Name": "entry_name", "Reviewed": "reviewed_raw",
    "Protein names": "protein_name", "Gene Names (primary)": "gene_primary",
    "Gene Names (synonym)": "gene_synonym", "Gene Names (ordered locus)": "gene_oln",
    "Gene Names (ORF)": "gene_orf", "Length": "length", "Sequence": "sequence",
    "Organism (ID)": "taxid", "RefSeq": "refseq", "GeneID": "geneid",
    "KEGG": "kegg", "STRING": "string", "EMBL": "embl", "eggNOG": "eggnog", "BioCyc": "biocyc",
    "InterPro": "interpro", "Pfam": "pfam", "PANTHER": "panther",
    "Gene Ontology IDs": "go_id", "EC number": "ec", "Protein families": "protein_families",
    "PDB": "pdb", "AlphaFoldDB": "alphafolddb",
}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def banner(args: argparse.Namespace) -> None:
    rule("=")
    say("STAGE 00 - reference proteomes, and one identified table per species")
    rule("=")
    say(f"  registry : {REGISTRY_PATH.relative_to(REPO_ROOT)}")
    say(f"  raw out  : {RAW_DIR.relative_to(REPO_ROOT)}/")
    say(f"  stage out: {OUT_DIR.relative_to(REPO_ROOT)}/  (+ accessory/)")
    say(f"  tiers    : {', '.join(args.tier)}")
    if args.only:
        say(f"  only     : {', '.join(args.only)}")
    say(f"  refresh  : {args.refresh}")
    say()

# ------------------------------------------------------------------ registry

def load_registry() -> pd.DataFrame:
    if not REGISTRY_PATH.exists():
        sys.exit(f"registry not found: {REGISTRY_PATH}")
    reg = pd.read_csv(REGISTRY_PATH, sep="\t", dtype=str, keep_default_na=False)
    missing = [c for c in ("label", "tier", "source", "id", "role", "species", "why") if c not in reg]
    if missing:
        sys.exit(f"registry is missing required columns: {missing}")
    blank = reg.loc[reg["why"].str.strip() == "", "label"].tolist()
    if blank:
        sys.exit("every registry row must carry a `why`; these do not: " + ", ".join(blank))
    dups = reg["label"][reg["label"].duplicated()].tolist()
    if dups:
        sys.exit(f"duplicate registry labels: {dups}")
    reg["expected_n"] = pd.to_numeric(reg["expected_n"], errors="coerce")
    return reg


def species_key(row: pd.Series) -> str:
    """Short species slug used for output filenames: kpneumoniae, ecoli, saureus, human."""
    return row["label"].split("__")[0]


# ------------------------------------------------------------------ HTTP

@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def _get(url: str, **kw) -> requests.Response:
    r = requests.get(url, timeout=kw.pop("timeout", 900), **kw)
    r.raise_for_status()
    return r


def uniprot_query(row: pd.Series) -> str:
    """The UniProt query for a registry row -- proteome for A/C/D, taxonomy for the B donor pools."""
    if row["source"] == "uniprot_species":
        return f"taxonomy_id:{row['id']} AND gene:*"
    q = f"proteome:{row['id']}"
    if row.get("reviewed_only", "no") == "yes":
        q += " AND reviewed:true"
    return q


def expected_from_uniprot(query: str) -> int:
    """Authoritative row count, read from X-Total-Results before streaming anything.

    This is the guard against v1's most common failure: a short or empty payload behind an HTTP
    200 (MobiDB's 61-byte 200, figshare's empty 202, a 200 from the wrong PMCID).
    """
    url = f"{UNIPROT_SEARCH}?query={quote(query, safe=':*')}&format=list&size=1"
    r = _get(url, timeout=120)
    n = r.headers.get("x-total-results")
    if n is None:
        raise RuntimeError(f"UniProt returned no X-Total-Results for: {query}")
    return int(n)


def uniprot_release(resp: requests.Response) -> tuple[str, str]:
    return (resp.headers.get("x-uniprot-release", "?"),
            resp.headers.get("x-uniprot-release-date", "?"))


def stream_tsv(query: str, fields: list[str]) -> tuple[pd.DataFrame, str, tuple[str, str]]:
    """Stream a UniProt TSV, decompressing explicitly and asserting the header parses.

    UniProt gzips when compressed=true and a naive save yields a file that reads as UTF-8 garbage,
    so the decompression is done here rather than trusted to Content-Encoding.
    """
    url = (f"{UNIPROT_STREAM}?compressed=true&format=tsv"
           f"&query={quote(query, safe=':*')}&fields={','.join(fields)}")
    r = _get(url)
    raw = r.content
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    text = raw.decode("utf-8")
    if not text.startswith("Entry"):
        raise RuntimeError(f"unexpected TSV header (first 80 bytes): {text[:80]!r}")
    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str, keep_default_na=False)
    return df, text, uniprot_release(r)


def stream_fasta(query: str) -> tuple[str, int]:
    url = f"{UNIPROT_STREAM}?compressed=true&format=fasta&query={quote(query, safe=':*')}"
    r = _get(url)
    raw = r.content
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    text = raw.decode("utf-8")
    return text, text.count("\n>") + (1 if text.startswith(">") else 0)


# ------------------------------------------------------------------ UniRef90 via ID mapping

def uniref90_map(accessions: list[str], chunk: int = 5000) -> dict[str, str]:
    """accession -> UniRef90 cluster id, via UniProt's ID-mapping service.

    This is what supplies the `species_uniref90` name tier: joining anchor and donors on their
    cluster gives >=90%-identity grouping with no alignment step and no DIAMOND dependency.

    The protocol has two sharp edges, both handled here:
      * `/idmapping/status/{job}` answers 200 `{"jobStatus":"RUNNING"}` while working, then **303
        with a Location header** naming the results URL. Following that redirect blindly yields a
        400, so redirects are disabled and the Location is read explicitly.
      * results paginate at 500 rows; the next page is in the `Link: <...>; rel="next"` header.
        Reading only the first response silently keeps 500 of N clusters.
    """
    out: dict[str, str] = {}
    total = len(accessions)
    for i in range(0, total, chunk):
        part = accessions[i:i + chunk]
        say(f"      UniRef90 map {i + 1}-{i + len(part)} of {total} ...")
        job = requests.post(IDMAP_RUN, timeout=120, data={
            "from": "UniProtKB_AC-ID", "to": "UniRef90", "ids": ",".join(part)}).json()
        jid = job.get("jobId")
        if not jid:
            raise RuntimeError(f"ID-mapping did not return a jobId: {job}")

        results_url = None
        for _ in range(300):
            st = requests.get(f"{IDMAP_STATUS}/{jid}", timeout=120, allow_redirects=False)
            if st.status_code in (303, 302):
                results_url = st.headers.get("Location") or f"{IDMAP_UNIREF}/{jid}"
                break
            if st.status_code != 200:
                raise RuntimeError(f"ID-mapping status {st.status_code}: {st.text[:200]}")
            body = st.json()
            if body.get("jobStatus") == "ERROR" or body.get("messages"):
                raise RuntimeError(f"ID-mapping job failed: {body}")
            if body.get("jobStatus") == "FINISHED" or "results" in body:
                results_url = f"{IDMAP_UNIREF}/{jid}"
                break
            time.sleep(2)
        if not results_url:
            raise RuntimeError(f"ID-mapping job {jid} did not finish in time")

        sep = "&" if "?" in results_url else "?"
        url = f"{results_url}{sep}format=tsv&fields=id&size=500"
        got, first = 0, True
        while url:
            r = _get(url, timeout=300)
            lines = r.text.splitlines()
            for line in (lines[1:] if first else lines):
                cols = line.split("\t")
                if len(cols) >= 2 and cols[1]:
                    out[cols[0]] = cols[1]
                    got += 1
            first = False
            m = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("Link", ""))
            url = m.group(1) if m else None
        say(f"        -> {got} clusters resolved")
    return out


# ------------------------------------------------------------------ name enrichment

# Automated annotation numbers duplicated loci as `phoE_2`, `manX_1`. That suffix is not part of
# the gene name -- underscore-digit endings are not valid bacterial nomenclature -- so it is
# stripped before voting. Without this, `phoE` vs `phoE_2` looks like a tie and abstains, which
# discards 23 Kp proteins for no reason. Genuinely different names (`maf` vs `yhdE`) still tie.
_GENE_SUFFIX = re.compile(r"_\d+$")


def _clean_gene(g: str) -> str:
    return _GENE_SUFFIX.sub("", (g or "").strip())


def _norm(df: pd.DataFrame) -> pd.DataFrame:
    df = df.rename(columns=RENAME)
    if "reviewed_raw" in df:
        df["is_reviewed"] = df["reviewed_raw"].str.strip().eq("reviewed")
        df = df.drop(columns=["reviewed_raw"])
    return df


PLACEHOLDER_RX = re.compile(r"^y[a-z]{2}[A-Z]?[0-9]*$")


def _is_placeholder(gene: str) -> bool:
    """`yhdE`, `yqiE`, `yjgP` are systematic placeholders assigned before a function was known.

    Enterobacterial convention: once a gene is characterised it gets a real name (`maf`, `lptF`),
    and the y-name survives only as a synonym. So a non-y name is always the better preferred name.
    """
    return bool(PLACEHOLDER_RX.match(gene))


def build_donor_index(donor: pd.DataFrame) -> tuple[dict, dict]:
    """sequence -> {gene: {count, reviewed, acc}}, and sequence -> set(synonyms).

    `acc` is an accession that actually SUPPLIES that gene name, preferring a reviewed one. It must
    not be "the first accession sharing this sequence" -- the anchor's own entry is usually in the
    pool (it matches `gene:*` on its locus tag) while carrying no primary name, so that shortcut
    attributes every name to the anchor itself and destroys the audit trail.
    """
    from collections import defaultdict
    by_seq: dict[str, dict] = defaultdict(dict)
    syn_by_seq: dict[str, set] = defaultdict(set)
    rev_col = donor["is_reviewed"] if "is_reviewed" in donor else [False] * len(donor)
    for acc, seq, gene, syn, rev in zip(donor["uniprot_ac"], donor["sequence"],
                                        donor["gene_primary"], donor["gene_synonym"], rev_col):
        seq = (seq or "").strip()
        gene = _clean_gene(gene)
        if not seq or not gene:
            continue
        slot = by_seq[seq].setdefault(gene, {"count": 0, "reviewed": False, "acc": acc})
        slot["count"] += 1
        if bool(rev) and not slot["reviewed"]:
            slot["acc"] = acc                    # a reviewed supporter outranks an unreviewed one
        slot["reviewed"] = slot["reviewed"] or bool(rev)
        for token in (syn or "").split():
            syn_by_seq[seq].add(token)
    return by_seq, syn_by_seq


def _pick(cands: dict) -> tuple[str, list[str], bool, str]:
    """Choose ONE preferred name under a total order, so a tie can never leave it empty.

    Preference, in order:
      1. a name supported by a REVIEWED (Swiss-Prot) donor beats an unreviewed one - curation wins
      2. a real name beats a `y###` systematic placeholder
      3. the more frequently attested name
      4. the shorter name, then alphabetical - a final deterministic tiebreak

    Returns (preferred, all_candidates_in_preference_order, was_contested, deciding_rule).
    """
    if not cands:
        return "", [], False, ""

    def key(g: str):
        d = cands[g]
        return (0 if d["reviewed"] else 1, 1 if _is_placeholder(g) else 0,
                -d["count"], len(g), g)

    ordered = sorted(cands, key=key)
    best, contested = ordered[0], len(cands) > 1
    rule = "sole candidate"
    if contested:
        runner = ordered[1]
        if cands[best]["reviewed"] != cands[runner]["reviewed"]:
            rule = "reviewed donor"
        elif _is_placeholder(runner) and not _is_placeholder(best):
            rule = "non-placeholder over y-name"
        elif cands[best]["count"] != cands[runner]["count"]:
            rule = f"more attested ({cands[best]['count']} vs {cands[runner]['count']})"
        else:
            rule = "shortest then alphabetical"
    return best, ordered, contested, rule


def enrich_names(ident: pd.DataFrame, donor: pd.DataFrame | None,
                 uniref: dict[str, str] | None) -> tuple[pd.DataFrame, list[dict]]:
    """Fill gene_name in three labelled tiers, all within the species. Never overwrites."""
    ident = ident.copy()
    ident["gene_name"] = ident["gene_primary"].fillna("").str.strip()
    ident["gene_synonyms"] = ident["gene_synonym"].fillna("").str.strip()
    ident["gene_name_source"] = ident["gene_name"].map(lambda g: "anchor" if g else "none")
    ident["gene_name_donor"] = ""
    ident["gene_name_candidates"] = ""

    audit: list[dict] = []
    if donor is None or donor.empty:
        return ident, audit

    by_seq, syn_by_seq = build_donor_index(donor)

    # --- tier 2: identical sequence, same species
    filled = 0
    for i, row in ident.iterrows():
        if row["gene_name"]:
            continue
        seq = (row["sequence"] or "").strip()
        if seq not in by_seq:
            continue
        chosen, cands, contested, rule = _pick(by_seq[seq])
        donor_acc = by_seq[seq][chosen]["acc"]
        ident.at[i, "gene_name_candidates"] = ";".join(cands)
        ident.at[i, "gene_name"] = chosen
        ident.at[i, "gene_name_source"] = "species_exact"
        ident.at[i, "gene_name_donor"] = donor_acc
        filled += 1
        audit.append(dict(uniprot_ac=row["uniprot_ac"], chosen=chosen, source="species_exact",
                          donor=donor_acc, candidates=";".join(cands),
                          contested=contested, rule=rule))
    say(f"      species_exact  filled {filled}")

    # synonyms travel with an identical sequence regardless of whether the name was already known
    gained = 0
    for i, row in ident.iterrows():
        seq = (row["sequence"] or "").strip()
        extra = syn_by_seq.get(seq)
        if not extra:
            continue
        have = set(filter(None, (row["gene_synonyms"] or "").split()))
        merged = sorted(have | set(extra))
        if merged and merged != sorted(have):
            ident.at[i, "gene_synonyms"] = " ".join(merged)
            gained += 1
    say(f"      synonyms       gained {gained}")

    # --- tier 3: same UniRef90 cluster, same species
    if uniref:
        from collections import defaultdict
        by_cluster: dict[str, dict] = defaultdict(dict)
        rev_col = donor["is_reviewed"] if "is_reviewed" in donor else [False] * len(donor)
        for acc, gene, rev in zip(donor["uniprot_ac"], donor["gene_primary"], rev_col):
            cid, gene = uniref.get(acc), _clean_gene(gene)
            if cid and gene:
                slot = by_cluster[cid].setdefault(gene, {"count": 0, "reviewed": False, "acc": acc})
                slot["count"] += 1
                if bool(rev) and not slot["reviewed"]:
                    slot["acc"] = acc
                slot["reviewed"] = slot["reviewed"] or bool(rev)
        filled = 0
        for i, row in ident.iterrows():
            if row["gene_name"]:
                continue
            cid = uniref.get(row["uniprot_ac"])
            if not cid or cid not in by_cluster:
                continue
            chosen, cands, contested, rule = _pick(by_cluster[cid])
            donor_acc = by_cluster[cid][chosen]["acc"]
            ident.at[i, "gene_name_candidates"] = ";".join(cands)
            ident.at[i, "gene_name"] = chosen
            ident.at[i, "gene_name_source"] = "species_uniref90"
            ident.at[i, "gene_name_donor"] = donor_acc
            filled += 1
            audit.append(dict(uniprot_ac=row["uniprot_ac"], chosen=chosen,
                              source="species_uniref90", donor=donor_acc,
                              candidates=";".join(cands), contested=contested, rule=rule))
        say(f"      species_uniref90 filled {filled}")

    return ident, audit


# ------------------------------------------------------------------ fetch one registry row

def source_md(row: pd.Series, query: str, expected: int, observed: int,
              release: tuple[str, str], urls: list[str]) -> str:
    return "\n".join([
        f"# {row['label']}", "",
        f"- species: {row['species']}  {row['strain']}".rstrip(),
        f"- source: {row['source']}  id: {row['id']}",
        f"- taxid: {row['taxid']}   assembly: {row['assembly'] or '-'}",
        f"- locus prefix: {row['locus_prefix'] or '-'}",
        f"- role/tier: {row['role']} / {row['tier']}",
        f"- reviewed_only: {row.get('reviewed_only', 'no')}",
        f"- expected: {expected}   observed: {observed}",
        f"- UniProt release: {release[0]} ({release[1]})",
        f"- fetched: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"- query: {query}", "",
        "## URLs", "", *[f"- {u}" for u in urls], "",
        "## Why this proteome is in the registry", "", row["why"], "",
    ])


def fetch_uniprot_row(row: pd.Series, refresh: bool) -> dict:
    label, is_pool = row["label"], row["source"] == "uniprot_species"
    out_dir = RAW_DIR / "uniprot"
    out_dir.mkdir(parents=True, exist_ok=True)
    tsv_path, fasta_path = out_dir / f"{label}.tsv", out_dir / f"{label}.fasta"
    src_path = out_dir / f"{label}.SOURCE.md"

    query = uniprot_query(row)
    fields = DONOR_FIELDS if is_pool else IDENTITY_FIELDS + ANNOTATION_FIELDS[1:]

    if tsv_path.exists() and src_path.exists() and not refresh:
        n = sum(1 for _ in tsv_path.open()) - 1
        if row["expected_n"] and n == int(row["expected_n"]):
            say(f"    [skip] {label}  ({n} rows already on disk)")
            return dict(label=label, skipped=True, n=n, path=tsv_path)
        say(f"    [refetch] {label}  on-disk {n} != expected {row['expected_n']}")

    expected = expected_from_uniprot(query)
    if row["expected_n"] and expected != int(row["expected_n"]):
        say(f"    [!] {label}: registry says {int(row['expected_n'])}, UniProt now reports "
            f"{expected} -- UniProt release drift, using the live count")

    say(f"    fetching {label}  (expect {expected} rows)")
    df, text, release = stream_tsv(query, fields)
    if len(df) != expected:
        raise RuntimeError(
            f"{label}: streamed {len(df)} rows but X-Total-Results said {expected}. "
            "Refusing a short payload -- this is the v1 silent-truncation failure mode.")
    tsv_path.write_text(text)
    say(f"      TSV   {len(df):>6} rows -> {tsv_path.name}")

    urls = [f"{UNIPROT_STREAM}?format=tsv&query={quote(query, safe=':*')}&fields={','.join(fields)}"]
    n_fa = 0
    if not is_pool:                      # a donor pool needs no FASTA; it is only a name source
        fa, n_fa = stream_fasta(query)
        if n_fa != expected:
            raise RuntimeError(f"{label}: FASTA has {n_fa} records, expected {expected}")
        fasta_path.write_text(fa)
        say(f"      FASTA {n_fa:>6} records -> {fasta_path.name}")
        urls.append(f"{UNIPROT_STREAM}?format=fasta&query={quote(query, safe=':*')}")

    src_path.write_text(source_md(row, query, expected, len(df), release, urls))
    return dict(label=label, skipped=False, n=len(df), n_fasta=n_fa, path=tsv_path,
                release=release[0], release_date=release[1], url=urls[0],
                sha256=hashlib.sha256(text.encode()).hexdigest())


def fetch_ncbi_row(row: pd.Series, refresh: bool) -> dict:
    """Tier D bridge strains: protein FASTA + GFF, for locus-tag bridging only."""
    label, acc = row["label"], row["id"]
    out_dir = RAW_DIR / "ncbi"
    out_dir.mkdir(parents=True, exist_ok=True)
    faa_path, gff_path = out_dir / f"{label}.faa", out_dir / f"{label}.gff"
    if faa_path.exists() and gff_path.exists() and not refresh:
        n = sum(1 for ln in faa_path.open() if ln.startswith(">"))
        say(f"    [skip] {label}  ({n} proteins already on disk)")
        return dict(label=label, skipped=True, n=n, path=faa_path)

    url = NCBI_DL.format(acc=acc) + "?include_annotation_type=PROT_FASTA,GENOME_GFF"
    say(f"    fetching {label}  ({acc})")
    r = _get(url)
    zf = zipfile.ZipFile(io.BytesIO(r.content))
    faa = [n for n in zf.namelist() if n.endswith("protein.faa")]
    gff = [n for n in zf.namelist() if n.endswith(".gff")]
    if not faa:
        raise RuntimeError(f"{label}: NCBI returned no protein.faa (members: {zf.namelist()[:6]})")
    faa_text = zf.read(faa[0]).decode()
    faa_path.write_text(faa_text)
    n = faa_text.count("\n>") + (1 if faa_text.startswith(">") else 0)
    say(f"      FAA {n:>6} proteins -> {faa_path.name}")
    if gff:
        gff_path.write_text(zf.read(gff[0]).decode())
        say(f"      GFF          -> {gff_path.name}")
    (out_dir / f"{label}.SOURCE.md").write_text(
        source_md(row, f"NCBI datasets {acc}", int(row["expected_n"] or 0), n, ("-", "-"), [url]))
    return dict(label=label, skipped=False, n=n, path=faa_path, url=url,
                sha256=hashlib.sha256(faa_text.encode()).hexdigest())


# ------------------------------------------------------------------ locus-tag bridge

def locus_bridge_from_gff(gff_path: Path) -> pd.DataFrame:
    """protein_id -> every locus-tag namespace NCBI knows for it.

    `old_locus_tag` is why this matters: for ECL8 it carries BOTH BN373_00001 and KPNEcl8_00001,
    which is the bridge v1 ran Prokka to reconstruct.
    """
    rows = []
    if not gff_path.exists():
        return pd.DataFrame(columns=["protein_id", "locus_tag", "old_locus_tags"])
    for line in gff_path.read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue
        attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        pid = attrs.get("protein_id")
        if not pid:
            continue
        rows.append(dict(
            protein_id=pid,
            locus_tag=attrs.get("locus_tag", ""),
            old_locus_tags=attrs.get("old_locus_tag", "").replace("%2C", ";"),
        ))
    return pd.DataFrame(rows).drop_duplicates("protein_id")


# ------------------------------------------------------------------ assemble the deliverable

# Nine columns: identity and nothing else. Everything dropped from the first run's 17 was measured
# rather than guessed:
#   taxid, species     -- exactly 1 distinct value per file, now that there is one table per species
#   sequence_md5       -- verified identical to md5(sequence), and `sequence` stays
#   length             -- len(sequence)
#   gene_name_donor, gene_name_candidates -> accessory/name_audit.tsv (provenance detail)
#   locus_tag, locus_tag_all             -> accessory/locus_tags_<species>.tsv
# The locus tags are the join key for published bacterial data, so they are kept in full next door
# rather than discarded -- including `locus_tag_all`, which carries the Keio JW ids for 4,252 of
# 4,403 E. coli rows. Join them back with src.proteomes.with_locus_tags().
IDENTITY_OUT = [
    "uniprot_ac", "is_reviewed", "gene_name", "gene_name_source", "gene_synonyms",
    "protein_name", "sequence", "refseq", "geneid",
]
LOCUS_OUT = ["uniprot_ac", "locus_tag", "locus_tag_all"]


def build_identity(anchor: pd.DataFrame, donor: pd.DataFrame | None,
                   uniref: dict | None) -> tuple[pd.DataFrame, list[dict], int, pd.DataFrame]:
    df = _norm(anchor)
    df, audit = enrich_names(df, donor, uniref)

    # locus_tag: the ordered-locus name is the real join key for published bacterial data.
    # UniProt returns ";" (or "; ; ; ;") for an organism with no ordered locus names -- human -- so
    # strip separator-only tokens rather than emitting punctuation as an identifier.
    oln = df.get("gene_oln", pd.Series([""] * len(df))).fillna("").astype(str)
    toks = oln.map(lambda v: [t for t in v.replace(";", " ").split() if t])
    df["locus_tag"] = toks.map(lambda t: t[0] if t else "")
    df["locus_tag_all"] = toks.map(" ".join)
    for c in ("refseq", "geneid"):
        if c not in df:
            df[c] = ""

    # gene_name is a PREFERRED name, not a maybe: any protein with a candidate must have one. This
    # is checked here rather than in source_tally() because gene_name_candidates does not survive
    # into the 11-column output -- it lives in accessory/name_audit.tsv.
    orphan = df[df["gene_name"].str.strip().eq("")
                & df["gene_name_candidates"].str.strip().ne("")]
    if len(orphan):
        raise RuntimeError(
            f"{len(orphan)} rows carry gene_name_candidates but no preferred gene_name "
            f"(e.g. {orphan['uniprot_ac'].head(3).tolist()}) - the total order in _pick() should "
            "make this impossible")
    n_contested = int(df["gene_name_candidates"].str.contains(";", na=False).sum())

    out = df.reindex(columns=IDENTITY_OUT)
    locus = df.reindex(columns=LOCUS_OUT)
    return out, audit, n_contested, locus


def coverage_table(df: pd.DataFrame, title: str) -> None:
    n = len(df)
    say(f"    {title}  (n = {n})")
    say(f"      {'column':<24} {'filled':>7} {'pct':>7}")
    for c in df.columns:
        if c == "sequence":
            continue
        s = df[c]
        filled = int(s.astype(str).str.strip().replace({"nan": "", "False": ""}).ne("").sum()) \
            if s.dtype == object else int(s.notna().sum())
        say(f"      {c:<24} {filled:>7} {100 * filled / n:>6.1f}%")


def source_tally(df: pd.DataFrame, n_contested: int = 0) -> None:
    say("      gene_name_source:")
    counts = df["gene_name_source"].value_counts()
    for k in ("anchor", "species_exact", "species_uniref90", "none"):
        v = int(counts.get(k, 0))
        say(f"        {k:<20} {v:>7} {100 * v / len(df):>6.1f}%")
    named = int(df["gene_name"].str.strip().ne("").sum())
    say(f"      == named total        {named:>7} {100 * named / len(df):>6.1f}%")
    bad = set(df["gene_name_source"]) - {"anchor", "species_exact", "species_uniref90", "none"}
    if bad:
        raise RuntimeError(f"unexpected gene_name_source values (cross-species leak?): {bad}")

    say(f"      contested (>1 candidate, all resolved to a preferred name): {n_contested}")


# ------------------------------------------------------------------ main

def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tier", action="append", choices=list("ABCD"),
                    help="tiers to fetch (repeatable; default A and B)")
    ap.add_argument("--only", nargs="+", metavar="LABEL", help="fetch only these registry labels")
    ap.add_argument("--refresh", action="store_true", help="refetch even if files are on disk")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the resolved plan and expected counts, fetch nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    args.tier = args.tier or ["A", "B"]

    reg = load_registry()
    banner(args)

    sel = reg[reg["tier"].isin(args.tier)]
    if args.only:
        sel = reg[reg["label"].isin(args.only)]
        unknown = set(args.only) - set(reg["label"])
        if unknown:
            sys.exit(f"unknown registry labels: {sorted(unknown)}")

    say(f"Registry: {len(reg)} rows total, {len(sel)} selected")
    for t in sorted(sel["tier"].unique()):
        part = sel[sel["tier"] == t]
        say(f"  tier {t}: {len(part)} rows, {int(part['expected_n'].fillna(0).sum()):,} expected records")
    say()

    if args.dry_run:
        rule()
        say(f"{'label':<48} {'tier':<5} {'source':<16} {'id':<18} {'expect':>8}")
        rule()
        for _, r in sel.iterrows():
            say(f"{r['label']:<48} {r['tier']:<5} {r['source']:<16} {r['id']:<18} "
                f"{'' if pd.isna(r['expected_n']) else int(r['expected_n']):>8}")
        rule()
        say("\ndry run - nothing fetched.")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ACC_DIR.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []

    # ---- fetch
    rule()
    say("FETCH")
    rule()
    for tier in sorted(sel["tier"].unique()):
        part = sel[sel["tier"] == tier]
        say(f"  tier {tier}  ({len(part)} rows)")
        for i, (_, r) in enumerate(part.iterrows(), 1):
            say(f"   [{i}/{len(part)}] {r['label']}")
            res = fetch_ncbi_row(r, args.refresh) if r["source"] == "ncbi" \
                else fetch_uniprot_row(r, args.refresh)
            res["tier"], res["source"], res["id"] = r["tier"], r["source"], r["id"]
            manifest.append(res)
    say()

    # ---- assemble, per species
    anchors = reg[(reg["tier"] == "A") & reg["label"].isin(sel["label"])]
    if len(anchors):
        rule()
        say("ASSEMBLE")
        rule()
    stacked, all_audit, bridge_rows = [], [], []
    for _, r in anchors.iterrows():
        sp = species_key(r)
        say(f"  {sp}  ({r['species']} {r['strain']})".rstrip())
        tsv = RAW_DIR / "uniprot" / f"{r['label']}.tsv"
        anchor = pd.read_csv(tsv, sep="\t", dtype=str, keep_default_na=False)

        pool = reg[(reg["tier"] == "B") & (reg["label"].str.startswith(sp + "__"))]
        donor, uniref = None, None
        if len(pool) and (RAW_DIR / "uniprot" / f"{pool.iloc[0]['label']}.tsv").exists():
            donor = _norm(pd.read_csv(RAW_DIR / "uniprot" / f"{pool.iloc[0]['label']}.tsv",
                                      sep="\t", dtype=str, keep_default_na=False))
            say(f"    donor pool: {len(donor)} entries with a gene field")
            unnamed = int(_norm(anchor)["gene_primary"].fillna("").str.strip().eq("").sum())
            say(f"    anchor proteins with no gene name: {unnamed}")
            if unnamed == 0:
                # Nothing for the UniRef90 tier to fill, so skip a clustering pass that would cost
                # ~22 s per 5,000 accessions for zero names. The donor pool is still used, for
                # synonyms. Note E. coli K-12 does NOT hit this branch: it has exactly one unnamed
                # protein, so it clusters once and then reads its cache.
                say("    -> anchor fully named; skipping UniRef90 (nothing to fill)")
                accs = []
            else:
                accs = sorted(set(_norm(anchor)["uniprot_ac"]) | set(donor["uniprot_ac"]))
            # Cache the clustering: it is ~22 s per 5,000 accessions, so a re-run of the assembly
            # step would otherwise cost ~10 min of pure re-mapping. Keyed on the accession set.
            cache = ACC_DIR / f".uniref90_{sp}.json"
            key = hashlib.md5("\n".join(accs).encode()).hexdigest()
            if not accs:
                pass
            elif cache.exists() and not args.refresh:
                blob = json.loads(cache.read_text())
                if blob.get("key") == key:
                    uniref = blob["map"]
                    say(f"    UniRef90 clusters from cache ({len(uniref)} accessions)")
            if uniref is None and accs:
                say(f"    UniRef90 clustering for {len(accs)} accessions ...")
                uniref = uniref90_map(accs)
                cache.write_text(json.dumps({"key": key, "map": uniref}))
                say(f"    cached -> {cache.name}")
        else:
            say("    no donor pool selected - names come from the anchor only")

        ident, audit, n_contested, locus = build_identity(anchor, donor, uniref)
        for a in audit:
            a["species"] = sp
        all_audit += audit

        ip = OUT_DIR / f"proteome_{sp}.tsv"
        ident.to_csv(ip, sep="\t", index=False)
        say(f"    wrote {ip.name}  ({len(ident)} rows x {ident.shape[1]} cols)")

        ann = _norm(anchor).reindex(
            columns=["uniprot_ac"] + [RENAME.get(c, c) for c in
                                      ["KEGG", "STRING", "EMBL", "eggNOG", "BioCyc", "InterPro",
                                       "Pfam", "PANTHER", "Gene Ontology IDs", "EC number",
                                       "Protein families", "PDB", "AlphaFoldDB"]])
        lp = ACC_DIR / f"locus_tags_{sp}.tsv"
        locus.to_csv(lp, sep="\t", index=False)
        n_lt = int(locus["locus_tag"].str.strip().ne("").sum())
        say(f"    wrote accessory/{lp.name}  ({n_lt}/{len(locus)} rows carry a locus tag)")

        ap_ = ACC_DIR / f"annotation_{sp}.tsv"
        ann.to_csv(ap_, sep="\t", index=False)
        say(f"    wrote accessory/{ap_.name}  ({ann.shape[1]} columns)")

        source_tally(ident, n_contested)
        stacked.append((sp, ident))
        say()

    # ---- tier D locus bridges
    for _, r in sel[sel["source"] == "ncbi"].iterrows():
        b = locus_bridge_from_gff(RAW_DIR / "ncbi" / f"{r['label']}.gff")
        if len(b):
            b["label"] = r["label"]
            bridge_rows.append(b)
            say(f"  bridge {r['label']}: {len(b)} CDS, "
                f"{int(b['old_locus_tags'].ne('').sum())} with legacy tags")

    # ---- outputs
    rule()
    say("OUTPUTS")
    rule()
    if stacked:
        say(f"  {len(stacked)} tables at the top level:")
        for sp, d in stacked:
            say(f"    proteome_{sp}.tsv{'':<{max(0, 14 - len(sp))}} {len(d):>6} rows x {d.shape[1]} cols")

        # No stacked parquet and no id_bridge: both were verified pure derivations of these tables
        # (19 MB and 8 MB for zero new information). Use src.proteomes.load_all() for a stacked view.
        if all_audit:
            adf = pd.DataFrame(all_audit)
            adf.to_csv(ACC_DIR / "name_audit.tsv", sep="\t", index=False)
            say(f"  accessory/name_audit.tsv  {len(adf)} fills, "
                f"{int(adf['contested'].sum())} contested (all resolved to a preferred name)")

    if bridge_rows:
        pd.concat(bridge_rows, ignore_index=True).to_csv(
            ACC_DIR / "locus_bridge_strains.tsv", sep="\t", index=False)
        say("  accessory/locus_bridge_strains.tsv written")

    sel.to_csv(ACC_DIR / "registry.tsv", sep="\t", index=False)
    pd.DataFrame(manifest).to_csv(ACC_DIR / "manifest.tsv", sep="\t", index=False)
    say(f"  accessory/registry.tsv + manifest.tsv  ({len(manifest)} fetched rows)")

    if stacked:
        say()
        rule()
        say("COVERAGE")
        rule()
        for sp, d in stacked:
            coverage_table(d, sp)
            say()
    rule("=")
    say("stage 00 complete.")
    rule("=")


if __name__ == "__main__":
    main()
