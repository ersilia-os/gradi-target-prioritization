"""Annotate a downloaded proteome with a "well-studied-ness" / popularity signal.

Pharos and Open Targets — the canonical target-popularity resources — are
human-only and don't cover bacterial pathogens. This script computes the next
best thing using two free, UniProt-AC-keyed signals:

1. UniProt's own annotation-depth fields, fetched in a single proteome stream
   query (``reviewed``, protein existence level, annotation score, count of
   cited PubMed references).
2. Europe PMC literature counts, queried as ``<gene_symbol> AND <organism>``
   (e.g. ``cpxR AND klebsiella``).

Why not query Europe PMC's ``ACCESSION_ID`` field directly? It is essentially
empty for bacterial protein accessions — e.g. ``Q9F663`` (KPC β-lactamase, a
heavily-studied carbapenemase) returns only ~4 hits via ``ACCESSION_ID``, while
``KPC AND klebsiella`` returns over 13,000. Strain-specific UniProt entries
like the ``A0A0H3G*`` entries for HS11286 are virtually never cited by AC in
the literature; gene symbols anchored to a genus keyword recover the real
signal while keeping homonyms in check (``asdf123notagene AND klebsiella``
returns 0).

From these we derive a 3-level ``popularity_tier`` (top-down; first match wins):

    well_studied: europepmc_n_articles >= 200, OR
                  (reviewed AND PE in {1,2} AND uniprot_n_references >= 5)
    studied:      europepmc_n_articles >= 5, OR reviewed, OR PE in {1,2}
    dark:         everything else (typically: no canonical gene symbol, not
                  reviewed, PE 3+, and < 5 organism-anchored hits)

Thresholds were calibrated against a 100-protein HS11286 sample. Two facts
shape the rules:

* In bacterial proteomes ~95% of UniProt entries are PE=3 (Inferred from
  homology) — the default. PE=3 is therefore NOT evidence of being studied
  and does not promote a row out of ``dark``. Only PE 1/2 do.
* For HS11286, only ~2% of entries are reviewed (Swiss-Prot) — the strain's
  proteins live mostly in TrEMBL even when the underlying biology is well
  characterized via E. coli orthologs. Literature counts therefore have to
  carry the bulk of the popularity signal; ``europepmc_n_articles`` 200 sits
  around the 75th percentile of querable rows (famous core genes like rho,
  groEL, pgi, lamB, KPC clear it comfortably).

Revisit by inspecting the end-of-run histograms; tier thresholds are local
constants in ``compute_tier``.

Writes ``data/processed/<organism>_popularity.tsv``. ``src/assemble.py`` and
the v1 annotation table are intentionally NOT touched here; merging into the
v1 schema can happen in a follow-up once the raw signals look sensible.

Europe PMC calls are cached on disk under
``tmp/europepmc/<organism>/<gene_symbol>.json``, one tiny JSON per unique gene
symbol. Re-runs and KeyboardInterrupt resumes are near-instant. Accessions
sharing a gene symbol (paralogs, isoforms) all hit the same cache entry.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"
EUROPEPMC_SEARCH_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# Keep in sync with ORGANISMS in 00_download_proteome.py / 01_annotate_panther.py.
# ``europepmc_token`` is appended as ``AND <token>`` to every Europe PMC query
# to anchor results to the relevant taxonomic context and filter homonyms.
ORGANISMS: dict[str, dict[str, str]] = {
    "kpneumoniae": {
        "proteome_id": "UP000007841",
        "label": "Klebsiella pneumoniae subsp. pneumoniae HS11286",
        "filename": "klebsiella_pneumoniae_popularity.tsv",
        "europepmc_token": "klebsiella",
    },
    "ecoli": {
        "proteome_id": "UP000000625",
        "label": "Escherichia coli K-12 MG1655",
        "filename": "escherichia_coli_popularity.tsv",
        "europepmc_token": '"escherichia coli"',
    },
}

DEFAULT_ORGANISM = "kpneumoniae"
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_CACHE_DIR = Path("tmp/europepmc")
DEFAULT_RATE_LIMIT = 0.2  # seconds between Europe PMC calls

PE_TEXT_TO_LEVEL = {
    "evidence at protein level": 1,
    "evidence at transcript level": 2,
    "inferred from homology": 3,
    "predicted": 4,
    "uncertain": 5,
}

# Patterns we use to *reject* gene-name tokens that are locus tags rather than
# canonical gene symbols. Matches `anchor.py` (KPHS_*, b####, JW####) plus a
# generic-symbol shape check.
LOCUS_TAG_RXES = [
    re.compile(r"^KPHS_[0-9p]+$"),
    re.compile(r"^b\d{4}$"),
    re.compile(r"^JW\d+(?:\.\d+)?$"),
]
GENE_SYMBOL_RX = re.compile(r"^[a-zA-Z][a-zA-Z0-9_]{1,6}$")

# Cache filename safety: only allow letters/digits/underscore/dash.
CACHE_FILENAME_SAFE_RX = re.compile(r"[^A-Za-z0-9_.-]+")

OUTPUT_COLUMNS = [
    "accession",
    "gene_names",
    "gene_symbol",
    "uniprot_reviewed",
    "uniprot_pe_level",
    "uniprot_annotation_score",
    "uniprot_n_references",
    "europepmc_n_articles",
    "popularity_tier",
]


def http_get_text(url: str, timeout: int = 600) -> str:
    request = urllib.request.Request(url, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def http_get_json(url: str, timeout: int = 60) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_uniprot_popularity_tsv(proteome_id: str) -> str:
    params = urllib.parse.urlencode(
        {
            "query": f"proteome:{proteome_id}",
            "format": "tsv",
            "fields": (
                "accession,gene_names,reviewed,"
                "protein_existence,annotation_score,lit_pubmed_id"
            ),
        }
    )
    return http_get_text(f"{UNIPROT_STREAM_URL}?{params}")


def _find_col(idx: dict[str, int], *candidates: str) -> int | None:
    """Return the index of the first matching column header (exact, then
    case-insensitive substring). Returns None if nothing matches — caller is
    responsible for treating that as 'missing data'."""
    for c in candidates:
        if c in idx:
            return idx[c]
    lower = {k.lower(): v for k, v in idx.items()}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    for k_lower, v in lower.items():
        for c in candidates:
            if c.lower() in k_lower:
                return v
    return None


def _parse_pe_level(raw: str) -> int | None:
    if not raw:
        return None
    return PE_TEXT_TO_LEVEL.get(raw.strip().lower())


_INT_RX = re.compile(r"\d+")


def _parse_annotation_score(raw: str) -> int | None:
    if not raw:
        return None
    m = _INT_RX.search(raw)
    return int(m.group(0)) if m else None


def _count_pmids(raw: str) -> int:
    if not raw:
        return 0
    return sum(1 for tok in raw.split(";") if tok.strip())


def _is_locus_tag(token: str) -> bool:
    return any(rx.match(token) for rx in LOCUS_TAG_RXES)


def extract_gene_symbol(gene_names: str) -> str | None:
    """Pick the first canonical gene-symbol-looking token from a UniProt
    ``Gene Names`` value. Returns None if every token is a locus tag or fails
    the symbol shape check."""
    if not gene_names:
        return None
    for tok in gene_names.split():
        tok = tok.strip()
        if not tok or _is_locus_tag(tok):
            continue
        if GENE_SYMBOL_RX.match(tok):
            return tok
    return None


def parse_uniprot_popularity_tsv(tsv: str) -> list[dict]:
    lines = tsv.splitlines()
    if not lines:
        return []
    header = lines[0].split("\t")
    idx = {col: i for i, col in enumerate(header)}

    i_entry = _find_col(idx, "Entry", "accession")
    i_names = _find_col(idx, "Gene Names", "gene_names")
    i_reviewed = _find_col(idx, "Reviewed", "reviewed")
    i_pe = _find_col(idx, "Protein existence", "protein_existence")
    i_score = _find_col(idx, "Annotation", "Annotation Score", "annotation_score")
    i_pmids = _find_col(idx, "PubMed ID", "lit_pubmed_id")

    if i_entry is None:
        raise RuntimeError(
            f"UniProt TSV missing Entry/accession column; header was: {header!r}"
        )

    rows: list[dict] = []
    for line in lines[1:]:
        if not line:
            continue
        cells = line.split("\t")

        def cell(i: int | None) -> str:
            if i is None or i >= len(cells):
                return ""
            return cells[i]

        accession = cell(i_entry).strip()
        if not accession:
            continue
        gene_names = cell(i_names).strip()
        rows.append(
            {
                "accession": accession,
                "gene_names": gene_names,
                "gene_symbol": extract_gene_symbol(gene_names) or "",
                "uniprot_reviewed": cell(i_reviewed).strip().lower() == "reviewed",
                "uniprot_pe_level": _parse_pe_level(cell(i_pe)),
                "uniprot_annotation_score": _parse_annotation_score(cell(i_score)),
                "uniprot_n_references": _count_pmids(cell(i_pmids)),
            }
        )
    return rows


def _safe_cache_name(name: str) -> str:
    """Make a token safe to use as a filename. Empty input → 'NONE'."""
    if not name:
        return "NONE"
    safe = CACHE_FILENAME_SAFE_RX.sub("_", name)
    return safe[:80] or "NONE"


def europepmc_count_for_gene(
    gene_symbol: str,
    organism_token: str,
    cache_dir: Path,
    rate_limit: float,
    timeout: int = 60,
    retries: int = 1,
) -> tuple[int | None, bool]:
    """Return (hit_count, from_cache) for ``<gene_symbol> AND <organism_token>``.

    Cached by gene symbol (the cache dir is per-organism), so paralogs sharing
    a gene symbol all share one network call.
    """
    cache_path = cache_dir / f"{_safe_cache_name(gene_symbol)}.json"
    if cache_path.exists():
        try:
            data = json.loads(cache_path.read_text(encoding="utf-8"))
            hit_count = data.get("hitCount")
            if isinstance(hit_count, int):
                return hit_count, True
        except (json.JSONDecodeError, OSError):
            pass

    query = f"{gene_symbol} AND {organism_token}"
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "json",
            "pageSize": "1",
            "resultType": "lite",
        }
    )
    url = f"{EUROPEPMC_SEARCH_URL}?{params}"

    last_err: Exception | None = None
    for attempt in range(retries + 1):
        try:
            data = http_get_json(url, timeout=timeout)
            hit_count = int(data.get("hitCount", 0))
            cache_dir.mkdir(parents=True, exist_ok=True)
            cache_path.write_text(
                json.dumps({"query": query, "hitCount": hit_count}),
                encoding="utf-8",
            )
            time.sleep(rate_limit)
            return hit_count, False
        except Exception as e:
            last_err = e
            if attempt < retries:
                time.sleep(min(2.0, rate_limit * 5))
            else:
                break

    print(f"  ! europepmc failed for {gene_symbol}: {last_err}", file=sys.stderr)
    return None, False


def compute_tier(
    reviewed: bool,
    pe_level: int | None,
    n_articles: int | None,
    n_refs: int,
) -> str:
    """3-level popularity tier, top-down; first match wins.

    Missing n_articles (no gene symbol or Europe PMC failure) is treated as 0;
    it cannot promote a row on its own but ``reviewed`` / PE 1-2 still can.
    PE=3 (the default annotation for ~95% of bacterial entries) is treated as
    no signal — it must be paired with reviewed or literature to promote."""
    pe_in_top2 = pe_level in {1, 2}
    n_art = n_articles or 0
    if n_art >= 200 or (reviewed and pe_in_top2 and n_refs >= 5):
        return "well_studied"
    if n_art >= 5 or reviewed or pe_in_top2:
        return "studied"
    return "dark"


def _fmt_int(v: int | None) -> str:
    return "" if v is None else str(v)


def _fmt_bool(v: bool) -> str:
    return "True" if v else "False"


def write_tsv(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for r in rows:
            f.write(
                "\t".join(
                    [
                        r["accession"],
                        r["gene_names"],
                        r["gene_symbol"],
                        _fmt_bool(r["uniprot_reviewed"]),
                        _fmt_int(r["uniprot_pe_level"]),
                        _fmt_int(r["uniprot_annotation_score"]),
                        _fmt_int(r["uniprot_n_references"]),
                        _fmt_int(r["europepmc_n_articles"]),
                        r["popularity_tier"],
                    ]
                )
                + "\n"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--organism",
        choices=sorted(ORGANISMS),
        default=DEFAULT_ORGANISM,
        help="Organism shortname (default: %(default)s).",
    )
    parser.add_argument(
        "--proteome-id",
        default=None,
        help="Override the UniProt proteome identifier for --organism.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Destination TSV. Defaults to "
            f"{DEFAULT_PROCESSED_DIR}/<organism>_popularity.tsv."
        ),
    )
    parser.add_argument(
        "--cache",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help=(
            "Base directory for cached Europe PMC responses. One JSON per "
            "unique gene symbol is written under <cache>/<organism>/. "
            f"Default: {DEFAULT_CACHE_DIR}."
        ),
    )
    parser.add_argument(
        "--rate-limit",
        type=float,
        default=DEFAULT_RATE_LIMIT,
        help=(
            "Seconds to sleep between Europe PMC requests (default: "
            f"{DEFAULT_RATE_LIMIT}). Europe PMC permits up to 10 req/s."
        ),
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=200,
        help="Print a progress line every N accessions (default: %(default)s).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on the number of accessions processed (for debugging).",
    )
    args = parser.parse_args()

    organism = ORGANISMS[args.organism]
    proteome_id = args.proteome_id or organism["proteome_id"]
    organism_token = organism["europepmc_token"]
    output_path = args.output or DEFAULT_PROCESSED_DIR / organism["filename"]
    cache_dir = args.cache / args.organism

    print(
        f"Fetching UniProt popularity fields for "
        f"{organism['label']} ({proteome_id})"
    )
    tsv = fetch_uniprot_popularity_tsv(proteome_id)
    rows = parse_uniprot_popularity_tsv(tsv)
    if args.limit is not None:
        rows = rows[: args.limit]
    print(f"  parsed {len(rows)} UniProt entries")

    reviewed_n = sum(1 for r in rows if r["uniprot_reviewed"])
    pe_hist = Counter(r["uniprot_pe_level"] for r in rows)
    with_symbol = sum(1 for r in rows if r["gene_symbol"])
    print(
        f"  reviewed: {reviewed_n}, unreviewed: {len(rows) - reviewed_n}; "
        f"with gene symbol: {with_symbol}/{len(rows)}; "
        f"PE histogram: "
        + ", ".join(
            f"{lvl or 'NA'}={pe_hist.get(lvl, 0)}"
            for lvl in [1, 2, 3, 4, 5, None]
        )
    )

    print(
        f"Querying Europe PMC ('<gene> AND {organism_token}') "
        f"for {with_symbol} accessions with a gene symbol "
        f"(cache: {cache_dir}, rate-limit: {args.rate_limit}s)"
    )
    t0 = time.time()
    cache_hits = 0
    api_calls = 0
    failures = 0
    no_symbol = 0
    processed = 0
    for i, r in enumerate(rows, start=1):
        gene_symbol = r["gene_symbol"]
        if not gene_symbol:
            r["europepmc_n_articles"] = None
            no_symbol += 1
        else:
            count, from_cache = europepmc_count_for_gene(
                gene_symbol, organism_token, cache_dir, args.rate_limit
            )
            r["europepmc_n_articles"] = count
            if count is None:
                failures += 1
            elif from_cache:
                cache_hits += 1
            else:
                api_calls += 1
            processed += 1

        if i % args.progress_every == 0 or i == len(rows):
            elapsed = time.time() - t0
            rate = processed / elapsed if elapsed > 0 else 0
            print(
                f"  [{i}/{len(rows)}] no_symbol={no_symbol} "
                f"cache={cache_hits} api={api_calls} fail={failures} "
                f"elapsed={elapsed:.0f}s ({rate:.1f}/s queried)"
            )

    for r in rows:
        r["popularity_tier"] = compute_tier(
            r["uniprot_reviewed"],
            r["uniprot_pe_level"],
            r["europepmc_n_articles"],
            r["uniprot_n_references"],
        )

    write_tsv(rows, output_path)

    tier_counts = Counter(r["popularity_tier"] for r in rows)
    articles_values = sorted(
        r["europepmc_n_articles"]
        for r in rows
        if isinstance(r["europepmc_n_articles"], int)
    )
    median_articles = (
        articles_values[len(articles_values) // 2] if articles_values else 0
    )
    max_articles = max(articles_values) if articles_values else 0

    print(f"\nWrote {len(rows)} rows to {output_path}")
    print(
        "  tier counts: "
        + ", ".join(
            f"{t}={tier_counts.get(t, 0)}"
            for t in ["dark", "studied", "well_studied"]
        )
    )
    print(
        f"  europepmc_n_articles (over {len(articles_values)} queryable rows): "
        f"median={median_articles}, max={max_articles}; "
        f"no_symbol={no_symbol}, failures={failures}"
    )
    if failures:
        print(
            f"  NOTE: {failures} gene symbols have europepmc_n_articles blank — "
            "re-run to retry (cached hits are not re-fetched)."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
