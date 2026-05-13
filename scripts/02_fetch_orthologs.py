"""Resolve OrthoDB ortholog sets for the K. pneumoniae HS11286 proteome.

Background
----------
The HS11286 reference proteome (UP000007841) is ~99.9% TrEMBL. UniProt only
exposes OrthoDB cross-references for SwissProt-reviewed entries, so a direct
``xref_orthodb`` query on the Kp proteome returns zero hits.

Pivot strategy
--------------
For each unique Kp gene symbol, look up a SwissProt entry of the same gene
symbol (any species) that has an OrthoDB cross-reference, take that
OrthoDB group id, and reverse-query UniProt with ``xref:orthodb-<group>`` to
enumerate every UniProt entry in the group across all species.

OrthoDB exposes one group per protein, typically at taxonomic level 2
(Bacteria) for bacterial query proteins, so the resulting ortholog sets
cover every bacterial species OrthoDB tracks. Cross-kingdom (e.g. human)
orthologs are NOT recovered through this path — that would require a
different ortholog source (eggNOG / KEGG OC).

Output
------
``data/processed/kp_orthodb_orthologs.tsv`` with one row per
``(kp_uniprot, ortholog_uniprot, ortholog_species, ortholog_taxid,
   orthodb_group_id, orthodb_level, truncated)``.

The Kp protein itself is emitted as the first "ortholog" row of its own group
so downstream joins can treat direct/ortholog evidence with the same shape.

Caches under ``tmp/``:
- ``tmp/orthodb_pivot/<symbol>.tsv`` — SwissProt lookup per Kp gene symbol
- ``tmp/orthodb_cache/<group_id>.tsv`` — full member list per OrthoDB group
"""

from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

from tenacity import retry, stop_after_attempt, wait_exponential

from _common import DEFAULT_ORGANISM, ORGANISMS


UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"

PROCESSED_DIR = Path("data/processed")
PIVOT_CACHE_DIR = Path("tmp/orthodb_pivot")
GROUP_CACHE_DIR = Path("tmp/orthodb_cache")

# OrthoDB group id format used by UniProt: "<digits>at<taxlevel>"
ORTHODB_TOKEN_RX = re.compile(r"^(\d+)at(\d+)$")

# Mirror the gene-symbol parsing from src/anchor.py without importing the
# package (scripts/ live outside it).
KP_LOCUS_RX = re.compile(r"(KPHS_[0-9p]+)")
GENE_SYMBOL_RX = re.compile(r"^[a-z][a-zA-Z0-9_]{2,5}$")


def parse_kp_proteome_tsv(body: str) -> list[dict[str, str]]:
    """Return one row per Kp protein with uniprot, gene_symbol, locus_tag."""
    rows: list[dict[str, str]] = []
    lines = body.splitlines()
    if not lines:
        return rows
    header = lines[0].split("\t")
    idx = {col: i for i, col in enumerate(header)}
    for line in lines[1:]:
        if not line:
            continue
        cells = line.split("\t")
        accession = cells[idx["Entry"]] if "Entry" in idx else ""
        names_raw = cells[idx["Gene Names"]] if "Gene Names" in idx else ""
        tokens = [tok for tok in (names_raw or "").split() if tok]
        locus = next(
            (m.group(1) for tok in tokens if (m := KP_LOCUS_RX.search(tok))), ""
        )
        symbol = ""
        for tok in tokens:
            if KP_LOCUS_RX.search(tok):
                continue
            if GENE_SYMBOL_RX.match(tok):
                symbol = tok
                break
        rows.append(
            {"uniprot": accession, "gene_symbol": symbol, "locus_tag": locus}
        )
    return rows


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, max=60))
def _http_get(url: str, timeout: int = 600) -> tuple[str, dict[str, str]]:
    req = urllib.request.Request(url, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        body = resp.read().decode("utf-8")
        headers = {k.lower(): v for k, v in resp.getheaders()}
    return body, headers


def _parse_next_link(link_header: str) -> Optional[str]:
    """Extract the ``next`` URL from a UniProt ``Link:`` header (RFC 5988)."""
    if not link_header:
        return None
    for part in link_header.split(","):
        seg = part.strip()
        if 'rel="next"' in seg:
            m = re.search(r"<([^>]+)>", seg)
            if m:
                return m.group(1)
    return None


def fetch_kp_anchor(proteome_id: str) -> list[dict[str, str]]:
    """Return one row per Kp protein with uniprot, gene_symbol, locus_tag."""
    params = urllib.parse.urlencode(
        {
            "query": f"proteome:{proteome_id}",
            "format": "tsv",
            "fields": "accession,gene_names",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/stream?{params}"
    body, _ = _http_get(url)
    return parse_kp_proteome_tsv(body)


def resolve_orthodb_group_for_symbol(symbol: str) -> Optional[tuple[str, str]]:
    """Find an OrthoDB group id for a gene symbol via a SwissProt lookup.

    Returns ``(group_id, level)`` (e.g. ``("9805706at2", "2")``) or None.
    Prefers bacterial hits (orthologs at level 2) when multiple are returned.
    """
    cache_path = PIVOT_CACHE_DIR / f"{_safe(symbol)}.tsv"
    if cache_path.exists() and cache_path.stat().st_size > 0:
        body = cache_path.read_text(encoding="utf-8")
    else:
        params = urllib.parse.urlencode(
            {
                "query": f"gene:{symbol} AND reviewed:true AND database:orthodb",
                "format": "tsv",
                "fields": "accession,organism_name,organism_id,xref_orthodb",
                "size": "10",
            }
        )
        body, _ = _http_get(f"{UNIPROT_SEARCH_URL}?{params}")
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(body, encoding="utf-8")
        time.sleep(0.1)

    lines = [ln for ln in body.splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    try:
        idx_ortho = header.index("OrthoDB")
    except ValueError:
        return None

    # Collect every (group_id, level) candidate from every returned row.
    candidates: list[tuple[str, str]] = []
    for line in lines[1:]:
        cells = line.split("\t")
        if idx_ortho >= len(cells):
            continue
        ortho_raw = cells[idx_ortho]
        for tok in ortho_raw.split(";"):
            tok = tok.strip()
            m = ORTHODB_TOKEN_RX.match(tok)
            if m:
                candidates.append((m.group(0), m.group(2)))

    if not candidates:
        return None

    # Prefer the broadest bacterial level (taxid 2) when multiple are present —
    # it'll yield the richest ortholog set across bacterial species.
    for group_id, level in candidates:
        if level == "2":
            return group_id, level
    return candidates[0]


def fetch_orthodb_group_members(
    group_id: str, max_members: int = 5000
) -> tuple[list[dict[str, str]], bool]:
    """Reverse-query UniProt for every entry in an OrthoDB group.

    Returns ``(rows, truncated)`` where each row has accession, organism_name,
    organism_id, length. Caches the final concatenated TSV under
    ``tmp/orthodb_cache/<group_id>.tsv``.
    """
    cache_path = GROUP_CACHE_DIR / f"{group_id}.tsv"
    truncated = False
    if cache_path.exists() and cache_path.stat().st_size > 0:
        body = cache_path.read_text(encoding="utf-8")
    else:
        params = urllib.parse.urlencode(
            {
                "query": f"xref:orthodb-{group_id}",
                "format": "tsv",
                "fields": "accession,organism_name,organism_id,length",
                "size": "500",
            }
        )
        url = f"{UNIPROT_SEARCH_URL}?{params}"
        all_lines: list[str] = []
        header_line: Optional[str] = None
        n_collected = 0
        while url:
            body, headers = _http_get(url)
            page_lines = body.splitlines()
            if not page_lines:
                break
            if header_line is None:
                header_line = page_lines[0]
            page_data = page_lines[1:]
            for ln in page_data:
                if n_collected >= max_members:
                    truncated = True
                    break
                all_lines.append(ln)
                n_collected += 1
            if truncated:
                break
            url = _parse_next_link(headers.get("link", ""))
            time.sleep(0.1)
        body = (header_line or "") + "\n" + "\n".join(all_lines) + "\n"
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        cache_path.write_text(body, encoding="utf-8")

    rows: list[dict[str, str]] = []
    lines = [ln for ln in body.splitlines() if ln.strip()]
    if len(lines) < 2:
        return rows, truncated
    header = lines[0].split("\t")
    idx = {col: i for i, col in enumerate(header)}
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) < len(header):
            cells = cells + [""] * (len(header) - len(cells))
        rows.append(
            {
                "uniprot": cells[idx.get("Entry", 0)],
                "organism_name": cells[idx["Organism"]] if "Organism" in idx else "",
                "organism_id": cells[idx["Organism (ID)"]]
                if "Organism (ID)" in idx
                else "",
                "length": cells[idx["Length"]] if "Length" in idx else "",
            }
        )
    return rows, truncated


def _safe(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", s)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
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
        help=f"Destination TSV. Defaults to {PROCESSED_DIR}/<organism>_orthodb_orthologs.tsv.",
    )
    parser.add_argument(
        "--max-members-per-group",
        type=int,
        default=5000,
        help="Cap each OrthoDB group at this many UniProt entries (default: %(default)s).",
    )
    parser.add_argument(
        "--limit-symbols",
        type=int,
        default=None,
        help="Only process the first N gene symbols (for smoke-testing).",
    )
    args = parser.parse_args()

    organism = ORGANISMS[args.organism]
    proteome_id = args.proteome_id or organism["proteome_id"]
    output_path = args.output or PROCESSED_DIR / f"{organism['slug']}_orthodb_orthologs.tsv"

    print(f"Fetching {organism['label']} proteome ({proteome_id}) ...")
    kp_rows = fetch_kp_anchor(proteome_id)
    print(f"  parsed {len(kp_rows)} proteins; {_n_with_symbol(kp_rows)} have a gene symbol")

    symbol_to_kp: dict[str, list[dict[str, str]]] = {}
    for row in kp_rows:
        sym = row["gene_symbol"]
        if not sym:
            continue
        symbol_to_kp.setdefault(sym, []).append(row)
    symbols = sorted(symbol_to_kp)
    if args.limit_symbols:
        symbols = symbols[: args.limit_symbols]
    print(f"  {len(symbols)} unique gene symbols to resolve")

    PIVOT_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    GROUP_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    symbol_to_group: dict[str, tuple[str, str]] = {}
    n_resolved = 0
    for i, symbol in enumerate(symbols, start=1):
        try:
            res = resolve_orthodb_group_for_symbol(symbol)
        except Exception as e:
            print(f"  [{i}/{len(symbols)}] {symbol}: ERROR {e!r}", file=sys.stderr)
            continue
        if res is not None:
            symbol_to_group[symbol] = res
            n_resolved += 1
        if i % 100 == 0 or i == len(symbols):
            print(f"  [{i}/{len(symbols)}] resolved {n_resolved} OrthoDB groups so far")
    print(
        f"Resolved {n_resolved}/{len(symbols)} gene symbols to OrthoDB groups "
        f"({100.0 * n_resolved / max(1, len(symbols)):.1f}%)"
    )

    unique_groups = sorted({g for g, _ in symbol_to_group.values()})
    print(f"Expanding {len(unique_groups)} unique OrthoDB groups ...")
    group_members: dict[str, list[dict[str, str]]] = {}
    group_truncated: dict[str, bool] = {}
    for i, group_id in enumerate(unique_groups, start=1):
        try:
            rows, truncated = fetch_orthodb_group_members(
                group_id, max_members=args.max_members_per_group
            )
        except Exception as e:
            print(f"  [{i}/{len(unique_groups)}] {group_id}: ERROR {e!r}", file=sys.stderr)
            continue
        group_members[group_id] = rows
        group_truncated[group_id] = truncated
        if i % 50 == 0 or i == len(unique_groups):
            total = sum(len(v) for v in group_members.values())
            n_trunc = sum(1 for t in group_truncated.values() if t)
            print(
                f"  [{i}/{len(unique_groups)}] {total} ortholog rows accumulated "
                f"({n_trunc} groups truncated)"
            )

    print(f"Writing {output_path} ...")
    columns = [
        "kp_uniprot",
        "kp_gene_symbol",
        "kp_locus_tag",
        "ortholog_uniprot",
        "ortholog_organism_name",
        "ortholog_organism_id",
        "orthodb_group_id",
        "orthodb_level",
        "truncated",
    ]
    n_rows = 0
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for symbol, kp_proteins in symbol_to_kp.items():
            grp = symbol_to_group.get(symbol)
            for kp in kp_proteins:
                if grp is None:
                    f.write(
                        "\t".join(
                            [
                                kp["uniprot"],
                                symbol,
                                kp["locus_tag"] or "",
                                kp["uniprot"],
                                "Klebsiella pneumoniae HS11286",
                                "1125630",
                                "",
                                "",
                                "False",
                            ]
                        )
                        + "\n"
                    )
                    n_rows += 1
                    continue
                group_id, level = grp
                truncated = "True" if group_truncated.get(group_id, False) else "False"
                members = group_members.get(group_id, [])
                # Always emit the Kp protein itself first so downstream
                # direct-vs-ortholog joins have a canonical anchor.
                f.write(
                    "\t".join(
                        [
                            kp["uniprot"],
                            symbol,
                            kp["locus_tag"] or "",
                            kp["uniprot"],
                            "Klebsiella pneumoniae HS11286",
                            "1125630",
                            group_id,
                            level,
                            truncated,
                        ]
                    )
                    + "\n"
                )
                n_rows += 1
                for m in members:
                    if m["uniprot"] == kp["uniprot"]:
                        continue
                    f.write(
                        "\t".join(
                            [
                                kp["uniprot"],
                                symbol,
                                kp["locus_tag"] or "",
                                m["uniprot"],
                                m["organism_name"],
                                m["organism_id"],
                                group_id,
                                level,
                                truncated,
                            ]
                        )
                        + "\n"
                    )
                    n_rows += 1
        # Kp proteins without any gene symbol: emit a self-only row so they
        # still appear in the downstream join (with empty group).
        for kp in kp_rows:
            if kp["gene_symbol"]:
                continue
            f.write(
                "\t".join(
                    [
                        kp["uniprot"],
                        "",
                        kp["locus_tag"] or "",
                        kp["uniprot"],
                        "Klebsiella pneumoniae HS11286",
                        "1125630",
                        "",
                        "",
                        "False",
                    ]
                )
                + "\n"
            )
            n_rows += 1

    print(f"Wrote {n_rows} rows to {output_path}")
    n_kp_with_group = sum(
        1 for row in kp_rows if row["gene_symbol"] and row["gene_symbol"] in symbol_to_group
    )
    print(
        f"Kp proteins linked to an OrthoDB group: {n_kp_with_group}/{len(kp_rows)} "
        f"({100.0 * n_kp_with_group / max(1, len(kp_rows)):.1f}%)"
    )
    unique_orthologs = {
        m["uniprot"]
        for members in group_members.values()
        for m in members
    }
    print(f"Unique ortholog UniProt accessions: {len(unique_orthologs)}")
    return 0


def _n_with_symbol(kp_rows: list[dict[str, str]]) -> int:
    return sum(1 for r in kp_rows if r["gene_symbol"])


if __name__ == "__main__":
    sys.exit(main())
