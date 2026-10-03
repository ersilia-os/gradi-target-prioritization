"""PubTator3 literature counts, by TWO routes, because for bacteria they are different datasets.

    data/processed/studiedness/scratch/pubtator_geneid_counts.tsv   geneid -> n_pubs  (route A)
    data/processed/studiedness/scratch/pubtator_symbol_counts.tsv   symbol -> n_pubs  (route B)
    data/processed/studiedness/evidence/pubtator_<species>.tsv      per anchor protein
    data/processed/studiedness/evidence/pubtator_route_comparison.tsv   the measured verdict

PubTator3 (Wei et al., *Nucleic Acids Res* 2024) text-mines 36M PubMed abstracts and 6.3M full
texts -- two orders of magnitude more literature than the curated routes see. The question this
script answers is whether any of that reaches *bacterial* proteins.

THE FINDING THE WHOLE DESIGN TURNS ON: GeneID AND SYMBOL ARE NOT THE SAME DATASET
-----------------------------------------------------------------------------------
PubTator3 addresses genes two ways, and for bacteria they disagree by three orders of magnitude.
Measured against the live API, 2026-10-03:

    gene             by GeneID (bulk / API)       by symbol (@GENE_FTSZ, API)
    ftsZ   E. coli                  13  /      2                          198
    clpP   E. coli                   6  /      1                        3,130
    clpP   S. aureus                10  /      1                            -
    lpxC   E. coli                 156  /     51                           42
    rpoB   (all species)             -                                   2,587
    ABCB1  human                     -  / 51,222                            -

The GeneID column shows **bulk / API**, which disagree because the bulk file is the complete
annotation set and `search/` returns a narrower ranked document set. Route A uses the bulk.

**PubTator3 holds the bacterial literature but normalises bacterial mentions to species-agnostic
SYMBOL concepts, not to strain GeneIDs.** A bacterial GeneID count is therefore near-empty, while
a symbol count is real but pools every organism that has a gene of that name.

**THE BULK FILE IS GeneID-ONLY** -- verified over 4,000,000 rows, zero non-numeric concept ids.
So route A streams the 756 MB download and route B must go to the API. They are not substitutes.

WHY AN EARLIER VERSION OF docs/studiedness.md REJECTED THIS DATASET WRONGLY
-----------------------------------------------------------------------------
It probed the bulk file and saw `P-glycoprotein`, `S100`, `mTOR`, `CD8` -- concluding "visibly
human-centric, and bacterial gene normalisation is weak". The second half is right for the wrong
reason and the first half is an artifact: **`gene2pubtator3.gz` is sorted by PMID**, so its head
is whatever was published most recently, which skews human biomedical. Reproduced exactly on
2026-10-03: the first four rows are those same four names. **Never characterise a sorted file
from its head.**

THE CEILING ROUTE A CANNOT BEAT, AND IT IS NOT PUBTATOR'S FAULT
-----------------------------------------------------------------
**Only 42.1% of the 374,731 prokaryotic SwissProt entries carry a GeneID at all** (157,756). Any
GeneID-keyed count is blind to the other 58% before PubTator3 is even consulted -- which is also
why gene2pubmed scored only 2,071 of the control's 2,945 proteins against UniProt's full 2,945.

WHAT THIS SCRIPT DOES NOT DO
------------------------------
**It does not change the deliverable.** Both counts land in `evidence/`, and `transfer.py` carries
them as extra donor columns so the held-out E. coli control can score all routes on identical
folds. Promotion is a decision for that measurement, not for this script -- the `interpro2go` and
gene2pubmed precedents. **Three definitions mean three columns and never a `max()` across them**:
mixing definitions per protein is exactly what killed the 0-1 composite on 2026-09-22.

**Route B's count is SPECIES-AGNOSTIC AND CANNOT BE DONOR-SCOPED.** `@GENE_CLPP` = 3,130 pools
E. coli, S. aureus, human mitochondrial CLPP and plant homologs indiscriminately. Promoting it
without saying so would silently reintroduce "this Kp protein is well studied because its human
homolog is" -- the failure `docs/studiedness.md` 2b rejects for donor scope. Its column name
carries `_symbol_anyspecies` so the limitation travels with the number.

**Route B keys on a GENE SYMBOL, which CLAUDE.md forbids as a join key** -- and rightly, since a
symbol join onto proteins lost ~27% of v1's essentials. It is legitimate here only because the
symbol is a lookup key into an external text-mining index, never a key joining two protein
tables: the protein identity is still carried by the donor accession DIAMOND chose.

Run with the `gradi` env, after `transfer.py` (route B reads its donor symbols).
  python scripts/studiedness/pubtator.py
  python scripts/studiedness/pubtator.py --route symbol --limit 50
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import gzip
import json
import io
import sys
import tarfile
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "studiedness"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
LIT_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "literature"
BULK = REPO_ROOT / "data" / "source" / "ncbi" / "pubtator" / "gene2pubtator3.gz"
GENEID_COUNTS = SCRATCH_DIR / "pubtator_geneid_counts.tsv"
SYMBOL_COUNTS = SCRATCH_DIR / "pubtator_symbol_counts.tsv"
SPECIES_COUNTS = SCRATCH_DIR / "pubtator_species_counts.tsv"
TAXDUMP = REPO_ROOT / "data" / "source" / "ncbi" / "taxonomy" / "taxdump.tar.gz"

SPECIES = ("kpneumoniae", "ecoli", "saureus")

API = "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/"
# The API answers in 1.3-2.7s, so SERIAL IS NOT VIABLE at this size: 5,076 symbols would take
# ~3 hours and the whole cost is latency, not throughput. Three workers hold NCBI's <=3 req/s
# guidance while cutting the wall clock to ~1 hour. Do NOT raise it without an API key.
API_WORKERS = 3
# NCBI answered 429 at ~4.5 req/s with 3 unpaced workers. The gate below is on the SHARED
# timeline (see `_pace`), so this is the real ceiling regardless of worker count.
MIN_INTERVAL = 0.4      # <=2.5 req/s
RATE_LIMIT_BACKOFF = 5  # seconds, multiplied by the attempt number, on a 429
API_RETRIES = 5
API_TIMEOUT = 30
# Checkpoint often: the cache IS the resume record, and a 250-symbol gap is ~10 minutes of work
# to redo. Writing a few thousand rows costs milliseconds.
CHECKPOINT_EVERY = 100

# The bulk file must reproduce these, or the parse is wrong: E. coli ftsZ/clpP/lpxC, Sa clpP.
#
# **THESE ARE THE BULK FILE'S NUMBERS, NOT THE API'S, AND THE TWO DISAGREE.** Measured
# 2026-10-03 -- bulk 13 / 6 / 10 / 156 against the API's 2 / 1 / 1 / 51. The bulk file is the
# complete annotation set while the `search/` endpoint returns a narrower, ranked document set,
# so the bulk is the right source for "how many papers mention this gene" and the two must not
# be compared as if they were one quantity. Asserting the API's answers here (the first version
# of this file did) flags a correct parse as broken.
GENEID_SPOT = {"947587": 13, "945538": 6, "3919354": 10, "948874": 156}
# Route B's equivalents. Text mining is not frozen, so these are checked with a wide tolerance --
# the assertion is "thousands, not one", which is the whole claim route B rests on.
SYMBOL_SPOT_MIN = {"clpP": 500, "rpoB": 500, "ftsZ": 50}

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def split_ids(cell: object) -> list[str]:
    """Split a `;`-joined list field. Empty, NaN and whitespace-only all give []."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    return [x.strip() for x in str(cell).split(";") if x.strip()]


# --------------------------------------------------------------------- route A: GeneID


def geneid_universe() -> set[str]:
    """The GeneIDs worth counting: our three anchors plus every SwissProt donor."""
    universe: set[str] = set()
    for sp in SPECIES:
        gids = P.load_full(sp)["geneid"].apply(split_ids)
        flat = {g for lst in gids for g in lst}
        universe |= flat
        say(f"  {sp:<14} {int(gids.apply(bool).sum()):,} proteins with a GeneID "
            f"-> {len(flat):,} distinct")
    meta = LIT_DIR / "swissprot_meta.tsv.gz"
    if not meta.exists():
        sys.exit(f"FATAL missing {meta.relative_to(REPO_ROOT)} -- run fetch.py first")
    n = n_with = 0
    for chunk in pd.read_csv(meta, sep="\t", dtype=str, keep_default_na=False,
                             usecols=["GeneID"], chunksize=100_000):
        n += len(chunk)
        for cell in chunk["GeneID"]:
            ids = split_ids(cell)
            n_with += bool(ids)
            universe.update(ids)
    say(f"  swissprot      {n:,} entries, {n_with:,} with a GeneID ({100 * n_with / n:.1f}%)")
    say(f"  universe       {len(universe):,} distinct GeneIDs")
    return universe


def stream_geneid_counts(universe: set[str], refresh: bool) -> dict[str, int]:
    """One pass over gene2pubtator3.gz, counting distinct PMIDs per GeneID in the universe."""
    if GENEID_COUNTS.exists() and not refresh:
        df = pd.read_csv(GENEID_COUNTS, sep="\t", dtype={"geneid": str, "n_pubs": int})
        say(f"  cached {GENEID_COUNTS.relative_to(REPO_ROOT)}: {len(df):,} GeneIDs")
        return dict(zip(df["geneid"], df["n_pubs"]))
    if not BULK.exists():
        sys.exit(f"FATAL missing {BULK.relative_to(REPO_ROOT)} -- see its SOURCE.md")

    # A PMID can mention the same gene many times; the quantity is PAPERS, so count pairs once.
    seen: set[tuple[str, str]] = set()
    counts: dict[str, int] = {}
    rows = kept = 0
    t0 = time.time()
    with gzip.open(BULK, "rt", errors="replace") as fh:
        for line in fh:
            rows += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[1] != "Gene":
                continue
            pmid = parts[0]
            for gid in parts[2].split(";"):
                gid = gid.strip()
                if not gid or gid not in universe:
                    continue
                key = (gid, pmid)
                if key in seen:
                    continue
                seen.add(key)
                counts[gid] = counts.get(gid, 0) + 1
                kept += 1
            if VERBOSE and rows % 20_000_000 == 0:
                say(f"    {rows:,} rows  {len(counts):,} genes  {time.time() - t0:.0f}s")
    say(f"  streamed {rows:,} rows in {time.time() - t0:.0f}s -> "
        f"{len(counts):,} GeneIDs, {kept:,} gene-paper pairs")
    if rows < 10_000_000:
        sys.exit(f"FATAL only {rows:,} rows in {BULK.name} -- the download is truncated. "
                 "An exit code of 0 is not evidence of a complete file; check Content-Length.")

    bad = {g: (counts.get(g, 0), want) for g, want in GENEID_SPOT.items()
           if counts.get(g, 0) != want}
    if bad:
        say(f"  WARN spot checks disagree with the API: {bad}")
        say("       (text mining is re-run monthly, so a small drift is expected; a large one "
            "means the parse or the universe is wrong)")
    else:
        say(f"  spot checks  {len(GENEID_SPOT)}/{len(GENEID_SPOT)} match the API exactly")

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    (pd.DataFrame({"geneid": list(counts), "n_pubs": list(counts.values())})
       .sort_values("geneid").to_csv(GENEID_COUNTS, sep="\t", index=False))
    say(f"  wrote {GENEID_COUNTS.relative_to(REPO_ROOT)}")
    return counts


# --------------------------------------------------------------------- route B: symbol


def human_symbols() -> set[str]:
    """Gene symbols that are ALSO human gene symbols -- the symbol route's one real failure mode.

    `@GENE_CRP` returns 345,630 papers, essentially all of them about human **C-reactive
    protein**, not the bacterial cAMP receptor protein. Likewise `relA` (NF-kB p65, not the ppGpp
    synthetase), `pth` (parathyroid hormone, not peptidyl-tRNA hydrolase), `gpt`, `rho`, `bax`.

    **Measured, and smaller than it looks: 3-4% of donors per species**, and excluding them does
    NOT change the route's held-out score (0.4428 against 0.4432 with them in) -- they are diluting
    noise, not the source of its advantage. Their median inflation over the curated count is
    **426x**, against 1.9x for everything else, so they are also detectable without this list.
    """
    meta = LIT_DIR / "swissprot_meta.tsv.gz"
    if not meta.exists():
        return set()
    out: set[str] = set()
    for chunk in pd.read_csv(meta, sep="\t", dtype=str, keep_default_na=False,
                             usecols=["Gene Names (primary)", "Organism"], chunksize=100_000):
        hit = chunk[chunk["Organism"].str.startswith("Homo sapiens")]
        out |= {g.strip().lower() for g in hit["Gene Names (primary)"] if g.strip()}
    return out


def donor_symbols() -> list[str]:
    """Every gene symbol DIAMOND actually selected as a donor, across the three species."""
    syms: set[str] = set()
    for sp in SPECIES:
        path = EVIDENCE_DIR / f"transfer_{sp}.tsv"
        if not path.exists():
            sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- run transfer.py first")
        t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        syms |= {g.strip() for g in t["donor_gene"] if g.strip()}
    return sorted(syms)


_THROTTLE = threading.Semaphore(1)
_LAST_CALL = [0.0]


def _pace() -> None:
    """One global gate: no two requests start closer together than MIN_INTERVAL.

    Worker count alone does NOT bound the request rate -- three workers against a 0.4s response
    spike to 7/s and NCBI answers 429. The limiter has to be on the SHARED timeline, not per
    thread, which is why this is a module-level lock rather than a sleep in the worker.
    """
    with _THROTTLE:
        wait = MIN_INTERVAL - (time.monotonic() - _LAST_CALL[0])
        if wait > 0:
            time.sleep(wait)
        _LAST_CALL[0] = time.monotonic()


def species_of(taxids: set[str]) -> dict[str, str]:
    """Roll each donor taxid up to its SPECIES ancestor, from the NCBI taxdump.

    **MANDATORY BEFORE ANY SPECIES-SCOPED QUERY.** SwissProt files donors under STRAIN taxids --
    3,004 of K. pneumoniae's donors sit under `83333` (E. coli K-12) -- and PubTator3 annotates
    at species level, so a strain-level query undercounts ~8x: `@GENE_FTSZ AND @SPECIES_562`
    returns 163 where `@SPECIES_83333` returns 21, and clpP goes 151 -> 2. Measured 2026-10-03,
    and it resolves 230 of 230 donor taxids to 180 distinct species.
    """
    if not TAXDUMP.exists():
        sys.exit(f"FATAL missing {TAXDUMP.relative_to(REPO_ROOT)} -- needed to roll strain "
                 "taxids up to species; see its SOURCE.md")
    par: dict[str, str] = {}
    rank: dict[str, str] = {}
    with tarfile.open(TAXDUMP, "r:gz") as tf:
        fh = tf.extractfile("nodes.dmp")
        for raw in io.TextIOWrapper(fh, encoding="utf-8"):
            f = [x.strip() for x in raw.split("|")]
            par[f[0]], rank[f[0]] = f[1], f[2]
    out: dict[str, str] = {}
    for t in taxids:
        cur, hops = t, 0
        while cur and cur != "1" and hops < 40:
            if rank.get(cur) == "species":
                out[t] = cur
                break
            cur = par.get(cur)
            hops += 1
    say(f"  taxonomy       {len(out):,} of {len(taxids):,} donor taxids -> "
        f"{len(set(out.values())):,} distinct species")
    return out


def donor_pairs() -> list[tuple[str, str]]:
    """Distinct (gene symbol, SPECIES taxid) pairs our donors actually need."""
    rows = []
    for sp in SPECIES:
        path = EVIDENCE_DIR / f"transfer_{sp}.tsv"
        if not path.exists():
            sys.exit(f"FATAL missing {path.relative_to(REPO_ROOT)} -- run transfer.py first")
        t = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
        t = t[t["donor_ac"].str.strip().ne("") & t["donor_gene"].str.contains("[A-Za-z]")]
        rows.append(t[["donor_gene", "donor_taxid"]])
    d = pd.concat(rows, ignore_index=True)
    d["donor_gene"] = d["donor_gene"].str.strip()
    tax = species_of(set(d["donor_taxid"]) - {""})
    d["sp_taxid"] = d["donor_taxid"].map(tax)
    d = d[d["sp_taxid"].notna()]
    pairs = {(g, t) for g, t in zip(d["donor_gene"], d["sp_taxid"])
             if len(g) >= MIN_FREETEXT_SYMBOL}
    dropped = len({g for g in d["donor_gene"] if len(g) < MIN_FREETEXT_SYMBOL})
    if dropped:
        say(f"  skipping      {dropped} symbols shorter than {MIN_FREETEXT_SYMBOL} chars -- "
            "as free text they match everything")
    return sorted(pairs)


# Free-text symbols shorter than this are not queried: as a SEARCH STRING rather than a
# normalised concept, "S" or "N" matches essentially every paper. This is not an arbitrary cut
# like the `len <= 2` flag it replaces -- free-text matching genuinely requires a minimum
# specificity, and 3 characters is the shortest real bacterial gene symbol (there are 610 of
# them, against 42 symbols of length 1-2 in the whole donor set).
MIN_FREETEXT_SYMBOL = 3


def api_count_species(symbol: str, taxid: str) -> int | None:
    """Papers mentioning this gene symbol AND this species. The scoped count.

    **THE SYMBOL IS FREE TEXT, NOT A `@GENE_` CONCEPT, AND THAT IS THE WHOLE POINT.**
    PubTator3's gene normaliser has no concept for most bacterial genes, so `@GENE_FNBA` returns
    **0** for a *S. aureus* adhesin with 1,158 papers, and `@GENE_SECM` returns 0 against 1,958.
    Measured 2026-10-03 over the donor set: **31.5% of primary symbols return 0 as a concept**,
    and 446 of those donors carry >= 5 CURATED papers -- false zeros, which for a novelty axis is
    the one direction that must not fail.

    Free text does not depend on the vocabulary, and keeps the species separation that motivated
    scoping: `crp` reads 15,210 for E. coli against 236,188 for human.

    | gene | `@GENE_x` + species | free text + species |
    |---|---|---|
    | fnbA (Sa) | 0 | 1,063 |
    | secM (Ec) | 0 | 294 |
    | clpP (Ec) | 495 | 1,938 |
    | ftsZ (Ec) | 163 | 3,919 |

    **The trade is recall for precision, deliberately.** This is co-occurrence of a string and a
    species concept in one paper, so a B. subtilis FtsZ paper that mentions E. coli counts. A
    false zero misranks a well-studied protein as novel; mild over-counting is a monotone
    distortion. For a ranking, recall wins -- but the column's definition must say this, which is
    why it is named `..._symbol_species_cooccurrence` rather than anything implying resolution.
    """
    if len(symbol) < MIN_FREETEXT_SYMBOL:
        return None
    return _api(f"{symbol} AND @SPECIES_{taxid}", f"{symbol}/{taxid}")


def api_count(symbol: str) -> int | None:
    """Papers PubTator3 has for this gene SYMBOL, across all species. None on a hard failure."""
    return _api(f"@GENE_{symbol.upper()}", symbol)


def _api(query: str, label: str) -> int | None:
    """One search call with pacing, 429 backoff and payload validation."""
    url = API + "?" + urllib.parse.urlencode({"text": query})
    for attempt in range(API_RETRIES):
        try:
            _pace()
            with urllib.request.urlopen(url, timeout=API_TIMEOUT) as r:
                payload = r.read()
            # An HTTP 200 is not evidence of data -- a short or non-JSON body is a failure.
            if len(payload) < 2:
                raise ValueError(f"{len(payload)}-byte body")
            data = json.loads(payload)
            if "count" not in data:
                raise ValueError(f"no 'count' key: {sorted(data)[:5]}")
            return int(data["count"])
        except urllib.error.HTTPError as exc:
            # 429 is the server asking for a slower rate, not a broken request -- back off hard
            # and keep the symbol. Treating it as a plain failure silently drops real data.
            if exc.code == 429:
                time.sleep(RATE_LIMIT_BACKOFF * (attempt + 1))
                continue
            if attempt == API_RETRIES - 1:
                say(f"    FAILED {label}: {exc}")
                return None
            time.sleep(2 ** attempt)
        except (urllib.error.URLError, ValueError, json.JSONDecodeError, OSError) as exc:
            if attempt == API_RETRIES - 1:
                say(f"    FAILED {label}: {exc}")
                return None
            time.sleep(2 ** attempt)
    say(f"    FAILED {label}: still rate-limited after {API_RETRIES} attempts")
    return None


def fetch_symbol_counts(symbols: list[str], refresh: bool, limit: int | None) -> dict[str, int]:
    """Query the API once per symbol, resumably -- the cache IS the progress record."""
    cached: dict[str, int] = {}
    if SYMBOL_COUNTS.exists() and not refresh:
        df = pd.read_csv(SYMBOL_COUNTS, sep="\t", dtype={"symbol": str, "n_pubs": int})
        cached = dict(zip(df["symbol"], df["n_pubs"]))
        say(f"  cached {SYMBOL_COUNTS.relative_to(REPO_ROOT)}: {len(cached):,} symbols")

    todo = [s for s in symbols if s not in cached]
    if limit:
        todo = todo[:limit]
    say(f"  {len(symbols):,} donor symbols, {len(todo):,} to fetch "
        f"on {API_WORKERS} workers (~{len(todo) * 2.0 / API_WORKERS / 60:.0f} min at ~2s/call)")

    t0 = time.time()
    failed = 0
    done = 0
    with cf.ThreadPoolExecutor(max_workers=API_WORKERS) as pool:
        futures = {pool.submit(api_count, s): s for s in todo}
        for fut in cf.as_completed(futures):
            sym = futures[fut]
            n = fut.result()
            if n is None:
                failed += 1
            else:
                cached[sym] = n
            done += 1
            if done % CHECKPOINT_EVERY == 0 or done == len(todo):
                _write_symbols(cached)
                rate = done / max(time.time() - t0, 1e-9)
                eta = (len(todo) - done) / rate / 60
                say(f"    {done:,}/{len(todo):,}  {time.time() - t0:.0f}s  "
                    f"{rate:.1f}/s  ETA {eta:.0f}m  {failed} failed")
    _write_symbols(cached)

    if failed > len(todo) * 0.1 and todo:
        sys.exit(f"FATAL {failed:,} of {len(todo):,} symbol lookups failed (>10%) -- "
                 "the API is unhealthy; re-run to resume from the cache.")
    for sym, floor in SYMBOL_SPOT_MIN.items():
        got = cached.get(sym)
        if got is not None and got < floor:
            say(f"  WARN spot check {sym} = {got:,}, below the {floor:,} floor -- route B is "
                "supposed to be the one that FINDS bacterial literature")
    say(f"  spot checks  " + "  ".join(f"{s}={cached.get(s, '?'):,}" if isinstance(
        cached.get(s), int) else f"{s}=?" for s in SYMBOL_SPOT_MIN))
    return cached


def fetch_species_counts(pairs: list[tuple[str, str]], refresh: bool,
                         limit: int | None) -> dict[tuple[str, str], int]:
    """One query per (symbol, species) pair -- the SCOPED count, resumable like route B."""
    cached: dict[tuple[str, str], int] = {}
    if SPECIES_COUNTS.exists() and not refresh:
        df = pd.read_csv(SPECIES_COUNTS, sep="\t",
                         dtype={"symbol": str, "sp_taxid": str, "n_pubs": int})
        cached = {(r.symbol, r.sp_taxid): r.n_pubs for r in df.itertuples()}
        say(f"  cached {SPECIES_COUNTS.relative_to(REPO_ROOT)}: {len(cached):,} pairs")

    todo = [p for p in pairs if p not in cached]
    if limit:
        todo = todo[:limit]
    say(f"  {len(pairs):,} (symbol, species) pairs, {len(todo):,} to fetch on "
        f"{API_WORKERS} workers (~{len(todo) * 0.5 / 60:.0f} min)")

    t0, failed, done = time.time(), 0, 0
    with cf.ThreadPoolExecutor(max_workers=API_WORKERS) as pool:
        futures = {pool.submit(api_count_species, g, t): (g, t) for g, t in todo}
        for fut in cf.as_completed(futures):
            key = futures[fut]
            n = fut.result()
            if n is None:
                failed += 1
            else:
                cached[key] = n
            done += 1
            if done % CHECKPOINT_EVERY == 0 or done == len(todo):
                _write_species(cached)
                rate = done / max(time.time() - t0, 1e-9)
                say(f"    {done:,}/{len(todo):,}  {time.time() - t0:.0f}s  {rate:.1f}/s  "
                    f"ETA {(len(todo) - done) / rate / 60:.0f}m  {failed} failed")
    _write_species(cached)
    if failed > len(todo) * 0.1 and todo:
        sys.exit(f"FATAL {failed:,} of {len(todo):,} scoped lookups failed (>10%) -- "
                 "re-run to resume from the cache.")
    return cached


def _write_species(counts: dict[tuple[str, str], int]) -> None:
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"symbol": g, "sp_taxid": t, "n_pubs": n}
                  for (g, t), n in counts.items()]).sort_values(
        ["symbol", "sp_taxid"]).to_csv(SPECIES_COUNTS, sep="\t", index=False)


def _write_symbols(counts: dict[str, int]) -> None:
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    (pd.DataFrame({"symbol": list(counts), "n_pubs": list(counts.values())})
       .sort_values("symbol").to_csv(SYMBOL_COUNTS, sep="\t", index=False))


# ------------------------------------------------------------------------- per species


def per_species(species: str, geneid: dict[str, int], symbol: dict[str, int],
                human: set[str] | None = None) -> pd.DataFrame:
    """One row per anchor protein: what each route says about the protein ITSELF."""
    human = human or set()
    prot = P.load_full(species)[["uniprot_ac", "gene_name", "geneid"]].copy()

    def by_geneid(cell: object) -> int:
        vals = [geneid.get(g, 0) for g in split_ids(cell)]
        return max(vals) if vals else 0

    prot["n_pubs_pubtator_geneid"] = prot["geneid"].apply(by_geneid)
    prot["n_pubs_pubtator_symbol_anyspecies"] = prot["gene_name"].apply(
        lambda g: symbol.get(str(g).strip(), 0) if str(g).strip() else 0)
    # THE FLAG TRAVELS WITH THE NUMBER. A symbol count on a human homonym is measuring a
    # different protein entirely, so the column saying so ships beside it rather than in a doc.
    prot["symbol_is_human_homonym"] = prot["gene_name"].map(
        lambda g: str(g).strip().lower() in human if str(g).strip() else False)
    return prot[["uniprot_ac", "gene_name", "geneid",
                 "n_pubs_pubtator_geneid", "n_pubs_pubtator_symbol_anyspecies",
                 "symbol_is_human_homonym"]]


def comparison(frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """The verdict table: what each route reaches, and how coarse it is."""
    rows = []
    for sp, df in frames.items():
        own = pd.read_csv(EVIDENCE_DIR / f"gene2pubmed_{sp}.tsv", sep="\t", dtype=str,
                          keep_default_na=False) if (
            EVIDENCE_DIR / f"gene2pubmed_{sp}.tsv").exists() else None
        for route, col in (("pubtator_geneid", "n_pubs_pubtator_geneid"),
                           ("pubtator_symbol_anyspecies",
                            "n_pubs_pubtator_symbol_anyspecies")):
            v = df[col]
            nz = v[v > 0]
            tie = nz.value_counts()
            row = {"species": sp, "route": route, "n": len(v),
                   "n_ge1": int((v > 0).sum()),
                   "pct_ge1": round(100 * float((v > 0).mean()), 1),
                   "distinct": int(v.nunique()), "max": int(v.max()),
                   "median_nonzero": int(nz.median()) if len(nz) else 0,
                   "largest_tie": int(tie.iloc[0]) if len(tie) else 0,
                   "largest_tie_pct": round(100 * float(tie.iloc[0] / len(v)), 1)
                   if len(tie) else 0.0}
            if own is not None:
                u = pd.to_numeric(own["n_pubs_uniprot"], errors="coerce")
                m = df.merge(own[["uniprot_ac"]].assign(u=u), on="uniprot_ac", how="left")
                row["spearman_vs_uniprot"] = round(
                    float(m[col].corr(m["u"], method="spearman")), 4)
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--route", choices=("geneid", "symbol", "species", "both"),
                    default="both")
    ap.add_argument("--limit", type=int, help="route B only: fetch at most N new symbols")
    ap.add_argument("--refresh", action="store_true", help="ignore the caches and refetch")
    ap.add_argument("--dry-run", action="store_true", help="measure and print, write nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("studiedness/pubtator.py -- PubTator3 literature counts, by GeneID and by symbol")
    say(f"   in   {BULK.relative_to(REPO_ROOT)}  +  {API}")
    say(f"   out  evidence/pubtator_<species>.tsv, evidence/pubtator_route_comparison.tsv")
    say(f"   route: {args.route}   species: {', '.join(args.species)}")
    rule("=")

    geneid: dict[str, int] = {}
    symbol: dict[str, int] = {}

    if args.route in ("geneid", "both"):
        say("\n[route A -- GeneID, from the bulk file]")
        geneid = stream_geneid_counts(geneid_universe(), args.refresh)

    if args.route in ("species", "both"):
        say("\n[route C -- gene symbol SCOPED TO THE DONOR'S SPECIES, from the API]")
        say("  free-text symbol AND @SPECIES_<taxid>, at SPECIES rank (strains rolled up).")
        say("  NOT @GENE_<symbol>: that concept is missing for most bacterial genes -- it reads")
        say("  0 for fnbA, which has 1,158 papers. See api_count_species.__doc__.")
        fetch_species_counts(donor_pairs(), args.refresh, args.limit)

    if args.route in ("symbol", "both"):
        say("\n[route B -- gene symbol, from the API]")
        say("  NOTE this count is SPECIES-AGNOSTIC and cannot be donor-scoped; see the docstring")
        symbol = fetch_symbol_counts(donor_symbols(), args.refresh, args.limit)

    if args.dry_run:
        rule()
        say("--dry-run: counts built, per-species tables not written")
        return

    say("")
    human = human_symbols() if symbol else set()
    if symbol:
        say(f"  {len(human):,} human gene symbols loaded, to flag homonym collisions")
    frames = {}
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    for sp in args.species:
        df = per_species(sp, geneid, symbol, human)
        frames[sp] = df
        path = EVIDENCE_DIR / f"pubtator_{sp}.tsv"
        df.to_csv(path, sep="\t", index=False)
        say(f"  {sp:<14} geneid >=1: {int((df.n_pubs_pubtator_geneid > 0).sum()):>5,}  "
            f"symbol >=1: {int((df.n_pubs_pubtator_symbol_anyspecies > 0).sum()):>5,}  "
            f"-> {path.name}")

    comp = comparison(frames)
    comp.to_csv(EVIDENCE_DIR / "pubtator_route_comparison.tsv", sep="\t", index=False)
    rule()
    say(comp.to_string(index=False))
    rule("=")


if __name__ == "__main__":
    main()
