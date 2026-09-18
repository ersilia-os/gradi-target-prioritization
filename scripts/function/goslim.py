"""GO slim -- the second broad functional classification.

Stage 02, scheme 2. COG (`function/cog.py`) answers "which orthologous group"; this answers
"what does it do", in Gene Ontology terms collapsed to **goslim_prokaryote** (97 terms, purpose-built
for bacteria).

Two tiers, both real annotation, every row labelled in `goslim_source`:

    curated   UniProt's own GO, from stage 00's annotation layer
    eggnog    eggNOG-mapper v2 (`function/eggnog.py`) -- orthology-based transfer

Why eggNOG-mapper and nothing else
----------------------------------
COG tops out at 79.1 / 84.5 / 73.1%, and that is COG's ceiling rather than a bug: NCBI's own
curators reach 81.6% on E. coli. The proteins it misses are short, uncharacterised, and for Kp/Sa
roughly 60% carry no Pfam and no InterPro either.

So a second scheme only helps if it uses a *different kind of evidence*. Measured, twice, that the
domain route does not:

    interpro2go as a GO tier               filled 5 proteins out of 13,020
    InterPro live API vs UniProt's xref    0 of 15 GO-less Kp proteins would gain a GO term

UniProt's electronic GO is already InterPro2GO-derived, so InterProScan, Pfam+pfam2go and the
InterPro API all re-derive what stage 00 already has. **Orthology is the only route left**, which
means eggNOG-mapper -- the field standard, and the tool v1 identified and never built.

Rejected: a nearest-neighbour transfer in ESM-C embedding space. It reached 100% coverage, but a
hand-rolled k-NN is not a well-established predictor, and its honest accuracy was far below its
headline -- leave-one-out on E. coli gave MF 86.8%, but reweighted to the donor distances that
actually occurred it was **MF 68% / BP 53% / CC 66%**. Replaced deliberately, accepting lower
coverage in exchange for a citable method. Do not reintroduce it.

Coverage is therefore NOT 100%, and that is the point. A protein no established tool can annotate
gets an empty row, not a guess.

One primary term per aspect
---------------------------
GO is multi-label across three aspects, so each aspect gets exactly one chosen term plus the full
set, the same shape COG has. The choice is an explicit total order, so a tie can never leave it
empty:

    deepest in the GO DAG (most specific)  ->  rarest in this proteome (most informative)  ->  lowest GO id

Run with the `gradi` env, after `function/eggnog.py`:
    python scripts/function/goslim.py
    python scripts/function/goslim.py --no-eggnog       # curated tier only
    python scripts/function/goslim.py --refresh
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

RAW_DIR = REPO_ROOT / "data" / "source" / "go"
OUT_DIR = REPO_ROOT / "data" / "processed" / "function"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
DEFAULT_SPECIES = ("kpneumoniae", "ecoli", "saureus")

GO_BASIC = "https://current.geneontology.org/ontology/go-basic.obo"
GO_SLIM = "https://current.geneontology.org/ontology/subsets/goslim_prokaryote.obo"
INTERPRO2GO = "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go"

ASPECTS = {"molecular_function": "mf", "biological_process": "bp", "cellular_component": "cc"}

OUT_COLUMNS = (
    ["uniprot_ac"]
    + [f"goslim_{a}{suf}" for a in ("mf", "bp", "cc") for suf in ("", "_name")]
    + [f"goslim_{a}_all" for a in ("mf", "bp", "cc")]
    + ["goslim_source", "eggnog_og", "eggnog_evalue", "eggnog_score"]
)

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------- resources

def fetch(url: str, refresh: bool) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    path = RAW_DIR / Path(url).name
    if path.exists() and not refresh:
        say(f"  [cache] {path.name:<26} {path.stat().st_size / 1e6:>7.1f} MB")
        return path
    r = requests.get(url, timeout=600)
    r.raise_for_status()
    if len(r.content) < 1000:
        sys.exit(f"refusing a {len(r.content)}-byte payload from {url} -- that is not the file")
    path.write_bytes(r.content)
    say(f"  [fetch] {path.name:<26} {path.stat().st_size / 1e6:>7.1f} MB")
    return path


def write_source_md(paths: dict[str, Path]) -> None:
    (RAW_DIR / "SOURCE.md").write_text(
        "# Gene Ontology resources for stage 02 (GO slim)\n\n"
        f"Fetched {datetime.now(timezone.utc).isoformat(timespec='seconds')} "
        "by scripts/function/goslim.py\n\n"
        "| file | bytes | url |\n|---|---|---|\n"
        + "".join(f"| `{p.name}` | {p.stat().st_size:,} | {u} |\n" for u, p in paths.items()),
        encoding="utf-8",
    )


def load_interpro2go(path: Path) -> dict[str, list[str]]:
    """`InterPro:IPR000003 ... > GO:DNA binding ; GO:0003677` -> {IPR: [GO, ...]}."""
    m: dict[str, list[str]] = defaultdict(list)
    pat = re.compile(r"^InterPro:(IPR\d+)\s.*;\s*(GO:\d{7})\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        hit = pat.match(line)
        if hit:
            m[hit.group(1)].append(hit.group(2))
    return dict(m)


# ---------------- slim mapping

def build_slim_mapper(obo: Path, slim_obo: Path):
    """Return (map_terms, term_meta): GO ids -> per-aspect slim ids, and id -> (name, aspect)."""
    from goatools.mapslim import mapslim
    from goatools.obo_parser import GODag

    dag = GODag(str(obo), prt=None)
    slim = GODag(str(slim_obo), prt=None)
    cache: dict[str, tuple[str, ...]] = {}

    def to_slim(go_id: str) -> tuple[str, ...]:
        if go_id in cache:
            return cache[go_id]
        try:
            direct, _covered = mapslim(go_id, dag, slim)
        except (KeyError, ValueError):
            direct = set()          # obsolete or retired id -- UniProt xrefs lag the GO release
        cache[go_id] = tuple(sorted(direct))
        return cache[go_id]

    def map_terms(go_ids) -> dict[str, list[str]]:
        out: dict[str, set[str]] = {a: set() for a in ASPECTS.values()}
        for g in go_ids:
            for s in to_slim(g):
                node = dag.get(s)
                if node is not None and node.namespace in ASPECTS:
                    out[ASPECTS[node.namespace]].add(s)
        return {a: sorted(v) for a, v in out.items()}

    meta = {t: (n.name, n.namespace, n.depth) for t, n in dag.items()}
    slim_ids = {n.id for n in slim.values()}   # GODag keys include alt_ids; .id dedupes
    return map_terms, meta, dag, slim_ids


def pick_primary(terms: list[str], meta: dict, freq: Counter) -> str:
    """Most specific wins; ties to the rarest term in this proteome; then the lowest GO id."""
    if not terms:
        return ""
    return sorted(terms, key=lambda t: (-meta.get(t, ("", "", 0))[2], freq[t], t))[0]


# ---------------- tiers

def split_go(field: str) -> list[str]:
    return [t for t in re.findall(r"GO:\d{7}", field or "")]


def build_species(species: str, map_terms) -> pd.DataFrame:
    """Tier 1. One row per protein, slim sets per aspect, with a source label."""
    ann = P.load_annotation(species)[["uniprot_ac", "go_id"]]
    base = P.load(species)[["uniprot_ac"]].merge(ann, on="uniprot_ac", how="left").fillna("")

    rows = []
    for ac, go_field in zip(base["uniprot_ac"], base["go_id"]):
        terms = split_go(go_field)
        source = "curated" if terms else ""
        slim = map_terms(sorted(set(terms))) if terms else {a: [] for a in ASPECTS.values()}
        # A protein can have GO that maps to no slim term at all -- that is not a fill.
        if source and not any(slim.values()):
            source = ""
        rows.append({"uniprot_ac": ac, "goslim_source": source,
                     **{f"_{a}": slim[a] for a in ASPECTS.values()}})
    return pd.DataFrame(rows)


def fill_from_eggnog(frames: dict[str, pd.DataFrame], map_terms) -> tuple[int, int]:
    """Tier 2: GO from eggNOG-mapper, for proteins UniProt does not curate.

    Reads `eggnog_<species>.tsv` written by `function/eggnog.py`. Nothing is computed here --
    this is a join onto a real tool's output, which is the whole point.
    """
    filled = missing = 0
    for sp, df in frames.items():
        need = df.index[df["goslim_source"] == ""].tolist()
        missing += len(need)
        path = EVIDENCE_DIR / f"eggnog_{sp}.tsv"
        if not path.exists():
            say(f"  {sp:<14} {len(need):>5} unfilled, but eggnog_{sp}.tsv is absent -- "
                "run scripts/function/eggnog.py")
            continue
        eg = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False).set_index("uniprot_ac")
        got = 0
        for i in need:
            ac = df.at[i, "uniprot_ac"]
            if ac not in eg.index:
                continue
            gos = split_go(eg.at[ac, "gos"])
            if not gos:
                continue
            slim = map_terms(sorted(set(gos)))
            if not any(slim.values()):
                continue          # eggNOG had GO, but none of it maps into the slim
            for a in ASPECTS.values():
                df.at[i, f"_{a}"] = slim[a]
            df.at[i, "goslim_source"] = "eggnog"
            df.at[i, "eggnog_og"] = eg.at[ac, "eggnog_ogs"].split(",")[0]
            df.at[i, "eggnog_evalue"] = eg.at[ac, "evalue"]
            df.at[i, "eggnog_score"] = eg.at[ac, "score"]
            got += 1
        filled += got
        say(f"  {sp:<14} {got:>5} of {len(need):,} filled from eggNOG-mapper")
    return filled, missing


def finalise(df: pd.DataFrame, meta: dict) -> pd.DataFrame:
    """Pick one primary term per aspect and flatten to the output columns."""
    freq = {a: Counter(t for lst in df[f"_{a}"] for t in lst) for a in ASPECTS.values()}
    for a in ASPECTS.values():
        prim = [pick_primary(lst, meta, freq[a]) for lst in df[f"_{a}"]]
        df[f"goslim_{a}"] = prim
        df[f"goslim_{a}_name"] = [meta.get(t, ("", "", 0))[0] if t else "" for t in prim]
        df[f"goslim_{a}_all"] = [";".join(lst) for lst in df[f"_{a}"]]
    for c in ("eggnog_og", "eggnog_evalue", "eggnog_score"):
        if c not in df.columns:
            df[c] = ""
    return df.fillna("")[OUT_COLUMNS]


# ---------------- controls

def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if (sa or sb) else 1.0


def compare_interpro2go(map_terms, ip2go, species: str = "ecoli") -> dict | None:
    """The rejected alternative, kept because the number is the reason it was rejected.

    E. coli has curated GO for 90.3% of its proteome, so for every protein that has BOTH curated GO
    and InterPro we can compute the tier-2 answer and compare it to the curated one. This is the
    same trick as the COG control: use the one species that has a real answer to price the inference.
    """
    ann = P.load_annotation(species)[["uniprot_ac", "go_id", "interpro"]]
    both = ann[(ann["go_id"] != "") & (ann["interpro"] != "")]
    if both.empty:
        return None

    stats = {a: {"n": 0, "exact": 0, "jac": 0.0} for a in ASPECTS.values()}
    for go_field, ip_field in zip(both["go_id"], both["interpro"]):
        truth = map_terms(sorted(set(split_go(go_field))))
        pred_go = sorted({g for ipr in re.findall(r"IPR\d{6}", ip_field)
                          for g in ip2go.get(ipr, [])})
        pred = map_terms(pred_go) if pred_go else {a: [] for a in ASPECTS.values()}
        for a in ASPECTS.values():
            if not truth[a]:
                continue
            stats[a]["n"] += 1
            stats[a]["exact"] += int(bool(set(pred[a]) & set(truth[a])))
            stats[a]["jac"] += jaccard(pred[a], truth[a])
    say(f"  REJECTED ALTERNATIVE - interpro2go as a second tier, scored on {species}:")
    say(f"    {'aspect':<8} {'n':>7} {'any slim term shared':>22} {'mean Jaccard':>15}")
    for a in ASPECTS.values():
        s = stats[a]
        if s["n"]:
            say(f"    {a:<8} {s['n']:>7,} {100 * s['exact'] / s['n']:>21.1f}% "
                f"{s['jac'] / s['n']:>15.3f}")
    return stats


def control_eggnog(map_terms, species: str = "ecoli") -> dict | None:
    """Price tier 2 where tier 1 already knows the answer.

    E. coli has curated GO for ~90% of its proteome. For every protein that has BOTH curated GO and
    an eggNOG-mapper annotation, compute the tier-2 answer and compare it to the curated one --
    scored exactly the way the rejected ESM-C tier was, so the numbers are directly comparable.
    """
    path = EVIDENCE_DIR / f"eggnog_{species}.tsv"
    if not path.exists():
        say(f"  skipped: eggnog_{species}.tsv absent")
        return None
    eg = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    ann = P.load_annotation(species)[["uniprot_ac", "go_id"]]
    both = ann.merge(eg[["uniprot_ac", "gos"]], on="uniprot_ac", how="inner")
    both = both[(both["go_id"] != "") & (both["gos"] != "")]
    if both.empty:
        say("  no overlap between curated GO and eggNOG-mapper output")
        return None

    stats = {a: {"n": 0, "exact": 0, "jac": 0.0} for a in ASPECTS.values()}
    for truth_f, pred_f in zip(both["go_id"], both["gos"]):
        truth = map_terms(sorted(set(split_go(truth_f))))
        pred = map_terms(sorted(set(split_go(pred_f))))
        for a in ASPECTS.values():
            if not truth[a]:
                continue
            stats[a]["n"] += 1
            stats[a]["exact"] += int(bool(set(pred[a]) & set(truth[a])))
            stats[a]["jac"] += jaccard(pred[a], truth[a])
    say(f"  tier 2 (eggNOG-mapper), scored on {species} where curated GO exists "
        f"(n={len(both):,} proteins):")
    say(f"    {'aspect':<8} {'n':>7} {'any slim term shared':>22} {'mean Jaccard':>15}")
    for a in ASPECTS.values():
        st = stats[a]
        if st["n"]:
            say(f"    {a:<8} {st['n']:>7,} {100 * st['exact'] / st['n']:>21.1f}% "
                f"{st['jac'] / st['n']:>15.3f}")
    say("\n    for reference, the REJECTED ESM-C k-NN tier scored, reweighted to its real donor")
    say("    distances: MF 68%  BP 53%  CC 66%   (its naive headline was MF 86.8%)")
    return stats


# ---------------- main

def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES),
                    choices=list(DEFAULT_SPECIES))
    ap.add_argument("--no-eggnog", action="store_true", help="curated tier only")
    ap.add_argument("--no-control", action="store_true")
    ap.add_argument("--refresh", action="store_true", help="re-download GO and InterPro2GO")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    rule("=")
    say("STAGE 02 - functional annotation :: GO slim (goslim_prokaryote)")
    rule("=")
    say("  in       : data/processed/proteomes/evidence/annotation_<species>.tsv  (go_id, interpro)")
    say("             data/processed/function/eggnog_<species>.tsv                (tier 2)")
    say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/goslim_<species>.tsv  (+ evidence/ + scratch/)")
    say("  tiers    : curated (UniProt GO) -> eggnog (eggNOG-mapper v2)")
    say(f"  species  : {', '.join(args.species)}")
    if args.dry_run:
        say("\n  --dry-run: nothing fetched, nothing written.")
        for sp in args.species:
            say(f"    {sp:<14} -> goslim_{sp}.tsv")
        return
    say()

    rule()
    say("RESOURCES")
    rule()
    obo = fetch(GO_BASIC, args.refresh)
    slim_obo = fetch(GO_SLIM, args.refresh)
    ip2go_path = fetch(INTERPRO2GO, args.refresh)
    write_source_md({GO_BASIC: obo, GO_SLIM: slim_obo, INTERPRO2GO: ip2go_path})
    ip2go = load_interpro2go(ip2go_path)
    say(f"  interpro2go : {len(ip2go):,} InterPro entries carry a GO mapping")
    map_terms, meta, dag, slim_ids = build_slim_mapper(obo, slim_obo)
    say(f"  ontology    : {len(dag):,} GO terms; slim has {len(slim_ids):,} terms")
    say()

    rule()
    say("TIER 1 - curated GO")
    rule()
    t0 = time.time()
    frames = {}
    for sp in args.species:
        df = build_species(sp, map_terms)
        frames[sp] = df
        vc = df["goslim_source"].value_counts()
        n = len(df)
        say(f"  {sp:<14} {n:>6}  curated {vc.get('curated', 0):>6,} "
            f"({100 * vc.get('curated', 0) / n:>5.1f}%)   "
            f"remaining {(df['goslim_source'] == '').sum():>5,}")
    say(f"  ({time.time() - t0:.0f}s)")
    say()

    if not args.no_eggnog:
        rule()
        say("TIER 2 - eggNOG-mapper")
        rule()
        filled, missing = fill_from_eggnog(frames, map_terms)
        say(f"  filled {filled:,} of {missing:,} remaining")
        say()

    rule()
    say("OUTPUTS")
    rule()
    summary = []
    for sp, df in frames.items():
        out = finalise(df, meta)
        path = EVIDENCE_DIR / f"goslim_{sp}.tsv"
        out.to_csv(path, sep="\t", index=False)
        vc = out["goslim_source"].value_counts()
        n = len(out)
        summary.append({"species": sp, "n": n,
                        **{k: int(vc.get(k, 0)) for k in ("curated", "eggnog")},
                        "unfilled": int((out["goslim_source"] == "").sum()),
                        "pct_covered": round(100 * (out["goslim_source"] != "").sum() / n, 2),
                        "mf": int((out["goslim_mf"] != "").sum()),
                        "bp": int((out["goslim_bp"] != "").sum()),
                        "cc": int((out["goslim_cc"] != "").sum())})
        say(f"  goslim_{sp}.tsv{'':<{max(0, 14 - len(sp))}} {path.stat().st_size / 1e3:>8.1f} kB")
    man = pd.DataFrame(summary)
    man.to_csv(EVIDENCE_DIR / "goslim_manifest.tsv", sep="\t", index=False)

    # The slim vocabulary itself, so a consumer can resolve any id in the table without the ontology.
    terms = pd.DataFrame(
        [(t, meta[t][0], meta[t][1], meta[t][2]) for t in sorted(slim_ids) if t in meta],
        columns=["go_id", "name", "aspect", "depth"])
    terms.to_csv(EVIDENCE_DIR / "goslim_terms.tsv", sep="\t", index=False)
    say(f"  evidence/ + scratch/    goslim_manifest.tsv, goslim_terms.tsv ({len(terms)} slim terms)")
    say()

    rule()
    say("COVERAGE")
    rule()
    say(f"  {'species':<14} {'n':>6} {'curated':>9} {'eggnog':>10} {'no annotation':>15} "
        f"{'covered':>9}")
    for r in summary:
        say(f"  {r['species']:<14} {r['n']:>6} {r['curated']:>9,} "
            f"{r['eggnog']:>10,} {r['unfilled']:>15,} {r['pct_covered']:>8.1f}%")
    say(f"\n  {'species':<14} {'has MF':>10} {'has BP':>10} {'has CC':>10}")
    for r in summary:
        say(f"  {r['species']:<14} {r['mf']:>10,} {r['bp']:>10,} {r['cc']:>10,}")
    say("\n  Not every protein gets all three aspects -- a slim term may exist for MF and not BP.")
    say()

    if not args.no_control:
        rule()
        say("CONTROLS - price the inference where E. coli already knows the answer")
        rule()
        control_eggnog(map_terms)
        say()
        compare_interpro2go(map_terms, ip2go)
        say()

    short = [r for r in summary if r["unfilled"]]
    if short:
        rule()
        say("NOT ANNOTATED")
        rule()
        for r in short:
            say(f"  {r['species']:<14} {r['unfilled']:,} proteins have no GO from either tier "
                f"({r['pct_covered']:.1f}% covered)")
        say("\n  This is the expected outcome, not a failure. No established tool annotates these:")
        say("  they are short, uncharacterised, and mostly carry no Pfam or InterPro either. An")
        say("  empty row is the honest answer -- see docs/function.md.")
        say()
    rule("=")
    say("stage 02 (GO slim) complete.")
    rule("=")


if __name__ == "__main__":
    main()
