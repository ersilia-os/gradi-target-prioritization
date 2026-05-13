"""Annotate a downloaded proteome with a Clp-protease degradability score.

There is no proteome-wide measurement of Clp substrates for *Klebsiella
pneumoniae*. The best available signal is a composite of (a) rule-based degron
features computed directly on each Kp protein sequence and (b) experimental
evidence transferred from *E. coli* K-12 via gene-symbol orthology.

Rule-based degron features (Flynn 2003, Mol Cell 11:671-683)
------------------------------------------------------------
ClpXP/ClpAP recognize five recurring motif classes; we score each protein for:

* ``cterm_ssra_like`` — C-terminal CM1 (-AA / -LAA / -YALAA family). The
  ssrA-tagged C-terminus is the archetypal ClpXP degron.
* ``cterm_mua_like``  — C-terminal CM2 (MuA-like; basic + small-aliphatic).
* ``nterm_destabilizing`` — bulky hydrophobic L/F/Y/W at position 2 after
  iMet cleavage (canonical bacterial N-end rule, ClpS → ClpAP).
* ``nterm_nm1`` / ``nm2`` / ``nm3`` — Flynn's three N-terminal motif families.

These are conservative regex matches on the raw UniProt sequence. ``cterm_*``
inspects the last 5 residues; ``nterm_*`` inspects positions 2-8 (after the
initiator Met, which is normally cleaved when residue 2 is small).

Orthology transfer
------------------
Two small curated reference TSVs (see ``data/raw/clp_substrates/SOURCE.md``)
carry:

* ``flynn2003_ecoli_clp_substrates.tsv`` — gene_symbol → (clp_class,
  degron_class). Hand-curated from Flynn 2003 + Sauer/Baker reviews.
* ``nagar2021_ecoli_halflives.tsv`` — gene_symbol → (halflife_class,
  halflife_min). Subset of Nagar 2021 pulsed-SILAC half-lives.

For each Kp protein the script picks a canonical gene symbol, looks it up in
the E. coli reference proteome (already in ``data/raw/escherichia_coli_proteome.tsv``)
to get an ortholog accession, and then joins to the two reference tables.

Composite score
---------------
``clp_degradability_score`` ∈ [0, 1]:

    score = degron_feature_score                  (sequence-only baseline)
          + 0.40 if E. coli ortholog trapped by ClpXP/ClpAP (Flynn 2003)
          + 0.20 if E. coli ortholog half-life 'fast'  (Nagar 2021)
          + 0.10 if E. coli ortholog half-life 'slow'
        (capped at 1.0)

``clp_degradability_tier``: ``high`` (>= 0.50), ``medium`` (>= 0.25), ``low``.

Writes ``data/processed/<organism>_clp_degradability.tsv``. Wire-up into
``src/assemble.py`` is handled there; this script touches no other files.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

from _common import DEFAULT_ORGANISM, ORGANISMS


DEFAULT_RAW_DIR = Path("data/raw")
DEFAULT_PROCESSED_DIR = Path("data/processed")
DEFAULT_REFERENCE_DIR = Path("data/raw/clp_substrates")

ECOLI_ORGANISM_KEY = "ecoli"

# Locus-tag patterns to reject when picking the canonical gene symbol. Matches
# the conventions used by ``src/anchor.py``.
LOCUS_TAG_RXES = [
    re.compile(r"^KPHS_[0-9p]+$"),
    re.compile(r"^b\d{4}$"),
    re.compile(r"^JW\d+(?:\.\d+)?$"),
]
# Lowercase-initial gene symbols, 3-6 chars (matches anchor.py so the join
# against the assembled annotation table stays consistent).
GENE_SYMBOL_RX = re.compile(r"^[a-z][a-zA-Z0-9_]{2,5}$")

# --- Degron motif patterns ---------------------------------------------------

# CM1 (ssrA-family): C-terminus ends with two small alanines (-AA), often
# preceded by L/V/I/A and one more hydrophobic. Flynn 2003 reports CM1 as a
# generalization of the ssrA tag's terminal -LAA / -YALAA.
CM1_RX = re.compile(r"[YAFWLIVM]?[ALV]A[ALV]A$")

# CM2 (MuA-family): basic + small-aliphatic patch at the extreme C-terminus.
# Flynn 2003 archetype is MuA -RRKKAI; we accept short basic+aliphatic tails.
CM2_RX = re.compile(r"[RK]{2,}[A-Z]{0,2}[AVILMG]{1,3}$")

# N-end rule destabilizing residues at position 2 (bacterial primary
# destabilizing residues recognized by ClpS-ClpAP).
NEND_DESTABILIZING = set("LFYW")

# Flynn 2003 N-terminal motif families. Patterns kept loose on purpose; treat
# as enrichment signals, not strict predictors. Anchored at position 2 (the
# residue after iMet cleavage), inspecting up to position 8.
NM1_RX = re.compile(r"^M?[AILVMFW]{2,}")             # hydrophobic stretch
NM2_RX = re.compile(r"^M?[KR][AILVMFW]")              # basic + hydrophobic
NM3_RX = re.compile(r"^M?[TS][AILVMFW]")              # polar-small + hydrophobic


# --- TSV helpers -------------------------------------------------------------

def _is_locus_tag(token: str) -> bool:
    return any(rx.match(token) for rx in LOCUS_TAG_RXES)


def extract_gene_symbol(gene_names: str) -> str | None:
    """Return the first canonical-looking gene symbol from a UniProt
    ``Gene Names`` cell, or None if no token qualifies."""
    if not gene_names:
        return None
    for tok in gene_names.split():
        tok = tok.strip()
        if not tok or _is_locus_tag(tok):
            continue
        if GENE_SYMBOL_RX.match(tok):
            return tok
    return None


def _find_col(header: list[str], *candidates: str) -> int | None:
    idx = {col: i for i, col in enumerate(header)}
    for c in candidates:
        if c in idx:
            return idx[c]
    lower = {k.lower(): v for k, v in idx.items()}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    return None


def load_proteome_tsv(path: Path) -> list[dict]:
    """Read a UniProt proteome TSV (accession, gene_names, sequence). Skip
    rows missing accession or sequence."""
    with open(path, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        i_acc = _find_col(header, "Entry", "accession")
        i_names = _find_col(header, "Gene Names", "gene_names")
        i_seq = _find_col(header, "Sequence", "sequence")
        if i_acc is None or i_seq is None:
            raise RuntimeError(
                f"Proteome TSV {path} missing Entry or Sequence column; "
                f"header was {header!r}"
            )
        rows: list[dict] = []
        for line in f:
            if not line.strip():
                continue
            cells = line.rstrip("\n").split("\t")

            def cell(i: int | None) -> str:
                if i is None or i >= len(cells):
                    return ""
                return cells[i]

            accession = cell(i_acc).strip()
            sequence = cell(i_seq).strip()
            if not accession or not sequence:
                continue
            gene_names = cell(i_names).strip()
            rows.append(
                {
                    "accession": accession,
                    "gene_names": gene_names,
                    "gene_symbol": extract_gene_symbol(gene_names) or "",
                    "sequence": sequence,
                }
            )
    return rows


def load_reference_tsv(path: Path) -> list[dict]:
    """Load a curated reference TSV with a header row."""
    if not path.exists():
        return []
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return [dict(r) for r in reader]


# --- Degron computation ------------------------------------------------------

def compute_degron_features(sequence: str) -> dict:
    """Compute rule-based degron features for one protein sequence."""
    seq = sequence.upper()
    last5 = seq[-5:] if len(seq) >= 5 else seq
    last3 = seq[-3:] if len(seq) >= 3 else seq

    cterm_ssra_like = bool(CM1_RX.search(last5))
    cterm_mua_like = bool(CM2_RX.search(seq[-8:] if len(seq) >= 8 else seq))

    # N-end rule: residue at position 2 (1-indexed). The initiator Met (pos 1)
    # is cleaved when pos 2 is small (A, C, G, P, S, T, V); for bulky residues
    # the Met stays, but for prediction we look at pos 2 as the destabilizer.
    pos2 = seq[1] if len(seq) >= 2 else ""
    nterm_destabilizing = pos2 in NEND_DESTABILIZING

    # N-terminal motifs: inspect the head of the sequence.
    head = seq[:8]
    nterm_nm1 = bool(NM1_RX.search(head))
    nterm_nm2 = bool(NM2_RX.search(head))
    nterm_nm3 = bool(NM3_RX.search(head))

    features = {
        "cterm_last3": last3,
        "cterm_last5": last5,
        "cterm_ssra_like": cterm_ssra_like,
        "cterm_mua_like": cterm_mua_like,
        "nterm_pos2_residue": pos2,
        "nterm_destabilizing": nterm_destabilizing,
        "nterm_nm1": nterm_nm1,
        "nterm_nm2": nterm_nm2,
        "nterm_nm3": nterm_nm3,
    }

    # Weighted feature score in [0, 1]. The two best-characterized signals
    # (ssrA-like C-terminus, N-end rule) carry the most weight; secondary
    # N-terminal motifs are weak enrichment markers and weighted lightly.
    w_ssra = 0.40 if cterm_ssra_like else 0.0
    w_nend = 0.25 if nterm_destabilizing else 0.0
    w_mua = 0.15 if cterm_mua_like else 0.0
    w_nm = 0.05 * sum([nterm_nm1, nterm_nm2, nterm_nm3])
    features["degron_feature_count"] = sum(
        [cterm_ssra_like, cterm_mua_like, nterm_destabilizing,
         nterm_nm1, nterm_nm2, nterm_nm3]
    )
    features["degron_feature_score"] = min(1.0, round(w_ssra + w_nend + w_mua + w_nm, 4))
    return features


# --- Score composition -------------------------------------------------------

def compose_score(
    degron_feature_score: float,
    ecoli_clp_trapped: bool,
    ecoli_halflife_class: str | None,
) -> float:
    s = float(degron_feature_score)
    if ecoli_clp_trapped:
        s += 0.40
    if ecoli_halflife_class == "fast":
        s += 0.20
    elif ecoli_halflife_class == "slow":
        s += 0.10
    return round(min(1.0, s), 4)


def tier_of(score: float) -> str:
    if score >= 0.50:
        return "high"
    if score >= 0.25:
        return "medium"
    return "low"


# --- Orthology ---------------------------------------------------------------

def build_ecoli_symbol_index(ec_rows: list[dict]) -> dict[str, str]:
    """Lowercased E. coli gene_symbol → accession. First wins on collision."""
    out: dict[str, str] = {}
    for r in ec_rows:
        sym = (r.get("gene_symbol") or "").lower()
        if sym and sym not in out:
            out[sym] = r["accession"]
    return out


def build_flynn_index(rows: list[dict]) -> dict[str, dict]:
    return {(r["gene_symbol"] or "").lower(): r for r in rows if r.get("gene_symbol")}


def build_nagar_index(rows: list[dict]) -> dict[str, dict]:
    return {(r["gene_symbol"] or "").lower(): r for r in rows if r.get("gene_symbol")}


# --- Output ------------------------------------------------------------------

OUTPUT_COLUMNS = [
    "accession",
    "gene_names",
    "gene_symbol",
    "cterm_last5",
    "cterm_ssra_like",
    "cterm_mua_like",
    "nterm_pos2_residue",
    "nterm_destabilizing",
    "nterm_nm1",
    "nterm_nm2",
    "nterm_nm3",
    "degron_feature_count",
    "degron_feature_score",
    "ecoli_ortholog_uniprot",
    "ortholog_status",
    "ecoli_clp_trapped",
    "ecoli_clp_class",
    "ecoli_degron_class",
    "ecoli_halflife_class",
    "ecoli_halflife_min",
    "clp_degradability_score",
    "clp_degradability_tier",
]


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "True" if v else "False"
    return str(v)


def write_tsv(rows: list[dict], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(_fmt(r.get(c, "")) for c in OUTPUT_COLUMNS) + "\n")


# --- Main --------------------------------------------------------------------

def annotate(
    kp_rows: list[dict],
    ec_rows: list[dict],
    flynn_rows: list[dict],
    nagar_rows: list[dict],
) -> list[dict]:
    ec_index = build_ecoli_symbol_index(ec_rows)
    flynn = build_flynn_index(flynn_rows)
    nagar = build_nagar_index(nagar_rows)

    out: list[dict] = []
    for r in kp_rows:
        feats = compute_degron_features(r["sequence"])
        sym = (r["gene_symbol"] or "").lower()

        ec_acc = ec_index.get(sym, "")
        if not sym:
            ortholog_status = "no_symbol"
        elif not ec_acc:
            ortholog_status = "no_ecoli_match"
        else:
            ortholog_status = "by_symbol"

        f = flynn.get(sym) if ec_acc else None
        n = nagar.get(sym) if ec_acc else None

        ecoli_clp_trapped = bool(f)
        ecoli_clp_class = (f.get("clp_class") if f else "") or ""
        ecoli_degron_class = (f.get("degron_class") if f else "") or ""
        ecoli_halflife_class = (n.get("halflife_class") if n else "") or ""
        ecoli_halflife_min = (n.get("halflife_min") if n else "") or ""

        score = compose_score(
            feats["degron_feature_score"],
            ecoli_clp_trapped,
            ecoli_halflife_class or None,
        )
        tier = tier_of(score)

        out.append(
            {
                "accession": r["accession"],
                "gene_names": r["gene_names"],
                "gene_symbol": r["gene_symbol"],
                **feats,
                "ecoli_ortholog_uniprot": ec_acc,
                "ortholog_status": ortholog_status,
                "ecoli_clp_trapped": ecoli_clp_trapped,
                "ecoli_clp_class": ecoli_clp_class,
                "ecoli_degron_class": ecoli_degron_class,
                "ecoli_halflife_class": ecoli_halflife_class,
                "ecoli_halflife_min": ecoli_halflife_min,
                "clp_degradability_score": score,
                "clp_degradability_tier": tier,
            }
        )
    return out


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--organism",
        choices=sorted(ORGANISMS),
        default=DEFAULT_ORGANISM,
        help="Target organism shortname (default: %(default)s).",
    )
    parser.add_argument(
        "--ecoli-proteome",
        type=Path,
        default=None,
        help=(
            "E. coli K-12 reference proteome TSV (Entry / Gene Names / Sequence). "
            f"Defaults to {DEFAULT_RAW_DIR}/escherichia_coli_proteome.tsv."
        ),
    )
    parser.add_argument(
        "--target-proteome",
        type=Path,
        default=None,
        help=(
            "Target proteome TSV. Defaults to "
            f"{DEFAULT_RAW_DIR}/<organism>_proteome.tsv."
        ),
    )
    parser.add_argument(
        "--flynn",
        type=Path,
        default=DEFAULT_REFERENCE_DIR / "flynn2003_ecoli_clp_substrates.tsv",
        help="Curated Flynn 2003 substrate reference TSV.",
    )
    parser.add_argument(
        "--nagar",
        type=Path,
        default=DEFAULT_REFERENCE_DIR / "nagar2021_ecoli_halflives.tsv",
        help="Curated Nagar 2021 half-life reference TSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help=(
            "Destination TSV. Defaults to "
            f"{DEFAULT_PROCESSED_DIR}/<organism>_clp_degradability.tsv."
        ),
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Optional cap on target proteins processed (for debugging).",
    )
    args = parser.parse_args()

    organism = ORGANISMS[args.organism]
    target_path = args.target_proteome or (
        DEFAULT_RAW_DIR / f"{organism['slug']}_proteome.tsv"
    )
    ecoli_path = args.ecoli_proteome or (
        DEFAULT_RAW_DIR / f"{ORGANISMS[ECOLI_ORGANISM_KEY]['slug']}_proteome.tsv"
    )
    output_path = args.output or (
        DEFAULT_PROCESSED_DIR / f"{organism['slug']}_clp_degradability.tsv"
    )

    if not target_path.exists():
        print(f"Target proteome not found: {target_path}", file=sys.stderr)
        print(
            f"Run: python scripts/00_download_proteome.py --organism {args.organism}",
            file=sys.stderr,
        )
        return 1
    if not ecoli_path.exists():
        print(f"E. coli reference proteome not found: {ecoli_path}", file=sys.stderr)
        print(
            "Run: python scripts/00_download_proteome.py --organism ecoli",
            file=sys.stderr,
        )
        return 1

    print(f"Loading target proteome: {target_path}")
    kp_rows = load_proteome_tsv(target_path)
    if args.limit is not None:
        kp_rows = kp_rows[: args.limit]
    print(f"  {len(kp_rows)} target proteins")

    print(f"Loading E. coli reference proteome: {ecoli_path}")
    ec_rows = load_proteome_tsv(ecoli_path)
    print(f"  {len(ec_rows)} E. coli proteins")

    print(f"Loading Flynn 2003 reference: {args.flynn}")
    flynn_rows = load_reference_tsv(args.flynn)
    print(f"  {len(flynn_rows)} Clp substrates")

    print(f"Loading Nagar 2021 reference: {args.nagar}")
    nagar_rows = load_reference_tsv(args.nagar)
    print(f"  {len(nagar_rows)} half-life entries")

    rows = annotate(kp_rows, ec_rows, flynn_rows, nagar_rows)
    write_tsv(rows, output_path)

    tier_counts = Counter(r["clp_degradability_tier"] for r in rows)
    ortholog_counts = Counter(r["ortholog_status"] for r in rows)
    ssra_n = sum(1 for r in rows if r["cterm_ssra_like"])
    nend_n = sum(1 for r in rows if r["nterm_destabilizing"])
    trapped_n = sum(1 for r in rows if r["ecoli_clp_trapped"])
    halflife_fast_n = sum(1 for r in rows if r["ecoli_halflife_class"] == "fast")
    scores = [r["clp_degradability_score"] for r in rows]
    median_score = sorted(scores)[len(scores) // 2] if scores else 0.0
    max_score = max(scores) if scores else 0.0

    print(f"\nWrote {len(rows)} rows to {output_path}")
    print(
        "  tier counts: "
        + ", ".join(f"{t}={tier_counts.get(t, 0)}" for t in ["low", "medium", "high"])
    )
    print(
        "  ortholog status: "
        + ", ".join(
            f"{k}={ortholog_counts.get(k, 0)}"
            for k in ["by_symbol", "no_ecoli_match", "no_symbol"]
        )
    )
    print(
        f"  rule-based: cterm_ssra_like={ssra_n}, nterm_destabilizing={nend_n}"
    )
    print(
        f"  ortholog evidence: ecoli_clp_trapped={trapped_n}, "
        f"ecoli_halflife_class=fast: {halflife_fast_n}"
    )
    print(f"  score: median={median_score}, max={max_score}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
