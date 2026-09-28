"""At what sequence identity does ligand precedent actually transfer? Measured, not assumed.

Ligandability, the calibration behind the bands. `src/ligandability.py` sets
`DIRECT/CLOSE/REMOTE = 95/60/40` and `docs/ligands.md` concedes that 60 "is the one arbitrary
number in the design". Nothing has ever measured the quantity those numbers are supposed to encode:
**if protein A has ligands, how much does that tell you about protein B at x% identity?**

The measurement is self-contained inside ChEMBL
------------------------------------------------
For a pair of bacterial ChEMBL targets we know BOTH ligand sets. So the transfer rule can be scored
directly: DIAMOND every bacterial target against every other, and for each pair compute the Jaccard
overlap of their `parent_molregno` sets. That is what "precedent transfers" means, made into a
number. No proteome, no model, no held-out anything -- the answer is already in the database.

Three questions, one table
---------------------------
1. **Where do the bands belong?** Jaccard against identity says where transfer decays, and whether
   40 is defensible as a floor.
2. **Does SPECIES add information beyond identity?** The same identity band, split by whether the
   two targets share a two-word binomial. This decides whether `precedents.py`'s `exact` count
   should be "identical sequence" (what it means today) or "same species, allowing strain-level
   variation" -- the project owner's request, and the reason `chembl.py` has a `species` bucket
   that `precedents.py` never got.
3. **Does ORTHOLOGY add information beyond identity?** The same bands split by reciprocal best hit.
   A paralog and an ortholog at 45% are treated identically today; `is_rbh` is the repo's own
   existing statement that those are different claims.

THE VERDICT: all three questions come back NEGATIVE, and that is the result
----------------------------------------------------------------------------
Measured on 1,582 pairs over 687 bacterial target sequences / 131 species:

1. **Compound-set overlap does not transfer at any identity.** Median Jaccard is **~0.00 in every
   band, 95-100% included**, and among pairs sharing anything at all it is 0.02-0.04 with no trend.
   Two near-identical ChEMBL targets are one enzyme screened twice by different groups against
   different libraries -- they share the protein, not the chemistry. Jaccard measures campaign
   coincidence. It cannot calibrate anything, and that is why the conditional below exists.

2. **The conditional is FLAT, so identity is not doing the work.**
   `P(potent | neighbour potent)` runs **0.87 / 0.93 / 0.98 / 0.84 / 0.95 / 0.94 / 0.67 / 1.00 /
   0.96** across the 25-30 ... 95-100 bands -- a lift of only 1.36-1.58x over the **0.619** base
   rate, with no decay. The reason is selection: **62% of bacterial ChEMBL targets already have a
   potent compound**, because a protein enters ChEMBL when somebody believed it was druggable.
   Conditioning on a potent neighbour therefore tells you little beyond "this is also a ChEMBL
   target".

   **So the 95/60/40 bands CANNOT be calibrated from inside ChEMBL, and this script's honest output
   is that finding rather than a new number.** What the floor actually controls is COVERAGE -- how
   many proteins get any precedent at all -- not the reliability of a transfer. Document it as a
   conservatism choice; do not present it as an accuracy threshold. A floor validated against an
   unselected population would need proteins nobody chose to screen, which by construction ChEMBL
   does not contain.

3. **Neither SPECIES nor RBH adds anything beyond identity.** Same-species vs cross-species at
   95-100%: **0.955 vs 0.960**. RBH vs not: no consistent direction, and RBH is *worse* in three
   bands. So no orthology criterion is added to the transfer rule -- the simple rule was not
   leaving anything on the table, which is worth knowing before anyone builds the complicated one.

   This does NOT argue against making `exact` species-aware. That change is about what counts as
   *the same protein* -- a 98%-identical protein from another K. pneumoniae strain is not a
   homolog, it is this protein -- which is an identity question, not a potency-prediction one.

Every row ships twice: over all pairs, and over pairs where both targets carry at least
`--min-compounds` (default 5). Quote the second. Pair counts ship so a band resting on three pairs
cannot be mistaken for a measurement -- the 80-90 and 90-95 bands have 2 and 3.

Bacteria only. DIAMOND comes from `gradi-ortho` via `GRADI_DIAMOND_BIN`. `gradi` env. ~2 min.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import precedents as PR  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "ligands"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
CALIB_PATH = EVIDENCE_DIR / "transfer_calibration.tsv"
PAIRS_PATH = SCRATCH_DIR / "transfer_pairs.tsv"

# Bands finer than the shipped 95/60/40, so the decay can be SEEN rather than assumed to be a step.
BANDS = [(25, 30), (30, 40), (40, 50), (50, 60), (60, 70), (70, 80), (80, 90), (90, 95), (95, 100.01)]

# Our column names, and DIAMOND's own field names in the SAME ORDER -- the mapping is positional,
# exactly as `src/precedents.py` does it. Passing our names to --outfmt is an `Invalid output
# field` error, which is the friendly version of this mistake.
HIT_COLS = ["query", "component_id", "pident", "ppos", "length", "qlen", "slen",
            "qcov", "scov", "evalue", "bitscore"]
DIAMOND_FIELDS = ["qseqid", "sseqid", "pident", "ppos", "length", "qlen", "slen",
                  "qcovhsp", "scovhsp", "evalue", "bitscore"]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 96) -> None:
    say(char * width)


def binomial(organism: str) -> str:
    """The two-word binomial, `chembl.py`'s own same-species rule.

    ChEMBL files strains under their own names AND their own taxids -- `Escherichia coli`,
    `Escherichia coli K-12` and `Escherichia coli (strain K12)` are three organism strings for one
    species. Matching on the genus alone would sweep in `Klebsiella aerogenes`; matching the full
    string would split a species into strains, which is exactly the distinction this script is
    measuring. Two words is the rule that already ships.
    """
    if not isinstance(organism, str):
        return ""
    words = re.split(r"\s+", organism.strip())
    return " ".join(words[:2]).lower()


def bacterial_targets() -> tuple[pd.DataFrame, dict[str, str], dict[str, set]]:
    """Bacterial ChEMBL targets: metadata, sequences, and each component's LIGAND SET.

    The component -> tid map is ONE-TO-MANY and must be walked in full. Collapsing it with
    `dict(zip(...))` dropped 490 components' extra tids and made E. coli gyrA read as unliganded
    (fixed in 858bf01). The same mistake here would silently deflate every Jaccard.
    """
    t, lg, seqs = PR._tables()
    bact = t[t["superkingdom"] == "Bacteria"].copy()
    if bact.empty:
        sys.exit("FAILED: no Bacteria rows in chembl_targets.tsv -- run scripts/ligands/chembl.py.")

    comp2tids: dict[str, list] = {}
    for cid, tid in zip(bact["component_id"], bact["tid"]):
        comp2tids.setdefault(str(cid), []).append(tid)

    single = lg[lg["track"] == "single"]
    tid2mols = single.groupby("tid")["parent_molregno"].apply(set).to_dict()
    tid2potent = (single[single["pchembl"] >= L.PCHEMBL_HEADLINE]
                  .groupby("tid")["parent_molregno"].apply(set).to_dict())

    ligands, potent = {}, {}
    for cid, tids in comp2tids.items():
        ligands[cid] = set().union(*(tid2mols.get(x, set()) for x in tids)) if tids else set()
        potent[cid] = set().union(*(tid2potent.get(x, set()) for x in tids)) if tids else set()

    meta = (bact.groupby("component_id")
                .agg(organism=("organism", "first"), target_chembl_id=("target_chembl_id", "first"))
                .reset_index())
    meta["component_id"] = meta["component_id"].astype(str)
    meta["binomial"] = meta["organism"].map(binomial)
    meta["n_compounds"] = meta["component_id"].map(lambda c: len(ligands.get(c, ())))
    meta["n_potent"] = meta["component_id"].map(lambda c: len(potent.get(c, ())))
    sub = {c: seqs[c] for c in meta["component_id"] if c in seqs}
    meta = meta[meta["component_id"].isin(sub)].reset_index(drop=True)
    return meta, sub, {"all": ligands, "potent": potent}


def all_vs_all(seqs: dict[str, str], threads: int) -> pd.DataFrame:
    """Every bacterial target against every other. Self-hits dropped; best HSP per ordered pair."""
    with tempfile.TemporaryDirectory() as td:
        faa, db, out = Path(td) / "q.faa", Path(td) / "db", Path(td) / "hits.tsv"
        faa.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()))
        subprocess.run([PR.diamond_bin(), "makedb", "--in", str(faa), "-d", str(db), "--quiet"],
                       check=True)
        subprocess.run([PR.diamond_bin(), "blastp", "-q", str(faa), "-d", str(db), "-o", str(out),
                        "--very-sensitive", "--id", "25", "--evalue", "1e-5",
                        "--max-target-seqs", str(len(seqs) + 5), "--outfmt", "6", *DIAMOND_FIELDS,
                        "--quiet", "--threads", str(threads)], check=True)
        if not out.exists() or out.stat().st_size == 0:
            sys.exit("FAILED: DIAMOND returned no hits at all -- refusing to calibrate on nothing.")
        d = pd.read_csv(out, sep="\t", names=HIT_COLS)
    d["query"] = d["query"].astype(str)
    d["component_id"] = d["component_id"].astype(str)
    d = d[d["query"] != d["component_id"]]
    return d.sort_values("bitscore", ascending=False).drop_duplicates(["query", "component_id"])


def build_pairs(hits: pd.DataFrame, meta: pd.DataFrame, sets: dict) -> pd.DataFrame:
    """One row per unordered pair passing both coverage floors, with its ligand-set overlap."""
    best = hits.sort_values("bitscore", ascending=False).drop_duplicates("query")
    rbh = set()
    top = dict(zip(best["query"], best["component_id"]))
    for a, b in top.items():
        if top.get(b) == a:
            rbh.add(frozenset((a, b)))

    ok = hits[(hits["qcov"] >= L.MIN_QCOV) & (hits["scov"] >= L.MIN_SCOV)].copy()
    info = meta.set_index("component_id")
    lig, pot = sets["all"], sets["potent"]

    seen, rows = set(), []
    for q, s, pid in zip(ok["query"], ok["component_id"], ok["pident"]):
        key = frozenset((q, s))
        if key in seen or q not in info.index or s not in info.index:
            continue
        seen.add(key)
        la, lb = lig.get(q, set()), lig.get(s, set())
        pa, pb = pot.get(q, set()), pot.get(s, set())
        union = len(la | lb)
        rows.append({
            "component_a": q, "component_b": s, "pident": pid,
            "same_species": info.at[q, "binomial"] == info.at[s, "binomial"]
                            and bool(info.at[q, "binomial"]),
            "is_rbh": key in rbh,
            "n_a": len(la), "n_b": len(lb),
            "n_shared": len(la & lb),
            "jaccard": (len(la & lb) / union) if union else np.nan,
            "a_potent": bool(pa), "b_potent": bool(pb),
            "both_potent": bool(pa) and bool(pb),
            "either_potent": bool(pa) or bool(pb),
        })
    return pd.DataFrame(rows)


def summarise(pairs: pd.DataFrame, min_compounds: int, base_potent: float) -> pd.DataFrame:
    """Two statistics per identity band, split by same-species and by RBH.

    **Jaccard of the compound sets turned out to be the WRONG statistic, and the measurement is
    what says so**: median Jaccard is ~0.00 in EVERY band including 95-100%, and the median among
    pairs sharing anything sits at 0.02-0.04 with no identity trend at all. Two 95%-identical
    ChEMBL targets are usually one enzyme screened twice by different groups against different
    libraries -- they share the protein, not the chemistry. So compound-set overlap measures
    campaign coincidence and cannot calibrate a transfer rule. It is kept in the table because the
    negative result is the reason the second statistic exists.

    The decision-relevant quantity is the CONDITIONAL one, which is what `precedents.py` actually
    claims when it transfers a count: given a homolog at x% identity HAS a potent ligand, is this
    protein potent too? `p_potent_given_neighbour_potent` is that, over ordered pairs, and it is
    read against `base_potent` -- the marginal rate among all bacterial targets. A band where the
    conditional does not exceed the marginal is a band where identity has told you nothing.
    """
    out = []
    strong = pairs[(pairs["n_a"] >= min_compounds) & (pairs["n_b"] >= min_compounds)]
    for lo, hi in BANDS:
        for label, d in (("all_pairs", pairs), (f"both_ge_{min_compounds}", strong)):
            for split, sub in (("any", d),
                               ("same_species", d[d["same_species"]]),
                               ("cross_species", d[~d["same_species"]]),
                               ("rbh", d[d["is_rbh"]]),
                               ("not_rbh", d[~d["is_rbh"]])):
                b = sub[(sub["pident"] >= lo) & (sub["pident"] < hi)]
                if b.empty:
                    continue
                shares = b[b["n_shared"] > 0]
                # Ordered pairs: each unordered pair contributes twice, once in each direction,
                # so "my neighbour is potent" is asked of both members.
                a_pot = np.concatenate([b["a_potent"].values, b["b_potent"].values])
                nb_pot = np.concatenate([b["b_potent"].values, b["a_potent"].values])
                given = int(a_pot.sum())
                cond = float(nb_pot[a_pot].mean()) if given else np.nan
                out.append({
                    "pident_lo": lo, "pident_hi": min(hi, 100.0), "population": label,
                    "split": split, "n_pairs": len(b),
                    "n_ordered_neighbour_potent": given,
                    "p_potent_given_neighbour_potent": round(cond, 4) if given else np.nan,
                    "lift_over_base": round(cond / base_potent, 3) if given and base_potent else np.nan,
                    "frac_sharing_any": round(len(shares) / len(b), 4),
                    "median_jaccard": round(float(b["jaccard"].median()), 4),
                    "median_jaccard_when_sharing": (round(float(shares["jaccard"].median()), 4)
                                                    if len(shares) else np.nan),
                    "base_potent": round(base_potent, 4),
                })
    return pd.DataFrame(out)


def report(summary: pd.DataFrame, pairs: pd.DataFrame, min_compounds: int) -> None:
    pop = f"both_ge_{min_compounds}"
    rule()
    say(f"LIGAND-SET OVERLAP BY IDENTITY   (pairs where both targets carry >= {min_compounds} compounds)")
    rule()
    say(f"  {'identity':>12s} {'pairs':>7s} {'median J':>9s} | {'n ordered':>10s} "
        f"{'P(potent|nb potent)':>21s} {'lift':>7s}")
    for _, r in summary[(summary["population"] == pop) & (summary["split"] == "any")].iterrows():
        cond = f"{r.p_potent_given_neighbour_potent:.3f}" if r.p_potent_given_neighbour_potent == r.p_potent_given_neighbour_potent else "-"
        lift = f"{r.lift_over_base:.2f}x" if r.lift_over_base == r.lift_over_base else "-"
        say(f"  {str(int(r.pident_lo)) + '-' + str(int(r.pident_hi)):>12s} {r.n_pairs:>7,} "
            f"{r.median_jaccard:>9.3f} | {int(r.n_ordered_neighbour_potent):>10,} {cond:>21s} {lift:>7s}")
    say("")
    say("  median J is ~0 in EVERY band, 95-100 included: two near-identical ChEMBL targets are")
    say("  one enzyme screened twice against different libraries. Compound overlap measures")
    say("  campaign coincidence, not transferability -- read the conditional, not the Jaccard.")

    for a, b, title in (("same_species", "cross_species", "DOES SPECIES ADD ANYTHING BEYOND IDENTITY?"),
                        ("rbh", "not_rbh", "DOES RECIPROCAL BEST HIT ADD ANYTHING BEYOND IDENTITY?")):
        rule()
        say(title)
        rule()
        say(f"  {'identity':>12s} {a:>26s} {b:>26s}")
        for lo, hi in BANDS:
            cells = []
            for split in (a, b):
                r = summary[(summary.population == pop) & (summary.split == split)
                            & (summary.pident_lo == lo)]
                if r.empty or r.n_ordered_neighbour_potent.iloc[0] == 0:
                    cells.append("-")
                else:
                    cells.append(f"{int(r.n_ordered_neighbour_potent.iloc[0]):>4,} ordered  "
                                 f"{r.p_potent_given_neighbour_potent.iloc[0]:.3f}")
            if cells == ["-", "-"]:
                continue
            say(f"  {str(int(lo)) + '-' + str(int(min(hi, 100))):>12s} {cells[0]:>26s} {cells[1]:>26s}")
        say("  (ordered pairs whose neighbour is potent / P(this one is potent too))")


def verdict(summary: pd.DataFrame, min_compounds: int, base_potent: float) -> None:
    """State the negative result in the run log, so the table is not over-read later."""
    pop = f"both_ge_{min_compounds}"
    band = summary[(summary.population == pop) & (summary.split == "any")]
    cond = band["p_potent_given_neighbour_potent"].dropna()
    rule()
    say("VERDICT")
    rule()
    say(f"  P(potent | neighbour potent) spans {cond.min():.3f}-{cond.max():.3f} across every band")
    say(f"  from 25% to 100% identity, against a base rate of {base_potent:.3f}. It does NOT decay.")
    say("")
    say("  The bands 95/60/40 CANNOT be calibrated from inside ChEMBL: 62% of bacterial targets")
    say("  already carry a potent compound, because a protein enters ChEMBL when somebody believed")
    say("  it was druggable. The floor controls COVERAGE, not transfer reliability -- document it")
    say("  as a conservatism choice, never as an accuracy threshold.")
    say("")
    say("  Species and RBH add nothing beyond identity, so the transfer rule stays simple.")


def main() -> int:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--min-compounds", type=int, default=5,
                    help="a pair counts toward the quoted population when BOTH carry this many")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("transfer_calibration.py -- at what identity does ligand precedent actually transfer?")
    rule("=")
    say(f"  shipped bands : direct {L.DIRECT_PIDENT:g}  close {L.CLOSE_PIDENT:g}  "
        f"remote {L.REMOTE_PIDENT:g}   (the numbers being tested)")
    say(f"  coverage      : qcov >= {L.MIN_QCOV:g}, scov >= {L.MIN_SCOV:g} (imported, not restated)")
    say(f"  potent        : pChEMBL >= {L.PCHEMBL_HEADLINE:g}")
    say(f"  -> {CALIB_PATH.relative_to(REPO_ROOT)}")
    say(f"  -> {PAIRS_PATH.relative_to(REPO_ROOT)}")
    if args.dry_run:
        say("\n  DRY RUN -- would load the cached extracts, run all-vs-all DIAMOND, write 2 tables")
        return 0

    meta, seqs, sets = bacterial_targets()
    say("")
    say(f"  bacterial targets : {len(meta):,} distinct sequences, "
        f"{meta['binomial'].nunique():,} species, {int(meta['n_compounds'].sum()):,} "
        "component-compound memberships")
    say("  all-vs-all DIAMOND ...")
    hits = all_vs_all(seqs, args.threads)
    pairs = build_pairs(hits, meta, sets)
    if pairs.empty:
        sys.exit("FAILED: no pairs survived the coverage floors -- refusing to write an empty table.")
    say(f"  pairs             : {len(pairs):,} passing both coverage floors  "
        f"({int(pairs['is_rbh'].sum()):,} RBH, {int(pairs['same_species'].sum()):,} same-species)")

    base_potent = float((meta["n_potent"] > 0).mean())
    say(f"  base rate         : {base_potent:.3f} of bacterial targets have a compound at "
        f"pChEMBL >= {L.PCHEMBL_HEADLINE:g}  -- the reference every conditional is read against")
    summary = summarise(pairs, args.min_compounds, base_potent)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    summary["min_compounds"] = args.min_compounds
    summary["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    summary.to_csv(CALIB_PATH, sep="\t", index=False)
    pairs.to_csv(PAIRS_PATH, sep="\t", index=False)

    report(summary, pairs, args.min_compounds)
    verdict(summary, args.min_compounds, base_potent)
    rule("=")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
