"""Screening EFFORT per ChEMBL bacterial target -- the denominator this axis never had.

Ligandability, part 1b. One read-only SQL pass over the ChEMBL 37 dump, writing per-target counts
of compounds ACTUALLY ASSAYED -- with no pChEMBL filter -- alongside how many of them reached a
real potency. `ligands/chembl.py` and `ligands/ligands.py` both require `pchembl_value IS NOT
NULL`, so everything a chemist tried that did not work is invisible to them.

Why the axis needs this
-----------------------
Without a denominator, "nobody ever screened this protein" and "people screened it and nothing
worked" are the same zero. They are opposite pieces of evidence for a target-prioritization
shortlist: the first is an open question, the second is a measured discouragement.

`chembl_<species>.tsv` ships `n_compounds_tested`, and `docs/ligands.md` claimed it separated those
two. It does not: `chembl.py:594` computes it from the already-pChEMBL-filtered ligand table, so it
only separates WEAKLY POTENT from NEVER MEASURED. This script is what that column was supposed to
be, and the doc claim is corrected in the same commit.

Measured on ChEMBL 37, bacterial targets, under this script's own predicate
---------------------------------------------------------------------------
  2,088  bacterial targets of any type
  1,073  bacterial SINGLE PROTEIN targets
    987  with B/F activity at confidence >= 8, any relation
    869  restricted to standard_relation '=' -- the population this script reports
    661  with any potency-measurable compound -- all the current axis can see
    430  with a compound at pChEMBL >= 6

    208  targets have NOTHING potency-measurable, i.e. they are invisible today
226,038  compound-target pairs were assayed with no measurable potency

And genuine negatives are plentiful, which is the useful part: of 453 bacterial targets with >= 10
compounds assayed, **148 (32.7%) never reached pChEMBL 6**; at >= 5 it is 251 of 603 (41.6%). Those
are proteins somebody tried and failed to drug -- the only real negatives anywhere in this axis.

What a zero means here, and what it does not
---------------------------------------------
`n_compounds_assayed = 0` still means "not in ChEMBL's bacterial target set", NOT "not
ligandable". This script narrows the unknown; it does not eliminate it. And ChEMBL's potency
fields only exist for `=` relations on IC50/EC50/Ki/Kd/Potency in nM, so MIC and %-inhibition
remain absent by construction -- a compound with a published MIC and no enzyme assay counts as
assayed here only if some B/F activity row exists for it.

Bacteria only. Human is deliberately out of scope for this part -- the human columns already
shipped by `chembl.py` are a SELECTIVITY LIABILITY signal, not a precedent claim, and they are
untouched.

Needs the dump restored (see `data/raw/other/chembl/SOURCE.md`); `gradi` env, sqlite3 only, no
rdkit. ~2 min.
"""

from __future__ import annotations

import argparse
import importlib.util
import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "ligands"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
EFFORT_PATH = SCRATCH_DIR / "chembl_effort.tsv"
FUNNEL_PATH = EVIDENCE_DIR / "chembl_effort_funnel.tsv"
RELATION_PATH = EVIDENCE_DIR / "chembl_effort_relations.tsv"
ASSAYED_PATH = SCRATCH_DIR / "chembl_assayed.tsv"
EXTRA_TSV = SCRATCH_DIR / "chembl_effort_targets.tsv"
EXTRA_FAA = SCRATCH_DIR / "chembl_effort_targets.faa"

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def stage() -> object:
    """`scripts/ligands/chembl.py`, loaded by path so its SQL predicate can be REUSED not restated.

    Two things are imported rather than copied, and both have already cost this project something:
    `assert_version` (the dump's `version` table holds 11 rows and `ChEMBL_37` is not the first, so
    a `fetchone()` reads `Bioassay Ontology 2.0`) and `_activity_where` (the target/assay predicate
    that defines which targets exist at all). If the effort denominator were computed over a
    different target population from the numerator, every ratio built on it would be wrong in a way
    no shape check could see. `scripts/` is not a package, hence the explicit spec load.
    """
    spec = importlib.util.spec_from_file_location(
        "ligands_chembl", REPO_ROOT / "scripts" / "ligands" / "chembl.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.VERBOSE = False              # it prints through its own `say`; we print our own
    for name in ("assert_version", "find_db", "_activity_where"):
        if not hasattr(mod, name):
            sys.exit(f"FATAL scripts/ligands/chembl.py has no {name} -- it was refactored; "
                     "this extract must be re-aligned rather than guessed at")
    return mod


KINGDOMS = ("bacteria", "human")


def taxa_for(kingdom: str) -> set[int]:
    """Taxa for one scope, from the superkingdom map `chembl.py` already wrote.

    `bacteria` is `l1 == "Bacteria"` -- 1,493 taxa, because ChEMBL files strains under their own
    taxids (the same reason `chembl.py` matches same-species by organism NAME, not tax_id).
    `human` is the single Homo sapiens taxon; it is a LIABILITY scope, so it needs no species
    prefix logic.
    """
    path = EVIDENCE_DIR / "organism_class.tsv"
    if not path.exists():
        sys.exit(f"FAILED: {path.relative_to(REPO_ROOT)} is missing -- run scripts/ligands/chembl.py first.")
    oc = pd.read_csv(path, sep="\t")
    if kingdom == "bacteria":
        taxa = {int(t) for t in oc.loc[oc["l1"] == "Bacteria", "tax_id"]}
    elif kingdom == "human":
        taxa = {9606}
    else:
        sys.exit(f"FAILED: unknown kingdom {kingdom!r}; choose from {KINGDOMS}")
    if not taxa:
        sys.exit(f"FAILED: organism_class.tsv named no {kingdom} taxa -- refusing to write an "
                 "empty extract.")
    return taxa


def target_predicate(mod: object) -> tuple[str, str]:
    """`chembl.py`'s single-protein predicate, twice: pChEMBL relaxed, and relation relaxed too.

    The removals are the entire point of this script, so they are done by explicit string surgery
    on the imported predicate and ASSERTED, rather than by writing a second predicate that would
    drift. Everything else -- target type, confidence, assay type, validity, duplicates -- stays
    byte-identical to the numerator's.

    Why TWO predicates, measured on ChEMBL 37 (bacterial SINGLE PROTEIN, B/F, confidence >= 8):

        no activity-level filter            1,007 targets   420,012 rows
        + data_validity_comment IS NULL       987           411,751
        + not potential_duplicate             987           410,662
        + standard_relation = '='             869           370,493

    That last clause costs **118 targets**, and with them (measured under the relaxed predicate, in
    `evidence/chembl_effort_relations.tsv`) **6,682 `>` rows covering 4,834 compounds**, plus 32,093
    rows whose relation is null. A result published as `IC50 > 100 uM` is the clearest statement in
    the whole database that a compound does NOT bind -- so for a DENOMINATOR, dropping it discards
    exactly the negatives this extract exists to recover. But for a RATIO the denominator must match
    the numerator's population, or the hit rate is computed over two different sets of rows.

    So both are kept and named, and nothing is silently resolved: `matched` is population-identical
    to `chembl.py` and is what `hit_rate` divides by; `any_relation` is the honest effort count.
    """
    where = mod._activity_where("single")
    drop = "AND act.pchembl_value IS NOT NULL"
    if drop not in where:
        sys.exit("FATAL scripts/ligands/chembl.py no longer filters on `act.pchembl_value IS NOT "
                 "NULL` in the single-protein predicate. This script exists to relax exactly that "
                 "clause; re-align it rather than guessing.")
    matched = where.replace(drop, "")
    rel = "AND act.standard_relation = '='"
    if rel not in matched:
        sys.exit("FATAL scripts/ligands/chembl.py no longer filters on `act.standard_relation = "
                 "'='`. The two denominators in this script are defined by that clause; re-align "
                 "rather than guessing.")
    return matched, matched.replace(rel, "")


def extract(db: Path, mod: object, kingdom: str = "bacteria") -> tuple[pd.DataFrame, ...]:
    """One pass over ONE scope: per-target effort counts, the funnel, relations, assayed pairs."""
    taxa = taxa_for(kingdom)
    ids = ",".join(str(t) for t in sorted(taxa))
    where, where_any = target_predicate(mod)
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        sql = f"""
            SELECT td.tid                                            AS tid,
                   td.chembl_id                                      AS target_chembl_id,
                   td.organism                                       AS organism,
                   COUNT(DISTINCT act.molregno)                      AS n_compounds_assayed,
                   COUNT(DISTINCT CASE WHEN act.pchembl_value IS NOT NULL
                                       THEN act.molregno END)        AS n_compounds_measurable,
                   COUNT(DISTINCT CASE WHEN act.pchembl_value >= {L.PCHEMBL_HEADLINE}
                                       THEN act.molregno END)        AS n_compounds_potent,
                   COUNT(DISTINCT act.doc_id)                        AS n_docs,
                   COUNT(DISTINCT act.assay_id)                      AS n_assays,
                   SUM(CASE WHEN a.assay_type='B' THEN 1 ELSE 0 END) AS n_binding,
                   SUM(CASE WHEN a.assay_type='F' THEN 1 ELSE 0 END) AS n_functional,
                   MAX(act.pchembl_value)                            AS max_pchembl,
                   MIN(d.year)                                       AS first_doc_year,
                   MAX(d.year)                                       AS last_doc_year
              FROM target_dictionary td
              JOIN assays a             ON a.tid       = td.tid
              JOIN activities act       ON act.assay_id = a.assay_id
              LEFT JOIN docs d          ON d.doc_id     = act.doc_id
             WHERE td.tax_id IN ({ids}) AND {where}
             GROUP BY 1, 2, 3
        """
        eff = pd.read_sql_query(sql, con)

        # The second denominator: every relation, so an `IC50 > 100 uM` counts as assayed. Joined
        # on, never merged into the matched count -- see `target_predicate`.
        any_sql = f"""
            SELECT td.tid                                        AS tid,
                   COUNT(DISTINCT act.molregno)                  AS n_compounds_assayed_any_relation,
                   COUNT(DISTINCT CASE WHEN act.standard_relation IN ('>','>=','>>')
                                       THEN act.molregno END)    AS n_compounds_reported_inactive
              FROM target_dictionary td
              JOIN assays a       ON a.tid       = td.tid
              JOIN activities act ON act.assay_id = a.assay_id
             WHERE td.tax_id IN ({ids}) AND {where_any}
             GROUP BY 1
        """
        eff = eff.merge(pd.read_sql_query(any_sql, con), on="tid", how="left")

        # The per-tid COUNTS above cannot be unioned: a compound assayed against three homologous
        # targets must count ONCE, and summing counts would inflate the denominator exactly where
        # the homology pool is widest. `precedents.py` counts distinct molecules over the union for
        # its numerator, so the denominator needs the same treatment -- which needs the ids, not
        # the totals. Parent molregno, matching chembl_ligands.tsv's key.
        assayed = pd.read_sql_query(f"""
            SELECT DISTINCT td.tid                                   AS tid,
                   COALESCE(mh.parent_molregno, act.molregno)        AS parent_molregno
              FROM target_dictionary td
              JOIN assays a       ON a.tid       = td.tid
              JOIN activities act ON act.assay_id = a.assay_id
              LEFT JOIN molecule_hierarchy mh ON mh.molregno = act.molregno
             WHERE td.tax_id IN ({ids}) AND {where_any}
        """, con)

        # SEQUENCES for the targets `chembl.py` never put in the FASTA. Without these the whole
        # `screened_clean` tier is UNREACHABLE and would ship as a permanently empty category:
        # `chembl.py` builds chembl_targets.faa from the tids that survived the pChEMBL filter, so
        # a target whose compounds were all assayed and none measurable has no sequence, DIAMOND
        # can never hit it, and no protein can ever be told "somebody screened your homolog and
        # nothing came out". Measured: 208 bacterial targets, 1,593 assayed compounds, 0 measurable.
        # Same column shape as chembl_targets.tsv so the two concatenate.
        extra_t = pd.read_sql_query(f"""
            SELECT DISTINCT tc.component_id        AS component_id,
                   td.tid                          AS tid,
                   td.chembl_id                    AS target_chembl_id,
                   td.target_type                  AS target_type,
                   td.tax_id                       AS tax_id,
                   td.organism                     AS organism,
                   td.pref_name                    AS pref_name,
                   cs.accession                    AS accession,
                   COALESCE(oc.l1, 'Unclassified') AS superkingdom,
                   cs.sequence                     AS sequence
              FROM target_dictionary td
              JOIN target_components   tc ON tc.tid = td.tid
              JOIN component_sequences cs ON cs.component_id = tc.component_id
              LEFT JOIN organism_class oc ON oc.tax_id = td.tax_id
             WHERE td.tax_id IN ({ids}) AND cs.sequence IS NOT NULL
               AND td.tid IN (
                   SELECT DISTINCT td2.tid
                     FROM target_dictionary td2
                     JOIN assays a2       ON a2.tid       = td2.tid
                     JOIN activities act2 ON act2.assay_id = a2.assay_id
                    WHERE td2.tax_id IN ({ids})
                      AND {where_any.replace('td.', 'td2.').replace('a.', 'a2.')
                                    .replace('act.', 'act2.')})
        """, con)

        rel = pd.read_sql_query(f"""
            SELECT COALESCE(act.standard_relation,'(null)') AS standard_relation,
                   COUNT(*) AS n_rows, COUNT(DISTINCT act.molregno) AS n_compounds
              FROM target_dictionary td
              JOIN assays a       ON a.tid       = td.tid
              JOIN activities act ON act.assay_id = a.assay_id
             WHERE td.tax_id IN ({ids}) AND {where_any}
             GROUP BY 1 ORDER BY 2 DESC""", con)

        # The funnel is EVIDENCE, not decoration: it is what makes "656 targets" a reproducible
        # number rather than a remembered one, and it says where each loss happens.
        base = f"""FROM target_dictionary td
                   JOIN assays a       ON a.tid       = td.tid
                   JOIN activities act ON act.assay_id = a.assay_id
                  WHERE td.tax_id IN ({ids}) AND {where}"""
        base_any = base.replace(where, where_any)
        rows = [
            (f"{kingdom} targets, any type",
             f"SELECT COUNT(*) FROM target_dictionary WHERE tax_id IN ({ids})"),
            (f"{kingdom} SINGLE PROTEIN targets",
             f"SELECT COUNT(*) FROM target_dictionary WHERE target_type='SINGLE PROTEIN' "
             f"AND tax_id IN ({ids})"),
            ("...with B/F activity, any relation",
             f"SELECT COUNT(DISTINCT td.tid) {base_any}"),
            ("...restricted to standard_relation '='",
             f"SELECT COUNT(DISTINCT td.tid) {base}"),
            ("...with any potency-measurable compound",
             f"SELECT COUNT(DISTINCT td.tid) {base} AND act.pchembl_value IS NOT NULL"),
            (f"...with a compound at pChEMBL >= {L.PCHEMBL_HEADLINE:g}",
             f"SELECT COUNT(DISTINCT td.tid) {base} AND act.pchembl_value >= {L.PCHEMBL_HEADLINE}"),
        ]
        funnel = pd.DataFrame(
            [{"stage": name, "n_targets": con.execute(q).fetchone()[0]} for name, q in rows])
    finally:
        con.close()

    eff["n_compounds_assayed"] = eff["n_compounds_assayed"].astype("Int64")
    for c in ("n_compounds_measurable", "n_compounds_potent", "n_docs", "n_assays",
              "n_binding", "n_functional", "first_doc_year", "last_doc_year",
              "n_compounds_assayed_any_relation", "n_compounds_reported_inactive"):
        eff[c] = pd.to_numeric(eff[c], errors="coerce").astype("Int64")
    eff["max_pchembl"] = pd.to_numeric(eff["max_pchembl"], errors="coerce").astype("Float64")

    # hit_rate is NULL, never 0, where nothing was assayed -- a zero there would claim a
    # measurement nobody made. Same rule as the studiedness axis's `no_hit` vs `below_floor`.
    denom = eff["n_compounds_assayed"]
    eff["hit_rate"] = (eff["n_compounds_potent"] / denom).where(denom > 0).astype("Float64")
    # Tag the scope on every frame: one file serves both, and a row whose kingdom is unstated is a
    # row that will eventually be counted against the wrong pool.
    for f in (eff, funnel, rel, assayed, extra_t):
        f.insert(0, "kingdom", kingdom)
    return eff.sort_values("tid").reset_index(drop=True), funnel, rel, assayed, extra_t


def report(eff: pd.DataFrame, funnel: pd.DataFrame, rel: pd.DataFrame) -> None:
    rule()
    say("FUNNEL  -- why the target set is the size it is")
    rule()
    for _, r in funnel.iterrows():
        say(f"  {r['stage']:44s} {r['n_targets']:>7,}")

    invisible = int((eff["n_compounds_measurable"].fillna(0) == 0).sum())
    pairs = int((eff["n_compounds_assayed"].fillna(0) - eff["n_compounds_measurable"].fillna(0)).sum())
    say("")
    say(f"  targets with NOTHING potency-measurable        {invisible:>7,}   "
        "(invisible to chembl.py and precedents.py)")
    say(f"  compound-target pairs assayed, no potency      {pairs:>7,}")

    rule()
    say("RELATION  -- what the '=' clause costs, and why both denominators ship")
    rule()
    for _, r in rel.head(6).iterrows():
        say(f"  {r['standard_relation']:>8s}  {r['n_rows']:>9,} rows  {r['n_compounds']:>8,} compounds")
    gt = int(eff["n_compounds_reported_inactive"].fillna(0).sum())
    say("")
    say(f"  compounds with a '>' result (explicit NON-binders): {gt:,}")
    say("  `n_compounds_assayed` matches chembl.py's population and is what hit_rate divides by;")
    say("  `n_compounds_assayed_any_relation` is the honest effort count. Neither is merged.")

    rule()
    say("REAL NEGATIVES  -- assayed, and nothing ever reached pChEMBL 6")
    rule()
    say(f"  {'compounds assayed':>18s} {'targets':>9s} {'never potent':>14s} {'':>7s}")
    for thr in (5, 10, 25, 50, 100):
        sub = eff[eff["n_compounds_assayed"].fillna(0) >= thr]
        neg = sub[sub["n_compounds_potent"].fillna(0) == 0]
        pct = 100 * len(neg) / len(sub) if len(sub) else 0.0
        say(f"  {'>= ' + str(thr):>18s} {len(sub):>9,} {len(neg):>14,} {pct:>6.1f}%")
    say("")
    say("  these are the only measured DISCOURAGEMENTS in the axis -- a protein somebody tried")
    say("  and failed to drug is a different object from one nobody has opened.")


def main() -> int:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--kingdoms", nargs="+", choices=list(KINGDOMS), default=list(KINGDOMS),
                    help="scopes to extract (default both). `human` is needed for the precedent "
                         "table's n_assayed_human; it is ~13x the bacterial volume.")
    ap.add_argument("--refresh", action="store_true", help="re-extract even if the cache exists")
    ap.add_argument("--dry-run", action="store_true", help="say what would be done, write nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("effort.py -- compounds ASSAYED per ChEMBL target, with no pChEMBL filter")
    rule("=")
    say(f"  source   : ChEMBL {L.CHEMBL_VERSION} SQLite dump")
    say(f"  scopes   : {', '.join(args.kingdoms)}")
    say(f"  targets  : SINGLE PROTEIN, confidence >= {L.MIN_CONFIDENCE_SINGLE}, "
        f"assay type {'/'.join(L.ASSAY_TYPES)}")
    say("  relaxed  : pchembl_value may be NULL -- that is the entire point of this extract")
    say(f"  potent   : pChEMBL >= {L.PCHEMBL_HEADLINE:g}")
    for pth in (EFFORT_PATH, FUNNEL_PATH, RELATION_PATH, ASSAYED_PATH, EXTRA_TSV):
        say(f"  -> {pth.relative_to(REPO_ROOT)}")

    if EFFORT_PATH.exists() and not args.refresh:
        say("")
        say(f"  cached   : {EFFORT_PATH.name} (pass --refresh to rebuild)")
        eff = L._read(EFFORT_PATH)
        say(f"             {len(eff):,} targets")
        return 0

    mod = stage()
    if args.dry_run:
        say("")
        say(f"  DRY RUN -- would restore-check the dump and run one SQL pass per scope "
            f"({len(args.kingdoms)})")
        return 0

    db = mod.find_db(required=True)           # exits with the tar recipe if the dump is absent
    mod.assert_version(db)                    # the dump must BE the release we claim
    say("")
    say(f"  db       : {db.relative_to(REPO_ROOT)}  (version asserted)")

    effs, funnels, rels, assays, extras = [], [], [], [], []
    for kingdom in args.kingdoms:
        say(f"  extracting {kingdom} ...")
        eff, funnel, rel, assayed, extra = extract(db, mod, kingdom)
        if eff.empty:
            sys.exit(f"FAILED: the {kingdom} extract is empty -- refusing to write. Check the "
                     "dump and the taxa map.")
        say(f"    {len(eff):>7,} targets   {len(assayed):>9,} (target, assayed compound) pairs")
        effs.append(eff); funnels.append(funnel); rels.append(rel)
        assays.append(assayed); extras.append(extra)

    eff = pd.concat(effs, ignore_index=True)
    assayed = pd.concat(assays, ignore_index=True)
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    eff.to_csv(EFFORT_PATH, sep="\t", index=False)
    assayed.to_csv(ASSAYED_PATH, sep="\t", index=False)
    for frames, path in ((funnels, FUNNEL_PATH), (rels, RELATION_PATH)):
        d = pd.concat(frames, ignore_index=True)
        d["built_utc"] = stamp
        d.to_csv(path, sep="\t", index=False)

    # The DIAMOND supplement is BACTERIA-ONLY on purpose: it exists so bacterial targets whose
    # compounds were all assayed and none measurable can be reached by the search at all. Human is
    # a liability scope where only potent matters, so it needs no equivalent -- and adding human
    # sequences to the subject database would change `n_ligands_human` for reasons unrelated to
    # this extract. Only targets the existing FASTA does NOT carry: it is a SUPPLEMENT that
    # `src/precedents.py` concatenates, never a replacement, and overlap would double-count.
    bact = [e for e, k in zip(extras, args.kingdoms) if k == "bacteria"]
    if bact:
        known = set(pd.read_csv(SCRATCH_DIR / "chembl_targets.tsv", sep="\t")["tid"])
        new = bact[0][~bact[0]["tid"].isin(known)].copy()
        new.drop(columns=["sequence"]).to_csv(EXTRA_TSV, sep="\t", index=False)
        seqs = new.drop_duplicates("component_id")
        EXTRA_FAA.write_text("".join(f">{c}\n{q}\n" for c, q in
                                     zip(seqs["component_id"], seqs["sequence"])))
        say(f"  supplement: {len(new):,} target rows / {len(seqs):,} sequences the pChEMBL filter "
            f"had hidden -> {EXTRA_FAA.name}")

    say(f"  wrote    : {len(eff):,} targets x {eff.shape[1]} columns, "
        f"{len(assayed):,} assayed pairs -- the denominator can be UNIONED, not summed")
    for kingdom in args.kingdoms:
        sub = eff[eff["kingdom"] == kingdom]
        rule()
        say(f"SCOPE: {kingdom}")
        report(sub, pd.concat([f for f, k in zip(funnels, args.kingdoms) if k == kingdom]),
               pd.concat([r for r, k in zip(rels, args.kingdoms) if k == kingdom]))
    rule("=")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
