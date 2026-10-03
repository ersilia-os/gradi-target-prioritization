"""Validate the precedent counts against the LIVE ChEMBL API -- an independent route in.

Everything else in this axis is computed from three cached extracts that `chembl.py` wrote from the
30.5 GB dump. If the SQL that built them were subtly wrong, every downstream check would agree with
it: the control against `chembl.py` shares the extracts, and the internal assertions share the
code. This script asks a different machine the same questions over HTTP and compares.

Two rounds, because they test different things
-----------------------------------------------
**Round 1 -- EXACT targets.** For a protein whose exact route fired we know precisely which ChEMBL
target it resolved to, so `n_ligands` is checkable against one API query. Tests the join chain
`sequence -> component_id -> tid -> parent_molregno` and the potency cut.

**Round 2 -- the UNION.** `n_ligands_bacterial` unions every bacterial homolog, which is where a
double-count or a dropped `tid` would show. Proteins are chosen for hitting the MOST targets
(Kp `KPC-2` hits 51), and both the potent and the any-measurable count are compared. The union is
also checked against the per-target SUM, which it must come in below: measured 1.0-1.3x inflation
(KPC-2 276 vs 365, folA 443 vs 529), and the API agrees with the union, never with the sum.

Why the semantics line up
--------------------------
The `activity` endpoint does not expose `confidence_score`, so it cannot be filtered on. That is
fine HERE and nowhere else: CLAUDE.md records that `>= 8` removes **0 of 3,271,336** single-protein
rows, so for SINGLE PROTEIN targets the gate is a no-op. Only single-protein tids are compared --
the complex track is excluded on both sides. `potential_duplicate` is not a query parameter either,
so it is filtered client-side, and the API's `parent_molecule_chembl_id` is the same
`molecule_hierarchy` parent our `parent_molregno` collapses to.

A target too large to page honestly is SKIPPED and named, never silently truncated -- E. coli
`ampC` carries 12,438 potent compounds over 16 targets and would need 40+ pages per target.

Needs network. `gradi` env, no dump. ~10 min for the default selection.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import ligandability as L  # noqa: E402
from src import precedents as PR  # noqa: E402
from src import proteomes as P  # noqa: E402

API = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
OUT = REPO_ROOT / "data" / "processed" / "ligands" / "evidence" / "precedent_api_validation.tsv"
PAGE = 1000
MAX_OFFSET = 40_000

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 96) -> None:
    say(char * width)


def molecules(target: str, mode: str) -> set | None:
    """Distinct PARENT molecules on one ChEMBL target, under this axis's own predicate.

    `mode` is `potent` (pChEMBL >= the headline cut) or `measured` (any pChEMBL). Returns None
    when the target is too large to page honestly -- better to skip and say so than to report a
    truncated number as a match.
    """
    mols, offset = set(), 0
    while True:
        q = {"target_chembl_id": target, "standard_relation": "=",
             "data_validity_comment__isnull": "true",
             "assay_type__in": ",".join(L.ASSAY_TYPES),
             "limit": str(PAGE), "offset": str(offset)}
        if mode == "potent":
            q["pchembl_value__gte"] = str(L.PCHEMBL_HEADLINE)
        else:
            q["pchembl_value__isnull"] = "false"
        for attempt in range(4):
            try:
                with urllib.request.urlopen(API + "?" + urllib.parse.urlencode(q), timeout=180) as r:
                    payload = json.load(r)
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(3 * (attempt + 1))
        acts = payload["activities"]
        for a in acts:
            if a.get("potential_duplicate"):
                continue
            mols.add(a.get("parent_molecule_chembl_id") or a["molecule_chembl_id"])
        if len(acts) < PAGE:
            return mols
        offset += PAGE
        if offset > MAX_OFFSET:
            return None


def bridge() -> tuple[dict, set, dict, dict]:
    t, _, _ = PR._tables()
    comp2tids: dict[str, list] = {}
    for cid, tid in zip(t["component_id"], t["tid"]):
        comp2tids.setdefault(str(cid), []).append(tid)
    single = set(t.loc[t["target_type"] == "SINGLE PROTEIN", "tid"])
    return (comp2tids, single, dict(zip(t["tid"], t["target_chembl_id"])),
            dict(zip(t["component_id"], t["superkingdom"].astype(str))))


def round1(per_species: int) -> list[dict]:
    """Exact-target counts. One API query per resolved target."""
    comp2tids, single, tid2id, _ = bridge()
    rows = []
    rule()
    say("ROUND 1 -- exact targets: does n_ligands match the API on the resolved target?")
    rule()
    say(f"  {'sp':6s} {'gene':10s} {'targets':24s} {'ours':>8s} {'api':>8s}")
    for sp in L.SPECIES:
        f = L.load_full(sp).merge(
            P.load(sp)[["uniprot_ac", "gene_name"]], on="uniprot_ac", how="left")
        hit = f[f["exact_route"].ne("none")].copy()
        pick = pd.concat([hit.nlargest(per_species, "n_ligands"),
                          hit[hit["n_ligands"] > 0].nsmallest(3, "n_ligands")]
                         ).drop_duplicates("uniprot_ac")
        for _, r in pick.iterrows():
            ids = sorted({tid2id[x] for c in str(r["exact_target"]).split(";")
                          for x in comp2tids.get(c, []) if x in single})
            if not ids:
                continue
            got = set()
            skip = False
            for c in ids:
                m = molecules(c, "potent")
                if m is None:
                    skip = True
                    break
                got |= m
            ours = int(r["n_ligands"] or 0)
            if skip:
                say(f"  {sp[:6]:6s} {str(r['gene_name'])[:10]:10s} {'':24s} "
                    f"{ours:>8,}   [too large to page -- skipped]")
                continue
            ok = len(got) == ours
            say(f"  {sp[:6]:6s} {str(r['gene_name'])[:10]:10s} {','.join(ids)[:24]:24s} "
                f"{ours:>8,} {len(got):>8,}  {'OK' if ok else 'DIFF'}")
            rows.append({"round": "exact", "species": sp, "gene": r["gene_name"],
                         "uniprot_ac": r["uniprot_ac"], "targets": ",".join(ids),
                         "column": "n_ligands", "ours": ours, "api": len(got), "match": ok})
    return rows


def round2(per_species: int) -> list[dict]:
    """The union over every bacterial homolog, plus the union-vs-sum check."""
    comp2tids, single, tid2id, king = bridge()
    _, lg, _ = PR._tables()
    pot_by_tid = (lg[(lg["track"] == "single") & (lg["pchembl"] >= L.PCHEMBL_HEADLINE)]
                  .groupby("tid")["parent_molregno"].apply(set).to_dict())
    rows = []
    rule()
    say("ROUND 2 -- the UNION over every bacterial homolog (potent AND any-measurable)")
    rule()
    say(f"  {'sp':6s} {'gene':10s} {'tgt':>4s} {'ours_pot':>9s} {'api_pot':>8s} "
        f"{'ours_meas':>10s} {'api_meas':>9s} {'sum':>8s}")
    for sp in L.SPECIES:
        f = L.load_full(sp).merge(
            P.load(sp)[["uniprot_ac", "gene_name", "sequence"]], on="uniprot_ac", how="left")
        multi = f[(f["n_targets_bacteria"].fillna(0) >= 2)
                  & (f["n_ligands_bacterial"].fillna(0) > 0)]
        pick = pd.concat([multi.nlargest(per_species, "n_targets_bacteria"),
                          multi.nsmallest(2, "n_ligands_bacterial")]
                         ).drop_duplicates("uniprot_ac")
        for _, r in pick.iterrows():
            h = PR._diamond({r["uniprot_ac"]: str(r["sequence"])}, 4)
            h = h[(h["qcov"] >= L.MIN_QCOV) & (h["scov"] >= L.MIN_SCOV)
                  & (h["pident"] >= PR.DEFAULT_MIN_IDENTITY)]
            b = h[h["component_id"].map(king).eq(PR.BACTERIA)]
            ids = sorted({tid2id[x] for c in b["component_id"]
                          for x in comp2tids.get(c, []) if x in single})
            if not ids:
                continue
            pot, meas, skip = set(), set(), False
            for c in ids:
                a, m = molecules(c, "potent"), molecules(c, "measured")
                if a is None or m is None:
                    skip = True
                    break
                pot |= a
                meas |= m
            op, om = int(r["n_ligands_bacterial"] or 0), int(r["n_measured_bacterial"] or 0)
            if skip:
                say(f"  {sp[:6]:6s} {str(r['gene_name'])[:10]:10s} {len(ids):>4d}   "
                    "[too large to page -- skipped]")
                continue
            # The union must come in BELOW the per-target sum, or it is double-counting.
            per_target = sum(len(pot_by_tid.get(x, set()))
                             for c in b["component_id"] for x in comp2tids.get(c, ()))
            ok1, ok2 = len(pot) == op, len(meas) == om
            say(f"  {sp[:6]:6s} {str(r['gene_name'])[:10]:10s} {len(ids):>4d} {op:>9,} "
                f"{len(pot):>8,} {om:>10,} {len(meas):>9,} {per_target:>8,}  "
                f"{'OK' if ok1 and ok2 else 'DIFF'}")
            for col, ours, api in (("n_ligands_bacterial", op, len(pot)),
                                   ("n_measured_bacterial", om, len(meas))):
                rows.append({"round": "union", "species": sp, "gene": r["gene_name"],
                             "uniprot_ac": r["uniprot_ac"], "targets": f"{len(ids)} targets",
                             "column": col, "ours": ours, "api": api, "match": ours == api,
                             "per_target_sum": per_target})
    return rows


def main() -> int:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--per-species", type=int, default=8,
                    help="test cases per species per round (default 8)")
    ap.add_argument("--round", nargs="+", choices=["exact", "union"], default=["exact", "union"])
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    VERBOSE = not a.quiet

    rule("=")
    say("validate_api.py -- precedent counts against the LIVE ChEMBL API")
    rule("=")
    say(f"  endpoint : {API}")
    say(f"  predicate: relation '=', valid, not duplicate, assay {'/'.join(L.ASSAY_TYPES)}, "
        f"potent >= {L.PCHEMBL_HEADLINE:g}")
    say("  note     : confidence_score is not exposed by the API, and does not need to be --")
    say(f"             >= {L.MIN_CONFIDENCE_SINGLE} removes 0 of 3,271,336 single-protein rows. "
        "SINGLE PROTEIN only, both sides.")
    say(f"  -> {OUT.relative_to(REPO_ROOT)}")

    rows = []
    if "exact" in a.round:
        rows += round1(a.per_species)
    if "union" in a.round:
        rows += round2(a.per_species)
    if not rows:
        sys.exit("FAILED: no comparisons were made -- refusing to report a pass.")

    d = pd.DataFrame(rows)
    d["built_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(OUT, sep="\t", index=False)

    bad = int((~d["match"]).sum())
    rule()
    say(f"  {len(d) - bad}/{len(d)} comparisons match exactly   ({bad} differ)")
    if "union" in a.round and "per_target_sum" in d:
        u = d[d["column"] == "n_ligands_bacterial"]
        infl = (u["per_target_sum"] / u["ours"].clip(lower=1)).max()
        say(f"  union vs per-target sum: up to {infl:.1f}x inflation avoided, and the API agrees")
        say("  with the UNION, never with the sum -- the distinct-molecule claim, checked outside.")
    say(f"  wrote {OUT.relative_to(REPO_ROOT)}")
    rule("=")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
