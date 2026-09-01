"""Stage 08b — validate the web-app export (app/data/*.json).

A lightweight schema/coverage smoke-check so front-end regressions or a broken 08a run are caught
before deploy. Exits non-zero on any failure. Run with any python3 (stdlib only):

    python scripts/08b_validate_export.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DATA = REPO_ROOT / "app" / "data"

# columns the front-end relies on (config.js views / cards / composite)
REQUIRED = [
    "uniprot_accession", "gene", "name", "protein_name",
    "comp_essentiality", "comp_ligandability", "comp_degradability", "comp_novelty",
    "essentiality_tier", "ligandability_tier", "degradability_tier",
    "clp_accessibility", "localization", "comp_human_selective", "selectivity",
    "conservation_score", "human_closeness", "tsne_x", "tsne_y", "family",
    "evidence_experimental", "evidence_predictor",
]
# columns that must lie in [0, 1] when present/non-null
UNIT_RANGE = [
    "comp_essentiality", "comp_ligandability", "comp_degradability", "comp_novelty",
    "comp_breadth", "clp_accessibility", "conservation_score", "human_closeness",
    "evidence_binding", "evidence_pocket", "evidence_experimental", "evidence_predictor",
    "evidence_transfer", "essentiality_score", "ligandability_score",
]
ORGS = ["kp", "ec"]


def fail(errs: list[str], msg: str) -> None:
    errs.append(msg)


def check_org(prefix: str, errs: list[str]) -> None:
    p = DATA / f"{prefix}.json"
    if not p.exists():
        fail(errs, f"{prefix}.json missing"); return
    d = json.loads(p.read_text())
    tag = f"[{prefix}]"

    if d.get("format") != "columnar":
        fail(errs, f"{tag} format is {d.get('format')!r}, expected 'columnar'")
    cols = d.get("columns") or []
    idx = {c: i for i, c in enumerate(cols)}
    rows = d.get("rows") or []

    for c in REQUIRED:
        if c not in idx:
            fail(errs, f"{tag} required column missing: {c}")
    if d.get("n") != len(rows):
        fail(errs, f"{tag} n={d.get('n')} != rows={len(rows)}")

    ncol = len(cols)
    bad_len = sum(1 for r in rows if len(r) != ncol)
    if bad_len:
        fail(errs, f"{tag} {bad_len} rows have != {ncol} cells")

    # NaN/Inf must not leak into JSON (08a maps them to None)
    nan = 0
    for r in rows:
        for v in r:
            if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
                nan += 1
    if nan:
        fail(errs, f"{tag} {nan} NaN/Inf cells (should be null)")

    # 0–1 range
    for c in UNIT_RANGE:
        if c not in idx:
            continue
        i = idx[c]
        oob = [r[i] for r in rows if isinstance(r[i], (int, float)) and not (-1e-6 <= r[i] <= 1 + 1e-6)]
        if oob:
            fail(errs, f"{tag} {c}: {len(oob)} values out of [0,1] (e.g. {oob[:3]})")

    # accessions unique + non-empty
    ai = idx.get("uniprot_accession")
    if ai is not None:
        accs = [r[ai] for r in rows]
        if any(not a for a in accs):
            fail(errs, f"{tag} empty uniprot_accession present")
        if len(set(accs)) != len(accs):
            fail(errs, f"{tag} duplicate accessions ({len(accs) - len(set(accs))})")

    # pockets sidecar
    pk = DATA / f"pockets_{prefix}.json"
    if pk.exists():
        pm = json.loads(pk.read_text())
        if not isinstance(pm, dict):
            fail(errs, f"{tag} pockets file is not an object")
        elif pm:
            k = next(iter(pm))
            if not isinstance(pm[k], list) or not all(isinstance(x, int) for x in pm[k]):
                fail(errs, f"{tag} pockets[{k}] is not a list[int]")

    # coverage report (informational)
    def cov(c):
        i = idx.get(c)
        if i is None:
            return "absent"
        nn = sum(1 for r in rows if r[i] not in (None, "", ))
        return f"{100 * nn / max(1, len(rows)):.0f}%"
    print(f"{tag} {len(rows)} rows, {ncol} cols · "
          f"ess {cov('comp_essentiality')} lig {cov('comp_ligandability')} "
          f"deg {cov('comp_degradability')} clp-access {cov('clp_accessibility')} "
          f"loc {cov('localization')} name {cov('protein_name')}")


def main() -> int:
    errs: list[str] = []
    for prefix in ORGS:
        check_org(prefix, errs)
    if errs:
        print("\nFAILED:")
        for e in errs:
            print("  ✗", e)
        return 1
    print("\nOK — export is valid.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
