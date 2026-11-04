"""Map the ASKA (-) clone library onto the E. coli K-12 MG1655 reference proteome, and CHECK it.

Reads   data/source/aska/ASKA_clone_minus_library_list.xlsx   (4,123 clones; NO UniProt accession)
Writes  data/processed/validation/aska_to_uniprot_ecoli.tsv      THE DELIVERABLE: a SUBSET of the proteome,
                one row per mapped clone -- `uniprot_ac` · `aska_id` · `confidence` (3/2/1), in the
                canonical order of `proteome_ecoli.tsv`. Clones with no protein in the reference
                proteome are NOT rows here (they are in the evidence file, tier U).
        data/processed/validation/evidence/aska_map_ecoli.tsv    every route and check, per clone
        data/processed/validation/evidence/aska_clones_ecoli.tsv ALL clones, with tier letters and flags
        data/processed/validation/evidence/aska_audit.tsv        the rate of every route and tier

THREE ROUTES propose an accession (the project's identifier-mapping rule: measure each, then union):

  jw       the clone's Keio JW id (or the JW id its Comment says it was "cloned as") against the JW
           tags in `locus_tag_all` -- the locus-tag route, `CLAUDE.md`'s preferred key
  gene     `Gene Name` against `gene_name` and `gene_synonyms`; a name hitting >1 protein is ambiguous
  primer   proteins whose termini match BOTH of the clone's own primers (see below)

THE CHECK THAT DOES NOT SHARE AN IDENTIFIER WITH THE ROUTES is the clone's own primer sequences,
which are the physical DNA the library was built from:

  N primer   `GCC` + 20 nt of the ORF STARTING AT ITS SECOND CODON (the start codon comes from the
             vector: ybeM is `MMTTILT...`, its primer encodes `MTTILT`). Translated, it must equal
             `protein[1:7]`.
  C primer   `CC` + the reverse complement of the last 21 nt of the ORF, stop codon excluded.
             Reverse-complemented and translated, it must equal the protein's last 6 residues.

Both are 6-residue exact matches, so a chance hit is ~20^-6 per end. **The test is shown to be able to
FAIL**: the same checks on a SHUFFLED mapping must pass for ~0 clones, and the script exits non-zero
if they do not.

The CHOSEN accession is the candidate (from the three routes) with the most primer ends verified,
ties broken jw > gene > primer -- so the JW id decides unless the primers positively prefer another
protein, and that override is recorded as `rule = primer_override`. Nothing is resolved silently.

SHIPPED `confidence` (the project's 1-3 evidence ladder): 3 = tier A, 2 = tier B, 1 = tier C or D.

CONFIDENCE TIERS
  A  both primer ends match the chosen protein
  B  one end matches (the other differs, usually a start/stop revised since the library was made)
  C  neither end matches, but the JW id AND the gene name agree on the protein
  D  neither end matches and the identifier routes do not agree -- needs a manual look
  U  not found in the proteome

Run with the `gradi` env:  python scripts/validation/aska.py [-q] [--floor 0.90]
"""

from __future__ import annotations

import argparse
import random
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from src import proteomes  # noqa: E402

SPECIES = "ecoli"
SRC = REPO_ROOT / "data/source/aska/ASKA_clone_minus_library_list.xlsx"
PROCESSED = REPO_ROOT / "data/processed/validation"
OUT = PROCESSED / "evidence"
JW_RE = re.compile(r"JW\d{4}(?:\.\d+)?")
K = 6  # residues compared at each end

_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON = {a + b + c: _AAS[i] for i, (a, b, c) in enumerate(
    (x, y, z) for x in _BASES for y in _BASES for z in _BASES)}
_COMP = str.maketrans("ACGT", "TGCA")


def log(msg: str = "", quiet: bool = False) -> None:
    if not quiet:
        print(msg)


def translate(nt: str) -> str:
    return "".join(CODON.get(nt[i:i + 3], "X") for i in range(0, len(nt) - len(nt) % 3, 3))


def n_peptide(primer: str) -> str | None:
    """First K residues encoded by the N primer, which begins at the protein's SECOND codon."""
    p = str(primer or "").strip().upper()
    if not p.startswith("GCC"):
        return None
    pep = translate(p[3:])
    return pep[:K] if len(pep) >= K and "*" not in pep[:K] and "X" not in pep[:K] else None


def c_peptide(primer: str) -> str | None:
    """Last K residues encoded by the C primer (reverse complement of the ORF's 3' end)."""
    p = str(primer or "").strip().upper()
    if not p.startswith("CC"):
        return None
    pep = translate(p[2:].translate(_COMP)[::-1])
    return pep[-K:] if len(pep) >= K and "*" not in pep[-K:] and "X" not in pep[-K:] else None


def load_library() -> pd.DataFrame:
    if not SRC.exists():
        sys.exit(f"missing {SRC}; stage the library there first (see data/source/aska/SOURCE.md)")
    d = pd.read_excel(SRC, dtype=str).fillna("")
    d.columns = ["jw_id_raw", "gene_name", "primer_n", "primer_c", "comment", "rrc_id"]
    d["jw_id"] = d["jw_id_raw"].str.replace(r"-AM$", "", regex=True)
    d["comment_jw"] = d["comment"].apply(lambda c: ";".join(JW_RE.findall(c)))
    d["n_pep"] = d["primer_n"].apply(n_peptide)
    d["c_pep"] = d["primer_c"].apply(c_peptide)
    return d


def build_indexes(p: pd.DataFrame, tags: pd.DataFrame):
    jw: dict[str, set[str]] = defaultdict(set)
    for ac, allt in zip(tags["uniprot_ac"], tags["locus_tag_all"].fillna("")):
        for t in JW_RE.findall(allt):
            jw[t].add(ac)
    gene: dict[str, set[str]] = defaultdict(set)
    for ac, g, syn in zip(p["uniprot_ac"], p["gene_name"].fillna(""), p["gene_synonyms"].fillna("")):
        for name in [g, *re.split(r"[;,\s]+", syn)]:
            if name:
                gene[name].add(ac)
    nt: dict[str, set[str]] = defaultdict(set)
    ct: dict[str, set[str]] = defaultdict(set)
    for ac, seq in zip(p["uniprot_ac"], p["sequence"]):
        nt[seq[1:1 + K]].add(ac)
        ct[seq[-K:]].add(ac)
    return jw, gene, nt, ct


def unique(s: set[str]) -> str:
    return next(iter(s)) if len(s) == 1 else ""


def n_ok(seq: str, pep: str | None) -> bool:
    return bool(pep) and seq[1:1 + K] == pep


def c_ok(seq: str, pep: str | None) -> bool:
    return bool(pep) and seq[-K:] == pep


def map_clones(lib: pd.DataFrame, p: pd.DataFrame, tags: pd.DataFrame) -> pd.DataFrame:
    jw_ix, gene_ix, nt_ix, ct_ix = build_indexes(p, tags)
    seq = p.set_index("uniprot_ac")["sequence"]
    names = {ac: {g, *re.split(r"[;,\s]+", s)} - {""} for ac, g, s in
             zip(p["uniprot_ac"], p["gene_name"].fillna(""), p["gene_synonyms"].fillna(""))}
    rows = []
    for r in lib.itertuples(index=False):
        jw_set = set(jw_ix.get(r.jw_id, set()))
        jw_via = "jw_id" if jw_set else ""
        if not jw_set and r.comment_jw:
            for t in r.comment_jw.split(";"):
                jw_set |= jw_ix.get(t, set())
            jw_via = "comment_jw" if jw_set else ""
        a_jw = unique(jw_set)
        a_gene = unique(gene_ix.get(r.gene_name, set()))
        joint = (nt_ix.get(r.n_pep, set()) & ct_ix.get(r.c_pep, set())) if r.n_pep and r.c_pep else set()
        a_pr = unique(joint)

        cand = {k: v for k, v in (("jw", a_jw), ("gene", a_gene), ("primer", a_pr)) if v}
        score = {ac: n_ok(seq[ac], r.n_pep) + c_ok(seq[ac], r.c_pep) for ac in set(cand.values())}
        chosen, rule = "", "unmapped"
        if cand:
            best = max(score.values())
            for k in ("jw", "gene", "primer"):  # tie-break order
                if k in cand and score[cand[k]] == best:
                    chosen, rule = cand[k], k
                    break
            if "jw" in cand and rule != "jw":
                rule = "primer_override"  # the primers positively prefer a protein the JW id did not name
        nn = c_ = False
        if chosen:
            nn, c_ = n_ok(seq[chosen], r.n_pep), c_ok(seq[chosen], r.c_pep)
        ends = int(nn) + int(c_)
        if not chosen:
            tier = "U"
        elif ends == 2:
            tier = "A"
        elif ends == 1:
            tier = "B"
        elif a_jw and a_jw == a_gene:
            tier = "C"
        else:
            tier = "D"

        flags = []
        if len(set(cand.values())) > 1:
            flags.append("routes_disagree")
        if chosen and len(joint) > 1 and chosen in joint:
            flags.append("primers_match_paralog")
        if chosen and not nn and r.n_pep and seq[chosen].find(r.n_pep) > 0:
            flags.append(f"n_primer_at_residue_{seq[chosen].find(r.n_pep) + 1}")
        if chosen and r.comment and not r.comment.startswith("Clone OK") or "<br" in r.comment:
            flags.append("library_comment")
        rows.append({
            "jw_id": r.jw_id, "aska_gene_name": r.gene_name, "comment": r.comment,
            "uniprot_jw": a_jw, "jw_via": jw_via, "uniprot_gene": a_gene,
            "uniprot_primer_joint": a_pr, "n_joint_primer_matches": len(joint),
            "n_primer_peptide": r.n_pep or "", "c_primer_peptide": r.c_pep or "",
            "uniprot_ac": chosen, "rule": rule, "n_end_ok": nn, "c_end_ok": c_,
            "n_ends_verified": ends, "confidence": tier,
            "name_agrees": bool(chosen) and r.gene_name in names.get(chosen, set()),
            "flags": ";".join(flags),
        })
    m = pd.DataFrame(rows)
    # a protein hit by several clones: Keio annotated them as separate ORFs; UniProt now has ONE entry
    sib = m[m["uniprot_ac"] != ""].groupby("uniprot_ac")["jw_id"].apply(list)
    m["shares_protein_with"] = [
        ";".join(j for j in sib.get(ac, []) if j != jw) if ac else "" for ac, jw in zip(m["uniprot_ac"], m["jw_id"])]
    m.loc[m["shares_protein_with"] != "", "flags"] = (
        m["flags"] + ";shares_protein_with_other_clone").str.lstrip(";")
    return m


def null_control(m: pd.DataFrame, lib: pd.DataFrame, p: pd.DataFrame, seed: int = 0) -> tuple[int, int]:
    """Re-run both end checks on a SHUFFLED mapping. A test that cannot fail proves nothing."""
    seq = p.set_index("uniprot_ac")["sequence"]
    ok = m[m["uniprot_ac"] != ""].merge(lib[["jw_id", "n_pep", "c_pep"]], on="jw_id")
    acs = list(ok["uniprot_ac"])
    random.Random(seed).shuffle(acs)
    return (sum(n_ok(seq[a], pep) for a, pep in zip(acs, ok["n_pep"])),
            sum(c_ok(seq[a], pep) for a, pep in zip(acs, ok["c_pep"])))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    ap.add_argument("--floor", type=float, default=0.90,
                    help="exit non-zero if fewer than this fraction of clones reach tier A or B")
    a = ap.parse_args()
    q = a.quiet

    log("aska.py  ASKA (-) clone library -> E. coli K-12 MG1655 proteome, with verification", q)
    log(f"  input  : {SRC.relative_to(REPO_ROOT)}", q)
    log(f"  outputs: {PROCESSED.relative_to(REPO_ROOT)}/aska_to_uniprot_ecoli.tsv (+ evidence/)", q)

    lib = load_library()
    full = proteomes.load_full(SPECIES)
    tags = proteomes.load_locus_tags(SPECIES)
    n_clones, n_proteome = len(lib), len(full)
    assert lib["jw_id"].is_unique, "duplicate JW ids in the library"
    log(f"\n  clones {n_clones:,} ({lib['jw_id'].nunique():,} distinct JW ids)   proteome {n_proteome:,}", q)
    log(f"  usable primers: N {int(lib['n_pep'].notna().sum()):,}  C {int(lib['c_pep'].notna().sum()):,}", q)

    m = map_clones(lib, full, tags)
    mapped = m["uniprot_ac"] != ""

    log("\n  ROUTES (proposals)", q)
    audit = []
    for name, col in (("jw", "uniprot_jw"), ("gene", "uniprot_gene"), ("primer_joint", "uniprot_primer_joint")):
        n = int((m[col] != "").sum())
        audit.append({"item": f"route:{name}", "clones": n, "pct_clones": round(100 * n / n_clones, 2)})
        log(f"    {name:<13}{n:>6,} clones  {100 * n / n_clones:5.1f}%", q)
    log("\n  CHOSEN BY: " + ", ".join(f"{k} {v:,}" for k, v in m["rule"].value_counts().items()), q)

    log("\n  VERIFICATION against the clone's own primers", q)
    for lab, col in (("N end", "n_end_ok"), ("C end", "c_end_ok")):
        n = int(m.loc[mapped, col].sum())
        log(f"    {lab} matches the chosen protein: {n:,} / {int(mapped.sum()):,}  ({100 * n / mapped.sum():.1f}%)", q)
    nn, cc = null_control(m, lib, full)
    log(f"    NULL CONTROL (accessions shuffled across clones): N {nn}  C {cc}   <- must be ~0", q)

    log("\n  CONFIDENCE", q)
    tiers = m["confidence"].value_counts().reindex(list("ABCDU"), fill_value=0)
    meaning = {"A": "both primer ends match", "B": "one end matches", "C": "no end; JW id = gene name",
               "D": "no end; routes disagree", "U": "not found"}
    for t, n in tiers.items():
        log(f"    {t}  {n:>6,}  {100 * n / n_clones:5.1f}%   {meaning[t]}", q)
        audit.append({"item": f"tier:{t}", "clones": int(n), "pct_clones": round(100 * n / n_clones, 2)})
    log(f"    gene name agrees with the chosen protein: {int(m['name_agrees'].sum()):,} / {int(mapped.sum()):,}", q)
    log(f"    routes disagree: {int(m['flags'].str.contains('routes_disagree').sum())}   "
        f"primer-override: {int((m['rule'] == 'primer_override').sum())}   "
        f"shares a protein with another clone: {int((m['shares_protein_with'] != '').sum())}", q)

    n_prot = m.loc[mapped, "uniprot_ac"].nunique()
    log("\n  COVERAGE", q)
    log(f"    clones found         {int(mapped.sum()):>6,} / {n_clones:,}   {100 * mapped.mean():.1f}%", q)
    log(f"    DISTINCT PROTEINS    {n_prot:>6,} / {n_proteome:,}   {100 * n_prot / n_proteome:.1f}% of the proteome", q)
    log(f"    clones not found     {int((~mapped).sum()):>6,}", q)
    audit.append({"item": "distinct_proteins", "clones": n_prot, "pct_clones": round(100 * n_prot / n_proteome, 2)})
    audit.append({"item": "null_control_N_pass", "clones": nn, "pct_clones": ""})
    audit.append({"item": "null_control_C_pass", "clones": cc, "pct_clones": ""})

    gname = full.set_index("uniprot_ac")["gene_name"]
    ltag = tags.set_index("uniprot_ac")["locus_tag"]
    gname = full.set_index("uniprot_ac")["gene_name"]
    ltag = tags.set_index("uniprot_ac")["locus_tag"]
    raw_id = lib.set_index("jw_id")["jw_id_raw"]
    clones = pd.DataFrame({
        "jw_id": m["jw_id"], "aska_id": m["jw_id"].map(raw_id), "aska_gene_name": m["aska_gene_name"],
        "uniprot_ac": m["uniprot_ac"], "uniprot_gene_name": m["uniprot_ac"].map(gname).fillna(""),
        "locus_tag": m["uniprot_ac"].map(ltag).fillna(""), "tier": m["confidence"],
        "n_ends_verified": m["n_ends_verified"], "rule": m["rule"],
        "shares_protein_with": m["shares_protein_with"], "flags": m["flags"],
    })
    level = {"A": 3, "B": 2, "C": 1, "D": 1}
    sub = clones[clones["uniprot_ac"] != ""].copy()
    sub["confidence"] = sub["tier"].map(level).astype(int)
    order = {ac: i for i, ac in enumerate(full["uniprot_ac"])}
    sub = sub.assign(_o=sub["uniprot_ac"].map(order)).sort_values(["_o", "jw_id"], kind="stable")
    deliver = sub[["uniprot_ac", "aska_id", "confidence"]].reset_index(drop=True)
    assert deliver["uniprot_ac"].isin(set(full["uniprot_ac"])).all(), "deliverable is not a subset of the proteome"
    assert deliver["aska_id"].is_unique and deliver["aska_id"].notna().all()
    log("\n  DELIVERABLE (subset of the proteome)", q)
    log(f"    rows (mapped clones)  {len(deliver):,}   distinct proteins {deliver['uniprot_ac'].nunique():,}"
        f"   all in proteome_{SPECIES}.tsv: yes", q)
    log("    confidence 3/2/1:     " + " / ".join(str(int((deliver['confidence'] == k).sum())) for k in (3, 2, 1)), q)
    OUT.mkdir(parents=True, exist_ok=True)
    deliver.to_csv(PROCESSED / "aska_to_uniprot_ecoli.tsv", sep="\t", index=False)
    clones.to_csv(OUT / "aska_clones_ecoli.tsv", sep="\t", index=False)
    m.to_csv(OUT / "aska_map_ecoli.tsv", sep="\t", index=False)
    pd.DataFrame(audit).to_csv(OUT / "aska_audit.tsv", sep="\t", index=False)

    rc = 0
    if (nn + cc) > 0.01 * mapped.sum():
        print(f"\nFAIL: the shuffled control passes {nn}+{cc} checks; the test cannot discriminate", file=sys.stderr)
        rc = 1
    ab = (tiers["A"] + tiers["B"]) / n_clones
    if ab < a.floor:
        print(f"\nFAIL: only {100 * ab:.1f}% of clones reach tier A/B < floor {100 * a.floor:.0f}%", file=sys.stderr)
        rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
