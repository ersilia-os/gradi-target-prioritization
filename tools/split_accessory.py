"""One-shot: split each task's `accessory/` into `evidence/` and `scratch/`.

`accessory/` became one name doing four jobs: controls, raw tool dumps, multi-GB caches and
secondary results all in one directory. The split asks one question per file --

    evidence/  would you CITE or CHECK it?   controls, audits, benchmarks, CV tables, manifests
    scratch/   would you DELETE it to reclaim space?  caches, shard dirs, raw tool output, smoke runs

-- which also makes `scratch/` safe to purge and pointless to upload to eosvc.

Classification is by explicit rule, printed for every file. Anything no rule matches is reported
as UNCLASSIFIED and left where it is rather than guessed at.

  python tools/split_accessory.py --dry-run
  python tools/split_accessory.py
"""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PROC = REPO / "data" / "processed"

# Ordered: the FIRST matching rule wins, so specific patterns precede general ones.
RULES: list[tuple[str, str, str]] = [
    # --- scratch: regenerable working data ------------------------------------------------
    (r"^\.",                        "scratch",  "hidden resume-key / fetch cache"),
    (r"^smoke_",                    "scratch",  "smoke-test output of a --limit run"),
    (r"^shards_",                   "scratch",  "resume cache"),
    (r"^(lazyqsar|geptop|deepgo|orthofinder|hits|deeplocpro_cache|tmbed_shards)$",
                                    "scratch",  "tool run/cache directory"),
    (r"_rpsblast_",                 "scratch",  "raw RPS-BLAST output"),
    (r"^eggnog_raw_",               "scratch",  "raw emapper output"),
    (r"^prott5_input_",             "scratch",  "input written for the worker"),
    (r"\.(faa|dmnd)$",              "scratch",  "sequence / DIAMOND database working file"),
    (r"^orthodb_(reps|gene_ogs|genes_ours|og2genes_ours|hits|groups_long)",
                                    "scratch",  "OrthoDB build intermediate"),
    (r"^(chembl|bindingdb)_(ligands|targets)$", "scratch", "per-ligand join intermediate"),

    # --- evidence: why the deliverable should be believed ---------------------------------
    (r"manifest",                   "evidence", "what ran, when, with what"),
    (r"(control|crosscheck|audit|decoy)", "evidence", "control / audit"),
    (r"^cv_|^cross_activator|^cutoff_sensitivity|^domain_bands|^oof_",
                                    "evidence", "cross-validation / benchmark table"),
    (r"^proteome_join$|^method_disagreement$|_pooling$",
                                    "evidence", "join or method-agreement control"),
    (r"^(cog_counts|deeplocpro_counts)", "evidence", "coverage summary"),
    (r"^model_",                    "evidence", "the fitted model actually shipped"),
    (r"^(labels_|deg_datasets|deg_genes|labeled_proteins|goslim_terms|registry|organism_class)",
                                    "evidence", "label / vocabulary table behind the deliverable"),
    (r"^(annotation_|locus_tags_)",  "evidence", "identity layer kept out of the 9-column table"),
    (r"^(orthogroups|scaffolds|chembl_hits|bindingdb_gain)$",
                                    "evidence", "secondary result quoted in the docs"),
    (r"^prott5_ecoli_control",      "evidence", "the UniProt validation embeddings"),
    (r"^(projection_pacmap|projection_umap)", "evidence", "rejected-alternative comparison"),
    (r"^orthodb_transfer_audit$|^orthodb_query$", "evidence", "transfer audit"),
    (r"^(chembl_targets|tmbed_labels|tmbed_topology|deeplocpro_probabilities)",
                                    "evidence", "per-protein detail behind the deliverable"),
]


# (file suffix, the exact f-string body) -> tier. For sites where the rules cannot be right:
# the name is entirely behind a variable, or a prefix means something different in this stage.
OVERRIDES: dict[tuple[str, str], str] = {
    ("essentiality/labels.py", "{label}.tsv"): "evidence",        # deg_datasets, deg_genes
    ("essentiality/deg_proteomes.py", "{label}.tsv"): "pre",      # {pre}proteome_join/labeled_*
    ("essentiality/geptop.py", "cv_{fasta.stem}_{tag}.pkl"): "scratch",   # COMPOSITION VECTOR
    ("essentiality/geptop.py", "db_{ref.stem}"): "scratch",       # BLAST database
    ("essentiality/geptop.py", "validate_{ds}"): "scratch",       # per-dataset run directory
    ("essentiality/geptop.py", "validation.tsv"): "evidence",     # the validation RESULT
    ("function/deepgo.py", "excluded_too_long.tsv"): "evidence",  # what was left out, and why
    ("function/deepgo.py", "run_{b:02d}"): "scratch",             # batch run directory
    ("function/deepgo.py", "yield_sweep.tsv"): "evidence",        # threshold sweep
    ("orthology/orthofinder.py", "fasta"): "scratch",             # inputs written for OrthoFinder
    ("proteomes/download.py", "locus_bridge_strains.tsv"): "evidence",
}


def classify(name: str) -> tuple[str | None, str]:
    stem = name[:-4] if name.endswith((".tsv", ".npz", ".faa", ".txt")) else name
    for pat, dest, why in RULES:
        if re.search(pat, name) or re.search(pat, stem):
            return dest, why
    return None, "no rule matched"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--code", action="store_true", help="rewrite the code, not the data")
    a = ap.parse_args()
    if a.code:
        rewrite_code(a.dry_run)
        return

    total = {"evidence": 0, "scratch": 0, "unclassified": 0}
    for acc in sorted(PROC.glob("*/accessory")):
        task = acc.parent.name
        print(f"\n### {task}")
        for p in sorted(acc.iterdir()):
            dest, why = classify(p.name)
            if dest is None:
                total["unclassified"] += 1
                print(f"  {'UNCLASSIFIED':>10}  {p.name}")
                continue
            total[dest] += 1
            print(f"  {dest:>10}  {p.name:<46} {why}")
            if not a.dry_run:
                d = acc.parent / dest
                d.mkdir(parents=True, exist_ok=True)
                p.rename(d / p.name)
        if not a.dry_run and not any(acc.iterdir()):
            acc.rmdir()
    print(f"\nevidence {total['evidence']}  scratch {total['scratch']}  "
          f"UNCLASSIFIED {total['unclassified']}")
    if total["unclassified"]:
        print("Unclassified files were NOT moved. Add a rule for each before re-running.")




# =============================================================================================
# CODE side: repoint every `ACC_DIR / "x"` at evidence/ or scratch/ using the SAME rules above,
# so the directory a file lands in and the directory the code looks in cannot disagree.
# =============================================================================================

DEF_RE = re.compile(r"^(\s*)(ACC_DIR|ACCESSORY_DIR)(\s*)=(\s*)(.+?)\s*/\s*\"accessory\"\s*$", re.M)
USE_RE = re.compile(r"\b(ACC_DIR|ACCESSORY_DIR)\s*/\s*(f?)\"([^\"]*)\"")


def _literal(fstring_body: str) -> tuple[str, bool]:
    """(the part we can classify on, is_pre_prefixed).

    `{pre}` is the smoke prefix: empty on a real run, "smoke_" under --limit. A file written
    through it is therefore evidence OR scratch depending on the run, so those sites get a runtime
    choice rather than a fixed directory.
    """
    pre = fstring_body.startswith("{pre}")
    body = fstring_body[5:] if pre else fstring_body
    # classify on the text up to the first interpolation, e.g. "cog_counts_" from "cog_counts_{sp}"
    lit = re.sub(r"\{[^}]*\}", "", body)
    return ("" if lit.lstrip(".").strip() in ("", "tsv", "npz", "faa", "pkl") else lit), pre


def rewrite_code(dry: bool) -> None:
    print("\n" + "=" * 92 + "\nCODE\n" + "=" * 92)
    files = [p for p in list(REPO.glob("scripts/**/*.py")) + list(REPO.glob("src/*.py"))
             if "legacy" not in p.parts]
    leftovers: list[str] = []
    for p in sorted(files):
        t = orig = p.read_text()
        if "ACC_DIR" not in t and "ACCESSORY_DIR" not in t:
            continue
        shown = False

        def head(name=p):
            nonlocal shown
            if not shown:
                print(f"\n### {name.relative_to(REPO)}")
                shown = True

        # 1. definitions -> both directories
        def def_sub(m):
            head()
            indent, _, _, _, base = m.groups()
            print(f"  def   {base} / accessory  ->  evidence + scratch")
            return (f"{indent}EVIDENCE_DIR = {base} / \"evidence\"\n"
                    f"{indent}SCRATCH_DIR = {base} / \"scratch\"")
        t = DEF_RE.sub(def_sub, t)

        # 2. uses -> the tier the filename belongs to
        def use_sub(m):
            head()
            _, f, body = m.groups()
            lit, pre = _literal(body)
            ov = next((v for (suf, b), v in OVERRIDES.items()
                       if body == b and str(p).endswith(suf)), None)
            if ov == "pre":
                dest, why, pre = "evidence", "override (smoke-aware)", True
            elif ov:
                dest, why = ov, "override"
            else:
                dest, why = classify(lit)
            if dest is None:
                leftovers.append(f"{p.relative_to(REPO)}: {m.group(0)}  ({lit!r})")
                return m.group(0)
            if pre and dest == "scratch":
                print(f"  use   {body!r}  ->  scratch")
                return f"SCRATCH_DIR / {f}\"{body}\""
            if pre:
                # evidence on a real run, scratch under --limit: decided where `pre` is
                print(f"  use   {body!r}  ->  (SCRATCH if pre else {dest.upper()})")
                return f"(SCRATCH_DIR if pre else {dest.upper()}_DIR) / {f}\"{body}\""
            print(f"  use   {body!r}  ->  {dest}   [{why}]")
            return f"{dest.upper()}_DIR / {f}\"{body}\""
        t = USE_RE.sub(use_sub, t)

        if t != orig and not dry:
            p.write_text(t)

    # 3. anything left: mkdir, glob, relative_to, parent, and any unmatched literal
    rest: list[str] = []
    for p in sorted(files):
        for i, line in enumerate(p.read_text().splitlines(), 1):
            if re.search(r"\b(ACC_DIR|ACCESSORY_DIR)\b", line):
                rest.append(f"  {p.relative_to(REPO)}:{i}  {line.strip()}")
    if leftovers:
        print("\nUNCLASSIFIED USES (left alone -- add a rule):")
        print("\n".join("  " + x for x in leftovers))
    if rest:
        print(f"\nREMAINING references to fix by hand ({len(rest)}):")
        print("\n".join(rest[:60]))

if __name__ == "__main__":
    main()
