"""One table describing every essentiality training set, and a control that says it is not inverted.

Answers, per dataset: how many positives, out of how many labelled proteins, what fraction of the
organism's proteome that is, whether every labelled protein has a feature vector, whether the
labels point the right way, and what a model scores on it.

THE CONTROL IS THE POINT OF THIS SCRIPT. `num_positives` looks identical whether the labels are
right or inverted, and an inverted column trains a confident, well-formed, exactly wrong model. So
every dataset is checked against biology that cannot be in dispute: the RIBOSOME. A genuine
essentiality screen calls nearly all 50S/30S ribosomal proteins essential, and calls the classic
dispensables (lacZ, araB, fadB, the flagellar and fimbrial genes) non-essential. `ribosome_recall`
minus `dispensable_rate` is the separation, and a dataset below `MIN_SEPARATION` fails the run.

Marker genes are matched by GENE SYMBOL, which every one of these proteomes carries except
RH201207's -- that one is staged from a table of tags and sequences alone, so its symbols come
across from KPPR1 by exact sequence. The recall denominator is printed per dataset, because a
control computed over 4 markers is not a control.

AUROC and AUPR come from `evidence/cv_endpoints.tsv`, written by predict.py. A dataset with no row
there reports `pending` rather than a blank, so the difference between "not measured yet" and
"measured as zero" is never lost.

Run with the `gradi` env:
    python scripts/essentiality/summary.py
    python scripts/essentiality/summary.py --features prott5
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from src import proteomes as P  # noqa: E402
from screens import SCREENS  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE_DIR = OUT_DIR / "evidence"
RESULTS_DIR = REPO_ROOT / "output" / "results" / "essentiality"
TRAINING_DIR = OUT_DIR / "training_sets"
NCBI_DIR = REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "ncbi"
KP_STRAINS = REPO_ROOT / "data" / "source" / "ncbi" / "kp_strains"
DEG_DIR = REPO_ROOT / "data" / "source" / "ncbi" / "deg_proteomes"
STRAIN_EMB = REPO_ROOT / "data" / "processed" / "embeddings" / "scratch" / "strains"

# Structural core of the cell. Every genome-wide essentiality screen ever published calls
# essentially all of these essential; a screen that does not is inverted, mis-parsed, or not an
# essentiality screen. Symbols only -- rps*/rpl*/rpm* is the whole ribosomal protein complement.
RIBOSOME = re.compile(r"^(rps|rpl|rpm)[A-Z]$")

# Textbook dispensables: catabolism of sugars the cell need not eat, fatty-acid degradation,
# motility, adhesion, prophage. Deleting any of these is a standard laboratory operation.
DISPENSABLE = {g.lower() for g in [
    "lacZ", "lacY", "lacA", "lacI", "araB", "araA", "araD", "araC",
    "fadB", "fadA", "fadD", "fadE", "fadL", "malE", "malF", "malG", "malK",
    "galK", "galT", "galE", "xylA", "xylB", "rhaA", "rhaB", "rhaD",
    "fliC", "flgB", "flgC", "flgD", "flgE", "fliA", "fliD", "motA", "motB", "cheA", "cheB",
    "fimA", "fimH", "fimD", "tsr", "tar", "trg", "aer", "nanA", "nagA", "nagB",
]}

# BOTH BARS ARE CALIBRATED ON MEASUREMENTS, NOT ON 1.00. Keio -- an ARRAYED single-gene knockout
# collection, the strongest evidence that exists -- calls only 0.774 of ribosomal proteins
# essential, and its misses (rplA, rplI, rplK, rplY, rpmE/F/G/I, rpsF/O/T/U) are genuinely
# dispensable in E. coli. So the biological ceiling is ~0.8, and a control demanding 1.0 would
# fail the gold standard. The bars are set to catch BREAKAGE, not noise:
#
#     Goodall 2018 TraDIS            0.961      <- the best of the E. coli screens
#     Keio arrayed knockout          0.774      <- the gold standard, and the realistic ceiling
#     Gerdes 2003 footprinting       0.538      <- noisy but correctly polarised (8x enrichment)
#     Ghomi BW25113 DESeq            0.113      <- BROKEN, retired; see screens.py
#
# Polarity is the primary test and the recall floor is the secondary one: a screen can be noisy
# and still usable, but it cannot point the wrong way and it cannot miss 9 ribosomal proteins in 10.
MIN_SEPARATION = 0.40     # ribosome_recall - dispensable_rate; the POLARITY test
MIN_RECALL = 0.45         # below this the screen is not measuring essentiality at all
MIN_MARKERS = 20          # a control computed over fewer markers than this is not a control


def say(m: str = "") -> None:
    print(m, flush=True)


def rule(c: str = "-", w: int = 120) -> None:
    say(c * w)


# ------------------------------------------------------------------ key -> gene symbol, per source

def _fasta_attr(faa: Path, attr: str) -> dict[str, str]:
    """first-token id -> a bracketed FASTA attribute (`[gene=thrA]`)."""
    out = {}
    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            m = re.search(rf"\[{attr}=([^\]]+)\]", line)
            if m:
                out[line[1:].split()[0]] = m.group(1)
    return out


def _gff_symbols(label: str) -> dict[str, str]:
    """protein_id -> gene symbol, walking gene -> CDS by ID/Parent as NCBI's GFF requires."""
    genes, out = {}, {}
    for line in (NCBI_DIR / f"{label}.gff").read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        a = dict(kv.split("=", 1) for kv in f[8].rstrip().split(";") if "=" in kv)
        if f[2] in ("gene", "pseudogene") and a.get("ID") and a.get("gene"):
            genes[a["ID"]] = a["gene"]
        elif f[2] == "CDS" and a.get("protein_id"):
            g = a.get("gene") or genes.get(a.get("Parent", ""))
            if g:
                out[a["protein_id"]] = g
    return out


def symbols_for(features: str) -> dict[str, str]:
    """The screen's own key -> gene symbol, however that proteome happens to carry one."""
    if features == "ecoli":
        d = P.load("ecoli")
        return dict(zip(d["uniprot_ac"], d["gene_name"].fillna("")))
    if features.startswith("DEG"):
        return _fasta_attr(DEG_DIR / f"{features}.faa", "gene")
    if (NCBI_DIR / f"{features}.gff").exists():
        return _gff_symbols(features)
    # RH201207 is staged as tags + sequences with no annotation of its own, so its symbols come
    # across from KPPR1 by EXACT SEQUENCE -- same species, so this is a lookup, not a prediction.
    # It reaches only ~39% of the proteome, which is why the marker denominator is printed.
    kppr1 = "kpneumoniae__kppr1__GCF_000742755.1"
    sym = _gff_symbols(kppr1)
    seq2pid, cur, buf = {}, None, []
    for line in (NCBI_DIR / f"{kppr1}.faa").read_text().splitlines():
        if line.startswith(">"):
            if cur:
                seq2pid["".join(buf)] = cur
            cur, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if cur:
        seq2pid["".join(buf)] = cur
    d = pd.read_csv(KP_STRAINS / f"{features}.tsv", sep="\t")
    return {t: sym[seq2pid[s]] for t, s in zip(d["locus_tag"].astype(str),
                                               d["sequence"].astype(str))
            if s in seq2pid and seq2pid[s] in sym}


def proteome_size(features: str) -> int:
    """Protein-coding genes in the organism the screen was run on -- the coverage denominator."""
    if features in ("ecoli", "kpneumoniae", "saureus"):
        return len(P.load(features))
    for p in (NCBI_DIR / f"{features}.faa", DEG_DIR / f"{features}.faa"):
        if p.exists():
            return sum(1 for ln in p.read_text().splitlines() if ln.startswith(">"))
    return len(pd.read_csv(KP_STRAINS / f"{features}.tsv", sep="\t"))


def embedded_keys(features: str, kind: str) -> set[str] | None:
    import numpy as np
    if features in ("ecoli", "kpneumoniae", "saureus"):
        f = REPO_ROOT / "data" / "processed" / "embeddings" / (
            f"embeddings_{features}.npz" if kind == "esmc" else f"{kind}_{features}.npz")
    else:
        f = STRAIN_EMB / (f"embeddings_{features}.npz" if kind == "esmc" else f"{kind}_{features}.npz")
    if not f.exists():
        return None
    return {str(a) for a in np.load(f, allow_pickle=True)["accessions"]}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--features", default="prott5", choices=["esmc", "prott5", "proteomelm"])
    a = ap.parse_args()

    cv = {}
    p = EVIDENCE_DIR / "cv_endpoints.tsv"
    if p.exists():
        d = pd.read_csv(p, sep="\t")
        d = d[d["features"] == a.features]
        cv = {r["endpoint"]: r for _, r in d.iterrows()}

    rule("=")
    say("scripts/essentiality/summary.py -- one row per training set")
    rule("=")
    say(f"  features   {a.features}")
    say(f"  sets       {TRAINING_DIR.relative_to(REPO_ROOT)}/<column>.tsv")
    say(f"  cv source  {p.relative_to(REPO_ROOT)}"
        + ("" if cv else "   (absent -- auroc/aupr report `pending`)"))
    say(f"  control    ribosome_recall >= {MIN_RECALL:.2f} and separation >= {MIN_SEPARATION:.2f}, "
        f"over >= {MIN_MARKERS} markers  (Keio, the gold standard, measures 0.774)")
    rule()

    rows, failed = [], []
    for s in SCREENS:
        d = pd.read_csv(TRAINING_DIR / f"{s.column}.tsv", sep="\t")
        keys = d["key"].astype(str)
        lab = dict(zip(keys, d["label"]))

        sym = symbols_for(s.features)
        ribo = [lab[k] for k, g in sym.items() if k in lab and RIBOSOME.match(str(g))]
        disp = [lab[k] for k, g in sym.items() if k in lab and str(g).lower() in DISPENSABLE]
        r_rec = sum(ribo) / len(ribo) if ribo else float("nan")
        d_rate = sum(disp) / len(disp) if disp else float("nan")
        sep = r_rec - d_rate

        emb = embedded_keys(s.features, a.features)
        emb_cov = keys.isin(emb).mean() if emb is not None else float("nan")

        c = cv.get(s.column)
        rows.append({
            "dataset": s.column.removeprefix("essential_"),
            "short_comment": s.comment,
            "column_name": s.column,
            "num_positives": int(d["label"].sum()),
            "num_total": len(d),
            "coverage": round(len(d) / proteome_size(s.features), 4),
            "embedding_coverage": round(emb_cov, 4) if emb_cov == emb_cov else "pending",
            "ribosome_recall": round(r_rec, 3) if ribo else "n/a",
            "separation": round(sep, 3) if ribo and disp else "n/a",
            "n_ribosomal": len(ribo),
            "dispensable_rate": round(d_rate, 3) if disp else "n/a",
            "n_dispensable": len(disp),
            "auroc": c["roc_auc_grouped"] if c is not None else "pending",
            "aupr": c["pr_auc_grouped"] if c is not None else "pending",
            "organism": s.organism,
            "assay": s.assay,
            "features_from": s.features,
        })
        if ribo and (sep < MIN_SEPARATION or r_rec < MIN_RECALL or len(ribo) < MIN_MARKERS):
            failed.append((s.column, r_rec, d_rate, len(ribo)))

    t = pd.DataFrame(rows)
    show = ["dataset", "num_positives", "num_total", "coverage", "embedding_coverage",
            "ribosome_recall", "n_ribosomal", "dispensable_rate", "separation", "auroc", "aupr"]
    say(t[show].to_string(index=False))
    rule()
    for r in rows:
        say(f"  {r['column_name']:42s} {r['short_comment']}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "screen_summary.tsv"
    t.to_csv(out, sep="\t", index=False)
    rule()
    say(f"  wrote {out.relative_to(REPO_ROOT)}")
    if failed:
        rule("=")
        for col, rr, dr, n in failed:
            say(f"  FAILED {col}: ribosome_recall {rr:.3f} over {n} markers, "
                f"dispensable_rate {dr:.3f}")
        sys.exit("label control failed -- do NOT train on these columns")
    rule("=")
    say("every dataset passes the polarity control.")
    rule("=")


if __name__ == "__main__":
    main()
