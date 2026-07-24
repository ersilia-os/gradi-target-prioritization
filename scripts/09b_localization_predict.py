"""Stage 09b — predict subcellular localization with PSORTb 3.0 (Gram-negative) to complete the
Clp-accessibility axis (docs §5.1b). UniProt-curated calls (09a) take precedence; PSORTb fills the rest,
so every protein gets a localization and a Clp-accessibility score.

Runs the brinkmanlab/psortb_commandline Docker image on each proteome FASTA (terse output), then merges
onto the 09a localization TSV. Requires Docker. Run with the `gradi` env:
    python scripts/09b_localization_predict.py
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.ligandability import ORGANISMS, REPO_ROOT  # noqa: E402

IMAGE = "brinkmanlab/psortb_commandline:1.0.2"
# PSORTb Gram-negative final localization -> (our category, Clp-accessibility 0-1 or None)
PSORTB_MAP = {
    "cytoplasmic": ("cytoplasm", 1.0),
    "cytoplasmicmembrane": ("inner_membrane", 0.5),
    "periplasmic": ("periplasm", 0.0),
    "outermembrane": ("outer_membrane", 0.0),
    "extracellular": ("extracellular", 0.0),
    "cellwall": ("extracellular", 0.0),   # (Gram+ only; harmless)
}
# accessibility for our own localization categories (matches scripts/09a_localization.classify)
ACCESS = {"cytoplasm": 1.0, "inner_membrane": 0.5, "membrane": 0.5,
          "periplasm": 0.0, "outer_membrane": 0.0, "secreted": 0.0, "extracellular": 0.0}


def parse_result(path: Path) -> dict:
    """Parse a PSORTb terse result -> {accession: (category, clp_accessibility_or_None, score)}."""
    out: dict[str, tuple] = {}
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        next(rdr, None)  # header
        for row in rdr:
            if len(row) < 3:
                continue
            seqid, loc_raw, score = row[0], row[1].strip(), row[2].strip()
            acc = seqid.split("|")[1] if "|" in seqid else seqid.split()[0]
            key = loc_raw.lower().replace(" ", "")
            # multiple/unknown localizations -> unknown (blank accessibility)
            cat, acc_score = PSORTB_MAP.get(key, ("unknown", None))
            try:
                sc = float(score)
            except ValueError:
                sc = None
            out[acc] = (cat, acc_score, sc)
    return out


def run_psortb(organism: str, reuse: bool = False) -> dict:
    """{accession: (localization_category, clp_accessibility_or_None, psortb_score)} via PSORTb.

    With reuse=True, a pre-existing non-empty *_psortb_gramneg.txt is parsed instead of re-running
    Docker — used to pick up a run already completed out-of-band (e.g. a long emulated background run).
    """
    pdir = REPO_ROOT / "data" / "raw" / organism / "proteome"
    fastas = sorted(pdir.glob("*.fasta"))
    if not fastas:
        print(f"  [{organism}] no proteome FASTA — skipping PSORTb")
        return {}
    fasta = fastas[0]
    outdir = REPO_ROOT / "data" / "raw" / organism / "localization" / "psortb"
    outdir.mkdir(parents=True, exist_ok=True)
    existing = [p for p in sorted(outdir.glob("*_psortb_gramneg.txt")) if p.stat().st_size > 0]
    if reuse:
        # never launch Docker in reuse mode: parse an existing result if present, else UniProt-only
        if existing:
            print(f"  [{organism}] reusing existing PSORTb result {existing[0].name}")
            out = parse_result(existing[0])
            print(f"  [{organism}] PSORTb localized {len(out)} proteins (reused)")
            return out
        print(f"  [{organism}] no PSORTb result to reuse — UniProt-only localization")
        return {}
    for old in outdir.glob("*_psortb_*.txt"):
        old.unlink()
    print(f"  [{organism}] running PSORTb (Gram-negative) on {fasta.name} — this can take a while under emulation…")
    cmd = [
        "docker", "run", "--rm", "--entrypoint", "bash",
        "-v", f"{fasta.parent}:/in", "-v", f"{outdir}:/tmp/results", IMAGE,
        "-c", f"/usr/local/psortb/bin/psort -n -o terse -i /in/{fasta.name} > /tmp/psort.log 2>&1",
    ]
    subprocess.run(cmd, check=True)
    results = sorted(outdir.glob("*_psortb_gramneg.txt"))
    if not results:
        print(f"  [{organism}] PSORTb produced no result file — skipping")
        return {}
    out = parse_result(results[0])
    print(f"  [{organism}] PSORTb localized {len(out)} proteins")
    return out


def merge(organism: str, psortb: dict) -> None:
    """Overwrite <prefix>_localization.tsv: UniProt-curated wins; PSORTb fills the rest."""
    _, prefix = ORGANISMS[organism]
    tsv = REPO_ROOT / "data" / "raw" / organism / "localization" / f"{prefix}_localization.tsv"
    if not tsv.exists():
        print(f"  [{organism}] no 09a localization TSV — run 09a first; skipping")
        return
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    n_uni = n_ps = n_none = 0
    out_rows = []
    for r in rows:
        acc = r["uniprot_accession"]
        uni_loc = (r.get("localization") or "").strip()
        source = r.get("localization_source", "")
        psortb_score = ""
        if uni_loc and uni_loc != "unknown":
            loc, acc_score, source = uni_loc, r.get("clp_accessibility", ""), "uniprot"
            n_uni += 1
        elif acc in psortb and psortb[acc][0] != "unknown":
            loc, acc_score, sc = psortb[acc]
            source, psortb_score = "psortb", ("" if sc is None else sc)
            acc_score = "" if acc_score is None else acc_score
            n_ps += 1
        else:
            loc, acc_score, source = "unknown", "", "none"
            n_none += 1
        out_rows.append({
            "uniprot_accession": acc, "localization": loc,
            "clp_accessibility": acc_score if acc_score != "" else ("" if loc not in ACCESS else ACCESS[loc]),
            "has_signal_peptide": r.get("has_signal_peptide", "0"),
            "n_transmembrane": r.get("n_transmembrane", "0"),
            "localization_source": source, "psortb_score": psortb_score,
            "subcellular_raw": r.get("subcellular_raw", ""),
        })
    cols = ["uniprot_accession", "localization", "clp_accessibility", "has_signal_peptide",
            "n_transmembrane", "localization_source", "psortb_score", "subcellular_raw"]
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)
    print(f"  [{organism}] merged -> {tsv.relative_to(REPO_ROOT)} "
          f"(uniprot {n_uni}, psortb {n_ps}, unknown {n_none} of {len(out_rows)})")


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reuse-existing", action="store_true",
                    help="parse a pre-existing non-empty PSORTb result instead of re-running Docker")
    args = ap.parse_args()
    for organism in ("kpneumoniae", "ecoli"):
        ps = run_psortb(organism, reuse=args.reuse_existing)
        merge(organism, ps)
    print("\nDone.")


if __name__ == "__main__":
    main()
