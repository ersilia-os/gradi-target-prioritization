"""Stage 09e — signal-peptide typing and lipoprotein sorting (docs §5.1e).

The point of this track is a single, specific correction: **outer-membrane lipoproteins**. They
carry a Sec/SPII signal peptide, are lipid-anchored to the inner leaflet of the outer membrane, and
are routinely mis-binned as `periplasmic` by rule-based predictors — which would wrongly give them
a non-zero Clp-accessibility. Resolving them needs two things:

  1. the signal-peptide **type** (Sec/SPI vs Sec/SPII vs Tat/SPI vs Tat/SPII vs Sec/SPIII), which
     is what SignalP 6.0 adds over every earlier version;
  2. the Lol pathway **"+2 rule"** — an Asp two residues after the lipidated cysteine retains the
     lipoprotein in the inner membrane; anything else routes it to the outer membrane
     (`src/localization.py:lipoprotein_sorting`).

SignalP 6.0 is licensed software: academic users must accept the terms at
https://services.healthtech.dtu.dk/services/SignalP-6.0/ and download the package themselves, so it
cannot be installed unattended. Point this script at the binary with `SIGNALP6_BIN` (the same
env-var pattern as `FPOCKET_BIN`/`P2RANK_DIR` elsewhere in the pipeline).

**Without the binary the script still runs.** It falls back to the lipobox motif
(`src/localization.py:find_lipobox` — Prosite PS51257 plus the charged-n-region and hydrophobic
h-region filters that make it usable) to call Sec/SPII lipoproteins, and records
`signalp_source=lipobox` so downstream code and the eventual write-up can tell the two apart.
Measured against UniProt's 99 lipid-anchor entries for E. coli K-12 the fallback runs at
precision 0.74 / recall 0.82. It finds lipoproteins only — it does not type Tat or SPIII — so it
is a stopgap, not a substitute for SignalP.

Output: output/results/<organism>/<prefix>_loc_signalp.csv
Run with the `gradi` conda env.
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

SIGNALP6_BIN = os.environ.get("SIGNALP6_BIN", shutil.which("signalp6") or "")

# SignalP 6.0 `prediction_results.txt` prediction labels -> our compact type.
SP_TYPES = {
    "OTHER": "none",
    "SP": "sec_spi",
    "LIPO": "sec_spii",
    "TAT": "tat_spi",
    "TATLIPO": "tat_spii",
    "PILIN": "sec_spiii",
}
LIPOPROTEIN_TYPES = {"sec_spii", "tat_spii"}


def run_signalp(fasta: Path, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        SIGNALP6_BIN,
        "--fastafile", str(fasta),
        "--organism", "other",          # bacteria (SignalP's non-eukaryote mode)
        "--output_dir", str(out_dir),
        "--format", "none",             # skip the per-protein plots
        "--mode", "fast",
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out_dir / "prediction_results.txt"


def parse_signalp(path: Path) -> dict[str, dict]:
    """Parse SignalP 6.0's tab-separated prediction_results.txt."""
    out: dict[str, dict] = {}
    with open(path) as fh:
        header: list[str] = []
        for line in fh:
            if line.startswith("# ID"):
                header = [c.strip() for c in line.lstrip("# ").rstrip("\n").split("\t")]
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts)) if header else {}
            acc = LOC.acc_from_header(parts[0])
            pred = (rec.get("Prediction") or parts[1]).strip()
            sp_type = SP_TYPES.get(pred.split("(")[0].strip(), "none")
            # "CS pos: 22-23. Pr: 0.87" -> mature protein starts at 23
            cs = rec.get("CS Position", "") or ""
            m = re.search(r"CS pos:\s*(\d+)-(\d+)", cs)
            prob = re.search(r"Pr:\s*([0-9.]+)", cs)
            out[acc] = {
                "sp_type": sp_type,
                "sp_prob": float(prob.group(1)) if prob else None,
                "sp_cleavage_pos": int(m.group(2)) if m else None,
            }
    return out


def lipobox_fallback(seqs: dict[str, str]) -> dict[str, dict]:
    """Call Sec/SPII lipoproteins from the lipobox motif alone."""
    out: dict[str, dict] = {}
    for acc, seq in seqs.items():
        cys = LOC.find_lipobox(seq)
        out[acc] = ({"sp_type": "sec_spii", "sp_prob": None, "sp_cleavage_pos": cys} if cys
                    else {"sp_type": "none", "sp_prob": None, "sp_cleavage_pos": None})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--limit", type=int, default=None, help="first N proteins (debug)")
    ap.add_argument("--signalp-bin", default=SIGNALP6_BIN,
                    help="path to the licensed signalp6 executable (or set SIGNALP6_BIN)")
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]

    accs = LOC.load_accessions(org)
    if args.limit:
        accs = accs[: args.limit]
    all_seqs = LOC.read_fasta(LOC.proteome_fasta(org))
    seqs = {a: all_seqs[a] for a in accs if all_seqs.get(a)}

    binary = args.signalp_bin
    if binary and Path(binary).exists():
        print(f"[{org}] running SignalP 6.0 ({binary}) on {len(seqs)} proteins ...", flush=True)
        with tempfile.TemporaryDirectory() as td:
            fa = Path(td) / "in.fasta"
            LOC.write_fasta(seqs, fa)
            res = run_signalp(fa, LOC.localization_processed_dir(org, "signalp6"))
            calls = parse_signalp(res)
        source = "signalp6"
    else:
        print(f"[{org}] SignalP 6.0 not available (SIGNALP6_BIN unset or missing) — "
              f"falling back to the lipobox motif; Tat/SPIII types will not be called.", flush=True)
        calls = lipobox_fallback(seqs)
        source = "lipobox"

    rows = []
    for acc in accs:
        rec = calls.get(acc, {"sp_type": "none", "sp_prob": None, "sp_cleavage_pos": None})
        is_lipo = rec["sp_type"] in LIPOPROTEIN_TYPES
        sorting = None
        if is_lipo and rec["sp_cleavage_pos"]:
            sorting = LOC.lipoprotein_sorting(seqs.get(acc, ""), rec["sp_cleavage_pos"])
        rows.append({
            "uniprot_accession": acc,
            "sp_type": rec["sp_type"],
            "sp_prob": rec["sp_prob"],
            "sp_cleavage_pos": rec["sp_cleavage_pos"],
            "is_lipoprotein": is_lipo,
            "lipoprotein_sorting": sorting,
            "signalp_source": source if acc in calls else "none",
        })

    df = pd.DataFrame(rows)
    out = LOC.results_dir(org) / f"{prefix}_loc_signalp.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(df)} proteins)", flush=True)
    print(f"[{org}] sp types: {df['sp_type'].value_counts().to_dict()}", flush=True)
    print(f"[{org}] lipoprotein sorting: {df['lipoprotein_sorting'].value_counts().to_dict()}",
          flush=True)


if __name__ == "__main__":
    main()
