"""Stage 09a — subcellular localization + Clp-accessibility (docs §5).

Fetches UniProt subcellular-location / signal-peptide / transmembrane annotation for the focal
proteomes and derives, per protein:
  - localization  : cytoplasm / inner_membrane / periplasm / outer_membrane / secreted /
                    extracellular / membrane / unknown  (parsed from cc_subcellular_location,
                    backed off to transmembrane/signal-peptide features when the text is silent)
  - clp_accessibility [0-1] : reachability by the cytoplasmic Clp protease machinery — the gating
                    requirement for BacPROTAC / targeted degradation. cytoplasm 1.0, inner-membrane
                    (cytoplasm-facing) 0.5, periplasm / OM / secreted / extracellular 0.0, unknown -> None.
  - has_signal_peptide (bool), n_transmembrane (int)

Writes data/raw/<organism>/localization/<prefix>_localization.tsv (eosvc-tracked, not Git).
Run with the `gradi` conda env. Network fetch (UniProt REST stream), same endpoint as 00a.
"""
from __future__ import annotations

import argparse
import csv
import io
from pathlib import Path
from urllib.parse import quote

import requests
from tenacity import retry, stop_after_attempt, wait_exponential

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,cc_subcellular_location,ft_signal,ft_transmem"

# label -> (proteome id, organism folder, prefix)
ORGS = {
    "HS11286": {"id": "UP000007841", "organism": "kpneumoniae", "prefix": "kp"},
    "EcoliK12": {"id": "UP000000625", "organism": "ecoli", "prefix": "ec"},
}
REPO_ROOT = Path(__file__).resolve().parents[1]


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def download(url: str) -> str:
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    return r.text


def classify(subcell: str, has_sig: bool, n_tm: int):
    """(localization category, clp_accessibility[0-1] or None)."""
    s = (subcell or "").lower()
    def has(*ks):
        return any(k in s for k in ks)
    if has("outer membrane", "outer-membrane"):
        return "outer_membrane", 0.0
    if has("secreted"):
        return "secreted", 0.0
    if has("periplasm"):
        return "periplasm", 0.0
    if has("extracellular", "cell surface", "fimbri", "pilus", "pili", "flagell", "capsule"):
        return "extracellular", 0.0
    if has("inner membrane", "cytoplasmic membrane", "plasma membrane", "cell membrane", "cell inner membrane"):
        return "inner_membrane", 0.5
    if has("cytoplasm", "cytosol"):
        return "cytoplasm", 1.0
    if has("membrane"):
        return "membrane", 0.5
    # no explicit CC text -> back off to sequence features
    if n_tm and n_tm > 0:
        return "inner_membrane", 0.5
    if has_sig:
        return "secreted", 0.0
    return "unknown", None


def fetch_one(label: str) -> None:
    spec = ORGS[label]
    out_dir = REPO_ROOT / "data" / "raw" / spec["organism"] / "localization"
    out_dir.mkdir(parents=True, exist_ok=True)
    query = quote("proteome:" + spec["id"], safe=":")
    url = f"{STREAM}?compressed=false&format=tsv&query={query}&fields={FIELDS}"
    print(f"Downloading {label} ({spec['id']}) localization from UniProt ...")
    tsv = download(url)
    rows = list(csv.DictReader(io.StringIO(tsv), delimiter="\t"))
    # UniProt tsv headers are human-readable; find columns positionally-safe by fuzzy match
    def col(row, *names):
        for n in names:
            for k in row:
                if k.lower().replace(" ", "").startswith(n):
                    return row[k]
        return ""
    out_path = out_dir / f"{spec['prefix']}_localization.tsv"
    n_acc = 0
    counts: dict[str, int] = {}
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["uniprot_accession", "localization", "clp_accessibility",
                    "has_signal_peptide", "n_transmembrane", "subcellular_raw"])
        for r in rows:
            acc = col(r, "entry").strip() or list(r.values())[0].strip()
            subcell = col(r, "subcellular")
            sig = col(r, "signal")
            tm = col(r, "transmembrane")
            has_sig = bool(sig and sig.strip())
            n_tm = tm.upper().count("TRANSMEM") if tm else 0
            loc, acc_score = classify(subcell, has_sig, n_tm)
            counts[loc] = counts.get(loc, 0) + 1
            w.writerow([acc, loc, "" if acc_score is None else acc_score,
                        int(has_sig), n_tm, (subcell or "").replace("\t", " ")])
            n_acc += 1
    print(f"  wrote {n_acc} rows -> {out_path.relative_to(REPO_ROOT)}")
    print(f"  localization: {dict(sorted(counts.items(), key=lambda kv: -kv[1]))}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", nargs="+", choices=list(ORGS), metavar="LABEL")
    args = ap.parse_args()
    for label in args.only or list(ORGS):
        fetch_one(label)


if __name__ == "__main__":
    main()
