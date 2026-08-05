"""Stage 09a — UniProt-curated subcellular localization (docs §5.1a).

Fetches the curated subcellular-location annotation for the focal proteomes and normalises it onto
the canonical class vocabulary in `src/localization.py`. This is the **highest-precedence** track in
the axis: a curated compartment beats any predictor, and an experimentally-supported one beats a
merely curated one — so the fetch also pulls the ECO evidence codes and splits the two apart.

Per protein it records:
  * `uniprot_localization` — canonical class, or `unknown` when UniProt says nothing;
  * `uniprot_evidence` — `experimental` (an ECO experimental code is cited) / `curated` / `none`;
  * `has_signal_peptide`, `n_transmembrane` — the sequence features UniProt annotates;
  * `has_lipid_anchor`, `lipid_cys_pos` — the LIPID feature. A diacylglycerol-cysteine anchor is
    the signature of a Sec/SPII lipoprotein, and its position feeds the Lol "+2 rule" in 09g,
    which is what separates outer- from inner-membrane lipoproteins.

Note this script no longer writes the merged axis table — that is 09g's job. It emits only its own
track so every source stays independently inspectable.

Output: output/results/<organism>/<prefix>_loc_uniprot.csv
Run with the `gradi` conda env. Network fetch (UniProt REST stream), same endpoint as 00a.
"""
from __future__ import annotations

import argparse
import csv
import io
import re
import sys
from pathlib import Path
from urllib.parse import quote

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import localization as LOC  # noqa: E402

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
FIELDS = "accession,cc_subcellular_location,ft_signal,ft_transmem,ft_lipid"

# organism -> UniProt proteome id (the same ids 00a fetches).
PROTEOME_ID = {"kpneumoniae": "UP000007841", "ecoli": "UP000000625"}

_LIPID_POS = re.compile(r"LIPID\s+(\d+)")


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def download(url: str) -> str:
    r = requests.get(url, timeout=300)
    r.raise_for_status()
    return r.text


def parse_lipid(field: str) -> tuple[bool, int | None]:
    """(is a Sec/SPII-style lipoprotein anchor, 1-based position of the lipidated Cys)."""
    if not field or not field.strip():
        return False, None
    is_lipoprotein = "diacylglycerol cysteine" in field.lower()
    m = _LIPID_POS.search(field)
    return is_lipoprotein, (int(m.group(1)) if m else None)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(LOC.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()

    org = args.organism
    _, prefix = LOC.ORGANISMS[org]
    pid = PROTEOME_ID[org]

    query = quote("proteome:" + pid, safe=":")
    url = f"{STREAM}?compressed=false&format=tsv&query={query}&fields={FIELDS}"
    print(f"[{org}] fetching UniProt localization for {pid} ...", flush=True)
    rows = list(csv.DictReader(io.StringIO(download(url)), delimiter="\t"))

    def col(row: dict, prefix_: str) -> str:
        for k, v in row.items():
            if k.lower().replace(" ", "").startswith(prefix_):
                return v or ""
        return ""

    records = []
    for r in rows:
        acc = (col(r, "entry") or list(r.values())[0]).strip()
        subcell = col(r, "subcellular")
        signal = col(r, "signal")
        transmem = col(r, "transmembrane")
        lipid = col(r, "lipidation")

        has_signal = bool(signal.strip())
        n_tm = transmem.upper().count("TRANSMEM")
        is_lipo, cys_pos = parse_lipid(lipid)

        records.append({
            "uniprot_accession": acc,
            "uniprot_localization": LOC.classify_uniprot(subcell, has_signal, n_tm),
            "uniprot_evidence": LOC.uniprot_evidence(subcell),
            "has_signal_peptide": has_signal,
            "n_transmembrane": n_tm,
            "has_lipid_anchor": is_lipo,
            "lipid_cys_pos": cys_pos,
            "subcellular_raw": subcell.replace("\t", " "),
        })

    df = pd.DataFrame(records)
    # Keep the canonical row set and file order (proteome TSV), as every other track does.
    accs = LOC.load_accessions(org)
    df = df.set_index("uniprot_accession").reindex(accs).reset_index()

    out = LOC.results_dir(org) / f"{prefix}_loc_uniprot.csv"
    df.to_csv(out, index=False)
    print(f"[{org}] wrote {out.relative_to(LOC.REPO_ROOT)} ({len(df)} proteins)", flush=True)
    print(f"[{org}] classes: {df['uniprot_localization'].value_counts().to_dict()}", flush=True)
    print(f"[{org}] evidence: {df['uniprot_evidence'].value_counts().to_dict()}", flush=True)
    print(f"[{org}] lipoprotein anchors: {int(df['has_lipid_anchor'].sum())}", flush=True)


if __name__ == "__main__":
    main()
