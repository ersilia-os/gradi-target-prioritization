"""Annotate a reference proteome with PANTHER families and subfamilies (step 1.1a).

PANTHER is one of InterPro's member databases, exposed by UniProt as a cross-reference, so
the mechanism mirrors the InterPro step (02a). Two cheap, cached downloads + a local join:
  1. UniProt REST stream -> per-protein PANTHER IDs (one call). Each protein gets a family id
       (PTHR#####) and, usually, a subfamily id (PTHR#####:SF##): "PTHR11066;PTHR11066:SF34;"
  2. PANTHER HMM_classifications -> {PTHR id: name} (one flat file; PANTHER's analogue of
       InterPro's entry.list). http://data.pantherdb.org/ftp/hmm_classifications/current_release/

Organism selected with --organism (kpneumoniae default, or ecoli). Raw downloads cached under
data/raw/<organism>/panther/ (pass --refresh to force). Output, keyed by UniProt accession:
  data/processed/<organism>/families/<prefix>_panther_annotation.csv
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter
from pathlib import Path

import matplotlib
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
XREF_FIELDS = "accession,xref_panther"
PANTHER_DIR = "http://data.pantherdb.org/ftp/hmm_classifications/current_release/"
CLASSIFICATIONS_RE = re.compile(r"PANTHER[\d.]+_HMM_classifications")

REPO_ROOT = Path(__file__).resolve().parents[1]
ORGANISMS = {
    "kpneumoniae": {
        "proteome": "UP000007841",
        "prefix": "kp",
        "name": "K. pneumoniae HS11286",
    },
    "ecoli": {"proteome": "UP000000625", "prefix": "ec", "name": "E. coli K-12 MG1655"},
}
TOP_N = 25  # families shown in the barplot


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def download(url: str) -> str:
    resp = requests.get(url, timeout=300)
    resp.raise_for_status()
    return resp.text


def fetch_cached(url: str, path: Path, label: str, refresh: bool) -> str:
    """Download `url` to `path`, or reuse the cached copy unless --refresh."""
    if path.exists() and not refresh:
        print(f"Using cached {label}: {path.relative_to(REPO_ROOT)}")
        return path.read_text()
    print(f"Downloading {label} ...")
    text = download(url)
    path.write_text(text)
    print(f"  wrote {path.relative_to(REPO_ROOT)}")
    return text


def resolve_classifications_url() -> tuple[str, str]:
    """Discover the current PANTHER HMM_classifications filename (robust to version bumps)."""
    listing = download(PANTHER_DIR)
    m = CLASSIFICATIONS_RE.search(listing)
    if not m:
        raise RuntimeError(f"no PANTHER*_HMM_classifications found in {PANTHER_DIR}")
    fname = m.group(0)
    version = re.search(r"PANTHER([\d.]+)_", fname).group(1)
    return PANTHER_DIR + fname, version


def parse_classifications(text: str) -> dict[str, str]:
    """tab-delimited; col0 = PANTHER id (PTHR##### or PTHR#####:SF##), col1 = name."""
    names: dict[str, str] = {}
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) >= 2 and cols[0].startswith("PTHR"):
            names[cols[0].strip()] = cols[1].strip()
    return names


def split_ids(value: str | None) -> list[str]:
    """Split a UniProt semicolon-delimited xref cell into clean IDs."""
    if not value:
        return []
    return [tok.strip() for tok in value.split(";") if tok.strip()]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    ap.add_argument(
        "--refresh",
        action="store_true",
        help="re-download UniProt xrefs and PANTHER classifications even if cached",
    )
    args = ap.parse_args()

    spec = ORGANISMS[args.organism]
    proteome_id, prefix, org_name = spec["proteome"], spec["prefix"], spec["name"]
    raw_dir = REPO_ROOT / "data" / "raw" / args.organism / "panther"
    xrefs_path = raw_dir / f"{proteome_id}_uniprot_xrefs.tsv"
    processed_dir = REPO_ROOT / "data" / "processed" / args.organism / "families"
    csv_path = processed_dir / f"{prefix}_panther_annotation.csv"
    plots_dir = REPO_ROOT / "output" / "plots" / args.organism
    plot_path = plots_dir / f"{prefix}_panther_families.png"
    log_path = REPO_ROOT / "docs" / f"{prefix}_panther_annotation_log.md"

    raw_dir.mkdir(parents=True, exist_ok=True)

    xrefs_tsv = fetch_cached(
        f"{STREAM}?compressed=false&format=tsv&query=proteome:{proteome_id}&fields={XREF_FIELDS}",
        xrefs_path,
        f"UniProt PANTHER cross-references ({org_name})",
        args.refresh,
    )

    class_url, version = resolve_classifications_url()
    class_path = raw_dir / Path(class_url).name
    class_text = fetch_cached(
        class_url, class_path, f"PANTHER {version} HMM_classifications", args.refresh
    )
    names = parse_classifications(class_text)
    print(f"Loaded {len(names)} PANTHER family/subfamily names (PANTHER {version}).")

    # ---- join + bucket into family (PTHR#####) vs subfamily (PTHR#####:SF##) ----
    records: list[dict] = []
    unknown_ids: set[str] = set()
    family_counts: Counter[str] = Counter()  # family_id -> n proteins assigned
    family_label: dict[str, str] = {}  # family_id -> name (for plot labels)
    reader = csv.DictReader(xrefs_tsv.splitlines(), delimiter="\t")
    for row in reader:
        acc = (row.get("Entry") or "").strip()
        if not acc:
            continue
        ptids = split_ids(row.get("PANTHER"))

        fam_ids, fam_names = [], []
        sub_ids, sub_names = [], []
        for pid in ptids:
            name = names.get(pid)
            if name is None:
                unknown_ids.add(pid)
                name = ""  # keep the id, name unresolved
            if ":SF" in pid:
                sub_ids.append(pid)
                sub_names.append(name)
            else:
                fam_ids.append(pid)
                fam_names.append(name)
                family_counts[pid] += 1
                family_label[pid] = name

        records.append(
            {
                "uniprot_accession": acc,
                "panther_family_ids": fam_ids,
                "panther_family_names": fam_names,
                "panther_subfamily_ids": sub_ids,
                "panther_subfamily_names": sub_names,
                "n_panther_entries": len(ptids),
            }
        )

    df = pd.DataFrame.from_records(
        records,
        columns=[
            "uniprot_accession",
            "panther_family_ids",
            "panther_family_names",
            "panther_subfamily_ids",
            "panther_subfamily_names",
            "n_panther_entries",
        ],
    )

    # ---- verification ----
    n = len(df)
    assert df["uniprot_accession"].is_unique, "duplicate accessions"
    n_any = int((df["n_panther_entries"] > 0).sum())
    no_hit = df.loc[df["n_panther_entries"] == 0, "uniprot_accession"].tolist()
    n_fam = int((df["panther_family_ids"].str.len() > 0).sum())
    n_sub = int((df["panther_subfamily_ids"].str.len() > 0).sum())

    print(
        f"\n{org_name}: annotated {n} proteins; {n_any} ({100 * n_any / n:.1f}%) have >=1 PANTHER entry."
    )
    print(f"  with family:     {n_fam} ({100 * n_fam / n:.1f}%)")
    print(f"  with subfamily:  {n_sub} ({100 * n_sub / n:.1f}%)")
    if no_hit:
        print(
            f"  NO PANTHER entry ({len(no_hit)}): {', '.join(no_hit[:20])}"
            + (" ..." if len(no_hit) > 20 else "")
        )
    if unknown_ids:
        print(
            f"  NOTE: {len(unknown_ids)} PANTHER IDs referenced by UniProt are absent from the "
            f"classifications file (id kept, name blank): {', '.join(sorted(unknown_ids)[:20])}"
            + (" ..." if len(unknown_ids) > 20 else "")
        )

    # ---- write output (CSV; list columns joined with ';') ----
    processed_dir.mkdir(parents=True, exist_ok=True)
    list_cols = [c for c in df.columns if c.endswith(("_ids", "_names"))]
    for c in list_cols:
        df[c] = df[c].apply(lambda xs: ";".join(xs))
    df.to_csv(csv_path, index=False)
    print(f"\nWrote:\n  {csv_path.relative_to(REPO_ROOT)}")

    # ---- barplot: top-N most populated PANTHER families ----
    plots_dir.mkdir(parents=True, exist_ok=True)
    top = family_counts.most_common(TOP_N)
    labels, counts = [], []
    for fam_id, cnt in top:
        name = family_label.get(fam_id) or ""
        name = name if len(name) <= 42 else name[:39] + "..."
        labels.append(f"{fam_id}  {name}")
        counts.append(cnt)
    labels.reverse()
    counts.reverse()
    fig, ax = plt.subplots(figsize=(11, 9))
    bars = ax.barh(range(len(counts)), counts, color="#2b6cb0")
    ax.set_yticks(range(len(counts)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("number of proteins")
    ax.set_title(
        f"Top {len(counts)} PANTHER families in {org_name}\n"
        f"({n_fam}/{n} proteins assigned a PANTHER family; PANTHER {version})",
        fontsize=11,
    )
    ax.bar_label(bars, padding=2, fontsize=8)
    ax.margins(x=0.08)
    fig.tight_layout()
    fig.savefig(plot_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {plot_path.relative_to(REPO_ROOT)}")

    # ---- coverage log (Git-tracked) ----
    log = (
        "# PANTHER annotation — coverage log\n\n"
        f"- Proteome: `{proteome_id}` ({org_name}), {n} proteins\n"
        f"- PANTHER {version} (HMM_classifications)\n"
        f"- Source: UniProt xref_panther (cached `{xrefs_path.relative_to(REPO_ROOT)}`) "
        f"joined to PANTHER HMM_classifications\n"
        f"- Scope: PANTHER IDs bucketed into family (PTHR#####) vs subfamily (PTHR#####:SF##)\n\n"
        "## Coverage\n\n"
        f"- >=1 PANTHER entry: {n_any}/{n} ({100 * n_any / n:.1f}%)\n"
        f"- with family: {n_fam} ({100 * n_fam / n:.1f}%)\n"
        f"- with subfamily: {n_sub} ({100 * n_sub / n:.1f}%)\n"
        f"- no PANTHER entry: {len(no_hit)}\n"
        f"- PANTHER IDs unresolved against classifications: {len(unknown_ids)}\n"
    )
    log_path.write_text(log)
    print(f"  {log_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
