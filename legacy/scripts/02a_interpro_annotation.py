"""Annotate a reference proteome with InterPro families, homologous superfamilies, and
domains (step 1.1b of docs/01_task_agnostic.md).

InterPro aggregates ~13 member databases (Pfam, SUPERFAMILY, Gene3D, SMART, PROSITE,
PANTHER, NCBIfam, ...) and folds overlapping signatures into a single InterPro entry
that carries a *type* (Family, Homologous_superfamily, Domain, Repeat, Active_site,
Binding_site, Conserved_site, PTM). We annotate against InterPro *entries only* (already
de-duplicated + typed) and bucket by entry type; Pfam is kept as one supplementary column.
No residue positions, no local InterProScan.

Two cheap, cached downloads + a local join:
  1. UniProt REST stream -> per-protein InterPro/Pfam cross-reference IDs (one call).
  2. InterPro entry.list -> IPR###### -> (type, name) dictionary (one flat file).
       https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list

Organism selected with --organism (kpneumoniae default, or ecoli). Raw downloads cached under
data/raw/<organism>/interpro/ (pass --refresh to force). Output, keyed by UniProt accession:
  data/processed/<organism>/families/<prefix>_interpro_annotation.csv
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
XREF_FIELDS = "accession,xref_interpro,xref_pfam"
ENTRY_LIST_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list"
)
RELEASE_NOTES_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/release_notes.txt"
)

# InterPro entry.list ENTRY_TYPE values -> output bucket
TYPE_FAMILY = "Family"
TYPE_SUPERFAMILY = "Homologous_superfamily"
TYPE_DOMAIN = "Domain"

REPO_ROOT = Path(__file__).resolve().parents[1]
ORGANISMS = {
    "kpneumoniae": {
        "proteome": "UP000007841",
        "prefix": "kp",
        "name": "K. pneumoniae HS11286",
    },
    "ecoli": {"proteome": "UP000000625", "prefix": "ec", "name": "E. coli K-12 MG1655"},
}
TOP_N = 15  # entries shown per panel in the barplot


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


def parse_entry_list(text: str) -> dict[str, tuple[str, str]]:
    """ENTRY_AC\\tENTRY_TYPE\\tENTRY_NAME -> {IPR######: (type, name)}."""
    entries: dict[str, tuple[str, str]] = {}
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    for row in reader:
        acc = (row.get("ENTRY_AC") or "").strip()
        if acc:
            entries[acc] = (
                (row.get("ENTRY_TYPE") or "").strip(),
                (row.get("ENTRY_NAME") or "").strip(),
            )
    return entries


def split_ids(value: str | None) -> list[str]:
    """Split a UniProt semicolon-delimited xref cell into clean IDs."""
    if not value:
        return []
    return [tok.strip() for tok in value.split(";") if tok.strip()]


def interpro_release(text: str) -> str:
    """Best-effort InterPro release version from release_notes.txt (else 'unknown')."""
    for line in text.splitlines():
        line = line.strip()
        if re.match(r"^Release\s+\d", line):
            return line
    return "unknown"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    ap.add_argument(
        "--refresh",
        action="store_true",
        help="re-download UniProt xrefs and entry.list even if cached",
    )
    args = ap.parse_args()

    spec = ORGANISMS[args.organism]
    proteome_id, prefix, org_name = spec["proteome"], spec["prefix"], spec["name"]
    raw_dir = REPO_ROOT / "data" / "raw" / args.organism / "interpro"
    xrefs_path = raw_dir / f"{proteome_id}_uniprot_xrefs.tsv"
    entry_list_path = raw_dir / "entry.list"
    processed_dir = REPO_ROOT / "data" / "processed" / args.organism / "families"
    csv_path = processed_dir / f"{prefix}_interpro_annotation.csv"
    plots_dir = REPO_ROOT / "output" / "plots" / args.organism
    plot_path = plots_dir / f"{prefix}_interpro_entries.png"
    log_path = REPO_ROOT / "docs" / f"{prefix}_interpro_annotation_log.md"

    raw_dir.mkdir(parents=True, exist_ok=True)

    xrefs_tsv = fetch_cached(
        f"{STREAM}?compressed=false&format=tsv&query=proteome:{proteome_id}&fields={XREF_FIELDS}",
        xrefs_path,
        f"UniProt InterPro/Pfam cross-references ({org_name})",
        args.refresh,
    )
    entry_text = fetch_cached(
        ENTRY_LIST_URL, entry_list_path, "InterPro entry.list", args.refresh
    )
    entries = parse_entry_list(entry_text)
    print(f"Loaded {len(entries)} InterPro entry definitions.")

    try:
        release = interpro_release(download(RELEASE_NOTES_URL))
    except Exception as exc:  # noqa: BLE001 - release version is best-effort metadata
        release = "unknown"
        print(f"  (could not read InterPro release notes: {exc})")
    print(f"InterPro release: {release}")

    # ---- join + bucket by entry type ----
    records: list[dict] = []
    unknown_iprs: set[str] = set()
    # per-bucket {ipr_id: n proteins} and {ipr_id: name}, for the barplot
    counts = {"family": Counter(), "superfamily": Counter(), "domain": Counter()}
    labels: dict[str, dict[str, str]] = {"family": {}, "superfamily": {}, "domain": {}}
    reader = csv.DictReader(xrefs_tsv.splitlines(), delimiter="\t")
    for row in reader:
        acc = (row.get("Entry") or "").strip()
        if not acc:
            continue
        ipr_ids = split_ids(row.get("InterPro"))
        pfam_ids = split_ids(row.get("Pfam"))

        fam_ids, fam_names = [], []
        sup_ids, sup_names = [], []
        dom_ids, dom_names = [], []
        for ipr in ipr_ids:
            etype, ename = entries.get(ipr, (None, None))
            if etype is None:
                unknown_iprs.add(ipr)
                continue
            if etype == TYPE_FAMILY:
                fam_ids.append(ipr)
                fam_names.append(ename)
                counts["family"][ipr] += 1
                labels["family"][ipr] = ename
            elif etype == TYPE_SUPERFAMILY:
                sup_ids.append(ipr)
                sup_names.append(ename)
                counts["superfamily"][ipr] += 1
                labels["superfamily"][ipr] = ename
            elif etype == TYPE_DOMAIN:
                dom_ids.append(ipr)
                dom_names.append(ename)
                counts["domain"][ipr] += 1
                labels["domain"][ipr] = ename
            # other types (Repeat, *_site, PTM) intentionally not bucketed here

        records.append(
            {
                "uniprot_accession": acc,
                "interpro_family_ids": fam_ids,
                "interpro_family_names": fam_names,
                "interpro_superfamily_ids": sup_ids,
                "interpro_superfamily_names": sup_names,
                "interpro_domain_ids": dom_ids,
                "interpro_domain_names": dom_names,
                "pfam_ids": pfam_ids,
                "n_interpro_entries": len(ipr_ids),
            }
        )

    df = pd.DataFrame.from_records(
        records,
        columns=[
            "uniprot_accession",
            "interpro_family_ids",
            "interpro_family_names",
            "interpro_superfamily_ids",
            "interpro_superfamily_names",
            "interpro_domain_ids",
            "interpro_domain_names",
            "pfam_ids",
            "n_interpro_entries",
        ],
    )

    # ---- verification ----
    n = len(df)
    assert df["uniprot_accession"].is_unique, "duplicate accessions"
    n_any = int((df["n_interpro_entries"] > 0).sum())
    no_hit = df.loc[df["n_interpro_entries"] == 0, "uniprot_accession"].tolist()
    n_fam = int((df["interpro_family_ids"].str.len() > 0).sum())
    n_sup = int((df["interpro_superfamily_ids"].str.len() > 0).sum())
    n_dom = int((df["interpro_domain_ids"].str.len() > 0).sum())
    n_pfam = int((df["pfam_ids"].str.len() > 0).sum())

    print(
        f"\n{org_name}: annotated {n} proteins; {n_any} ({100 * n_any / n:.1f}%) have >=1 InterPro entry."
    )
    print(f"  with family:        {n_fam} ({100 * n_fam / n:.1f}%)")
    print(f"  with superfamily:   {n_sup} ({100 * n_sup / n:.1f}%)")
    print(f"  with domain:        {n_dom} ({100 * n_dom / n:.1f}%)")
    print(f"  with Pfam:          {n_pfam} ({100 * n_pfam / n:.1f}%)")
    if no_hit:
        print(
            f"  NO InterPro entry ({len(no_hit)}): {', '.join(no_hit[:20])}"
            + (" ..." if len(no_hit) > 20 else "")
        )
    if unknown_iprs:
        print(
            f"  NOTE: {len(unknown_iprs)} InterPro IDs referenced by UniProt are absent from "
            f"entry.list (UniProt xrefs lag the InterPro release; these are dropped from "
            f"buckets): {', '.join(sorted(unknown_iprs)[:20])}"
            + (" ..." if len(unknown_iprs) > 20 else "")
        )

    # ---- write output (CSV; list columns joined with ';') ----
    processed_dir.mkdir(parents=True, exist_ok=True)
    list_cols = [c for c in df.columns if c.endswith(("_ids", "_names"))]
    for c in list_cols:
        df[c] = df[c].apply(lambda xs: ";".join(xs))
    df.to_csv(csv_path, index=False)
    print(f"\nWrote:\n  {csv_path.relative_to(REPO_ROOT)}")

    # ---- barplot: top-N entries per type (Family / Homologous superfamily / Domain) ----
    plots_dir.mkdir(parents=True, exist_ok=True)
    panels = [
        ("family", "Family", "#2b6cb0"),
        ("superfamily", "Homologous superfamily", "#2f855a"),
        ("domain", "Domain", "#b7791f"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(20, 8))
    for ax, (bucket, title, color) in zip(axes, panels):
        top = counts[bucket].most_common(TOP_N)
        ylabels, vals = [], []
        for ipr, cnt in top:
            name = labels[bucket].get(ipr) or ""
            name = name if len(name) <= 38 else name[:35] + "..."
            ylabels.append(f"{ipr}  {name}")
            vals.append(cnt)
        ylabels.reverse()
        vals.reverse()
        bars = ax.barh(range(len(vals)), vals, color=color)
        ax.set_yticks(range(len(vals)))
        ax.set_yticklabels(ylabels, fontsize=7)
        ax.set_xlabel("number of proteins")
        ax.set_title(f"{title}\n(top {len(vals)})", fontsize=11)
        ax.bar_label(bars, padding=2, fontsize=7)
        ax.margins(x=0.10)
    fig.suptitle(
        f"Top InterPro entries by type in {org_name}  ({release})", fontsize=13
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(plot_path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {plot_path.relative_to(REPO_ROOT)}")

    # ---- coverage log (Git-tracked) ----
    log = (
        "# InterPro annotation — coverage log\n\n"
        f"- Proteome: `{proteome_id}` ({org_name}), {n} proteins\n"
        f"- InterPro {release}\n"
        f"- Source: UniProt xref_interpro/xref_pfam (cached `{xrefs_path.relative_to(REPO_ROOT)}`) "
        f"joined to `entry.list`\n"
        f"- Scope: InterPro entries bucketed by type (Family / Homologous_superfamily / Domain) "
        f"+ Pfam supplementary\n\n"
        "## Coverage\n\n"
        f"- >=1 InterPro entry: {n_any}/{n} ({100 * n_any / n:.1f}%)\n"
        f"- with Family: {n_fam} ({100 * n_fam / n:.1f}%)\n"
        f"- with Homologous_superfamily: {n_sup} ({100 * n_sup / n:.1f}%)\n"
        f"- with Domain: {n_dom} ({100 * n_dom / n:.1f}%)\n"
        f"- with Pfam: {n_pfam} ({100 * n_pfam / n:.1f}%)\n"
        f"- no InterPro entry: {len(no_hit)}\n"
        f"- InterPro IDs unresolved against entry.list: {len(unknown_iprs)}\n"
    )
    log_path.write_text(log)
    print(f"  {log_path.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
