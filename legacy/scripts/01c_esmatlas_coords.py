"""Download ABSOLUTE Biohub ESM Atlas coordinates for the K. pneumoniae HS11286 proteome.

The latest Biohub ESM Atlas (ESMC + ESMFold2, https://biohub.ai/esm/protein/atlas) places
every protein at a fixed (x, y) in its global UMAP of the protein universe. The atlas web app
serves those exact coordinates from a public, unauthenticated REST API; this script pulls them
for all 5,728 HS11286 proteins so we get the atlas's OWN absolute coordinates, not a local
re-projection.

How it works (reverse-engineered from the live atlas + its JS bundle, all verified):
  - The atlas keys every protein by `protein_hash = md5(sequence)` (verified exact).
  - POST {BASE}/umap/coords/by-hash  body {"protein_hashes": [...<=10...]}
       -> {"by_hash": {"<hash>": [x, y], ...}}  in absolute world space (~0-1024;
          cf. {BASE}/umap/manifest `transform`: space_size=1024, center≈(573, 514)).
    Max 10 hashes per request (HTTP 413 above that), so we batch by 10.

Outputs (keyed by UniProt accession, per CLAUDE.md), under the organism folder:
  output/results/<organism>/<prefix>_esmatlas_coords.csv
  output/plots/<organism>/<prefix>_esmatlas_projection.png

Organism selected with --organism (kpneumoniae default, or ecoli).
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[1]
# organism -> (proteome file stem, output prefix, display name)
ORGANISMS = {
    "kpneumoniae": ("UP000007841_HS11286", "kp", "K. pneumoniae HS11286"),
    "ecoli": ("UP000000625_EcoliK12", "ec", "E. coli K-12 MG1655"),
}

BASE = "https://biohub.ai/esm/protein/api/v1alpha1"
COORDS_BY_HASH = f"{BASE}/umap/coords/by-hash"
BATCH = 10  # API hard limit: max 10 hashes per request
SOURCE_TAG = "biohub_esm_atlas_absolute"


def read_proteome(tsv_path: Path) -> list[tuple[str, str]]:
    """Return [(uniprot_accession, sequence), ...] from the UniProt proteome TSV."""
    rows: list[tuple[str, str]] = []
    with open(tsv_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            acc, seq = row.get("Entry"), row.get("Sequence")
            if acc and seq:
                rows.append((acc, seq))
    return rows


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=1, min=1, max=20))
def fetch_coords(hashes: list[str]) -> dict[str, list[float]]:
    # No keep-alive: occasional hung sockets on this API; a fresh connection + short
    # timeout lets a stuck request fail fast and be retried rather than blocking.
    resp = requests.post(
        COORDS_BY_HASH,
        json={"protein_hashes": hashes},
        headers={"Content-Type": "application/json", "Connection": "close"},
        timeout=(10, 30),
    )
    resp.raise_for_status()
    return resp.json().get("by_hash", {})


def chunked(seq: list, n: int):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ORGANISMS), default="kpneumoniae")
    ap.add_argument("--workers", type=int, default=8, help="concurrent requests")
    ap.add_argument(
        "--max-passes",
        type=int,
        default=5,
        help="re-query passes for on-demand-folded proteins",
    )
    args = ap.parse_args()

    stem, prefix, display = ORGANISMS[args.organism]
    tsv_path = REPO_ROOT / "data" / "raw" / args.organism / "proteome" / f"{stem}.tsv"
    results_dir = REPO_ROOT / "output" / "results" / args.organism
    plots_dir = REPO_ROOT / "output" / "plots" / args.organism
    csv_path = results_dir / f"{prefix}_esmatlas_coords.csv"
    plot_path = plots_dir / f"{prefix}_esmatlas_projection.png"

    proteome = read_proteome(tsv_path)
    n = len(proteome)
    print(f"Loaded {n} proteins from {tsv_path.relative_to(REPO_ROOT)}")

    # protein_hash == md5(sequence) (verified against the atlas)
    records = [
        {
            "uniprot_accession": acc,
            "protein_hash": hashlib.md5(seq.encode()).hexdigest(),
        }
        for acc, seq in proteome
    ]
    hashes = [r["protein_hash"] for r in records]

    def fetch_all(to_fetch: list[str], workers: int) -> dict[str, list[float]]:
        out: dict[str, list[float]] = {}
        batches = list(chunked(to_fetch, BATCH))
        done = 0
        with ThreadPoolExecutor(max_workers=workers) as pool:
            for result in pool.map(fetch_coords, batches):
                out.update(result)
                done += 1
                if done % 50 == 0 or done == len(batches):
                    print(
                        f"    {done}/{len(batches)} batches, {len(out)} coords",
                        flush=True,
                    )
        return out

    # The atlas folds some proteins ON DEMAND: a hash whose structure isn't cached yet is
    # silently omitted from the first response, then returns a coord once folded. So we make
    # an initial pass, then re-query the still-missing hashes until the set stops shrinking.
    coords: dict[str, list[float]] = {}
    pending = hashes
    for attempt in range(1, args.max_passes + 1):
        print(
            f"Pass {attempt}: fetching {len(pending)} hashes "
            f"({args.workers} workers, batch <= {BATCH}) ..."
        )
        coords.update(fetch_all(pending, args.workers))
        pending = [h for h in hashes if h not in coords]
        print(f"  -> {len(coords)} resolved, {len(pending)} still missing")
        if not pending or attempt == args.max_passes:
            break

    for r in records:
        xy = coords.get(r["protein_hash"])
        r["matched"] = xy is not None
        r["atlas_x"] = float(xy[0]) if xy else float("nan")
        r["atlas_y"] = float(xy[1]) if xy else float("nan")
        r["source"] = SOURCE_TAG
        r["atlas_api"] = COORDS_BY_HASH

    df = pd.DataFrame.from_records(
        records,
        columns=[
            "uniprot_accession",
            "protein_hash",
            "atlas_x",
            "atlas_y",
            "matched",
            "source",
            "atlas_api",
        ],
    )

    # ---- verification ----
    n_matched = int(df["matched"].sum())
    unmatched = df.loc[~df["matched"], "uniprot_accession"].tolist()
    assert len(df) == n, f"expected {n} rows, got {len(df)}"
    assert df["uniprot_accession"].is_unique, "duplicate accessions"
    assert not df.loc[df["matched"], ["atlas_x", "atlas_y"]].isna().any().any(), (
        "NaN in matched"
    )
    print(f"Matched {n_matched}/{n} ({100 * n_matched / n:.1f}%).")
    if unmatched:
        print(
            f"  UNMATCHED ({len(unmatched)}): {', '.join(unmatched[:20])}"
            + (" ..." if len(unmatched) > 20 else "")
        )
    m = df[df["matched"]]
    print(
        f"  x range [{m.atlas_x.min():.1f}, {m.atlas_x.max():.1f}]  "
        f"y range [{m.atlas_y.min():.1f}, {m.atlas_y.max():.1f}]"
    )

    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(11, 11))
    ax.scatter(m.atlas_x, m.atlas_y, s=6, c="#1f77b4", linewidths=0, alpha=0.6)
    ax.set_aspect("equal")
    ax.set_xlabel("ESM Atlas UMAP-1 (absolute world space)")
    ax.set_ylabel("ESM Atlas UMAP-2 (absolute world space)")
    ax.set_title(
        f"{display} proteome in the Biohub ESM Atlas (absolute coordinates)\n"
        f"n={n_matched}/{n} proteins located via umap/coords/by-hash (hash = md5 of sequence)",
        fontsize=10,
    )
    ax.invert_yaxis()  # atlas screen-space convention (y grows downward)
    fig.savefig(plot_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(
        f"Wrote:\n  {csv_path.relative_to(REPO_ROOT)}\n  {plot_path.relative_to(REPO_ROOT)}"
    )


if __name__ == "__main__":
    main()
