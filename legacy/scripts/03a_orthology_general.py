"""Build ortholog "synonym" tables mapping other proteomes onto a focal anchor proteome.

Generalized over the project's focal organisms (`--organism kpneumoniae|ecoli`): for the chosen
anchor it maps proteins from other proteomes onto it, keyed by UniProt accession, in four tiers:
  (a) close relatives, (b) other Gram-negatives, (c) other bacteria  -> OrthoFinder (DIAMOND)
  (d) human                                                          -> DIAMOND similarity

The close-relative tier is anchor-specific: other Klebsiella for K. pneumoniae, Escherichia +
Shigella for E. coli. The two focal organisms cross-include each other in the Gram-negative tier
(the anchor itself is always excluded). Why two methods: OrthoFinder (community-standard
phylogenetic orthology) gives proper orthologs within bacteria; across the bacteria->human
distance orthology is unreliable, so the human tier uses sequence-similarity search
(subtractive-genomics standard) to flag a human homolog (selectivity / off-target safety).

Output: one long tidy table per anchor at
  data/processed/other/orthology/<prefix>_orthologs_long.tsv   (kp_… or ec_…)
columns: anchor_uniprot, tier, species, target_uniprot, method, orthogroup, pident, coverage,
evalue, bitscore, is_human_homolog.

Run with the `gradi-ortho` conda env (orthofinder + diamond on PATH):
    conda activate gradi-ortho
    python scripts/03a_orthology_general.py --organism kpneumoniae --tiers a,b,c,d
    python scripts/03a_orthology_general.py --organism ecoli       --tiers a,b,c,d
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import requests
from tenacity import retry, stop_after_attempt, wait_exponential

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = (
    REPO_ROOT / "data" / "raw" / "other" / "orthology"
)  # shared bare-AC panel FASTAs
PROC_DIR = (
    REPO_ROOT / "data" / "processed" / "other" / "orthology"
)  # per-anchor results

STREAM = "https://rest.uniprot.org/uniprotkb/stream"
PROTEOMES = "https://rest.uniprot.org/proteomes/search"

# OrthoFinder/DIAMOND live in this interpreter's env bin; the orthofinder launcher uses a
# `#!/usr/bin/env python3` shebang, so that dir must be on PATH for subprocesses to resolve
# the right python (with ete4) even when the env isn't "activated".
ENV_BIN = Path(sys.executable).parent


def _anchor_proteome(organism: str, stem: str) -> Path:
    return REPO_ROOT / "data" / "raw" / organism / "proteome" / stem


# Focal anchors. close = which CLOSE tier to use; proteome_id is used to exclude the anchor
# from the shared panel tiers (so it is never compared to itself).
ANCHORS: dict[str, dict] = {
    "kpneumoniae": {
        "label": "HS11286",
        "proteome_id": "UP000007841",
        "prefix": "kp",
        "close": "klebsiella",
        "fasta": _anchor_proteome("kpneumoniae", "UP000007841_HS11286.fasta"),
        "tsv": _anchor_proteome("kpneumoniae", "UP000007841_HS11286.tsv"),
    },
    "ecoli": {
        "label": "EcoliK12",
        "proteome_id": "UP000000625",
        "prefix": "ec",
        "close": "escherichia",
        "fasta": _anchor_proteome("ecoli", "UP000000625_EcoliK12.fasta"),
        "tsv": _anchor_proteome("ecoli", "UP000000625_EcoliK12.tsv"),
    },
}

# Curated reference-proteome panel. Entries: (label, taxon_id, proteome_id_or_None).
# A pinned proteome_id is used directly; otherwise the taxon's reference proteome is resolved.
# Anchor-specific "close relatives" tier:
CLOSE: dict[str, list[tuple[str, int, str | None]]] = {
    # other Klebsiella (NTUH-K2044 / KPNIH1 / KPPR1 omitted: redundant proteomes serve no sequences)
    "klebsiella": [
        ("Kpneumoniae_MGH78578", 272620, "UP000000265"),
        ("Kvaricola", 244366, None),
        ("Koxytoca", 571, None),
        ("Kaerogenes", 548, None),
        ("Kquasipneumoniae", 1463165, None),
    ],
    # Escherichia + Shigella (Shigella is genomically E. coli)
    "escherichia": [
        ("Shigella_flexneri", 623, None),
        ("Shigella_sonnei", 624, None),
        ("Ecoli_fergusonii", 564, None),
        ("Ecoli_albertii", 208962, None),
    ],
}
# Shared tiers. gram_negative includes BOTH focal organisms; the anchor is dropped at runtime.
GRAM_NEGATIVE: list[tuple[str, int, str | None]] = [
    ("Kpneumoniae_HS11286", 1125630, "UP000007841"),
    ("Ecoli_K12_MG1655", 83333, "UP000000625"),
    ("Paeruginosa_PAO1", 208964, "UP000002438"),
    ("Abaumannii", 470, None),
    ("Senterica_Typhimurium_LT2", 99287, "UP000001014"),
    ("Enterobacter_hormaechei", 158836, None),
    ("Serratia_marcescens", 615, None),
    ("Proteus_mirabilis", 584, None),
    ("Haemophilus_influenzae", 71421, None),
    ("Vibrio_cholerae", 666, None),
    ("Neisseria_gonorrhoeae", 485, None),
]
BACTERIA: list[tuple[str, int, str | None]] = [
    ("Bsubtilis_168", 224308, "UP000001570"),
    ("Saureus_NCTC8325", 93061, "UP000008816"),
    ("Mtuberculosis_H37Rv", 83332, "UP000001584"),
    ("Caulobacter_vibrioides", 565050, "UP000001364"),
    ("Spneumoniae", 1313, None),
    ("Efaecium", 1352, None),
    ("Lmonocytogenes", 1639, None),
    ("Cdifficile", 1496, None),
    ("Hpylori", 85962, None),
]
HUMAN = ("Homo_sapiens", 9606, "UP000005640")


def _run(cmd: list[str]) -> None:
    env = os.environ.copy()
    env["PATH"] = f"{ENV_BIN}{os.pathsep}{env.get('PATH', '')}"
    subprocess.run(cmd, check=True, env=env)


def _require(tool: str) -> None:
    if not (ENV_BIN / tool).exists() and not shutil.which(tool):
        sys.exit(
            f"`{tool}` not found. Run this script with the gradi-ortho env "
            f"(its python), e.g. ~/miniconda3/envs/gradi-ortho/bin/python (see install.sh)."
        )


# --------------------------- fetching proteomes ---------------------------
@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def _get(url: str, **params) -> requests.Response:
    r = requests.get(url, params=params or None, timeout=300)
    r.raise_for_status()
    return r


def resolve_reference_proteome(taxon_id: int) -> str:
    """Return the UniProt reference (or representative) proteome UPID for a taxon."""
    r = _get(PROTEOMES, query=f"taxonomy_id:{taxon_id}", format="json", size=50)
    results = r.json().get("results", [])
    order = {
        "Reference and representative proteome": 0,
        "Reference proteome": 0,
        "Representative proteome": 1,
        "Other proteome": 2,
        "Redundant proteome": 3,
    }
    results.sort(
        key=lambda p: (
            order.get(p.get("proteomeType", ""), 5),
            -(p.get("proteinCount") or 0),
        )
    )
    if not results:
        raise RuntimeError(f"No proteome found for taxon {taxon_id}")
    return results[0]["id"]


def _fetch_fasta(upid: str) -> list[str]:
    """Fetch a proteome's FASTA and rewrite headers to the bare UniProt accession."""
    text = _get(
        STREAM, query=f"proteome:{upid}", format="fasta", compressed="false"
    ).text
    lines = []
    for line in text.splitlines():
        if line.startswith(">"):
            parts = line[1:].split("|")
            lines.append(">" + (parts[1] if len(parts) >= 2 else line[1:].split()[0]))
        else:
            lines.append(line)
    return lines


def fetch_proteome(
    label: str, taxon_id: int, proteome_id: str | None, force: bool
) -> Path | None:
    """Download a proteome FASTA (bare-AC headers). Returns None if no sequences are served."""
    out = RAW_DIR / f"{label}.fasta"
    if out.exists() and not force and out.stat().st_size > 0:
        return out
    upid = proteome_id or resolve_reference_proteome(taxon_id)
    try:
        lines = _fetch_fasta(upid)
    except requests.HTTPError:
        upid = resolve_reference_proteome(taxon_id)
        lines = _fetch_fasta(upid)
    n = sum(1 for ln in lines if ln.startswith(">"))
    if n == 0:  # e.g. "redundant" proteomes serve no sequences via the stream
        print(f"  SKIP {label} ({upid}): 0 proteins served")
        return None
    out.write_text("\n".join(lines) + "\n")
    print(f"  fetched {label} ({upid}): {n} proteins")
    return out


def anchor_accession_fasta(anchor: dict, force: bool) -> Path:
    """Copy the anchor proteome FASTA with bare-accession headers into the shared panel dir."""
    out = RAW_DIR / f"{anchor['label']}.fasta"
    if out.exists() and not force:
        return out
    lines = []
    for line in anchor["fasta"].read_text().splitlines():
        if line.startswith(">"):
            lines.append(">" + line[1:].split("|")[1])
        else:
            lines.append(line)
    out.write_text("\n".join(lines) + "\n")
    return out


# --------------------------- OrthoFinder (bacterial tiers) ---------------------------
def run_orthofinder(
    anchor: dict,
    species: list[tuple[str, int, str | None]],
    of_input: Path,
    threads: int,
    force: bool,
) -> Path:
    """Stage anchor + bacterial proteomes, run OrthoFinder, return its Results dir."""
    if of_input.exists() and force:
        shutil.rmtree(of_input)
    of_input.mkdir(parents=True, exist_ok=True)
    shutil.copy(
        anchor_accession_fasta(anchor, force), of_input / f"{anchor['label']}.fasta"
    )
    for label, taxon, pid in species:
        fasta = fetch_proteome(label, taxon, pid, force)
        if fasta is not None:
            shutil.copy(fasta, of_input / f"{label}.fasta")

    existing = sorted(of_input.glob("OrthoFinder/Results_*"))
    if existing and not force:
        print(f"  reusing OrthoFinder results {existing[-1].name}")
        return existing[-1]
    _require("orthofinder")
    # Full run (NOT -og): the default produces the Orthologues/ pairwise ortholog files.
    print("  running OrthoFinder (DIAMOND) ...")
    _run(["orthofinder", "-f", str(of_input), "-t", str(threads)])
    return sorted(of_input.glob("OrthoFinder/Results_*"))[-1]


def parse_orthofinder(
    results_dir: Path, anchor_label: str, tier_of_species: dict[str, str]
) -> pd.DataFrame:
    """Parse Orthologues/ pairwise files (anchor vs each species) into long rows."""
    orth_dir = results_dir / "Orthologues" / f"Orthologues_{anchor_label}"
    rows = []
    for f in sorted(orth_dir.glob(f"{anchor_label}__v__*.tsv")):
        species = f.stem.split("__v__")[1]
        df = pd.read_csv(f, sep="\t")  # cols: Orthogroup, <anchor>, <species>
        for _, r in df.iterrows():
            anchors = [a.strip() for a in str(r[anchor_label]).split(",")]
            targets = [t.strip() for t in str(r[species]).split(",")]
            for a in anchors:
                for t in targets:
                    rows.append(
                        {
                            "anchor_uniprot": a,
                            "tier": tier_of_species[species],
                            "species": species,
                            "target_uniprot": t,
                            "method": "orthofinder",
                            "orthogroup": r["Orthogroup"],
                        }
                    )
    return pd.DataFrame(rows)


# --------------------------- DIAMOND (human tier) ---------------------------
def run_diamond_human(
    human_fasta: Path,
    anchor_fasta: Path,
    prefix: str,
    threads: int,
    evalue: float,
    min_pident: float,
    min_cov: float,
) -> pd.DataFrame:
    _require("diamond")
    db = PROC_DIR / f"human_db_{prefix}.dmnd"
    hits = PROC_DIR / f"human_hits_{prefix}.tsv"
    _run(["diamond", "makedb", "--in", str(human_fasta), "-d", str(db.with_suffix(""))])
    _run(
        [
            "diamond",
            "blastp",
            "-q",
            str(anchor_fasta),
            "-d",
            str(db),
            "-o",
            str(hits),
            "-p",
            str(threads),
            "--evalue",
            str(evalue),
            "--max-target-seqs",
            "5",
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "qcovhsp",
            "evalue",
            "bitscore",
        ]
    )
    cols = [
        "anchor_uniprot",
        "target_uniprot",
        "pident",
        "length",
        "coverage",
        "evalue",
        "bitscore",
    ]
    df = pd.read_csv(hits, sep="\t", names=cols)
    # keep the best hit per anchor protein (highest bitscore)
    df = df.sort_values("bitscore", ascending=False).drop_duplicates("anchor_uniprot")
    df["tier"] = "human"
    df["species"] = "Homo_sapiens"
    df["method"] = "diamond"
    df["is_human_homolog"] = (
        (df["pident"] >= min_pident)
        & (df["coverage"] >= min_cov)
        & (df["evalue"] <= evalue)
    )
    return df.drop(columns=["length"])


# --------------------------- main ---------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(ANCHORS), default="kpneumoniae")
    ap.add_argument("--tiers", default="a,b,c,d")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--evalue", type=float, default=1e-4)
    ap.add_argument("--min-pident", type=float, default=30.0)
    ap.add_argument("--min-cov", type=float, default=50.0)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    anchor = ANCHORS[args.organism]
    prefix, label, anchor_pid = anchor["prefix"], anchor["label"], anchor["proteome_id"]
    of_input = PROC_DIR / f"of_input_{prefix}"
    out_tsv = PROC_DIR / f"{prefix}_orthologs_long.tsv"
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    PROC_DIR.mkdir(parents=True, exist_ok=True)

    tiers = [t.strip() for t in args.tiers.split(",") if t.strip()]
    # tier letter -> (tier_label, species list); the anchor is excluded by proteome_id
    candidate = {}
    if "a" in tiers:
        candidate[anchor["close"]] = CLOSE[anchor["close"]]
    if "b" in tiers:
        candidate["gram_negative"] = GRAM_NEGATIVE
    if "c" in tiers:
        candidate["bacteria"] = BACTERIA

    species, tier_of_species = [], {}
    for tier_label, entries in candidate.items():
        for lbl, taxon, pid in entries:
            if pid == anchor_pid:  # never compare the anchor to itself
                continue
            species.append((lbl, taxon, pid))
            tier_of_species[lbl] = tier_label

    frames: list[pd.DataFrame] = []
    if species:
        print(
            f"Anchor {label} ({args.organism}); bacterial tiers: "
            f"{ {t: sum(1 for s in tier_of_species.values() if s == t) for t in candidate} }"
        )
        results_dir = run_orthofinder(
            anchor, species, of_input, args.threads, args.force
        )
        frames.append(parse_orthofinder(results_dir, label, tier_of_species))

    if "d" in tiers:
        human_fasta = fetch_proteome(*HUMAN, force=args.force)
        anchor_fasta = anchor_accession_fasta(anchor, args.force)
        frames.append(
            run_diamond_human(
                human_fasta,
                anchor_fasta,
                prefix,
                args.threads,
                args.evalue,
                args.min_pident,
                args.min_cov,
            )
        )

    long = pd.concat(frames, ignore_index=True)
    for col in [
        "orthogroup",
        "pident",
        "coverage",
        "evalue",
        "bitscore",
        "is_human_homolog",
    ]:
        if col not in long.columns:
            long[col] = pd.NA
    long = long[
        [
            "anchor_uniprot",
            "tier",
            "species",
            "target_uniprot",
            "method",
            "orthogroup",
            "pident",
            "coverage",
            "evalue",
            "bitscore",
            "is_human_homolog",
        ]
    ]
    long.to_csv(out_tsv, sep="\t", index=False)

    anchor_ids = set(pd.read_csv(anchor["tsv"], sep="\t")["Entry"])
    print(f"\nWrote {len(long)} rows to {out_tsv.relative_to(REPO_ROOT)}")
    for tier_name in long["tier"].unique():
        covered = long.loc[long.tier == tier_name, "anchor_uniprot"].nunique()
        print(
            f"  tier {tier_name}: {covered}/{len(anchor_ids)} anchor proteins have >=1 mapping"
        )


if __name__ == "__main__":
    main()
