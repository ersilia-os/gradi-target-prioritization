"""Find a NON-REGISTRY proteome by label, whatever directory and format it happens to be in.

A "strain" here is any proteome that is not one of the four registry species: a tier-D screen
strain (ECL8, KPPR1), a Kp strain staged by `essentiality/strain_proteomes.py`, a DEG dataset
proteome, or an OGEE taxon. They share three properties that separate them from the anchors:

  * no UniProt proteome, so no `uniprot_ac` and no stage-00 table
  * no canonical row order, so `src.matrices` has no reference to check them against
  * keyed on whatever identifier their own FASTA uses -- a RefSeq `WP_*`, a locus tag, a DEG
    `lcl|...` header. That is deliberate: label, sequence and embedding then share one namespace
    and nothing has to be joined across annotations.

WHY THIS MODULE EXISTS. The same 12-line resolver had been copy-pasted into
`embeddings/prott5.py`, `embeddings/proteomelm.py` and `essentiality/paralog_clusters.py`, and
`embeddings/esmc.py` was about to become the fourth. Four copies means four chances for the search
path to drift, and a strain that resolves under one script but not another looks like a missing
file rather than a missing directory entry.

`uniprot_ac` keeps its name even though these ids are not UniProt accessions: the column is the row
key, and renaming it per source would fork every consumer downstream.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]

# Ordered: the first location holding `<label><ext>` wins. Add new providers here and every script
# that resolves a strain gains them at once.
SOURCES: tuple[tuple[Path, str], ...] = (
    (REPO_ROOT / "data" / "source" / "uniprot" / "proteomes" / "ncbi", ".faa"),   # tier-D registry
    (REPO_ROOT / "data" / "source" / "ncbi" / "kp_strains", ".tsv"),              # strain_proteomes
    (REPO_ROOT / "data" / "source" / "ncbi" / "deg_proteomes", ".faa"),           # DEG datasets
    (REPO_ROOT / "data" / "source" / "ncbi" / "ogee_proteomes", ".faa"),          # OGEE taxa
)


def resolve(label: str) -> tuple[Path, str]:
    """(path, extension) for `label`. Exits with every location searched if it is absent."""
    for d, ext in SOURCES:
        p = d / f"{label}{ext}"
        if p.exists():
            return p, ext
    searched = "\n    ".join(str((d / f"{label}{e}").relative_to(REPO_ROOT)) for d, e in SOURCES)
    raise FileNotFoundError(f"no proteome for strain {label!r}. Looked in:\n    {searched}")


def load_frame(label_or_path: str | Path, ext: str | None = None) -> pd.DataFrame:
    """(uniprot_ac, sequence) for a strain, from either a protein FASTA or a two-column TSV.

    Accepts a label or an already-resolved path. Rows are de-duplicated on the key and SORTED by
    it: sorting makes shard boundaries deterministic across runs and independent of whatever order
    the source file happens to be in.
    """
    if ext is None:
        path, ext = resolve(str(label_or_path))
    else:
        path = Path(label_or_path)

    if ext == ".tsv":
        d = pd.read_csv(path, sep="\t").rename(columns={"locus_tag": "uniprot_ac"})
        d = d[["uniprot_ac", "sequence"]]
    else:
        ids, seqs, cur, buf = [], [], None, []
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                if cur is not None:
                    ids.append(cur)
                    seqs.append("".join(buf))
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
        if cur is not None:
            ids.append(cur)
            seqs.append("".join(buf))
        d = pd.DataFrame({"uniprot_ac": ids, "sequence": seqs})

    d["uniprot_ac"] = d["uniprot_ac"].astype(str)
    d["sequence"] = d["sequence"].astype(str)
    return (d[d["sequence"].str.strip().ne("")]
            .drop_duplicates("uniprot_ac")
            .sort_values("uniprot_ac")
            .reset_index(drop=True))
