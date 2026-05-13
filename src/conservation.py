"""Conservation annotation from BV-BRC PATtyFam protein-family assignments.

The BV-BRC `genome_feature` table for HS11286 (genome_id 1125630.4) has been
staged at `data/raw/bvbrc/hs11286_features.tsv` with columns:
    patric_id, refseq_locus_tag (KPHS_*), gene, product, plfam_id, pgfam_id

For a proper conservation_class call, we need the number of K. pneumoniae genomes
in which each PLFam appears. That aggregation requires either a separate BV-BRC
bulk query or downloading and counting against the species-level pan-genome —
deferred to a v2 pass. For v1 we record:
    - plfam_id (present / absent)
    - pgfam_id (present / absent)
    - has_plfam (boolean — quick proxy for "this protein is part of the Kp
      pan-protein-family system at all"; ~99% True for chromosomal genes)

Downstream, the user can join this to a `kp_plfam_counts.tsv` (PLFam -> n_genomes)
to compute conservation_class properly.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

BVBRC_FEATURES = Path("data/raw/bvbrc/hs11286_features.tsv")


def load_bvbrc_features() -> pd.DataFrame:
    if not BVBRC_FEATURES.exists():
        return pd.DataFrame(columns=["kp_locus_tag", "plfam_id", "pgfam_id", "product"])
    df = pd.read_csv(BVBRC_FEATURES, sep="\t", dtype=str)
    df = df.rename(columns={"refseq_locus_tag": "kp_locus_tag"})
    # BV-BRC TSV quotes string fields; strip them
    for c in ["kp_locus_tag", "plfam_id", "pgfam_id", "product", "gene"]:
        if c in df.columns:
            df[c] = df[c].str.strip('"').replace({"": None})
    df = df.drop_duplicates(subset=["kp_locus_tag"], keep="first")
    df["has_plfam"] = df["plfam_id"].notna()
    return df[["kp_locus_tag", "plfam_id", "pgfam_id", "product", "has_plfam"]]
