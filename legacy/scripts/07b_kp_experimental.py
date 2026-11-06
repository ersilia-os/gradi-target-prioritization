"""Direct K. pneumoniae experimental essentiality (docs §4.1a/b/c).

Ingests the primary Kp transposon-insertion / CRISPRi screens fetched by 07a and lands each per-gene
call on the HS11286 anchor. The screens are keyed by strain-specific locus tags (ECL8 `ecl8_`, KPPR1
`VK055_`) plus a gene symbol; since HS11286's essential genes are the well-conserved, named ones, we
map by **gene symbol** (normalised: lowercased, trailing `_<n>` paralog suffixes stripped), which
captures the essentiality signal cleanly. Coverage (named genes matched) is reported per track.

Tracks:
  * 4.1a in-vitro essential — Eichelberger/Short 2024 ECL8 (fig1-data1 `Essential` bimodal call).
  * 4.1a conditional (niche) — ECL8 urine (fig4) & serum (fig6): required in-niche if logFC<0 & q<0.05.
  * 4.1b in-vivo fitness — Mike & Bachman 2023 KPPR1 TnSeq (s011): defect if log2FC<0 & p<0.05.
  * 4.1c vulnerability — the hand-curated Jana/Zhu 2023 Mobile-CRISPRi highlights (gene symbols).

This is a Kp-specific track; `--organism ecoli` writes an all-NA table so the 07h merge schema is
uniform (E. coli's own essentiality is the direct 07c call). Output:
  output/results/<org>/<prefix>_ess_kp.csv   (keyed by uniprot_accession)
Run with the `gradi` conda env interpreter. No network (07a already fetched the data).
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

RAW = lambda *p: E.essentiality_raw_dir("kpneumoniae", *p)  # noqa: E731
ECL8 = "eichelberger2024_ECL8"
MB2023 = "mikebachman2023_KPPR1"
COND_Q = 0.05          # niche/in-vivo significance cutoff
COND_LFC = 0.0         # depletion direction (required -> negative logFC)
SUFFIX_RE = re.compile(r"_\d+$")

OUT_COLS = [
    "uniprot_accession",
    "kp_ess_in_vitro_call", "kp_ess_in_vitro_score",
    "kp_ess_urine_call", "kp_ess_urine_score",
    "kp_ess_serum_call", "kp_ess_serum_score",
    "kp_ess_in_vivo_call", "kp_ess_in_vivo_score",
    "kp_ess_vulnerability_call",
    "kp_ess_sources",
]


def norm_gene(g: str) -> str:
    g = str(g).strip().lower()
    return SUFFIX_RE.sub("", g)


def gene_map() -> dict[str, str]:
    """Normalised gene symbol -> HS11286 accession, from ALL Gene-Names tokens (not just the first)."""
    locus_like = re.compile(r"^(kphs_|b\d{4}$|jw)", re.IGNORECASE)
    tsv = E.proteome_tsv("kpneumoniae")
    d = pd.read_csv(tsv, sep="\t", usecols=["Entry", "Gene Names"])
    out: dict[str, str] = {}
    for acc, gn in zip(d["Entry"], d["Gene Names"].fillna("")):
        for tok in str(gn).split():
            if tok and not locus_like.match(tok):
                out.setdefault(norm_gene(tok), str(acc))
    return out


def _read_offset(path: Path, sheet: str, header_row: int) -> pd.DataFrame:
    return pd.read_excel(path, sheet_name=sheet, header=header_row)


def ecl8_in_vitro(g2a: dict[str, str]) -> pd.DataFrame:
    """ECL8 in-vitro essential call -> {essential, non_essential, unclear} per HS11286 accession."""
    f = RAW(ECL8) / "elife-88971-fig1-data1-v1.xlsx"
    if not f.exists():
        return pd.DataFrame(columns=["uniprot_accession", "call", "score"])
    d = _read_offset(f, "Table S1_Essential gene table", 1)
    rows = []
    for _, r in d.iterrows():
        acc = g2a.get(norm_gene(r.get("Gene_name", "")))
        if not acc:
            continue
        ess = str(r.get("Essential")).strip().lower() == "true"
        unclear = str(r.get("Unclear")).strip().lower() == "true"
        call = "essential" if ess else ("unclear" if unclear else "non_essential")
        rows.append({"uniprot_accession": acc, "call": call,
                     "score": 1.0 if ess else (0.4 if unclear else 0.0)})
    return pd.DataFrame(rows).sort_values("score", ascending=False).drop_duplicates("uniprot_accession")


def ecl8_condition(g2a: dict[str, str], fname: str, sheet: str) -> pd.DataFrame:
    """ECL8 niche fitness (urine/serum): required-in-niche if logFC<0 & q<0.05; score = |logFC| scaled."""
    f = RAW(ECL8) / fname
    if not f.exists():
        return pd.DataFrame(columns=["uniprot_accession", "call", "score"])
    d = _read_offset(f, sheet, 1)
    lfc = pd.to_numeric(d.get("logFC"), errors="coerce")
    q = pd.to_numeric(d.get("q.value"), errors="coerce")
    rows = []
    for i, r in d.iterrows():
        acc = g2a.get(norm_gene(r.get("gene_name", "")))
        if not acc or pd.isna(lfc.iloc[i]):
            continue
        required = bool(lfc.iloc[i] < COND_LFC and (pd.notna(q.iloc[i]) and q.iloc[i] < COND_Q))
        rows.append({"uniprot_accession": acc,
                     "call": "conditional" if required else "not_required",
                     "score": float(min(1.0, abs(lfc.iloc[i]) / 5.0)) if required else 0.0})
    if not rows:
        return pd.DataFrame(columns=["uniprot_accession", "call", "score"])
    return pd.DataFrame(rows).sort_values("score", ascending=False).drop_duplicates("uniprot_accession")


def mikebachman_invivo(g2a: dict[str, str]) -> pd.DataFrame:
    """KPPR1 in-vivo TnSeq (Mike & Bachman 2023 s011): fitness defect if log2FC<0 & p<0.05."""
    f = RAW(MB2023) / "ppat.1011233.s011.xlsx"
    if not f.exists():
        return pd.DataFrame(columns=["uniprot_accession", "call", "score"])
    d = _read_offset(f, "S1 Table. TnSeq Results", 0)
    lfc = pd.to_numeric(d.get("log2FC(Output/Input)"), errors="coerce")
    p = pd.to_numeric(d.get("pvalue"), errors="coerce")
    rows = []
    for i, r in d.iterrows():
        acc = g2a.get(norm_gene(r.get("Gene_name", "")))
        if not acc or pd.isna(lfc.iloc[i]):
            continue
        defect = bool(lfc.iloc[i] < 0 and pd.notna(p.iloc[i]) and p.iloc[i] < COND_Q)
        rows.append({"uniprot_accession": acc,
                     "call": "in_vivo_defect" if defect else "no_defect",
                     "score": float(min(1.0, abs(lfc.iloc[i]) / 5.0)) if defect else 0.0})
    if not rows:
        return pd.DataFrame(columns=["uniprot_accession", "call", "score"])
    return pd.DataFrame(rows).sort_values("score", ascending=False).drop_duplicates("uniprot_accession")


def crispri_vulnerability(g2a: dict[str, str]) -> pd.DataFrame:
    """Hand-curated Jana/Zhu 2023 Mobile-CRISPRi highlights (gene symbols)."""
    f = E.legacy_essentiality_dir("literature", "zhu2023_kp_crispri", "curated_highlights.tsv")
    if not f.exists():
        return pd.DataFrame(columns=["uniprot_accession", "call"])
    d = pd.read_csv(f, sep="\t")
    rows = []
    for _, r in d.iterrows():
        acc = g2a.get(norm_gene(r.get("gene_symbol", "")))
        if acc:
            rows.append({"uniprot_accession": acc, "call": "vulnerable"})
    return pd.DataFrame(rows).drop_duplicates("uniprot_accession")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    accs = E.load_accessions(org)

    if org != "kpneumoniae":
        df = pd.DataFrame({"uniprot_accession": accs})
        for c in OUT_COLS[1:]:
            df[c] = pd.NA
        out = E.results_dir(org) / f"{prefix}_ess_kp.csv"
        df[OUT_COLS].to_csv(out, index=False)
        print(f"[{org}] Kp-experimental track is Kp-specific; wrote all-NA {out.relative_to(E.REPO_ROOT)}", flush=True)
        return

    g2a = gene_map()
    print(f"gene-symbol bridge: {len(g2a)} normalised symbols -> HS11286 accession", flush=True)

    invitro = ecl8_in_vitro(g2a)
    urine = ecl8_condition(g2a, "elife-88971-fig4-data1-v1.xlsx", "Table S5_ECL8_Urinome")
    serum = ecl8_condition(g2a, "elife-88971-fig6-data1-v1.xlsx", "Table S5_ECL8_Serum_Resistome")
    invivo = mikebachman_invivo(g2a)
    vuln = crispri_vulnerability(g2a)

    df = pd.DataFrame({"uniprot_accession": accs})

    def attach(sub, call_col, score_col=None, src=None):
        nonlocal df
        if sub.empty:
            df[call_col] = pd.NA
            if score_col:
                df[score_col] = pd.NA
            return 0
        m = sub.set_index("uniprot_accession")
        df[call_col] = df["uniprot_accession"].map(m["call"])
        if score_col and "score" in m:
            df[score_col] = df["uniprot_accession"].map(m["score"])
        return int(m.index.isin(accs).sum())

    n_iv = attach(invitro, "kp_ess_in_vitro_call", "kp_ess_in_vitro_score")
    n_ur = attach(urine, "kp_ess_urine_call", "kp_ess_urine_score")
    n_se = attach(serum, "kp_ess_serum_call", "kp_ess_serum_score")
    n_vv = attach(invivo, "kp_ess_in_vivo_call", "kp_ess_in_vivo_score")
    n_vu = attach(vuln, "kp_ess_vulnerability_call")

    # provenance string
    def sources(r):
        s = []
        if pd.notna(r.get("kp_ess_in_vitro_call")): s.append("ECL8")
        if r.get("kp_ess_urine_call") == "conditional": s.append("urine")
        if r.get("kp_ess_serum_call") == "conditional": s.append("serum")
        if r.get("kp_ess_in_vivo_call") == "in_vivo_defect": s.append("KPPR1_invivo")
        if pd.notna(r.get("kp_ess_vulnerability_call")): s.append("CRISPRi")
        return ";".join(s)
    df["kp_ess_sources"] = df.apply(sources, axis=1)

    out = E.results_dir(org) / f"{prefix}_ess_kp.csv"
    df[OUT_COLS].to_csv(out, index=False)
    n_ess = int((df["kp_ess_in_vitro_call"] == "essential").sum())
    print(f"[{org}] ECL8 in-vitro mapped: {n_iv} (essential={n_ess}); urine {n_ur}; serum {n_se}; "
          f"in-vivo {n_vv}; CRISPRi {n_vu}", flush=True)
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)} ({len(df)} proteins)", flush=True)


if __name__ == "__main__":
    main()
