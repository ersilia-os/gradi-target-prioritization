"""Publication-based (experimental) essentiality consolidation (docs §4.1/§4.2, no predictions).

Gathers the DIRECT experimental essentiality evidence — every transposon / CRISPRi screen we have —
into one per-protein table, keyed by HS11286 UniProt accession. This is deliberately prediction-free:
only measured data from published screens. The headline is the Kp Mobile-CRISPRi-seq screen.

Sources & mapping onto HS11286:
  * **Jana/Zhu 2023 CRISPRi** (aem.00956-23 s0001):
      - `870 selected essential genes` — the Mobile-CRISPRi conditionally-essential library
        (`KPN_` MGH78578 tags -> HS11286 by DIAMOND, 868/870 map);
      - `in vitro screening-KPNIH1` — Ramage 2017 KPNIH1 essential set (424; gene-symbol) + its
        KPNIH1-Tnseq / E. coli-essential flags;
      - `in vivo screening-KPPR1` — Bachman 2015 KPPR1 in-vivo depletion (`VK055_`; ratio + p-value,
        mapped VK055_ -> gene symbol via Mike & Bachman 2023).
  * **Eichelberger/Short 2024 ECL8** in-vitro essential + urine/serum niche (from 07b, gene-symbol).
  * **Mike & Bachman 2023 KPPR1** in-vivo TnSeq defect (from 07b).
  * **Enterobacteriaceae-TraDIS compendium** (Goodall/Gardner) — the per-genome experimental essential
    call across 12 genomes (2 Klebsiella, 3 E. coli, Salmonella, Citrobacter), decoded from the
    `giant-tab` numeric TraDIS log-ratio (essential ≤ -0.5, calibrated against the curated EcoGene set),
    mapped via the E. coli Keio b-number ortholog. Gives a cross-species conservation-of-essentiality count.

Output: output/results/<org>/<prefix>_ess_publications.csv (keyed by uniprot_accession).
Run with the `gradi` env (DIAMOND from gradi-ortho). No network (07a/Chrome already fetched the data).
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

TR = E.REPO_ROOT / "data" / "raw" / "other" / "essentiality" / "enterobacteriaceae_tradis"
JANA = E.essentiality_raw_dir("kpneumoniae", "jana2023_crispri") / "aem.00956-23-s0001.xlsx"
MB = E.essentiality_raw_dir("kpneumoniae", "mikebachman2023_KPPR1") / "ppat.1011233.s011.xlsx"
TRADIS_ESS_CUTOFF = -0.5   # decoded: giant-tab numeric TraDIS log-ratio <= -0.5 == essential
CORE_FRAC = 0.80           # essential in >= this fraction of covered genomes -> "core essential"
BNUM_RE = re.compile(r"^b\d{4}$")
KPN_RE = re.compile(r"KPN_\d+")


# --------------------------------------------------------------------------- mapping helpers
def kpn_to_acc(org: str) -> dict[str, str]:
    """MGH78578 KPN_ locus tag -> HS11286 accession (best DIAMOND hit), reusing the 07f map."""
    faa = E.essentiality_raw_dir(org, "strains") / "mgh78578_KPN.faa"
    tsv = E.essentiality_processed_dir(org, "fba") / "hs11286_vs_mgh78578.tsv"
    if not faa.exists():
        return {}
    m = E.map_strain_by_sequence(org, faa, tsv)  # uniprot_accession, strain_locus(KPN_), bitscore
    m = m.sort_values("bitscore", ascending=False).drop_duplicates("strain_locus")
    return dict(zip(m["strain_locus"], m["uniprot_accession"]))


def vk_to_symbol() -> dict[str, str]:
    """VK055_ locus -> gene symbol, from Mike & Bachman 2023 s011 (ID + Gene_name)."""
    if not MB.exists():
        return {}
    d = pd.read_excel(MB, sheet_name="S1 Table. TnSeq Results")
    out = {}
    for i, g in zip(d["ID"].astype(str), d["Gene_name"].astype(str)):
        if i and g and g.lower() != "nan":
            out[i.strip()] = g.strip()
    return out


# --------------------------------------------------------------------------- Jana CRISPRi
def jana_crispri(org: str, g2a: dict[str, str]) -> pd.DataFrame:
    base = pd.DataFrame({"uniprot_accession": E.load_accessions(org)})
    if not JANA.exists():
        for c in ["crispri_ce_library", "crispri_ce_group", "crispri_invivo_ratio",
                  "crispri_invivo_log2ratio", "crispri_invivo_hit", "kpnih1_essential"]:
            base[c] = pd.NA
        return base

    k2a = kpn_to_acc(org)

    # (1) 870-gene CRISPRi conditionally-essential library
    lib = pd.read_excel(JANA, sheet_name="870 selected essential genes", header=1)
    lib_map = {}   # acc -> group
    for kpn_raw, grp in zip(lib["Locus_tag (MGH 78578)"].astype(str), lib["Group "].astype(str)):
        for kpn in KPN_RE.findall(kpn_raw):
            acc = k2a.get(kpn)
            if acc:
                lib_map.setdefault(acc, grp)
    base["crispri_ce_library"] = base["uniprot_accession"].isin(lib_map)
    base["crispri_ce_group"] = base["uniprot_accession"].map(lib_map)

    # (2) in-vivo KPPR1 depletion (VK055_ ratio + Ceder p) -> map VK055_ -> symbol -> acc
    v2s = vk_to_symbol()
    iv = pd.read_excel(JANA, sheet_name="in vivo screening-KPPR1", header=1)
    inv = {}   # acc -> (ratio, pvalue)
    for locus, ratio, pval in zip(iv.iloc[:, 0].astype(str), iv["ratio"], iv["Ceder p-value"]):
        tag = locus if locus.startswith("VK055_") else f"VK055_{locus.strip()}"
        sym = v2s.get(tag)
        acc = g2a.get((sym or "").lower())
        if acc and pd.notna(ratio):
            r = float(ratio)
            if acc not in inv or r > inv[acc][0]:
                inv[acc] = (r, float(pval) if pd.notna(pval) else np.nan)
    base["crispri_invivo_ratio"] = base["uniprot_accession"].map(lambda a: inv.get(a, (np.nan,))[0])
    base["crispri_invivo_log2ratio"] = np.log2(base["crispri_invivo_ratio"].clip(lower=1e-3))
    base["crispri_invivo_hit"] = base["uniprot_accession"].map(
        lambda a: (a in inv and inv[a][0] > 2 and (np.isnan(inv[a][1]) or inv[a][1] < 0.05)))

    # (3) KPNIH1 in-vitro essential (Ramage 2017) via gene symbol
    kp = pd.read_excel(JANA, sheet_name="in vitro screening-KPNIH1", header=1)
    kpnih1 = {}
    for gene in kp["Gene"].astype(str):
        acc = g2a.get(gene.strip().lower())
        if acc:
            kpnih1[acc] = True
    base["kpnih1_essential"] = base["uniprot_accession"].isin(kpnih1)
    return base


# --------------------------------------------------------------------------- compendium cross-species
def cross_species(org: str) -> pd.DataFrame:
    accs = E.load_accessions(org)
    base = pd.DataFrame({"uniprot_accession": accs})
    gt = TR / "giant-tab_final.tsv"
    if not gt.exists():
        return base
    ess_cols = {  # column -> short genome label
        "TraDIS Essentiality: Klebsiella pneumoniae Ecl8": "K. pneumoniae ECL8",
        "TraDIS Essentiality: Klebsiella pneumoniae RH201207": "K. pneumoniae RH201207",
        "TraDIS Essentiality: Escherichia coli BW25113": "E. coli BW25113",
        "TraDIS Essentiality: Escherichia coli ST131 EC958": "E. coli EC958",
        "TraDIS Essentiality: Escherichia coli UPEC ST131 NCTC13441": "E. coli NCTC13441",
        "TraDIS Essentiality: Citrobacter rodentium ICC168": "C. rodentium ICC168",
        "TraDIS Essentiality: Salmonella Typhi Ty2": "S. Typhi Ty2",
        "TraDIS Essentiality: Salmonella Typhimurium A130": "S. Tm A130",
        "TraDIS Essentiality: Salmonella Typhimurium D23580": "S. Tm D23580",
        "TraDIS Essentiality: Salmonella Typhimurium SL3261": "S. Tm SL3261",
        "TraDIS Essentiality: Salmonella Typhimurium SL1344": "S. Tm SL1344",
        "TraDIS Essentiality: Salmonella Enteritidis P125109": "S. Enteritidis P125109",
    }
    bcol = "Locus: Escherichia coli BW25113 (Keio)"
    g = pd.read_csv(gt, sep="\t", low_memory=False, usecols=[bcol] + list(ess_cols))
    g["bnumber"] = g[bcol].astype(str).str.strip()
    g = g[g["bnumber"].str.match(BNUM_RE)]
    for col in ess_cols:
        g[col] = pd.to_numeric(g[col], errors="coerce") <= TRADIS_ESS_CUTOFF
    gb = g.groupby("bnumber")[list(ess_cols)].max()

    # anchor accession -> best E. coli b-number (via orthology), reusing 07c bridge
    acc2b_ec = _ecoli_acc_to_bnumber()
    if org == "ecoli":
        amap = {a: acc2b_ec.get(a) for a in accs}
    else:
        orth = E.load_orthologs(org)
        ec = orth[orth["species"] == "Ecoli_K12_MG1655"].copy()
        ec["bnumber"] = ec["target_uniprot"].map(acc2b_ec)
        ec = ec.dropna(subset=["bnumber"]).drop_duplicates("anchor_uniprot")
        amap = dict(zip(ec["anchor_uniprot"], ec["bnumber"]))

    short_labels = list(ess_cols.values())
    rows = []
    for a in accs:
        b = amap.get(a)
        r = gb.loc[b] if (b is not None and b in gb.index) else None
        rec = {"uniprot_accession": a, "pub_bnumber": b or ""}
        n = 0
        for col, lab in ess_cols.items():
            v = bool(r[col]) if r is not None else False
            rec[f"pub_ess__{lab}"] = v
            n += int(v)
        rec["pub_n_species_essential"] = n if r is not None else 0
        rec["pub_covered"] = r is not None
        rec["pub_frac_species_essential"] = (n / len(ess_cols)) if r is not None else np.nan
        rows.append(rec)
    out = pd.DataFrame(rows)
    out["pub_core_essential"] = out["pub_frac_species_essential"] >= CORE_FRAC
    out.attrs["genomes"] = short_labels
    return out


def _ecoli_acc_to_bnumber() -> dict[str, str]:
    d = pd.read_csv(E.proteome_tsv("ecoli"), sep="\t", usecols=["Entry", "Gene Names"])
    out = {}
    for acc, gn in zip(d["Entry"], d["Gene Names"].fillna("")):
        for tok in str(gn).split():
            if BNUM_RE.match(tok):
                out[str(acc)] = tok
                break
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = E.ORGANISMS[org]
    g2a = {g: a for g, a in _kp_gene_map(org).items()}

    df = cross_species(org)
    genomes = df.attrs.get("genomes", [])

    # Kp direct screens
    if org == "kpneumoniae":
        cr = jana_crispri(org, g2a)
        df = df.merge(cr, on="uniprot_accession", how="left")
        # ECL8 + KPPR1 in-vivo from 07b
        kpb = E.results_dir(org) / f"{prefix}_ess_kp.csv"
        if kpb.exists():
            b = pd.read_csv(kpb)
            df["ecl8_essential"] = df["uniprot_accession"].map(
                b.set_index("uniprot_accession")["kp_ess_in_vitro_call"].eq("essential"))
            df["ecl8_urine_req"] = df["uniprot_accession"].map(
                b.set_index("uniprot_accession")["kp_ess_urine_call"].eq("conditional"))
            df["ecl8_serum_req"] = df["uniprot_accession"].map(
                b.set_index("uniprot_accession")["kp_ess_serum_call"].eq("conditional"))
            df["kppr1_invivo_defect"] = df["uniprot_accession"].map(
                b.set_index("uniprot_accession")["kp_ess_in_vivo_call"].eq("in_vivo_defect"))
        screen_cols = ["crispri_ce_library", "kpnih1_essential", "ecl8_essential"]
    else:
        # E. coli — first-class direct screens from 07n (Keio, Goodall, Rousset18/Wang18 CRISPRi, ...)
        for c in ["crispri_ce_library", "crispri_invivo_hit", "kpnih1_essential",
                  "ecl8_essential", "ecl8_urine_req", "ecl8_serum_req", "kppr1_invivo_defect"]:
            df[c] = pd.NA  # Kp-specific columns kept NA for schema uniformity
        exp = E.results_dir(org) / "ec_ess_experimental.csv"
        if exp.exists():
            e = pd.read_csv(exp).set_index("uniprot_accession")
            for c in ["ecoli_keio_essential", "ecoli_goodall_essential",
                      "ecoli_crispri_rousset18_essential", "ecoli_crispri_wang18_essential",
                      "ecoli_experimental_essential", "ecoli_vulnerability_score",
                      "ecoli_crispri_rousset18_log2fc", "ecoli_crispri_wang18_fitness",
                      "ecoli_crispri_multistrain_frac", "ecoli_rbtnseq_min_fitness"]:
                if c in e.columns:
                    df[c] = df["uniprot_accession"].map(e[c])
        screen_cols = ["ecoli_keio_essential", "ecoli_goodall_essential",
                       "ecoli_crispri_rousset18_essential", "ecoli_crispri_wang18_essential"]

    # count of independent direct experimental screens calling this protein essential (per organism)
    present = [c for c in screen_cols if c in df]
    if present:
        bmat = pd.concat([df[c].map(lambda v: v is True) for c in present], axis=1)
        df["n_experimental_screens"] = bmat.sum(axis=1).astype(int)
    else:
        df["n_experimental_screens"] = 0
    df["n_kp_experimental_screens"] = df["n_experimental_screens"]  # back-compat alias
    if org == "ecoli" and "ecoli_experimental_essential" in df:
        df["experimental_essential"] = df["ecoli_experimental_essential"].fillna(False).astype(bool)
    else:
        ess_any = pd.Series(False, index=df.index)
        for c in screen_cols:
            if c in df:
                ess_any = ess_any | (df[c] == True)  # noqa: E712 — NA-safe
        df["experimental_essential"] = ess_any

    out = E.results_dir(org) / f"{prefix}_ess_publications.csv"
    df.to_csv(out, index=False)

    print(f"[{org}] {len(df)} proteins; cross-species covered {int(df['pub_covered'].sum())}; "
          f"core-essential (≥{int(CORE_FRAC*100)}% of {len(genomes)} genomes) "
          f"{int(df['pub_core_essential'].sum())}", flush=True)
    if org == "kpneumoniae":
        print(f"  CRISPRi 870-library mapped: {int(df.crispri_ce_library.sum())}; "
              f"KPNIH1 essential: {int(df.kpnih1_essential.sum())}; "
              f"CRISPRi in-vivo hits: {int(df.crispri_invivo_hit.sum())}; "
              f"ECL8 essential: {int(df.ecl8_essential.sum())}", flush=True)
    else:
        print(f"  experimentally essential (Keio/Goodall/CRISPRi): {int(df.experimental_essential.sum())}; "
              f">=2 screens: {int((df.n_experimental_screens>=2).sum())}", flush=True)
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)}", flush=True)


def _kp_gene_map(org: str) -> dict[str, str]:
    """Lowercased gene symbol -> accession (all Gene-Names tokens; paralog suffixes stripped)."""
    locus_like = re.compile(r"^(kphs_|kpn_|vk055|b\d{4}$|jw)", re.IGNORECASE)
    suf = re.compile(r"_\d+$")
    d = pd.read_csv(E.proteome_tsv(org), sep="\t", usecols=["Entry", "Gene Names"])
    out = {}
    for acc, gn in zip(d["Entry"], d["Gene Names"].fillna("")):
        for tok in str(gn).split():
            if tok and not locus_like.match(tok):
                out.setdefault(suf.sub("", tok.lower()), str(acc))
    return out


if __name__ == "__main__":
    main()
