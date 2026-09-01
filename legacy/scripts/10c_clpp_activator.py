"""Activated-ClpP (partnerless) substrate evidence, transferred onto Kp/Ec (docs §3.3c).

WHY THIS TRACK EXISTS — it is the only one that measures the criterion Gr-ADI actually gates on
------------------------------------------------------------------------------------------------
The Gr-ADI research vision (`Gr-ADI_ResearchVision_TPD-vs-Gram-Negs_Draft-v5(Final)`) validates every
candidate target with an in-vitro assay against "the activated ClpP of the target organism **only
(i.e. in the absence of an unfoldase partner)**", and states that "only targets that show degradation
in this assay will progress". That is a *different mechanism* from every other evidence track in this
stage:

  * 10a's trap papers (Flynn, Neher, Bhat, Feng, Graham, Lunge, Ziemski) are **unfoldase-dependent** —
    ClpXP/ClpAP/ClpCP. There an unstructured terminus is enough, because the AAA+ ring grips it and
    does the unfolding work itself.
  * Activated ClpP **has no unfoldase**. ADEP-class activators bind the hydrophobic clefts where the
    ClpX/ClpA (L/I)GF loops dock, flipping ClpP compressed->extended: the axial pore opens and the
    catalytic triad orders, but nothing pulls. So a substrate must ALREADY be unstructured enough to
    diffuse into an open pore.

Consequence: a good ClpXP substrate is not automatically an activated-ClpP substrate, and the two
bars can disagree. This script supplies the evidence for the partnerless bar; `10b_degrons.py` and the
trap pool supply the unfoldase-dependent one. They are kept as SEPARATE columns and never merged —
their disagreement is a deliverable, not noise.

The two datasets
----------------
Both are S. aureus (no equivalent exists for E. coli or K. pneumoniae — a gap the Gr-ADI SoW itself
notes and hopes to fill through the consortium), and each ships TWO readouts:

  Conlon 2013 (Nature 503:365), ADEP4 vs untreated MRSA, iTRAQ:
    'Table S1' 1,712 proteins  — fully tryptic peptides    -> ABUNDANCE, `Average` = log2(ADEP4/ctrl)
    'Table S2' 2,382 peptides  — partially tryptic peptides -> CLEAVAGE (one terminus made by an
                                 endogenous protease, i.e. a direct proteolysis product)
  Jacques 2020 (Genetics 214:1103), 30 uM ONC212 vs untreated S. aureus, label-free LC-MS/MS-FAIMS:
    'Table S3' 1,620 proteins  — `24H_log2_fold-change`                       -> ABUNDANCE
                                 `10-40_minutes_non-tryptic_peptides_log2...` -> CLEAVAGE

**The cleavage readouts are the better evidence and are weighted above abundance.** A 24 h abundance
drop also moves with growth arrest, regulon change and resynthesis; an early rise in non-/partially
tryptic peptides is a direct cleavage product. The clearest demonstration is AcpP, whose ONC212
abundance change is +0.03 (nothing at all) while its cleavage signal is +3.46 — actively cut, steady
state held up by resynthesis. Scoring abundance alone would have called that protein a non-substrate.

Identifier route (all four hops are reported with counts; none is assumed)
-------------------------------------------------------------------------
Both papers predate the current RefSeq annotation, and the bridge is not obvious:

  * They key on SACOL locus tags (S. aureus COL, 2013). The CURRENT RefSeq annotation of
    GCF_000012045.1 has re-tagged everything `SACOL_RS#####` and kept **no** `old_locus_tag`, so the
    SACOL tags cannot be recovered from the modern assembly at all.
  * UniProt has demoted the COL and Mu50 proteomes as redundant — only 939 and 1,033 entries survive,
    so a UniProt-first bridge would silently drop ~45% of Conlon.
  * What does work: each paper carries its own protein accessions (Conlon `Entrez` = `YP_*`, Jacques
    `Gene` = `ODV*`), and NCBI efetch still serves both. So we go to SEQUENCE and let sequence do the
    mapping, which is this repo's convention for external resources anyway.

  paper accession -> (efetch) sequence -> (DIAMOND RBH) Kp/Ec UniProt accession

Cross-phylum transfer (Firmicutes -> Gammaproteobacteria) is the weak link and is labelled as such:
reciprocal best hits only, with %identity and coverage retained on every row so a reader can raise the
bar. This is evidence, not ground truth.

Outputs
-------
  output/results/<org>/<prefix>_clpp_activator.csv   per-protein transferred activator evidence
  data/processed/other/degradability/activator/      cached sequences + DIAMOND hits

Run with the `gradi` conda env interpreter. DIAMOND comes from `gradi-ortho` (see src/ligandability).
"""

from __future__ import annotations

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
EFETCH_BATCH = 200
UA = {"User-Agent": "gradi-target-prioritization (miquel@ersilia.io)"}

# --- evidence thresholds ------------------------------------------------------------------------
# "Reduced by at least 2-fold" is the criterion the Gr-ADI proposal uses when it reads these tables,
# so log2 <= -1 is kept as the abundance call for comparability with the proposal's own wording.
ABUNDANCE_LOG2 = -1.0
# Cleavage is a RISE in endogenous-protease-generated peptides. log2 >= +1 is the mirror threshold.
CLEAVAGE_LOG2 = 1.0
PADJ = 0.05
# Cleavage outranks abundance (see the module docstring). Both are [0,1] sub-scores; a protein
# measured on only one readout keeps that one and is not penalised for the missing one.
W_CLEAVAGE = 0.65
W_ABUNDANCE = 0.35
# Reciprocal-best-hit floor for the cross-phylum hop. DIAMOND defaults (25% id / 50% query cover) are
# already permissive; we keep them and record pident so the operating point can be raised downstream.
RBH_MIN_PIDENT = 25.0


# --------------------------------------------------------------------------- source parsing
def _norm_sacol(v) -> str:
    m = re.match(r"^\s*(SACOL\d+)", str(v).strip(), re.I)
    return m.group(1).upper() if m else ""


def parse_conlon() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Conlon 2013 Table S1 (abundance) and Table S2 (cleavage), per SACOL locus tag.

    Table S2 is PEPTIDE-level: one row per partially tryptic peptide, so a protein can appear many
    times (DnaK has 25). Aggregated to per-protein as the count of significant cleavage peptides plus
    the maximum log2 ratio — a protein cut at many positions is stronger evidence than one cut once.
    """
    d = D.xbac_raw_dir("conlon2013_adep4_saureus")
    files = sorted(d.glob("*.xlsx"))
    if not files:
        return pd.DataFrame(), pd.DataFrame()
    src = files[0]
    # Header rows differ per sheet (multi-row title/caption preamble): S1 row 3, S2 row 5.
    s1 = pd.read_excel(src, sheet_name="Table S1", header=3)
    s2 = pd.read_excel(src, sheet_name="Table S2", header=5)

    s1 = s1.assign(sacol=s1["Protein"].map(_norm_sacol))
    s1 = s1[s1["sacol"] != ""]
    ab = (s1.groupby("sacol")
            .agg(refseq=("Entrez", "first"), gene=("ShortName", "first"),
                 description=("Description", "first"),
                 adep4_abundance_log2fc=("Average", "mean"),
                 adep4_abundance_padj=("Adjusted p-value", "min"))
            .reset_index())

    s2 = s2.assign(sacol=s2["Protein"].map(_norm_sacol))
    s2 = s2[s2["sacol"] != ""]
    sig = s2[(s2["Average"] >= CLEAVAGE_LOG2) & (s2["Adjusted p-value"] < PADJ)]
    cl = (s2.groupby("sacol")
            .agg(adep4_cleavage_log2fc_max=("Average", "max"),
                 adep4_cleavage_peptides=("Average", "size"))
            .reset_index())
    cl["adep4_cleavage_peptides_sig"] = (cl["sacol"]
                                         .map(sig.groupby("sacol").size())
                                         .fillna(0).astype(int))
    return ab, cl


def parse_jacques() -> pd.DataFrame:
    """Jacques 2020 Table S3 (ONC212, S. aureus), per SACOL via its own `COL_homolog` column."""
    d = D.xbac_raw_dir("jacques2020_onc212_saureus")
    files = sorted(d.glob("*.xlsx"))
    if not files:
        return pd.DataFrame()
    j = pd.read_excel(files[0])
    j = j.rename(columns={
        "Gene": "genbank",
        "24H_log2_fold-change": "onc212_abundance_log2fc",
        "10-40_minutes_non-tryptic_peptides_log2_fold-change": "onc212_cleavage_log2fc",
        "Mu50_homolog": "saureus_mu50_uniprot",
        "Description": "description_onc212",
    })
    j["sacol"] = j["COL_homolog"].map(_norm_sacol)
    keep = ["sacol", "genbank", "saureus_mu50_uniprot", "description_onc212",
            "onc212_abundance_log2fc", "onc212_cleavage_log2fc"]
    j = j[keep]
    # Jacques carries no p-values, so the abundance/cleavage calls for this paper rest on effect size
    # alone; recorded here rather than silently treated as equivalent to Conlon's tested values.
    return j[j["sacol"] != ""].drop_duplicates("sacol").reset_index(drop=True)


# --------------------------------------------------------------------------- sequences
def _cache_dir() -> Path:
    return D.xbac_processed_dir("activator")


def fetch_sequences(accessions: list[str], refresh: bool = False) -> dict[str, str]:
    """RefSeq/GenBank protein accession -> sequence, via NCBI efetch, cached as one FASTA.

    The 2013 accessions are superseded in the current annotation but efetch still serves them, which
    is the only complete route to these proteins' sequences (see the module docstring).
    """
    out = _cache_dir() / "saureus_activator_proteins.faa"
    have: dict[str, str] = {}
    if out.exists() and not refresh:
        acc = None
        for line in out.read_text().splitlines():
            if line.startswith(">"):
                acc = line[1:].split()[0]
                have[acc] = ""
            elif acc:
                have[acc] += line.strip()
    missing = [a for a in accessions if a and a not in have]
    if missing:
        print(f"  [efetch] {len(missing)} sequence(s) to fetch "
              f"({len(have)} cached)", flush=True)
        chunks: list[str] = []
        for i in range(0, len(missing), EFETCH_BATCH):
            batch = missing[i:i + EFETCH_BATCH]
            r = requests.get(EFETCH, params={"db": "protein", "id": ",".join(batch),
                                             "rettype": "fasta", "retmode": "text"},
                             headers=UA, timeout=180)
            r.raise_for_status()
            chunks.append(r.text)
            print(f"    batch {i // EFETCH_BATCH + 1}: {len(batch)} ids, "
                  f"{len(r.text)} bytes", flush=True)
            time.sleep(0.4)          # NCBI courtesy rate limit (no API key)
        acc = None
        for line in "\n".join(chunks).splitlines():
            if line.startswith(">"):
                acc = line[1:].split()[0]
                have.setdefault(acc, "")
            elif acc:
                have[acc] += line.strip()
        with open(out, "w") as fh:
            for a, s in sorted(have.items()):
                if s:
                    fh.write(f">{a}\n{s}\n")
        print(f"  [efetch] cache -> {out.relative_to(D.REPO_ROOT)} ({len(have)} seqs)", flush=True)
    return {a: s for a, s in have.items() if s}


# --------------------------------------------------------------------------- transfer
def rbh_map(organism: str, query_fasta: Path) -> pd.DataFrame:
    """Reciprocal best hits between the activator-dataset proteins and an anchor proteome.

    Forward: S. aureus -> anchor. Reverse: anchor -> S. aureus. A pair survives only when each is the
    other's top bitscore hit, which is the standard guard against transferring onto a paralogue.
    """
    tgt = D.proteome_fasta(organism)
    fwd = D.run_diamond_blastp(query_fasta, tgt, _cache_dir() / f"hits_sa_to_{organism}.tsv")
    rev = D.run_diamond_blastp(tgt, query_fasta, _cache_dir() / f"hits_{organism}_to_sa.tsv")
    f = D.load_diamond_hits(fwd)
    r = D.load_diamond_hits(rev)
    if f.empty or r.empty:
        return pd.DataFrame(columns=["saureus_acc", "uniprot_accession", "pident", "coverage"])
    f["sseqid"] = f["sseqid"].map(D.acc_from_header)
    r["qseqid"] = r["qseqid"].map(D.acc_from_header)
    best_f = f.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    best_r = r.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    pairs = best_f.merge(best_r, left_on=["qseqid", "sseqid"], right_on=["sseqid", "qseqid"],
                         suffixes=("_f", "_r"))
    pairs = pairs[pairs["pident_f"] >= RBH_MIN_PIDENT]
    return (pairs.rename(columns={"qseqid_f": "saureus_acc", "sseqid_f": "uniprot_accession",
                                  "pident_f": "pident", "qcovhsp_f": "coverage"})
                 [["saureus_acc", "uniprot_accession", "pident", "coverage"]]
                 .reset_index(drop=True))


def _sub(v, thresh: float, direction: str) -> float | None:
    """Effect size -> [0,1] sub-score, saturating two thresholds out. None when unmeasured."""
    if v is None or pd.isna(v):
        return None
    v = float(v)
    x = (-v / -thresh) if direction == "down" else (v / thresh)
    return round(min(1.0, max(0.0, x / 2.0)), 4)


def build() -> pd.DataFrame:
    ab, cl = parse_conlon()
    j = parse_jacques()
    if ab.empty and j.empty:
        raise SystemExit("no activator source tables found — run scripts/10a_fetch_degradability.py "
                         "--only conlon2013_adep4_saureus / jacques2020_onc212_saureus first")
    print(f"[parse] Conlon abundance {len(ab)} proteins | Conlon cleavage {len(cl)} proteins "
          f"| Jacques {len(j)} proteins", flush=True)

    sa = ab.merge(cl, on="sacol", how="outer").merge(j, on="sacol", how="outer")
    print(f"[parse] union on SACOL: {len(sa)} S. aureus proteins", flush=True)

    # Prefer Conlon's YP_ accession; fall back to Jacques' GenBank id for SACOLs only Jacques saw.
    sa["saureus_acc"] = sa["refseq"].fillna("").astype(str).str.strip()
    sa.loc[sa["saureus_acc"].isin(["", "nan"]), "saureus_acc"] = (
        sa["genbank"].fillna("").astype(str).str.strip())
    sa = sa[sa["saureus_acc"] != ""]
    sa["gene"] = sa["gene"].fillna("")
    sa["description"] = sa["description"].fillna(sa["description_onc212"]).fillna("")

    # per-paper sub-scores
    sa["adep4_abundance_sub"] = sa["adep4_abundance_log2fc"].map(
        lambda v: _sub(v, ABUNDANCE_LOG2, "down"))
    sa["adep4_cleavage_sub"] = sa["adep4_cleavage_log2fc_max"].map(
        lambda v: _sub(v, CLEAVAGE_LOG2, "up"))
    sa["onc212_abundance_sub"] = sa["onc212_abundance_log2fc"].map(
        lambda v: _sub(v, ABUNDANCE_LOG2, "down"))
    sa["onc212_cleavage_sub"] = sa["onc212_cleavage_log2fc"].map(
        lambda v: _sub(v, CLEAVAGE_LOG2, "up"))

    # Pool across papers by taking the STRONGER of the two activators per readout: they are two
    # chemically distinct activators of the same enzyme, so evidence from either supports the same
    # mechanistic claim, and neither is a replicate of the other.
    sa["cleavage_sub"] = sa[["adep4_cleavage_sub", "onc212_cleavage_sub"]].max(axis=1)
    sa["abundance_sub"] = sa[["adep4_abundance_sub", "onc212_abundance_sub"]].max(axis=1)
    scored = sa.apply(lambda r: D.renormalised_score(
        {"cleavage": r["cleavage_sub"], "abundance": r["abundance_sub"]},
        {"cleavage": W_CLEAVAGE, "abundance": W_ABUNDANCE}), axis=1)
    sa["activator_evidence"] = [s for s, _ in scored]
    sa["activator_confidence"] = [c for _, c in scored]

    sa["activator_depleted"] = (
        ((sa["adep4_abundance_log2fc"] <= ABUNDANCE_LOG2)
         & (sa["adep4_abundance_padj"] < PADJ))
        | (sa["onc212_abundance_log2fc"] <= ABUNDANCE_LOG2))
    sa["activator_cleaved"] = (
        (sa["adep4_cleavage_peptides_sig"].fillna(0) > 0)
        | (sa["onc212_cleavage_log2fc"] >= CLEAVAGE_LOG2))
    print(f"[score] depleted (>=2x, either activator): {int(sa['activator_depleted'].sum())} | "
          f"cleaved: {int(sa['activator_cleaved'].sum())}", flush=True)
    return sa


def transfer(sa: pd.DataFrame, organism: str, refresh: bool = False) -> pd.DataFrame:
    seqs = fetch_sequences(sorted(sa["saureus_acc"].unique().tolist()), refresh=refresh)
    q = _cache_dir() / "saureus_activator_proteins.faa"
    n_seq = len(seqs)
    m = rbh_map(organism, q)
    print(f"[{organism}] sequences {n_seq}/{sa['saureus_acc'].nunique()} | "
          f"RBH pairs {len(m)}", flush=True)

    out = (m.merge(sa, on="saureus_acc", how="inner")
            .sort_values("activator_evidence", ascending=False)
            .drop_duplicates("uniprot_accession"))
    cols = ["uniprot_accession", "activator_evidence", "activator_confidence",
            "activator_depleted", "activator_cleaved",
            "adep4_abundance_log2fc", "adep4_abundance_padj",
            "adep4_cleavage_log2fc_max", "adep4_cleavage_peptides",
            "adep4_cleavage_peptides_sig",
            "onc212_abundance_log2fc", "onc212_cleavage_log2fc",
            "saureus_acc", "sacol", "saureus_mu50_uniprot", "gene", "description",
            "pident", "coverage"]
    out = out.reindex(columns=cols)
    anchors = set(D.load_accessions(organism))
    out = out[out["uniprot_accession"].isin(anchors)]
    print(f"[{organism}] transferred onto {len(out)} / {len(anchors)} anchor proteins "
          f"({100 * len(out) / max(1, len(anchors)):.1f}% of the proteome); "
          f"with cleavage evidence: {int(out['activator_cleaved'].fillna(False).sum())}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", default="both", choices=["kpneumoniae", "ecoli", "both"])
    ap.add_argument("--refresh", action="store_true", help="re-fetch sequences from NCBI")
    args = ap.parse_args()

    sa = build()
    orgs = ["kpneumoniae", "ecoli"] if args.organism == "both" else [args.organism]
    for org in orgs:
        out = transfer(sa, org, refresh=args.refresh)
        _, prefix = D.ORGANISMS[org]
        dest = D.results_dir(org) / f"{prefix}_clpp_activator.csv"
        out.to_csv(dest, index=False)
        print(f"[{org}] wrote {dest.relative_to(D.REPO_ROOT)}", flush=True)
        top = out.head(12)[["uniprot_accession", "gene", "activator_evidence",
                            "activator_cleaved", "pident"]]
        print(top.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
