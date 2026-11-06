"""Degron motifs + AlphaFold terminal exposure (docs §3.1a / §3.1b).

This is the track the spec is built around — "downweighted when the C-terminus sits in a high-pLDDT
(folded) region rather than an exposed disordered tail" — and the one the legacy prototype never
implemented: its output TSV has no pLDDT column at all.

Two things happen here that did not before.

**1. The motifs are corrected.** The legacy set had defects a weight change cannot fix (see the long
note in `src/degradability.py`): its CM1 regex did not match the ssrA tag it claims to encode, its
N-end-rule feature had the biology inverted (and by itself produced 680 of the 698 `medium` calls),
its CM2 pattern shows no signal at all, and its NM1 fired on 30.7% of the proteome. Every motif here
is asserted against its published archetype before the pipeline runs, and every motif is reported
with a proteome hit rate *and* an enrichment odds ratio against real measured turnover data.

**2. Terminal structure carries the track.** For each protein we read the AlphaFold model from the
existing 06e pdb cache (B-factor == pLDDT) and compute, per terminus, the mean pLDDT and — more
importantly — the length of the contiguous low-pLDDT run, i.e. the candidate pore-initiation region.
ClpXP/ClpAP cannot begin unfolding without an unstructured segment to thread, so this is the
mechanistically load-bearing quantity, and unlike the motifs it is dense (>99% coverage) and
correctly signed against measured half-lives.

`degron_score` therefore weights structural exposure above motifs, and lets a motif contribute *only*
through the exposure of the terminus it sits on.

Validation (printed, and copied into docs/degradability_log.md): every feature is scored against the
Nagar 2021 measured turnover classes — directly for E. coli, and via the 03a sequence orthology for
K. pneumoniae. A feature whose odds ratio is ~1 is reported as such rather than quietly weighted.

Output: output/results/<org>/<prefix>_deg_degrons.csv
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402

TERM_WINDOWS = (10, 20)


def _proteome(organism: str) -> pd.DataFrame:
    df = pd.read_csv(D.proteome_tsv(organism), sep="\t",
                     usecols=["Entry", "Gene Names", "Sequence"])
    df = df.rename(columns={"Entry": "uniprot_accession", "Gene Names": "gene_names",
                            "Sequence": "sequence"})
    df["sequence"] = df["sequence"].fillna("").astype(str).str.upper()
    return df[df["sequence"] != ""].reset_index(drop=True)


def _row_features(acc: str, seq: str, organism: str) -> dict:
    """Motif + structural degron features for one protein."""
    r: dict = {"uniprot_accession": acc, "seq_len": len(seq),
               "cterm_last5": seq[-5:], "nterm_first5": seq[:5]}

    # --- sequence motifs (anchored at the relevant terminus by their own patterns)
    for name in D.MOTIFS:
        r[name] = D.match_motif(name, seq)
    # N-end-rule provenance. BOTH columns are emitted and BOTH are weighted 0 — keeping them makes
    # the legacy sign error auditable rather than merely deleted.
    pos2 = seq[1] if len(seq) >= 2 else ""
    r["nterm_pos2_residue"] = pos2
    r["nend_primary_destabilizing"] = pos2 in D.NEND_DESTABILIZING
    r["nend_imet_cleaved"] = pos2 in D.MAP_CLEAVED

    # --- structural exposure from the AlphaFold model
    plddt = D.plddt_series(organism, acc)
    r["af_available"] = bool(plddt)
    for w in TERM_WINDOWS:
        r[f"cterm_plddt_{w}"] = D.terminal_plddt(plddt, w, "c")
        r[f"nterm_plddt_{w}"] = D.terminal_plddt(plddt, w, "n")
    r["cterm_init_region_len"] = D.init_region_len(plddt, "c")
    r["nterm_init_region_len"] = D.init_region_len(plddt, "n")
    r["cterm_exposure"] = D.exposure(r["cterm_plddt_10"])
    r["nterm_exposure"] = D.exposure(r["nterm_plddt_10"])

    # --- composite degron prior
    # A motif contributes only through the exposure of ITS OWN terminus: a perfect ssrA-like tail
    # buried in a folded domain is not a usable degron. When no model exists we fall back to the
    # unmodulated motif weight and flag it, rather than silently assuming a buried terminus.
    r["degron_structure_missing"] = not bool(plddt)
    comps = []
    for name, weight in D.MOTIF_WEIGHTS.items():
        if weight <= 0 or not r.get(name):
            continue
        exp = r["cterm_exposure"] if name.startswith("cterm") else r["nterm_exposure"]
        comps.append(weight * exp if exp is not None else weight)
    r["degron_motif_score"] = round(max(comps), 4) if comps else 0.0

    exps = [e for e in (r["cterm_exposure"], r["nterm_exposure"]) if e is not None]
    r["degron_exposure_score"] = round(max(exps), 4) if exps else None

    score, conf = D.renormalised_score(
        {"motif": r["degron_motif_score"], "exposure": r["degron_exposure_score"]},
        {"motif": D.W_DEGRON_MOTIF, "exposure": D.W_DEGRON_EXPOSURE},
    )
    r["degron_score"] = score
    r["degron_confidence"] = conf
    r["degron_feature_count"] = int(sum(bool(r.get(m)) for m in D.MOTIFS))

    # Legacy column names the webapp already keys on (app/config.js). Semantics changed with the
    # motif fix, so the help text in app/config.js and docs/degradability_log.md must say so:
    # `degron_nterm` is now a real N-terminal motif (Flynn NM2), not the inverted N-end-rule bool.
    r["degron_cterm"] = bool(r["cterm_cm1_broad"] or r["cterm_cm1_strict"] or r["cterm_cm2"])
    r["degron_nterm"] = bool(r["nterm_nm2"])
    return r


# --------------------------------------------------------------------------- validation
def _nagar_labels_for(organism: str) -> pd.DataFrame:
    """Measured turnover class per accession of `organism` — direct for E. coli, orthology for Kp."""
    nagar = D.load_nagar()
    if nagar.empty:
        return pd.DataFrame(columns=["uniprot_accession", "halflife_class"])
    if organism == "ecoli":
        return nagar[["uniprot_accession", "halflife_class"]]
    orth = D.load_orthologs(organism)
    sub = orth[orth["species"] == D.SPECIES_ECOLI][["anchor_uniprot", "target_uniprot"]]
    ec = nagar[["uniprot_accession", "halflife_class"]].rename(
        columns={"uniprot_accession": "ec_uniprot"})
    m = sub.merge(ec, left_on="target_uniprot", right_on="ec_uniprot", how="inner")
    # One Kp anchor can have several E. coli orthologs; keep the fastest class seen (most
    # conservative for a *validation* set — it maximises the positives a motif could explain).
    rank = {"fast": 0, "intermediate": 1, "stable": 2}
    m["r"] = m["halflife_class"].map(rank)
    m = m.sort_values("r").drop_duplicates("anchor_uniprot")
    return m.rename(columns={"anchor_uniprot": "uniprot_accession"})[
        ["uniprot_accession", "halflife_class"]]


def _odds_ratio(pos_with: int, n_with: int, pos_without: int, n_without: int) -> float | None:
    """Odds ratio of `fast` given the feature, with a Haldane-Anscombe 0.5 correction."""
    a, b = pos_with, n_with - pos_with
    c, d = pos_without, n_without - pos_without
    if min(n_with, n_without) == 0:
        return None
    a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return round((a / b) / (c / d), 2)


def _validate(df: pd.DataFrame, organism: str) -> pd.DataFrame:
    """Score every feature against measured turnover. This is the axis's honesty check."""
    lab = _nagar_labels_for(organism)
    if lab.empty:
        print("  [validate] no Nagar turnover data available — skipping "
              "(run scripts/10a_fetch_degradability.py)", flush=True)
        return pd.DataFrame()
    j = df.merge(lab, on="uniprot_accession", how="inner")
    j["is_fast"] = j["halflife_class"] == "fast"
    n, pos = len(j), int(j["is_fast"].sum())
    base = pos / n if n else 0.0
    print(f"\n  [validate] {n} {organism} proteins with a measured turnover class "
          f"({pos} fast, baseline P(fast)={base:.3f})", flush=True)

    rows = []
    bool_feats = list(D.MOTIFS) + ["nend_primary_destabilizing", "nend_imet_cleaved",
                                   "degron_cterm", "degron_nterm"]
    for f in bool_feats:
        with_f = j[j[f] == True]  # noqa: E712
        without = j[j[f] != True]  # noqa: E712
        hit_rate = float((df[f] == True).mean())  # noqa: E712  proteome-wide prevalence
        p_with = float(with_f["is_fast"].mean()) if len(with_f) else float("nan")
        rows.append({
            "feature": f, "weight": D.MOTIF_WEIGHTS.get(f, 0.0),
            "proteome_hit_rate": round(hit_rate, 4), "n_labelled_hits": len(with_f),
            "p_fast_with": round(p_with, 4) if len(with_f) else None,
            "p_fast_without": round(float(without["is_fast"].mean()), 4) if len(without) else None,
            "odds_ratio": _odds_ratio(int(with_f["is_fast"].sum()), len(with_f),
                                      int(without["is_fast"].sum()), len(without)),
        })
    tab = pd.DataFrame(rows).sort_values("odds_ratio", ascending=False, na_position="last")
    print(tab.to_string(index=False), flush=True)

    # Continuous features: mean per measured class. Correctly signed = fast has the LOWEST pLDDT and
    # the LONGEST initiation region.
    cont = ["cterm_plddt_10", "nterm_plddt_10", "cterm_init_region_len", "nterm_init_region_len",
            "cterm_exposure", "nterm_exposure", "degron_motif_score", "degron_score"]
    means = j.groupby("halflife_class")[cont].mean().round(3)
    order = [c for c in ("fast", "intermediate", "stable") if c in means.index]
    print("\n  [validate] mean by measured turnover class (fast should have LOW pLDDT, "
          "LONG init region, HIGH exposure):", flush=True)
    print(means.loc[order].to_string(), flush=True)

    over = tab[(tab["proteome_hit_rate"] > 0.10) & (tab["weight"] > 0)]
    if not over.empty:
        print(f"\n  [WARN] weighted feature(s) firing on >10% of the proteome — that is a prior, "
              f"not a motif: {', '.join(over['feature'])}", flush=True)
    return tab


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", choices=list(D.ORGANISMS), default="kpneumoniae")
    args = ap.parse_args()
    org = args.organism
    _, prefix = D.ORGANISMS[org]

    # Fail fast if any motif cannot match the published sequence it claims to encode.
    D.selftest_motifs()
    print(f"[{org}] motif archetype self-tests passed ({len(D.MOTIFS)} motifs)", flush=True)

    prot = _proteome(org)
    print(f"[{org}] {len(prot)} proteins from {D.proteome_tsv(org).name}", flush=True)

    rows = [_row_features(a, s, org) for a, s in zip(prot["uniprot_accession"], prot["sequence"])]
    df = pd.DataFrame(rows)
    df = df.merge(prot[["uniprot_accession", "gene_names"]], on="uniprot_accession", how="left")

    n_af = int(df["af_available"].sum())
    print(f"[{org}] AlphaFold models: {n_af}/{len(df)} ({n_af / len(df):.1%}); "
          f"median C-term pLDDT(10) = {df['cterm_plddt_10'].median():.1f}", flush=True)
    for t in ("cterm", "nterm"):
        long_init = int((df[f"{t}_init_region_len"] >= D.MIN_INIT_REGION).sum())
        print(f"[{org}]   {t}: {long_init} proteins ({long_init / len(df):.1%}) with a "
              f">={D.MIN_INIT_REGION}-residue initiation region", flush=True)

    tab = _validate(df, org)

    out_dir = D.results_dir(org)
    out = out_dir / f"{prefix}_deg_degrons.csv"
    df.to_csv(out, index=False)
    print(f"\n[{org}] wrote {out.relative_to(D.REPO_ROOT)}  ({len(df)} rows, {df.shape[1]} cols)",
          flush=True)
    if not tab.empty:
        vout = D.degradability_processed_dir(org) / f"{prefix}_degron_validation.csv"
        tab.to_csv(vout, index=False)
        print(f"[{org}] wrote {vout.relative_to(D.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
