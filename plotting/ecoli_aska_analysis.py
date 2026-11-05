"""E. coli x ASKA: how much does each shortlist criterion cost, and is a cytosolic rule justified?

Phase A of the ASKA-restricted shortlist. This builds NO shortlist: it measures what each criterion
removes, so the gates in `shortlist_ecoli.py` can be chosen with their cost in view. Ligandability,
pockets and selectivity are decorators and do not appear here.

Reads    data/processed/validation/aska_to_uniprot_ecoli.tsv   the ASKA (-) clones mapped onto the proteome
         src loaders: localization, degradability, essentiality, orthology (via `filters.load_joined`)
         data/processed/essentiality/deg_ecoli.tsv             the MEASURED essentiality calls (DEG/Keio)
         data/processed/degradability/evidence/labels_saureus.tsv   the measured S. aureus labels
         data/processed/localization/evidence/deeplocpro_ecoli.tsv  DeepLocPro's confidence
Writes   output/results/shortlist/ecoli_aska_*.tsv             one table per question

**Definitions, all stated here so they are arguable.** `cytosolic` = DeepLocPro class `cytoplasm` OR
TMbed `cytoplasmic_fraction == 1`: the two predictors are independent opinions and DeepLocPro labels
many E. coli ribosomal proteins `outer_membrane`. `degradable` = a probability at or above the
re-derived base-rate cut (`D.BASE_RATE_THRESHOLD`), per activator. `essential` is shown several ways
on purpose (measured in any / all DEG screens, `essentiality_consensus` at the base-rate cut) because
the choice is the open question. **E. coli degradability is extrapolated from S. aureus labels and
the E. coli essentiality consensus is partly circular (built from the same screens as the measured
calls).**

Run with the `gradi` env:  python plotting/ecoli_aska_analysis.py [-q]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from plotting import filters as F  # noqa: E402
from src import degradability as D  # noqa: E402
from src import essentiality as E  # noqa: E402

SPECIES = "ecoli"
OUT = REPO_ROOT / "output" / "results" / "shortlist"
ASKA = REPO_ROOT / "data/processed/validation/aska_to_uniprot_ecoli.tsv"
LABELS = REPO_ROOT / "data/processed/degradability/evidence/labels_saureus.tsv"
DLP = REPO_ROOT / "data/processed/localization/evidence/deeplocpro_ecoli.tsv"
FRACTION_BINS = [(0, 0.2), (0.2, 0.5), (0.5, 0.8), (0.8, 1.0), (1.0, 1.0001)]
RIBOSOMAL = r"^(?:rps|rpl|rpm)[A-Z]"
QUIET = False


def say(msg: str = "") -> None:
    if not QUIET:
        print(msg, flush=True)


def emit(name: str, df: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / f"ecoli_aska_{name}.tsv", sep="\t", index=False)
    say(f"\n== {name}  ({len(df)} rows -> output/results/shortlist/ecoli_aska_{name}.tsv)")
    say(df.to_string(index=False))


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def rate_row(g: pd.DataFrame, **keys) -> dict:
    row = dict(keys)
    for act in ("adep4", "onc212"):
        y = g[act].dropna()
        lo, hi = wilson(int(y.sum()), len(y))
        row.update({f"{act}_n": len(y), f"{act}_rate": round(float(y.mean()), 3) if len(y) else np.nan,
                    f"{act}_ci": f"{lo:.3f}-{hi:.3f}" if len(y) else ""})
    return row


def load_ecoli() -> pd.DataFrame:
    """The whole proteome, joined, with every flag the analysis uses. ASKA is a flag, not a filter."""
    d = F.load_joined(SPECIES)
    n0 = len(d)
    deg = E.load_deg(SPECIES)[["uniprot_ac", "deg_essential_any", "deg_essential_all", "deg_n_datasets"]]
    d = d.merge(deg, on="uniprot_ac", how="left")
    d = d.merge(E.load(SPECIES)[["uniprot_ac", "essentiality_consensus", "essentiality_evidence"]], on="uniprot_ac")
    dg = D.load(SPECIES)[["uniprot_ac", "degradability_consensus", "adep4_ortholog_hit", "onc212_ortholog_hit"]]
    d = d.merge(dg, on="uniprot_ac")
    d = d.merge(pd.read_csv(DLP, sep="\t").rename(columns={"confidence": "dlp_confidence"})[["uniprot_ac", "dlp_confidence"]],
                on="uniprot_ac", how="left")
    assert len(d) == n0, "a join changed the row count"
    # ASKA: one row per PROTEIN. 34 proteins have two clones; keep the best-confidence clone.
    a = pd.read_csv(ASKA, sep="\t").sort_values(["uniprot_ac", "confidence"], ascending=[True, False])
    n_clones = a.groupby("uniprot_ac").size().rename("n_clones")
    a = a.drop_duplicates("uniprot_ac").rename(columns={"confidence": "aska_confidence"}).merge(n_clones, on="uniprot_ac")
    d = d.merge(a, on="uniprot_ac", how="left")
    d["in_aska"] = d["aska_id"].notna()
    d["length"] = pd.read_csv(REPO_ROOT / f"data/processed/proteomes/proteome_{SPECIES}.tsv", sep="\t")["sequence"].str.len().values
    d["ribosomal"] = d["gene_name"].astype(str).str.match(RIBOSOMAL) | d["protein_name"].astype(str).str.contains("ribosomal protein", case=False)
    d["cls"] = d["localization"] == "cytoplasm"
    d["f1"] = d["cytoplasmic_fraction"] == 1
    d["cytosolic"] = d["cls"] | d["f1"]
    d["cyto_route"] = np.select([d.cls & d.f1, d.f1, d.cls], ["both", "tmbed_only", "deeploc_only"], "none")
    d["ess_any"] = d["deg_essential_any"] == 1
    d["ess_all"] = d["deg_essential_all"] == 1
    d["no_deg_call"] = d["deg_essential_any"].isna()
    base = float(d["ess_any"].mean())
    d.attrs["base_rate"] = base
    d.attrs["cons_cut"] = float(d["essentiality_consensus"].quantile(1 - base))
    d["ess_cons"] = d["essentiality_consensus"] >= d.attrs["cons_cut"]
    ca, co = D.BASE_RATE_THRESHOLD["adep4"], D.BASE_RATE_THRESHOLD["onc212"]
    d["deg_a"] = d["adep4_prob"] >= ca
    d["deg_o"] = d["onc212_prob"] >= co
    d["deg_any"] = d.deg_a | d.deg_o
    d["deg_both"] = d.deg_a & d.deg_o
    # The consensus has no natural threshold (a mean of percentiles). CUT is DERIVED: the quantile at
    # which the consensus passes the same share of the proteome as "either activator above its
    # base-rate cut", so the two routes are comparable. Never hard-coded.
    d.attrs["pass_rate_either"] = float(d["deg_any"].mean())
    d.attrs["cons_deg_cut"] = float(d["degradability_consensus"].quantile(1 - d.attrs["pass_rate_either"]))
    d["cons_hi"] = d["degradability_consensus"] >= d.attrs["cons_deg_cut"]
    # Ortholog route: a measured S. aureus ortholog is a degrader in EITHER screen. NA = nothing to
    # say, so it is never a 0 and never a veto.
    sup = D.saureus_ortholog_support(SPECIES).reset_index()[["uniprot_ac", "n_orthologs", "adep4_support", "onc212_support"]]
    attrs = dict(d.attrs)  # merge() returns a new frame and drops .attrs
    d = d.merge(sup, on="uniprot_ac", how="left")
    d.attrs.update(attrs)
    d["orth_hit"] = (d["adep4_ortholog_hit"] == 1).fillna(False).astype(bool) | (d["onc212_ortholog_hit"] == 1).fillna(False).astype(bool)
    d["orth_measured"] = d["adep4_ortholog_hit"].notna() | d["onc212_ortholog_hit"].notna()
    d["orth_contradicts"] = d["orth_measured"] & ~d["orth_hit"]
    # The ortholog route is the weaker one (its label was training data, and the 90 proteins it alone
    # admitted scored like a random draw), so it needs a FLOOR on the consensus: the proteome median,
    # "at least an average score". Derived each run, never hard-coded. Measured on the E. coli
    # proteins with a measured ortholog, the ortholog-degrader rate is 0.04 below consensus 0.3,
    # 0.16-0.19 from 0.3 to 0.6, 0.29 at 0.6-0.7 and 0.79 above 0.9.
    d.attrs["orth_floor"] = float(d["degradability_consensus"].median())
    d["orth_ok"] = d["orth_hit"] & (d["degradability_consensus"] >= d.attrs["orth_floor"])
    d["orth_dropped_by_floor"] = d["orth_hit"] & ~d["orth_ok"] & ~d["cons_hi"]
    d["degradable"] = d["cons_hi"] | d["orth_ok"]
    d["deg_route"] = np.select([d.cons_hi & d.orth_ok, d.cons_hi, d.orth_ok], ["both", "consensus_only", "ortholog_only"], "none")
    assert len(d) == n0
    return d


# ---------------------------------------------------------------- 1. cytosolic
def cytosolic(d: pd.DataFrame) -> None:
    say("\n" + "=" * 92 + "\n1. CYTOSOLIC: does it make sense, and how should it be defined?\n" + "=" * 92)
    lab = pd.read_csv(LABELS, sep="\t")
    from src import localization as Lc

    sa = lab.merge(Lc.load("saureus"), on="uniprot_ac")
    rows = []
    for c, g in sa.groupby("localization"):
        rows.append(rate_row(g, block="S. aureus by compartment", group=c))
    for lo, hi in FRACTION_BINS:
        g = sa[(sa.cytoplasmic_fraction >= lo) & (sa.cytoplasmic_fraction < hi)]
        rows.append(rate_row(g, block="S. aureus by cytoplasmic_fraction", group=f"[{lo:g},{min(hi, 1):g}{']' if hi > 1 else ')'}"))
    sa["cls"], sa["f1"] = sa.localization == "cytoplasm", sa.cytoplasmic_fraction == 1
    for name, m in (("class AND fraction==1", sa.cls & sa.f1), ("class only", sa.cls & ~sa.f1),
                    ("fraction==1 only", ~sa.cls & sa.f1), ("neither", ~sa.cls & ~sa.f1),
                    ("UNION (class OR fraction==1)", sa.cls | sa.f1)):
        rows.append(rate_row(sa[m], block="S. aureus, the two predictors", group=name))
    emit("cytosolic_saureus_rates", pd.DataFrame(rows))

    k = d[d.in_aska]
    rows = []
    for c in (True, False):
        for f in (True, False):
            g = k[(k.cls == c) & (k.f1 == f)]
            rows.append({"deeploc_cytoplasm": c, "tmbed_fraction_1": f, "n": len(g), "essential_any": int(g.ess_any.sum()),
                         "degradable_either": int(g.deg_any.sum()), "essential_and_degradable": int((g.ess_any & g.deg_any).sum()),
                         "ribosomal": int(g.ribosomal.sum()), "mean_adep4": round(float(g.adep4_prob.mean()), 3)})
    emit("cytosolic_ecoli_agreement", pd.DataFrame(rows))

    variants = {"no gate": pd.Series(True, index=k.index), "class == cytoplasm": k.cls, "fraction == 1": k.f1,
                "class AND fraction == 1": k.cls & k.f1, "class OR fraction == 1  (proposed)": k.cytosolic,
                "class OR fraction >= 0.8": k.cls | (k.cytoplasmic_fraction >= 0.8), "not membrane": k.localization != "cytoplasmic_membrane"}
    rows = [{"gate": n, "passes": int(m.sum()), "pct_of_aska": round(100 * m.mean(), 1), "essential_any": int((m & k.ess_any).sum()),
             "essential_and_degradable": int((m & k.ess_any & k.deg_any).sum()),
             "lost_essential_and_degradable": int((k.ess_any & k.deg_any).sum() - (m & k.ess_any & k.deg_any).sum())}
            for n, m in variants.items()]
    emit("cytosolic_gate_variants", pd.DataFrame(rows))

    r = k[(k.cyto_route.isin(["tmbed_only", "deeploc_only"])) & (k.ess_any | k.deg_any)].copy()
    r = r[["gene_name", "cyto_route", "localization", "dlp_confidence", "cytoplasmic_fraction", "ribosomal", "ess_any",
           "adep4_prob", "onc212_prob"]].sort_values(["cyto_route", "ribosomal"], ascending=[True, False])
    say(f"\n  the OR admits {int((k.cyto_route == 'tmbed_only').sum())} proteins DeepLocPro does not call cytoplasm; "
        f"{int(((k.cyto_route == 'tmbed_only') & k.ribosomal).sum())} are ribosomal; their DeepLocPro classes: "
        f"{k[k.cyto_route == 'tmbed_only'].localization.value_counts().to_dict()}")
    emit("cytosolic_or_discordant", r[(r.ess_any) & (r.adep4_prob >= D.BASE_RATE_THRESHOLD['adep4']) | (r.ess_any) & (r.onc212_prob >= D.BASE_RATE_THRESHOLD['onc212'])])

    side = k[(~k.cytosolic) & (k.cytoplasmic_fraction >= 0.5) & k.ess_any & k.deg_any]
    say(f"\n  membrane-anchored side tier (fraction 0.5-<1, class not cytoplasm), essential and degradable: "
        f"{', '.join(side.gene_name.astype(str))}")
    final = k[k.cytosolic & k.ess_any & k.deg_any]
    say(f"  ribosomal share of the cytosolic essential-and-degradable set: {int(final.ribosomal.sum())} of {len(final)}")

    rows = []
    for c, g in k.groupby("localization"):
        rows.append({"localization": c, "n": len(g), "adep4_above_cut": round(float(g.deg_a.mean()), 3),
                     "onc212_above_cut": round(float(g.deg_o.mean()), 3), "essential_any": round(float((g.ess_any).mean()), 3),
                     "essential_and_degradable": int((g.ess_any & g.deg_any).sum())})
    emit("cytosolic_ecoli_by_compartment", pd.DataFrame(rows))


# ---------------------------------------------------------------- 2. dropout
def dropout(d: pd.DataFrame) -> None:
    say("\n" + "=" * 92 + "\n0. DROPOUT: what each criterion removes\n" + "=" * 92)
    cut = d.attrs["cons_cut"]
    crit = {
        "ASKA": d.in_aska,
        "cytosolic (class OR fraction==1)": d.cytosolic,
        "cytosolic: class == cytoplasm": d.cls,
        "cytosolic: fraction == 1": d.f1,
        "essential: measured any": d.ess_any,
        "essential: measured all": d.ess_all,
        f"essential: consensus >= {cut:.3f} (top {100 * d.attrs['base_rate']:.1f}%)": d.ess_cons,
        "essential: measured any OR consensus": d.ess_any | d.ess_cons,
        "essential: measured any AND consensus": d.ess_any & d.ess_cons,
        "degradable: ADEP4": d.deg_a, "degradable: ONC212": d.deg_o,
        "degradable: either": d.deg_any, "degradable: both": d.deg_both,
    }
    ak = d.in_aska
    rows = [{"criterion": n, "passes_proteome": int(m.sum()), "pct_proteome": round(100 * m.mean(), 1),
             "passes_within_aska": int((m & ak).sum()), "pct_of_aska": round(100 * (m & ak).sum() / ak.sum(), 1)}
            for n, m in crit.items()]
    emit("dropout_standalone", pd.DataFrame(rows))

    ess_defs = {"measured any": d.ess_any, "measured all": d.ess_all, "consensus cut": d.ess_cons,
                "measured any OR consensus": d.ess_any | d.ess_cons, "measured any AND consensus": d.ess_any & d.ess_cons}
    deg_defs = {"either": d.deg_any, "ADEP4": d.deg_a, "both": d.deg_both}
    rows = []
    for en, em in ess_defs.items():
        for dn, dm in deg_defs.items():
            cur, row = d.in_aska.copy(), {"essential": en, "degradable": dn, "aska": int(d.in_aska.sum())}
            for lab, m in (("cytosolic", d.cytosolic), ("essential", em), ("degradable", dm)):
                cur &= m
                row[f"after_{lab}"] = int(cur.sum())
            parts = {"aska": d.in_aska, "cytosolic": d.cytosolic, "essential": em, "degradable": dm}
            row["final"] = int(cur.sum())
            for k_ in parts:
                oth = np.logical_and.reduce([v for kk, v in parts.items() if kk != k_])
                row[f"marginal_dropout_pct_{k_}"] = round(100 * (1 - cur.sum() / oth.sum()), 1)
            rows.append(row)
    emit("dropout_funnels", pd.DataFrame(rows))


# ---------------------------------------------------------------- 3. degradability
def degradability(d: pd.DataFrame) -> None:
    say("\n" + "=" * 92 + "\n2. DEGRADABILITY considerations (E. coli rows are extrapolated from S. aureus labels)\n" + "=" * 92)
    k = d[d.in_aska]
    rho_len = k["adep4_prob"].corr(k["length"], method="spearman")
    rho_ao = k["adep4_prob"].corr(k["onc212_prob"], method="spearman")
    jac = (k.deg_a & k.deg_o).sum() / (k.deg_a | k.deg_o).sum()
    say(f"  ASKA: rho(ADEP4, ONC212) {rho_ao:.3f}; Jaccard of the two passing sets {jac:.3f}; "
        f"rho(length, ADEP4) {rho_len:.3f}; rho(length, ONC212) {k['onc212_prob'].corr(k['length'], method='spearman'):.3f}")
    rows = []
    for name, m in (("ADEP4", k.deg_a), ("ONC212", k.deg_o), ("either", k.deg_any), ("both", k.deg_both)):
        g = k[m]
        rows.append({"rule": name, "passes": int(m.sum()), "median_length": int(g.length.median()),
                     "pct_ribosomal": round(100 * g.ribosomal.mean(), 1), "pct_cytosolic": round(100 * g.cytosolic.mean(), 1),
                     "pct_essential_any": round(100 * g.ess_any.mean(), 1)})
    emit("degradability_rules", pd.DataFrame(rows))

    rows = []
    q = pd.qcut(k["length"], 4, labels=["Q1 shortest", "Q2", "Q3", "Q4 longest"])
    for lab, g in k.groupby(q, observed=True):
        rows.append({"length_quartile": lab, "median_length": int(g.length.median()), "n": len(g),
                     "adep4_above_cut": round(float(g.deg_a.mean()), 3), "onc212_above_cut": round(float(g.deg_o.mean()), 3),
                     "essential_any": round(float(g.ess_any.mean()), 3)})
    emit("degradability_length", pd.DataFrame(rows))

    fin = k[k.cytosolic & k.ess_any & k.deg_any]
    rows = []
    for act in ("adep4", "onc212"):
        h = fin[f"{act}_ortholog_hit"]
        rows.append({"activator": act, "candidates": len(fin), "with_measured_saureus_ortholog": int(h.notna().sum()),
                     "ortholog_is_degrader": int((h == 1).sum()), "ortholog_not_degrader": int((h == 0).sum()),
                     "no_ortholog_or_unmeasured": int(h.isna().sum())})
    emit("degradability_ortholog_support", pd.DataFrame(rows))


def degradability_gate(d: pd.DataFrame) -> None:
    say("\n" + "=" * 92 + "\n2b. THE DEGRADABILITY GATE:  consensus >= CUT  OR  a measured S. aureus ortholog is a degrader\n" + "=" * 92)
    cut, rate = d.attrs["cons_deg_cut"], d.attrs["pass_rate_either"]
    say(f"  CUT = consensus.quantile(1 - {rate:.4f}) = {cut:.4f}   (the pass rate of 'either activator above its base-rate cut' "
        f"{D.BASE_RATE_THRESHOLD['adep4']} / {D.BASE_RATE_THRESHOLD['onc212']}, applied to the whole proteome)")
    g = d[d.in_aska & d.cytosolic]
    say(f"  ORTHOLOG FLOOR = proteome median consensus = {d.attrs['orth_floor']:.4f}: an ortholog degrader enters only at or above it")
    say(f"  gate set (ASKA and cytosolic): {len(g)}")
    say(f"  consensus >= CUT {int(g.cons_hi.sum())} | probability rule (either) {int(g.deg_any.sum())} | both {int((g.cons_hi & g.deg_any).sum())} | "
        f"consensus only {int((g.cons_hi & ~g.deg_any).sum())} | probability only {int((~g.cons_hi & g.deg_any).sum())}")
    say(f"  S. aureus ortholog: has one {int((g.n_orthologs > 0).sum())} | measured {int(g.orth_measured.sum())} | degrader (either) "
        f"{int(g.orth_hit.sum())} (ADEP4 {int((g.adep4_ortholog_hit == 1).sum())}, ONC212 {int((g.onc212_ortholog_hit == 1).sum())}) | "
        f"measured non-degrader {int(g.orth_contradicts.sum())}")
    y = g[g.orth_measured]
    say(f"  AUROC consensus vs ortholog hit, among the {len(y)} with a measured ortholog: "
        f"{roc_auc_score(y.orth_hit, y.degradability_consensus):.3f}  [consistency, not validation: the ortholog label was training data]")
    say(f"  P(consensus high | ortholog degrader) {g[g.orth_hit].cons_hi.mean():.2f} | P(consensus high | measured non-degrader) "
        f"{g[g.orth_contradicts].cons_hi.mean():.2f}")
    tab = pd.crosstab(g.orth_hit.map({True: "ortholog degrader", False: "not"}), g.cons_hi.map({True: "consensus high", False: "consensus low"}))
    say(tab.to_string())
    rows = []
    ess = {"no essential gate": pd.Series(True, index=g.index), "measured any": g.ess_any, "measured all": g.ess_all,
           "measured any OR consensus": g.ess_any | g.ess_cons}
    for en, em in ess.items():
        row = {"essential": en, "essential_n": int(em.sum())}
        for dn, dm in (("consensus_high", g.cons_hi), ("ortholog_degrader", g.orth_hit), ("consensus_OR_ortholog_above_floor (gate)", g.degradable),
                       ("consensus_AND_ortholog", g.cons_hi & g.orth_hit), ("consensus_high_no_contradiction", g.cons_hi & ~g.orth_contradicts)):
            row[dn] = int((em & dm).sum())
        row["rescued_by_ortholog"] = int((em & g.orth_ok & ~g.cons_hi).sum())
        row["ortholog_dropped_by_floor"] = int((em & g.orth_dropped_by_floor).sum())
        rows.append(row)
    emit("degradability_gate", pd.DataFrame(rows))
    r = g[g.degradable]
    say(f"\n  cautions: median length gate set {int(g.length.median())} aa | consensus-high {int(g[g.cons_hi].length.median())} | "
        f"ortholog degraders {int(g[g.orth_hit].length.median())}; ribosomal: ortholog degraders {int(g[g.orth_hit].ribosomal.sum())} of "
        f"{int(g.orth_hit.sum())}, consensus-high {int(g[g.cons_hi].ribosomal.sum())} of {int(g.cons_hi.sum())}")
    say(f"  ortholog degraders below the consensus cut: {int((g.orth_hit & ~g.cons_hi).sum())}; kept by the floor {int((g.orth_ok & ~g.cons_hi).sum())}, "
        f"dropped {int(g.orth_dropped_by_floor.sum())}")
    say(f"  deg_route among {len(r)} passing: {r.deg_route.value_counts().to_dict()}; with a contradicting measured ortholog: {int(r.orth_contradicts.sum())}")


# ---------------------------------------------------------------- 4. essentiality
def essentiality(d: pd.DataFrame) -> None:
    say("\n" + "=" * 92 + "\n3. ESSENTIALITY: measured vs consensus (consensus is partly circular on E. coli)\n" + "=" * 92)
    k = d[d.in_aska]
    called = k[~k.no_deg_call]
    say(f"  ASKA: measured any {int(k.ess_any.sum())}, all {int(k.ess_all.sum())}, no DEG call {int(k.no_deg_call.sum())}; "
        f"consensus cut {d.attrs['cons_cut']:.4f} (proteome base rate {d.attrs['base_rate']:.3f})")
    for lab, y in (("any", called.ess_any), ("all", called.ess_all)):
        say(f"  AUROC of essentiality_consensus vs measured {lab} (called, n={len(called)}): "
            f"{roc_auc_score(y, called.essentiality_consensus):.3f}  [circular by construction]")
    tiers = {
        "T1 measured in ALL screens": k.ess_all,
        "T2 measured any, not all, consensus >= cut": k.ess_any & ~k.ess_all & k.ess_cons,
        "T3 measured any, not all, consensus < cut": k.ess_any & ~k.ess_all & ~k.ess_cons,
        "T4 no DEG call, consensus >= cut": k.no_deg_call & k.ess_cons,
        "T5 DEG says NOT essential, consensus >= cut": (k.deg_essential_any == 0) & k.ess_cons,
        "T6 none of the above": ~(k.ess_any | (k.no_deg_call & k.ess_cons) | ((k.deg_essential_any == 0) & k.ess_cons)),
    }
    rows = [{"tier": n, "aska": int(m.sum()), "cytosolic": int((m & k.cytosolic).sum()),
             "cytosolic_and_degradable": int((m & k.cytosolic & k.deg_any).sum()),
             "cytosolic_and_degradable_adep4": int((m & k.cytosolic & k.deg_a).sum()),
             "mean_consensus": round(float(k.loc[m, "essentiality_consensus"].mean()), 3)} for n, m in tiers.items()]
    emit("essentiality_tiers", pd.DataFrame(rows))
    ev = pd.crosstab(k.essentiality_evidence, [k.ess_any.map({True: "measured essential", False: "not measured essential"})])
    say("\n  essentiality_evidence (1-3) by measured status, ASKA:")
    say(ev.to_string())
    ev.reset_index().to_csv(OUT / "ecoli_aska_essentiality_evidence.tsv", sep="\t", index=False)
    disc = k[(k.ess_any & ~k.ess_cons)][["gene_name", "localization", "essentiality_consensus", "deg_n_datasets", "adep4_prob"]]
    say(f"\n  measured essential but consensus below the cut: {len(disc)} (top by lowest consensus below)")
    say(disc.sort_values("essentiality_consensus").head(12).to_string(index=False))


# ---------------------------------------------------------------- 5. ASKA itself
def aska(d: pd.DataFrame) -> None:
    say("\n" + "=" * 92 + "\n4. ASKA itself: what the library does not cover\n" + "=" * 92)
    miss = d[~d.in_aska]
    rows = [{"item": "proteome", "n": len(d)}, {"item": "in ASKA (distinct proteins)", "n": int(d.in_aska.sum())},
            {"item": "NOT in ASKA", "n": len(miss)},
            {"item": "  of which essential (measured any)", "n": int(miss.ess_any.sum())},
            {"item": "  of which cytosolic and essential and degradable", "n": int((miss.cytosolic & miss.ess_any & miss.deg_any).sum())},
            {"item": "ASKA confidence 3 / 2 / 1", "n": " / ".join(str(int((d.aska_confidence == c).sum())) for c in (3, 2, 1))},
            {"item": "proteins with two clones", "n": int((d.n_clones == 2).sum())}]
    emit("aska_coverage", pd.DataFrame(rows))
    g = miss[miss.ess_any]
    say(f"  essential genes with NO clone ({len(g)}): {', '.join(g.gene_name.astype(str).head(25))}")
    fin = d[d.in_aska & d.cytosolic & d.ess_any & d.deg_any]
    say(f"  confidence-1 clones among the cytosolic essential-and-degradable set: {int((fin.aska_confidence == 1).sum())} of {len(fin)}")


def main() -> int:
    global QUIET
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    QUIET = ap.parse_args().quiet
    d = load_ecoli()
    say(f"ecoli_aska_analysis.py  proteome {len(d)}, ASKA {int(d.in_aska.sum())} distinct proteins "
        f"(clone-table rows 4,110); base rate {d.attrs['base_rate']:.3f}, consensus cut {d.attrs['cons_cut']:.4f}")
    dropout(d)
    cytosolic(d)
    degradability(d)
    degradability_gate(d)
    essentiality(d)
    aska(d)
    return 0


if __name__ == "__main__":
    sys.exit(main())
