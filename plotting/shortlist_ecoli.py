"""E. coli K-12 shortlist restricted to the ASKA library: four stated gates, then decorators. NO composite score.

Replaces the earlier Kp-style cascade. Each gate is a single predicate with a written rationale, applied
in order, and every step reports what it removed; `leave_one_out` prices each gate by dropping it.

- `aska`        the protein has a clone in the ASKA (-) library (4,076 distinct proteins; 34 have two clones)
- `cytosolic`   DeepLocPro class `cytoplasm` OR TMbed `cytoplasmic_fraction == 1` (two independent opinions;
                DeepLocPro labels many ribosomal proteins `outer_membrane`)
- `degradable`  `degradability_consensus` >= CUT, OR a measured S. aureus ortholog is a degrader in either
                screen AND the consensus >= the proteome median. CUT is DERIVED (the quantile at which the
                consensus passes the share of the proteome that "either activator above its base-rate cut"
                passes), the floor is the proteome median; both are computed each run, never hard-coded
- `essential`   MEASURED in any DEG/Keio screen (`deg_essential_any == 1`) AND `essentiality_consensus` >=
                its base-rate cut (the quantile at 1 - the measured base rate, top 15.8%). Swappable with
                `--essential` (any, all, any_and_consensus, any_or_consensus, consensus)

Decorators are shown and never gated or scored: selectivity (`has_human_ortholog`), ligandability, pockets,
studiedness, the degradability route and probabilities, ASKA clone ids. Survivors are ordered by
`degradability_consensus` as a STATED TIEBREAK, which is not a score.

**Read it for what it is.** E. coli degradability is EXTRAPOLATED from S. aureus labels (ranking
hypotheses, not measurements), and the E. coli essentiality consensus is built from the same screens as the
measured calls, so `any AND consensus` reads "a screen says essential and the screens' own models do not
contradict it", not independent confirmation. Only about 55% of measured-any essentials clear the consensus
cut; the ones that do not are written to a second table, never dropped silently.

Writes   output/results/shortlist/shortlist_ecoli_aska.tsv                  the survivors
         output/results/shortlist/shortlist_ecoli_aska_measured_below_cut.tsv  passes the other gates, measured
                                                                          essential, consensus below the cut
         output/results/shortlist/shortlist_ecoli_aska_gates.tsv             funnel + leave-one-out
         output/results/shortlist/shortlist_ecoli_aska_top200.tsv            THE LIST: absolute top 50 by degradability
                                                                          among the first three gates + the essential-gate
                                                                          survivors by degradability, to 200; `list_source`
         output/results/shortlist/shortlist_ecoli_aska_top200_overflow.tsv   essential-gate proteins that did not fit
         output/results/shortlist/shortlist_ecoli_aska_control50.tsv         random control: 50 cytosolic ASKA proteins below
                                                                          the median degradability, not on the list, measured
                                                                          essential / not in the list's proportion

Run with the `gradi` env:  python plotting/shortlist_ecoli.py [--essential any_and_consensus] [-q]
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402

from plotting import filters as F  # noqa: E402

SPECIES = "ecoli"
OUT = REPO_ROOT / "output" / "results" / "shortlist"
ESSENTIAL_DEFINITIONS = {
    "any": "measured essential in any DEG screen",
    "all": "measured essential in every DEG screen that covers it",
    "any_and_consensus": "measured essential in any DEG screen AND essentiality_consensus >= base-rate cut",
    "any_or_consensus": "measured essential in any DEG screen OR essentiality_consensus >= base-rate cut",
    "consensus": "essentiality_consensus >= base-rate cut",
}
DECORATORS = [
    ("proteomes", ["gene_name", "protein_name"]),
    ("selectivity", ["has_human_ortholog", "bacterial_panel_orthologs"]),
    ("ligands", ["ligands_consensus", "ligands_evidence", "n_ligands_own", "n_ligands_bacterial", "n_drugs"]),
    ("pockets", ["pockets_evidence", "n_ligands_pdb", "n_ligands_alphafill", "n_pdb_structures", "p2rank_score"]),
]
COLUMNS = ["uniprot_ac", "gene_name", "protein_name", "aska_id", "aska_confidence", "n_clones",
           "localization", "cytoplasmic_fraction", "cyto_route",
           "degradability_consensus", "adep4_prob", "onc212_prob", "deg_route", "ortholog_support",
           "essentiality_consensus", "essentiality_evidence", "deg_essential_any", "deg_essential_all", "deg_n_datasets",
           "has_human_ortholog", "bacterial_panel_orthologs",
           "ligands_consensus", "ligands_evidence", "n_ligands_own", "n_ligands_bacterial", "n_drugs",
           "pockets_evidence", "n_ligands_pdb", "n_ligands_alphafill", "n_pdb_structures", "p2rank_score",
           "n_papers_uniprot_prokaryotic", "length", "ribosomal", "on_panel"]
QUIET = False


def say(msg: str = "") -> None:
    if not QUIET:
        print(msg, flush=True)


def load() -> pd.DataFrame:
    spec = importlib.util.spec_from_file_location("ecoli_aska_analysis", REPO_ROOT / "plotting/ecoli_aska_analysis.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    d = m.load_ecoli()
    attrs = dict(d.attrs)  # merge() returns a new frame and drops .attrs
    from src import interest, ligandability, orthology, pockets

    for mod, cols in ((orthology, ["has_human_ortholog", "bacterial_panel_orthologs"]),
                      (ligandability, ["ligands_consensus", "ligands_evidence", "n_ligands_own", "n_ligands_bacterial", "n_drugs"]),
                      (pockets, ["pockets_evidence", "n_ligands_pdb", "n_ligands_alphafill", "n_pdb_structures", "p2rank_score"])):
        part = mod.load(SPECIES)[["uniprot_ac", *[c for c in cols if c not in d.columns]]]
        d = d.merge(part, on="uniprot_ac", how="left", validate="1:1")
    d.attrs.update(attrs)
    d["on_panel"] = d["gene_name"].isin(set(interest.genes()))
    d["ortholog_support"] = pd.Series(
        pd.NA, index=d.index, dtype="object").mask(d.n_orthologs.fillna(0) == 0, "no ortholog")
    meas = d.n_orthologs.fillna(0) > 0
    d.loc[meas & ~d.orth_measured, "ortholog_support"] = "unmeasured"
    d.loc[meas & d.orth_measured & d.orth_hit, "ortholog_support"] = "degrader"
    d.loc[meas & d.orth_measured & ~d.orth_hit, "ortholog_support"] = "not a degrader"
    return d


def build_gates(d: pd.DataFrame, essential: str) -> list[F.Filter]:
    a = d.attrs
    ess = {"any": d.ess_any, "all": d.ess_all, "any_and_consensus": d.ess_any & d.ess_cons,
           "any_or_consensus": d.ess_any | d.ess_cons, "consensus": d.ess_cons}[essential]
    return [
        F.Filter("aska", "aska_id", "A clone exists, so the protein can be tested in the lab.",
                 lambda x: x.in_aska, "in the ASKA (-) library"),
        F.Filter("cytosolic", "cytosolic",
                 "Activated ClpP meets soluble cytoplasmic chains; the OR of two independent predictors recovers "
                 "ribosomal proteins DeepLocPro mislabels.",
                 lambda x: x.cytosolic, "DeepLocPro cytoplasm OR TMbed cytoplasmic_fraction == 1"),
        F.Filter("degradable", "degradability_consensus",
                 "Extrapolated from S. aureus labels; the ortholog route needs a measured S. aureus degrader and an "
                 "at-least-average consensus.",
                 lambda x: x.degradable,
                 f"consensus >= {a['cons_deg_cut']:.4f} OR (measured S. aureus ortholog degrader AND consensus >= {a['orth_floor']:.4f})"),
        F.Filter("essential", "deg_essential_any",
                 "Measured in rich-medium screens, optionally corroborated by the consensus (partly circular on E. coli).",
                 lambda x: ess.reindex(x.index), f"{ESSENTIAL_DEFINITIONS[essential]} (cut {a['cons_cut']:.4f}, base rate {a['base_rate']:.3f})"),
    ]


def run(d: pd.DataFrame, gates: list[F.Filter], skip: str | None = None):
    cur, steps = d, []
    for g in gates:
        if g.name == skip:
            continue
        before = len(cur)
        cur = cur[g.apply(cur)]
        steps.append(F.Step(g.name, g.rule, g.rationale, before, len(cur)))
    return cur, steps


def main() -> int:
    global QUIET
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--essential", choices=list(ESSENTIAL_DEFINITIONS), default="any_and_consensus")
    ap.add_argument("--top", type=int, default=50, help="absolute top-N by degradability_consensus among the proteins passing the first three gates")
    ap.add_argument("--n", type=int, default=200, help="size of the final list")
    ap.add_argument("--control", type=int, default=50, help="size of the random below-median control")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    QUIET = a.quiet

    say("shortlist_ecoli.py  ASKA -> cytosolic -> degradable -> essential, decorators, no composite")
    say(f"  essential rule: {a.essential} ({ESSENTIAL_DEFINITIONS[a.essential]})")
    d = load()
    say(f"  proteome {len(d):,} | CUT {d.attrs['cons_deg_cut']:.4f} | ortholog floor {d.attrs['orth_floor']:.4f} | "
        f"essentiality cut {d.attrs['cons_cut']:.4f} (base rate {d.attrs['base_rate']:.3f})")
    gates = build_gates(d, a.essential)
    g, steps = run(d, gates)

    # Independent recomputation with plain masks, so a bug in Filter/run cannot agree with itself.
    ess_mask = {"any": d.ess_any, "all": d.ess_all, "any_and_consensus": d.ess_any & d.ess_cons,
                "any_or_consensus": d.ess_any | d.ess_cons, "consensus": d.ess_cons}[a.essential]
    pool = d.in_aska & (d.cls | d.f1)
    check = pool & ((d.degradability_consensus >= d.attrs["cons_deg_cut"])
                    | (d.orth_hit & (d.degradability_consensus >= d.attrs["orth_floor"])))
    assert int((check & ess_mask).sum()) == len(g), "gate recomputation disagrees with the cascade"

    say(f"\n  {'gate':<12}{'before':>8}{'after':>8}{'removed':>9}   rule")
    gate_rows = []
    for s in steps:
        say(f"  {s.name:<12}{s.n_before:>8,}{s.n_after:>8,}{s.n_removed:>9,}   {s.rule}")
        gate_rows.append({"gate": s.name, "n_before": s.n_before, "n_after": s.n_after, "n_removed": s.n_removed, "rule": s.rule})
    say("\n  leave-one-out (survivors with that gate dropped):")
    for gate in gates:
        n = len(run(d, gates, skip=gate.name)[0])
        say(f"    drop {gate.name:<11}{n:>6,}  (+{n - len(g):,})")
        for r in gate_rows:
            if r["gate"] == gate.name:
                r["n_without_gate"] = n

    g = g.sort_values(["degradability_consensus", "essentiality_consensus"], ascending=False)
    gated_ex = run(d, gates, skip="essential")[0]
    below = gated_ex[gated_ex.ess_any & ~gated_ex.ess_cons].sort_values("degradability_consensus", ascending=False)

    # The final list: the absolute top-N by degradability among the first three gates, topped up with the
    # essential-gate survivors in degradability order until --n. A two-part rule, not a score.
    pool3 = gated_ex.sort_values("degradability_consensus", ascending=False)
    top = pool3.head(a.top)
    rest = g[~g.uniprot_ac.isin(top.uniprot_ac)].sort_values("degradability_consensus", ascending=False)
    room = max(a.n - len(top), 0)
    final = pd.concat([top, rest.head(room)])
    overflow = rest.iloc[room:]
    in_top, in_ess = final.uniprot_ac.isin(top.uniprot_ac), final.uniprot_ac.isin(g.uniprot_ac)
    final = final.assign(list_source=["both" if t and e else "top_degradability" if t else "essential_gate"
                                      for t, e in zip(in_top, in_ess)])
    assert len(final) == final.uniprot_ac.nunique() and len(final) <= a.n

    OUT.mkdir(parents=True, exist_ok=True)
    final[COLUMNS + ["list_source"]].to_csv(OUT / f"shortlist_ecoli_aska_top{a.n}.tsv", sep="\t", index=False)
    overflow[COLUMNS].to_csv(OUT / f"shortlist_ecoli_aska_top{a.n}_overflow.tsv", sep="\t", index=False)
    g[COLUMNS].to_csv(OUT / "shortlist_ecoli_aska.tsv", sep="\t", index=False)
    below[COLUMNS].to_csv(OUT / "shortlist_ecoli_aska_measured_below_cut.tsv", sep="\t", index=False)
    pd.DataFrame(gate_rows).to_csv(OUT / "shortlist_ecoli_aska_gates.tsv", sep="\t", index=False)

    say(f"\n  survivors {len(g):,} -> output/results/shortlist/shortlist_ecoli_aska.tsv")
    say(f"  measured essential, other gates pass, consensus below the cut: {len(below):,} "
        f"-> shortlist_ecoli_aska_measured_below_cut.tsv")
    # Random control: cytosolic ASKA proteins NOT on the list with degradability below the proteome median
    # (the ortholog floor), split measured-essential / not in the SAME proportion as the final list.
    cpool = d[d.in_aska & d.cytosolic & (d.degradability_consensus < d.attrs["orth_floor"]) & ~d.uniprot_ac.isin(final.uniprot_ac)]
    share = float(final.ess_any.mean())
    n_ess = round(a.control * share)
    ce, cn = cpool[cpool.ess_any], cpool[~cpool.ess_any]
    assert len(ce) >= n_ess and len(cn) >= a.control - n_ess, "pool too small for the requested split"
    ctrl = pd.concat([ce.sample(n_ess, random_state=a.seed), cn.sample(a.control - n_ess, random_state=a.seed)])
    ctrl = ctrl.assign(list_source="random_control_below_median")
    assert not set(ctrl.uniprot_ac) & set(final.uniprot_ac)
    ctrl[COLUMNS + ["list_source"]].to_csv(OUT / f"shortlist_ecoli_aska_control{a.control}.tsv", sep="\t", index=False)
    say(f"\n  random control {len(ctrl)} from {len(cpool):,} cytosolic ASKA proteins below the median ({d.attrs['orth_floor']:.4f}) not on the list; "
        f"measured essential {n_ess} / not {a.control - n_ess} (list share {share:.2f}; pool share {cpool.ess_any.mean():.3f}, seed {a.seed})")
    say(f"  control median length {ctrl.length.median():.0f} vs list {final.length.median():.0f} | median consensus {ctrl.degradability_consensus.median():.3f} "
        f"vs list {final.degradability_consensus.median():.3f} | ribosomal {int(ctrl.ribosomal.sum())}")

    say(f"\n  final list {len(final)}: top {a.top} by degradability {len(top)} (of which also pass the essential gate {int(top.uniprot_ac.isin(g.uniprot_ac).sum())}) "
        f"+ essential gate {len(final) - len(top)}; {len(overflow)} essential-gate proteins left out (lowest degradability, in the overflow table)")
    say(f"  list_source {final.list_source.value_counts().to_dict()} | lowest consensus kept from the essential gate "
        f"{final[final.list_source == 'essential_gate'].degradability_consensus.min():.4f}")
    say(f"  deg_route          {g.deg_route.value_counts().to_dict()}")
    say(f"  cyto_route         {g.cyto_route.value_counts().to_dict()}")
    say(f"  essentiality_evidence {g.essentiality_evidence.value_counts().sort_index().to_dict()}")
    say(f"  measured in all screens: {int(g.ess_all.sum())} | ribosomal {int(g.ribosomal.sum())} | on interest panel {int(g.on_panel.sum())}")
    say(f"  has a human ortholog (decorator, not a gate): {int(g.has_human_ortholog.astype(bool).sum())} of {len(g)}")
    say(f"  ligands_evidence   {g.ligands_evidence.value_counts().sort_index().to_dict()} | n_drugs > 0: {int((g.n_drugs > 0).sum())}")
    say(f"  median length survivors {g.length.median():.0f} | ASKA pool {d[d.in_aska].length.median():.0f} | proteome {d.length.median():.0f}")
    say(f"  ortholog_support   {g.ortholog_support.value_counts(dropna=False).to_dict()}")
    if not QUIET:
        say(g[["gene_name", "deg_route", "degradability_consensus", "essentiality_consensus", "deg_n_datasets",
               "has_human_ortholog", "n_drugs"]].head(15).to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
