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

THE LIST THAT IS DELIVERED is the FUNCTION- AND LENGTH-REPRESENTATIVE list of 200 (`balanced_list`), not the cascade above, which is kept
as the cost analysis (funnel, leave-one-out, `shortlist_ecoli_aska*.tsv`). Slots follow the cytosolic ASKA proteome on COG function and on
four length bands whose target is the equal-weight mean of three distributions (`length_targets`). Eligible = degradability consensus in
the top 30% of the protein's own length decile (`--len-adj-min`, CHOSEN: the consensus is length-driven, rho -0.73); inside a cell the
order is measured essential, then the essential rule (measured any OR `essentiality_consensus` >= `--ess-or-cut`, CHOSEN), then ligand evidence (never for
ribosomal proteins), then within-decile degradability. One control of 50, raw consensus below the proteome median, matched on function and
measured-essential share but not on length. The raw degradability gate and essentiality are FLAGS here.

Writes   output/results/shortlist/shortlist_ecoli_aska_balanced200.tsv           THE LIST
         output/results/shortlist/shortlist_ecoli_aska_balanced_control50.tsv    THE CONTROL
         output/results/shortlist/shortlist_ecoli_aska_length_targets.tsv        the three length distributions and their mean
         output/results/shortlist/shortlist_ecoli_aska{,_measured_below_cut,_gates}.tsv   the four-gate cascade analysis
         output/results/shortlist/shortlist_ecoli_aska_top200{,_overflow}.tsv, ..._control50.tsv   the EARLIER top-50-plus-essential-gate list (superseded)

Run with the `gradi` env:  python plotting/shortlist_ecoli.py [--essential any_and_consensus] [-q]
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
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
           "n_papers_uniprot_prokaryotic", "length", "ribosomal", "on_panel", "cog"]
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
    from src import function, interest, ligandability, orthology, pockets

    for mod, cols in ((orthology, ["has_human_ortholog", "bacterial_panel_orthologs"]),
                      (ligandability, ["ligands_consensus", "ligands_evidence", "n_ligands_own", "n_ligands_bacterial", "n_drugs"]),
                      (pockets, ["pockets_evidence", "n_ligands_pdb", "n_ligands_alphafill", "n_pdb_structures", "p2rank_score"])):
        part = mod.load(SPECIES)[["uniprot_ac", *[c for c in cols if c not in d.columns]]]
        d = d.merge(part, on="uniprot_ac", how="left", validate="1:1")
    fn = function.load(SPECIES)[["uniprot_ac", "cog_categories"]]
    d = d.merge(fn, on="uniprot_ac", how="left", validate="1:1")
    # The FIRST COG letter is the protein's function (9-16% carry several); no letter -> "none" (unclassified)
    d["cog"] = d.cog_categories.fillna("").astype(str).str.replace(r"[^A-Z]", "", regex=True).str[:1].replace("", "none")
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


LENGTH_EDGES = [0, 100, 150, 300, 100000]
LENGTH_LABELS = ["<=100", "100-150", "150-300", ">300"]


def band_labels(edges: list) -> list:
    """Labels for right-closed length bands, e.g. [0, 100, 300, 1e5] -> ["<=100", "100-300", ">300"]."""
    return [f"<={edges[1]}"] + [f"{lo}-{hi}" for lo, hi in zip(edges[1:-2], edges[2:-1])] + [f">{edges[-2]}"]


def length_targets(d: pd.DataFrame) -> pd.DataFrame:
    """Target length mix: the equal-weight average of three distributions, one row per band.

    E. coli cytosolic ASKA proteome (representativeness), the measured-essential proteins in it (the target
    space of the essential rule) and the S. aureus proteome (where the degradability labels come from). No
    single distribution is right: degradation references run shorter and ligandable proteins longer.
    """
    from src import proteomes

    cy = d[d.in_aska & d.cytosolic]
    sa = proteomes.load("saureus")
    sets = {"ecoli_cytosolic": cy.length, "ecoli_measured_essential": cy[cy.deg_essential_any == 1].length,
            "saureus_proteome": sa.sequence.str.len()}
    t = pd.DataFrame({k: pd.cut(v, LENGTH_EDGES, labels=LENGTH_LABELS).value_counts(normalize=True).reindex(LENGTH_LABELS)
                      for k, v in sets.items()})
    t["target"] = t.mean(axis=1)
    assert abs(t.target.sum() - 1) < 1e-9
    return t


def margin_greedy(cells: dict, n: int, target_f: pd.Series, target_b: pd.Series, tiebreak: dict, pick=None) -> list:
    """Pick `n` proteins so BOTH margins (function letter, length band) track their targets.

    `cells[(f, b)]` is the ordered candidate list of one cell. Each step adds one protein to the cell whose
    function and length band are the furthest below their targets (the squared-deviation cost of both
    margins), ties broken by `tiebreak` (smaller is better). A cell with no candidates left is skipped, so
    availability limits are absorbed by the other cells. Exact joint cells are not targeted: they are
    too small for 200 proteins.
    """
    nf, nb = {f: 0 for f in target_f.index}, {b: 0 for b in target_b.index}
    cells = {k: list(v) for k, v in cells.items() if v and k[0] in nf and k[1] in nb}  # a function/band with no target gets no slot
    chosen = []
    for _ in range(n):
        best = None
        for (f, b), lst in cells.items():
            if not lst:
                continue
            cost = (2 * (nf[f] - target_f[f]) + 1) + (2 * (nb[b] - target_b[b]) + 1)
            key = (round(cost, 9), tiebreak[lst[0]])
            if best is None or key < best[0]:
                best = (key, (f, b))
        if best is None:
            break
        f, b = best[1]
        i = pick(cells[(f, b)], chosen) if pick else 0
        chosen.append(cells[(f, b)].pop(i))
        nf[f] += 1
        nb[b] += 1
    return chosen


def balanced_list(d: pd.DataFrame, n: int, or_cut: float, len_adj_min: float, n_control: int, seed: int,
                  band_edges: list | None = None, band_shares: dict | None = None,
                  order: list | None = None, lig_len_adj_min: float | None = None, n_len_bins: int = 10):
    """Function- AND length-representative list by degradability. Returns (list, control, notes).

    Slots follow the cytosolic ASKA proteome on two margins: COG function letter and length band. A
    protein is ELIGIBLE if its degradability consensus is in the top `1 - len_adj_min` of its own length decile
    (the raw gate cannot be met by long proteins: the consensus is strongly length-driven). Inside a cell the
    order is, as stated predicates and never a score: (1) measured essential in any screen, (2) meets the
    essential rule (measured OR `essentiality_consensus` >= `or_cut`, a CHOSEN cut), (3) has ligand evidence, then (4) the
    percentile within its length decile. Ribosomal proteins never get the ligand preference: their drug
    records belong to the whole ribosome.
    The control matches the list's function mix and measured-essential share, and its length bands as far as the pool
    allows, drawn at random from proteins with RAW degradability below the proteome median that are not on the list.
    """
    pool = d[d.in_aska & d.cytosolic].copy()
    edges = band_edges or LENGTH_EDGES
    labels = band_labels(edges)
    pool["lband"] = pd.cut(pool.length, edges, labels=labels).astype(str)
    pool["len_decile"] = pd.qcut(pool.length, n_len_bins, labels=False)
    pool["deg_pct_in_len"] = pool.groupby("len_decile").degradability_consensus.rank(pct=True)
    pool["ligand_evidence"] = ((pool.n_ligands_own.fillna(0) > 0) | (pool.n_ligands_bacterial.fillna(0) > 0)
                               | (pool.n_ligands_pdb.fillna(0) > 0) | (pool.n_drugs.fillna(0) > 0))
    pool["meas_ess"] = pool.deg_essential_any == 1  # measured is better evidence than a prediction; a missing DEG call is not a measured negative
    pool["ess_rule"] = pool.meas_ess | (pool.essentiality_consensus >= or_cut)
    pool["lig_rank"] = pool.ligand_evidence & ~pool.ribosomal
    share_f = pool.cog.value_counts(normalize=True)
    share_b = pd.Series(band_shares) if band_shares else pool.lband.value_counts(normalize=True)
    share_b = share_b.reindex(labels).fillna(0)

    ok = pool.deg_pct_in_len >= len_adj_min
    if lig_len_adj_min is not None:  # ligand evidence may lower the degradability bar (still never ribosomal)
        ok = ok | (pool.lig_rank & (pool.deg_pct_in_len >= lig_len_adj_min))
    elig = pool[ok].sort_values(order or ["meas_ess", "ess_rule", "lig_rank", "deg_pct_in_len"], ascending=False)
    pos = {u: i for i, u in enumerate(elig.uniprot_ac)}
    cells = {k: g.uniprot_ac.tolist() for k, g in elig.groupby(["cog", "lband"])}
    chosen = margin_greedy(cells, n, share_f * n, share_b * n, pos)
    assert len(chosen) == n, f"only {len(chosen)} eligible proteins fit the quotas"
    top50 = set(pool.nlargest(50, "degradability_consensus").uniprot_ac)
    final = pool.set_index("uniprot_ac").loc[chosen].reset_index()
    final = final.assign(list_source=["essential_rule" if e else "degradability_only" for e in final.ess_rule],
                         top50=final.uniprot_ac.isin(top50), passes_degradability_gate=final.degradable.astype(bool),
                         essential_strict=(final.deg_essential_any == 1) & (final.essentiality_consensus >= d.attrs["cons_cut"]),
                         length_band=final.lband)
    meas, cons = final.deg_essential_any == 1, final.essentiality_consensus >= or_cut
    final["essential_route"] = ["measured_and_consensus" if m and c else "measured_only" if m else "consensus_only" if c else "neither"
                                for m, c in zip(meas, cons)]
    final = final.sort_values("degradability_consensus", ascending=False)

    # ONE control, drawn with the LIST's two margins and its measured-essential share, random inside each cell, from proteins
    # whose RAW degradability is below the proteome median. No protein up to 100 aa and one of 286 at 100-150 aa scores that
    # low, so the short bands cannot be filled and the control is longer than the list. A control below the median of its own
    # length decile would fill them but hold raw scores above the proteome median, which is not a low-degradability control.
    free = pool[~pool.uniprot_ac.isin(final.uniprot_ac)]
    want = float((final.deg_essential_any == 1).mean())
    lf, lb = final.cog.value_counts(normalize=True), final.lband.value_counts(normalize=True).reindex(labels).fillna(0)

    def draw(cp: pd.DataFrame, name: str) -> pd.DataFrame:
        cp = cp.copy()
        cp["ess"] = cp.deg_essential_any == 1
        cp["rnd"] = np.random.default_rng(seed).permutation(len(cp))
        cells_c = {k: g.sort_values("rnd").uniprot_ac.tolist() for k, g in cp.groupby(["cog", "lband"])}
        ess_of, rnd_of = cp.set_index("uniprot_ac").ess.to_dict(), cp.set_index("uniprot_ac").rnd.to_dict()

        def pick(lst, chosen_so_far):
            # steer the running measured-essential share toward the list's
            prefer = sum(ess_of[u] for u in chosen_so_far) / max(len(chosen_so_far), 1) < want
            return next((i for i, u in enumerate(lst) if ess_of[u] == prefer), 0)

        chosen = margin_greedy(cells_c, n_control, lf * n_control, lb * n_control, rnd_of, pick=pick)
        out = cp.set_index("uniprot_ac").loc[chosen].reset_index().assign(list_source=name, length_band=lambda x: x.lband)
        assert len(out) == n_control and not set(out.uniprot_ac) & set(final.uniprot_ac)
        return out

    cp_raw = free[free.degradability_consensus < d.attrs["orth_floor"]]
    ctrl = draw(cp_raw, "random_control_below_median")
    assert (ctrl.degradability_consensus < d.attrs["orth_floor"]).all()
    band = lambda x: x.lband.value_counts().reindex(labels).fillna(0).astype(int).to_dict()
    notes = {"eligible": len(elig), "control_pool": len(cp_raw), "share_b": share_b.to_dict(),
             "list_b": band(final), "ctrl_b": band(ctrl), "labels": labels}
    return final, ctrl, notes


def main() -> int:
    global QUIET
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--essential", choices=list(ESSENTIAL_DEFINITIONS), default="any_and_consensus")
    ap.add_argument("--top", type=int, default=50, help="absolute top-N by degradability_consensus among the proteins passing the first three gates")
    ap.add_argument("--n", type=int, default=200, help="size of the final list")
    ap.add_argument("--control", type=int, default=50, help="size of the random below-median control")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--balanced-n", type=int, default=200, help="size of the function-balanced list")
    ap.add_argument("--len-adj-min", type=float, default=0.7, help="CHOSEN: minimum percentile of degradability_consensus within the protein's own length decile")
    ap.add_argument("--ess-or-cut", type=float, default=0.7, help="CHOSEN, not derived: essentiality_consensus cut of the OR rule")
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

    targets = length_targets(d)
    targets.to_csv(OUT / "shortlist_ecoli_aska_length_targets.tsv", sep="\t")
    bal, bctrl, bn = balanced_list(d, a.balanced_n, a.ess_or_cut, a.len_adj_min, a.control, a.seed, band_shares=targets.target.to_dict())
    bal[COLUMNS + ["list_source", "essential_route", "essential_strict", "top50", "passes_degradability_gate", "ligand_evidence",
                   "length_band", "deg_pct_in_len"]].to_csv(OUT / f"shortlist_ecoli_aska_balanced{a.balanced_n}.tsv", sep="\t", index=False)
    bctrl[COLUMNS + ["list_source", "length_band", "deg_pct_in_len"]].to_csv(OUT / f"shortlist_ecoli_aska_balanced_control{a.control}.tsv", sep="\t", index=False)
    cyto = d[d.in_aska & d.cytosolic]
    sf, lf = cyto.cog.value_counts(normalize=True), bal.cog.value_counts(normalize=True).reindex(cyto.cog.value_counts().index).fillna(0)
    say(f"\n  BALANCED LIST {len(bal)}: slots on function AND length band (cytosolic ASKA proteome), eligible = top {100 * (1 - a.len_adj_min):.0f}% of "
        f"degradability inside its own length decile ({bn['eligible']:,} proteins); order inside a cell: measured essential, essential rule (measured any OR consensus >= "
        f"{a.ess_or_cut}, chosen), ligand evidence, within-decile degradability")
    say("    length bands (target % = mean of E. coli cytosolic, E. coli measured essential, S. aureus proteome -> list n): "
        + " | ".join(f"{k} {100 * bn['share_b'][k]:.0f}% -> {bn['list_b'][k]}" for k in LENGTH_LABELS))
    say(f"    largest function deviation {100 * (lf - sf).abs().max():.1f} points ({(lf - sf).abs().idxmax()}) | functions beyond 2 points: "
        f"{{{', '.join(f'{k} {100 * v:+.1f}' for k, v in (lf - sf).items() if abs(v) > 0.02)}}}")
    say(f"    measured essential {int((bal.deg_essential_any == 1).sum())} | essential rule {int((bal.list_source == 'essential_rule').sum())} | strict "
        f"{int(bal.essential_strict.sum())} | top 50 inside {int(bal.top50.sum())} | ribosomal {int(bal.ribosomal.sum())} | "
        f"fail the degradability gate {int((~bal.passes_degradability_gate).sum())} | median length {bal.length.median():.0f} | "
        f"ligand evidence {int(bal.ligand_evidence.sum())} | median degradability {bal.degradability_consensus.median():.2f}")
    say(f"    control: {len(bctrl)} from {bn['control_pool']:,} (raw consensus below the proteome median) | length bands {bn['ctrl_b']} | measured essential "
        f"{int((bctrl.deg_essential_any == 1).sum())} | raw consensus median {bctrl.degradability_consensus.median():.2f}, max {bctrl.degradability_consensus.max():.3f} | "
        f"median length {bctrl.length.median():.0f}")

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
