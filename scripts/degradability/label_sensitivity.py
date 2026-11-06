"""Does the way the S. aureus screens are BINARIZED change which E. coli / K. pneumoniae proteins rank?

A sensitivity analysis on the DOWNSTREAM ranking, not a label search. `docs/degradability.md`: the
log2FC cutoff must not be tuned, because AUROC is not comparable across different labels -- a
stricter label inflates it mechanically. So nothing here is chosen by its AUROC. Each variant is
refitted and the whole of *E. coli* and *K. pneumoniae* is rescored; the question is how far the
ranking moves against the shipped one, and whether the shortlist moves with it.

Variants (per activator; every one uses the SAME two embeddings, fused late as shipped):

  shipped         the shipped label: ADEP4 log2FC <= -1 and padj < 0.05, ONC212 log2FC <= -1.
                  CONTROL: must reproduce the shipped `<act>_prob` exactly (free, from the cache)
  grey_dropped    proteins with -1 < log2FC <= -0.5 (called negative today) are REMOVED, not relabelled
  strong          positive = log2FC <= -2 (4-fold), everything else measured is negative
  replicate       ADEP4 only. Conlon's table has two replicates. Positive = BOTH <= -1; negative = both
                  > -1; the ~70 proteins that straddle the cutoff between replicates are removed
  onc212_matched  ONC212 only. Cutoff set so its positive rate on the proteins BOTH screens measured
                  equals ADEP4's (about -2.4), removing the stringency mismatch

Reads   output/results/other/activator_features.csv, the labels table, both embeddings
Writes  output/results/degradability/label_sensitivity.tsv
        data/processed/degradability/scratch/label_sensitivity/pred_<act>_<variant>_<species>.tsv

SPENDS hosted TabPFN credits (about 10,000 per uncached call). It counts the uncached calls first,
prints the cost, and refuses to run above `--max-credits`. `GRADI_TABPFN_CACHE_ONLY=1` makes it free
and is how the shipped control is checked. Needs `TABPFN_TOKEN` for anything uncached.

Run with the `gradi` env:  python scripts/degradability/label_sensitivity.py [--max-credits N] [-q]
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from scipy.stats import spearmanr  # noqa: E402

from src import degradability as D  # noqa: E402
from src import essentiality, localization, orthology  # noqa: E402
from src import proteomes as P  # noqa: E402
from src import tabpfn as T  # noqa: E402

OUT = REPO_ROOT / "output/results/degradability/label_sensitivity.tsv"
SCRATCH = REPO_ROOT / "data/processed/degradability/scratch/label_sensitivity"
CONLON = (REPO_ROOT / "data/raw/other/degradability/conlon2013_adep4_saureus/"
          "conlon2013_ADEP4_MRSA_TableS1-S2.xlsx")
EVAL_SPECIES = ("ecoli", "kpneumoniae")
VARIANTS = {
    "adep4": ("shipped", "grey_dropped", "strong", "replicate"),
    "onc212": ("shipped", "grey_dropped", "strong", "onc212_matched"),
}
STRONG = -2.0
GREY = -0.5
# The shortlist rule from `plotting/shortlist_excel.py` (ESS_MIN, DEG_MIN); restated, not imported,
# because that module pulls in stylia, which deletes the matplotlib cache at import.
ESS_MIN, DEG_MIN = 0.6, 0.7


def say(msg: str = "", quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def stage04():
    """`predict.py` as a module, for its `feature_frame` -- the one definition of the features."""
    spec = importlib.util.spec_from_file_location(
        "stage04", REPO_ROOT / "scripts/degradability/predict.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.VERBOSE = False
    return mod


def conlon_replicates() -> pd.DataFrame:
    """Per-replicate log2FC from Conlon 2013 Table S1, keyed on SACOL."""
    import openpyxl

    ws = openpyxl.load_workbook(CONLON, read_only=True, data_only=True)["Table S1"]
    rows = list(ws.iter_rows(values_only=True))
    hdr = [str(c).replace(" ", "") if c else "" for c in rows[3]]
    ip, ia, ib = hdr.index("Protein"), hdr.index("Sample_A"), hdr.index("Sample_B")
    out = []
    for r in rows[4:]:
        try:
            out.append((str(r[ip]).strip(), float(r[ia]), float(r[ib])))
        except (TypeError, ValueError):
            continue
    return pd.DataFrame(out, columns=["sacol", "rep_a", "rep_b"]).drop_duplicates("sacol")


def variant_labels(lab: pd.DataFrame, rep: pd.DataFrame, act: str, name: str) -> pd.Series:
    """0 / 1 / NaN aligned to `lab`. NaN = not measured or deliberately removed, never negative."""
    y = lab[act].astype(float)
    fc = lab[f"{act}_log2fc"]
    nan = np.nan
    if name == "shipped":
        return y
    if name == "grey_dropped":
        return y.mask((y == 0) & (fc <= GREY))
    if name == "strong":
        return pd.Series(np.where(y.notna(), ((fc <= STRONG) & (y == 1)).astype(float), nan),
                         index=lab.index)
    if name == "replicate":
        r = lab["sacol"].map(rep.set_index("sacol")["rep_a"]), lab["sacol"].map(rep.set_index("sacol")["rep_b"])
        a, b = r
        pos, neg = (a <= -1) & (b <= -1), (a > -1) & (b > -1)
        return pd.Series(np.where(pos, 1.0, np.where(neg, 0.0, nan)), index=lab.index).where(y.notna())
    if name == "onc212_matched":
        both = lab[lab["adep4"].notna() & lab["onc212_log2fc"].notna()]
        cut = float(both["onc212_log2fc"].quantile(float((both["adep4"] == 1).mean())))
        return pd.Series(np.where(fc.notna(), (fc <= cut).astype(float), nan), index=lab.index)
    raise ValueError(name)


def top_overlap(a: pd.Series, b: pd.Series, frac: float = 0.10) -> float:
    k = int(round(frac * len(a)))
    return len(set(a.nlargest(k).index) & set(b.nlargest(k).index)) / k


def shortlist(cons: pd.Series, gates: pd.DataFrame) -> set[str]:
    ok = ((gates["localization"] == "cytoplasm") & ~gates["has_human_ortholog"].astype(bool)
          & (gates["essentiality_consensus"] >= ESS_MIN) & (cons >= DEG_MIN))
    return set(gates.index[ok])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--max-credits", type=int, default=1_000_000,
                    help="refuse to run if the uncached calls would cost more than this")
    ap.add_argument("-q", "--quiet", action="store_true")
    q = ap.parse_args()
    maxc, q = q.max_credits, q.quiet

    say("label_sensitivity.py  does the binarization of the two screens change the E. coli / Kp ranking?", q)
    s4 = stage04()
    lab = D.load_labels()
    rep = conlon_replicates()
    feats = s4.feature_frame("saureus")
    frames = {sp: s4.feature_frame(sp) for sp in EVAL_SPECIES}
    shipped = {sp: D.load(sp).set_index("uniprot_ac") for sp in EVAL_SPECIES}
    say(f"  labels {len(lab):,} S. aureus proteins; Conlon replicates for "
        f"{int(lab['sacol'].isin(rep['sacol']).sum()):,} of them", q)

    # ---- the labels each variant trains on
    train: dict[tuple[str, str], tuple[pd.DataFrame, np.ndarray]] = {}
    say(f"\n  {'activator':<8}{'variant':<16}{'n':>6}{'pos':>6}{'rate':>8}   note", q)
    for act, names in VARIANTS.items():
        for v in names:
            y = variant_labels(lab, rep, act, v)
            sub = lab[y.notna() & lab["uniprot_ac"].isin(feats.index)]
            yy = y[sub.index].to_numpy(dtype=float).astype(int)
            train[(act, v)] = (sub, yy)
            note = {"grey_dropped": f"removes {int(lab[act].notna().sum() - y.notna().sum())} grey-zone proteins",
                    "replicate": f"removes {int(lab[act].notna().sum() - y.notna().sum())} discordant / unmatched",
                    "onc212_matched": "cutoff from ADEP4's rate on the shared proteins"}.get(v, "")
            say(f"  {act:<8}{v:<16}{len(sub):>6}{int(yy.sum()):>6}{yy.mean():>8.3f}   {note}", q)

    # ---- cost, BEFORE spending
    def calls():
        for (act, v), (sub, yy) in train.items():
            for b in D.FEATURE_BLOCKS:
                cols = D.block_columns(b, feats)
                Xtr = feats.loc[sub["uniprot_ac"], cols].to_numpy(dtype=float)
                for sp in EVAL_SPECIES:
                    yield (Xtr, yy, frames[sp][cols].to_numpy(dtype=float))

    n_calls = sum(1 for _ in calls())
    n_unc, credits = T.would_cost(calls())
    say(f"\n  CALLS {n_calls}   uncached {n_unc}   cost {credits:,} credits   (budget {maxc:,})", q)
    if credits > maxc:
        print(f"FAIL: {credits:,} credits exceeds --max-credits {maxc:,}", file=sys.stderr)
        return 1

    # ---- score every variant
    SCRATCH.mkdir(parents=True, exist_ok=True)
    pred: dict[tuple[str, str, str], pd.Series] = {}
    for (act, v), (sub, yy) in train.items():
        for sp in EVAL_SPECIES:
            per = []
            for b in D.FEATURE_BLOCKS:
                cols = D.block_columns(b, feats)
                Xtr = feats.loc[sub["uniprot_ac"], cols].to_numpy(dtype=float)
                p = T.predict_fold(Xtr, yy, frames[sp][cols].to_numpy(dtype=float),
                                   tag=f"label_sens/{act}/{v}/{sp}/{b}")
                per.append(pd.Series(p, index=frames[sp].index))
            fused = pd.concat(per, axis=1).mean(axis=1)
            pred[(act, v, sp)] = fused
            fused.rename(f"{act}_prob").round(4).to_csv(
                SCRATCH / f"pred_{act}_{v}_{sp}.tsv", sep="\t", index_label="uniprot_ac")
        say(f"  scored {act}/{v}", q)

    # ---- guard: the control must reproduce the shipped table
    for act in VARIANTS:
        for sp in EVAL_SPECIES:
            d = (pred[(act, "shipped", sp)] - shipped[sp][f"{act}_prob"]).abs().max()
            say(f"  control {act:<7}{sp:<14} max |fresh - shipped| = {d:.5f}", q)
            if d > 2e-4:
                print(f"FAIL: the shipped control does not reproduce for {act}/{sp} ({d:.5f}); the "
                      "variants below cannot be compared to it", file=sys.stderr)
                return 1

    # ---- assess
    rows = []
    for sp in EVAL_SPECIES:
        gates = (localization.load(sp).set_index("uniprot_ac")[["localization"]]
                 .join(orthology.load(sp).set_index("uniprot_ac")[["has_human_ortholog"]])
                 .join(essentiality.load(sp).set_index("uniprot_ac")[["essentiality_consensus"]]))
        length = P.load(sp).set_index("uniprot_ac")["sequence"].str.len()
        base = shipped[sp][["adep4_prob", "onc212_prob"]].copy()
        base_cons = D.consensus(base.reset_index()).set_axis(base.index)
        base_list = shortlist(base_cons, gates)
        for act, names in VARIANTS.items():
            ref = base[f"{act}_prob"]
            for v in names:
                p = pred[(act, v, sp)].reindex(ref.index)
                alt = base.copy()
                alt[f"{act}_prob"] = p
                cons = D.consensus(alt.reset_index()).set_axis(alt.index)
                lst = shortlist(cons, gates)
                sub, yy = train[(act, v)]
                rows.append({
                    "species": sp, "activator": act, "variant": v, "n_train": len(sub),
                    "n_pos": int(yy.sum()), "base_rate": round(float(yy.mean()), 4),
                    "spearman_vs_shipped": round(float(spearmanr(p, ref)[0]), 4),
                    "top10pct_overlap": round(top_overlap(p, ref), 4),
                    "spearman_with_length": round(float(spearmanr(p, length.reindex(p.index))[0]), 4),
                    "consensus_spearman_vs_shipped": round(float(spearmanr(cons, base_cons)[0]), 4),
                    "shortlist_n": len(lst),
                    "shortlist_overlap_with_shipped": round(len(lst & base_list) / max(len(base_list), 1), 4),
                    "shipped_shortlist_n": len(base_list)})
    res = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT, sep="\t", index=False)

    pd.set_option("display.width", 250)
    for sp in EVAL_SPECIES:
        say(f"\n  {sp.upper()}  (shortlist rule: cytoplasm, no human ortholog, essentiality >= {ESS_MIN}, "
            f"degradability consensus >= {DEG_MIN})", q)
        say(res[res.species == sp].drop(columns=["species", "shipped_shortlist_n"]).to_string(index=False), q)
    say(f"\n  wrote {OUT.relative_to(REPO_ROOT)}", q)
    return 0


if __name__ == "__main__":
    sys.exit(main())
