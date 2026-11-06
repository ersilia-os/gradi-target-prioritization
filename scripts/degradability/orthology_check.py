"""Do the proteins the model UPRANKS in E. coli / K. pneumoniae make sense against what is MEASURED in S. aureus?

The model is trained on the two S. aureus activator screens and applied elsewhere by embedding
similarity. Orthology is the independent lens: if a Gram-negative protein ranks high, its S. aureus
ortholog -- where one exists and was measured -- should usually have been degradable too.

Reads   data/processed/degradability/degradability_<species>.tsv     the fused scores
        data/processed/orthology/orthologs.tsv                       OrthoFinder + reciprocal best hit
        data/processed/degradability/evidence/labels_saureus.tsv     the MEASURED S. aureus calls
Writes  output/results/degradability/orthology_check_<species>.tsv   one row per protein, audit trail

**READ THIS BEFORE QUOTING ANY AGREEMENT NUMBER: IT IS NOT INDEPENDENT VALIDATION.** The ortholog's
own measured label was in the training set, so a high score for its E. coli ortholog is partly the
model copying that label through sequence similarity. What the check CAN show:
  1. how much of the upranked set has any S. aureus support at all, and how much of it is
     extrapolation with NO ortholog (the part the model is guessing, mostly from size and embedding),
  2. the DISCORDANT cases -- upranked, but the measured S. aureus ortholog was NOT degradable,
  3. the reverse direction: do S. aureus degraders get upranked in their E. coli orthologs.

Ortholog evidence is the union of OrthoFinder and reciprocal best hit; `--rbh-only` restricts it to
RBH. Many-to-many groups are kept: a protein is 'ortholog_degradable' if ANY measured ortholog is a
hit. An empty ortholog set is `no_ortholog`, never 'not degradable'.

Run with the `gradi` env:  python scripts/degradability/orthology_check.py [--species ecoli] [--rbh-only] [-q]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from scipy.stats import fisher_exact  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score  # noqa: E402

from src import degradability as D  # noqa: E402
from src import orthology as O  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output/results/degradability"
CATS = ("ortholog_degradable", "ortholog_not_degradable", "ortholog_unmeasured", "no_ortholog")
DEG_MIN = 0.7  # the shortlist's degradability-consensus floor (plotting/shortlist_excel.py)


def say(msg: str = "", quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def pct(n: int, d: int) -> str:
    return f"{n:>5} {100 * n / d:5.1f}%" if d else f"{n:>5}     -"


def build(sp: str, rbh_only: bool):
    deg = D.load(sp).set_index("uniprot_ac")
    prot = P.load(sp).set_index("uniprot_ac")
    lab = D.load_labels().set_index("uniprot_ac")
    o = O.load_orthologs()
    o = o[(o.query_species == sp) & (o.target_species == "saureus")]
    o = o[o.is_rbh] if rbh_only else o[o.is_rbh | o.is_ortholog_orthofinder]
    df = deg.join(prot[["gene_name", "sequence"]])
    df["length"] = df["sequence"].str.len()
    df = df.drop(columns="sequence")
    grp = o.groupby("query_ac")
    df = df.drop(columns=list(D.ORTHOLOG_HIT_COLUMNS), errors="ignore")
    sup = D.saureus_ortholog_support(sp, accessions=df.index, rbh_only=rbh_only)
    df = df.join(sup)
    df["ortholog_tier"] = np.where(df.index.isin(o.loc[o.is_rbh, "query_ac"]), "rbh",
                                   np.where(df["n_orthologs"] > 0, "orthofinder", ""))
    df["ortholog_genes"] = grp["target_ac"].apply(
        lambda s: ";".join(str(lab["gene"].get(a, a)) for a in sorted(set(s)))).reindex(df.index).fillna("")
    return df, o, lab


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", default="ecoli", choices=("ecoli", "kpneumoniae"))
    ap.add_argument("--rbh-only", action="store_true", help="count only reciprocal best hits as orthologs")
    ap.add_argument("-q", "--quiet", action="store_true")
    a = ap.parse_args()
    q, sp = a.quiet, a.species

    say(f"orthology_check.py  {sp}: does the upranked set make sense against measured S. aureus?", q)
    say(f"  ortholog evidence: {'reciprocal best hit only' if a.rbh_only else 'OrthoFinder or RBH'}", q)
    df, o, lab = build(sp, a.rbh_only)
    n = len(df)
    say(f"  {n:,} proteins; {int((df.n_orthologs > 0).sum()):,} have >=1 S. aureus ortholog "
        f"({int((df.ortholog_tier == 'rbh').sum()):,} by RBH)", q)

    sets = {}
    for act in D.ACTIVATORS:
        p = df[f"{act}_prob"]
        sets[f"{act} top 10%"] = p >= p.quantile(0.90)
    sets[f"consensus >= {DEG_MIN}"] = df["degradability_consensus"] >= DEG_MIN
    sets["consensus >= 0.9"] = df["degradability_consensus"] >= 0.9

    rows = []
    for act in D.ACTIVATORS:
        sup = df[f"{act}_support"]
        say(f"\n  {act.upper()}  S. aureus support of the ortholog set  (counts and share of the set)", q)
        say(f"    {'set':<22}{'n':>6}   " + "".join(f"{c:>26}" for c in CATS), q)
        for name, mask in {"whole proteome": pd.Series(True, index=df.index), **sets}.items():
            m = mask & True
            cnt = [int((sup[m] == c).sum()) for c in CATS]
            say(f"    {name:<22}{int(m.sum()):>6}   " + "".join(f"{pct(c, int(m.sum())):>26}" for c in cnt), q)
            rows.append({"activator": act, "set": name, "n": int(m.sum()),
                         **{c: v for c, v in zip(CATS, cnt)}})
        # hit rate among proteins whose ortholog was measured, top set vs the rest
        top = sets[f"{act} top 10%"]
        meas = sup.isin(["ortholog_degradable", "ortholog_not_degradable"])
        a1, b1 = int(((sup == "ortholog_degradable") & top).sum()), int(((sup == "ortholog_not_degradable") & top).sum())
        a0, b0 = int(((sup == "ortholog_degradable") & ~top).sum()), int(((sup == "ortholog_not_degradable") & ~top).sum())
        orr = fisher_exact([[a1, b1], [a0, b0]])
        say(f"    measured-ortholog hit rate: top 10% {a1}/{a1 + b1} = {100 * a1 / max(a1 + b1, 1):.0f}%   "
            f"rest {a0}/{a0 + b0} = {100 * a0 / max(a0 + b0, 1):.0f}%   Fisher OR {orr.statistic:.1f} (p {orr.pvalue:.1e})", q)
        y = (df.loc[meas, f"{act}_ortholog_hit"] == 1).astype(int)
        s = df.loc[meas, f"{act}_prob"]
        ln = -np.log10(df.loc[meas, "length"])
        say(f"    on the {int(meas.sum())} proteins with a measured ortholog: score vs ortholog call  "
            f"AUROC {roc_auc_score(y, s):.3f} PR {average_precision_score(y, s):.3f} "
            f"(base {y.mean():.3f}) | shorter-first baseline AUROC {roc_auc_score(y, ln):.3f}", q)

    # ---- what share of the upranked set is pure extrapolation
    say(f"\n  UPRANKED WITH NO S. AUREUS SUPPORT  (consensus >= {DEG_MIN}, n={int(sets[f'consensus >= {DEG_MIN}'].sum())})", q)
    up = df[sets[f"consensus >= {DEG_MIN}"]]
    for lab_, g in (("has an ortholog", up[up.n_orthologs > 0]), ("NO ortholog", up[up.n_orthologs == 0]),
                    ("whole proteome", df)):
        say(f"    {lab_:<17} n={len(g):>5}  median length {g['length'].median():>4.0f} aa  "
            f"< 100 aa {100 * (g['length'] < 100).mean():4.1f}%  median nn_similarity {g['nn_similarity'].median():.3f}", q)

    # ---- the discordant cases
    say("\n  DISCORDANT: upranked, but the MEASURED S. aureus ortholog was not degradable", q)
    for act in D.ACTIVATORS:
        d = df[sets[f"{act} top 10%"] & (df[f"{act}_support"] == "ortholog_not_degradable")]
        d = d.sort_values(f"{act}_prob", ascending=False)
        say(f"    {act.upper()}: {len(d)} of {int(sets[f'{act} top 10%'].sum())} top-decile proteins"
            f" ({100 * len(d) / int(sets[f'{act} top 10%'].sum()):.0f}%)", q)
        for ac, r in d.head(12).iterrows():
            say(f"      {str(r['gene_name'])[:10]:<11}{ac:<11} prob {r[f'{act}_prob']:.3f}  "
                f"S. aureus {r['ortholog_genes'][:26]:<27} min log2FC {r[f'{act}_ortholog_log2fc_min']:>6.2f}"
                f"  len {r['length']}", q)

    # ---- reverse: do the S. aureus degraders get upranked in the target species
    say("\n  REVERSE: measured S. aureus proteins with an ortholog here -- is the ortholog upranked?", q)
    for act in D.ACTIVATORS:
        rev = o.copy()
        rev["sa_hit"] = rev["target_ac"].map(lab[act])
        rev["score"] = rev["query_ac"].map(df["degradability_consensus"])
        rev = rev.dropna(subset=["sa_hit"]).groupby("target_ac").agg(sa_hit=("sa_hit", "max"), score=("score", "max"))
        h, nh = rev[rev.sa_hit == 1], rev[rev.sa_hit == 0]
        say(f"    {act.upper()}: S. aureus degraders {len(h)}: ortholog consensus >= {DEG_MIN} for "
            f"{100 * (h.score >= DEG_MIN).mean():.0f}%  | non-degraders {len(nh)}: "
            f"{100 * (nh.score >= DEG_MIN).mean():.0f}%   (median consensus {h.score.median():.2f} vs {nh.score.median():.2f})", q)

    # ---- the top of the list, by eye
    say("\n  TOP 25 by degradability consensus", q)
    say(f"    {'gene':<10}{'cons':>6}{'adep4':>7}{'onc212':>7}  {'S. aureus ortholog':<22}{'ADEP4':<24}{'ONC212':<24}{'len':>5}", q)
    for ac, r in df.sort_values("degradability_consensus", ascending=False).head(25).iterrows():
        say(f"    {str(r['gene_name'])[:9]:<10}{r['degradability_consensus']:>6.3f}{r['adep4_prob']:>7.2f}"
            f"{r['onc212_prob']:>7.2f}  {r['ortholog_genes'][:20] or '-':<22}"
            f"{r['adep4_support'].replace('ortholog_', '')[:22]:<24}{r['onc212_support'].replace('ortholog_', '')[:22]:<24}{r['length']:>5}", q)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    keep = ["gene_name", "length", "adep4_prob", "onc212_prob", "degradability_consensus", "nn_similarity",
            "n_orthologs", "ortholog_tier", "ortholog_genes", "adep4_support", "adep4_ortholog_hit",
            "adep4_ortholog_log2fc_min", "onc212_support", "onc212_ortholog_hit", "onc212_ortholog_log2fc_min"]
    path = OUT_DIR / f"orthology_check_{sp}{'_rbh' if a.rbh_only else ''}.tsv"
    df[keep].to_csv(path, sep="\t", index_label="uniprot_ac")
    say(f"\n  wrote {path.relative_to(REPO_ROOT)}", q)
    return 0


if __name__ == "__main__":
    sys.exit(main())
