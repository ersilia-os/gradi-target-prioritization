"""Figures for stage 00: how the proteomes got their gene names, and which identifiers to join on.

Stage 00's real job is the naming gap. Gene names are what literature and external databases key on,
and on the anchor organism UniProt supplies them for **18.4%** of proteins. This shows what the
three-tier fill bought, how the ties were broken, and — the operational point — which identifier
column is actually safe to join on.

    gene_names.png   A  coverage by fill tier, per species
                     B  how the contested fills were decided
                     C  identifier coverage, i.e. what you can join on

Everything is read from stage 00's own outputs: `gene_name_source` in `proteome_<species>.tsv` and
`evidence/name_audit.tsv`. Nothing is recomputed.

Run with the `gradi` env:
    python scripts/plots/proteomes.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "proteomes"
AUDIT = REPO_ROOT / "data" / "processed" / "proteomes" / "evidence" / "name_audit.tsv"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli",
          "saureus": "S. aureus", "human": "H. sapiens"}
# Short codes for the narrow grouped-bar panel, where four italic binomials do not fit.
SHORT = {"kpneumoniae": "Kp", "ecoli": "Ec", "saureus": "Sa", "human": "Hs"}

# The fill tiers, in the order stage 00 applies them: what UniProt already had, then two
# within-species donor tiers. `none` is the honest remainder -- no cross-species guessing.
TIERS = ("anchor", "species_exact", "species_uniref90", "none")
TIER_LABEL = {
    "anchor": "anchor (UniProt's own)",
    "species_exact": "species_exact (identical sequence)",
    "species_uniref90": "species_uniref90 (same cluster)",
    "none": "still unnamed",
}
TIER_COLOR = {
    "anchor": NC.plum,
    "species_exact": NC.orange,
    "species_uniref90": NC.mint,
    "none": NC.gray,
}

# The identifier columns a downstream join might use, and where each comes from.
ID_COLS = ("gene_name", "locus_tag", "refseq", "geneid")
ID_COLOR = {"gene_name": NC.plum, "locus_tag": NC.orange,
            "refseq": NC.mint, "geneid": NC.blue}


def nonempty(s: pd.Series) -> float:
    """Percent of rows carrying a real value. Stage 00 writes '' for absent, never NaN."""
    return 100 * (s.fillna("").astype(str).str.strip() != "").mean()


def tier_table(species: tuple[str, ...]) -> pd.DataFrame:
    """Per species, the share of the proteome in each fill tier."""
    rows = []
    for sp in species:
        src = P.load_full(sp)["gene_name_source"].fillna("none").astype(str)
        n = len(src)
        row = {"species": sp, "n": n}
        for t in TIERS:
            row[t] = 100 * (src == t).sum() / n
        row["named"] = 100 - row["none"]
        rows.append(row)
    return pd.DataFrame(rows)


def id_table(species: tuple[str, ...]) -> pd.DataFrame:
    """Per species, coverage of each identifier column."""
    rows = []
    for sp in species:
        d = P.load_full(sp)   # ID_COLS needs refseq + geneid, now in evidence/
        lt = P.load_locus_tags(sp)[["uniprot_ac", "locus_tag"]]
        m = d.merge(lt, on="uniprot_ac", how="left")
        rows.append({"species": sp, **{c: nonempty(m[c]) for c in ID_COLS}})
    return pd.DataFrame(rows)


def rule_table() -> pd.DataFrame:
    """How each audited fill was decided, with the `more attested (n vs m)` variants pooled."""
    a = pd.read_csv(AUDIT, sep="\t")
    a["rule_group"] = np.where(a.rule.str.startswith("more attested"), "more attested", a.rule)
    out = (a.groupby(["rule_group", "contested"]).size()
             .rename("n").reset_index().sort_values("n"))
    return out


# ---------------------------------------------------------------- panels


def plot_tiers(ax, tiers: pd.DataFrame, abc: str) -> None:
    """Stacked coverage by tier. The first segment is what UniProt gave us unaided."""
    species = list(tiers.species)
    left = np.zeros(len(species))
    for t in TIERS:
        vals = tiers[t].to_numpy()
        ax.barh([LABELS[s] for s in species], vals, left=left,
                color=TIER_COLOR[t], label=TIER_LABEL[t])
        left += vals
    # No separate marker for "where the anchor tier ended": the first segment's right edge already
    # is that boundary, and an extra tick on top of it needs a legend entry to earn its place.
    ax.set_xlim(0, 100)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    stylia.label(ax, xlabel="% of proteome", ylabel="",
                 title="Gene names by fill tier", abc=abc)


def plot_rules(ax, rules: pd.DataFrame, abc: str) -> None:
    """What decided each fill. Most were uncontested; the rest went to the total order."""
    colors = [NC.gray if not c else NC.plum for c in rules.contested]
    ax.barh(rules.rule_group, rules.n, color=colors)
    for y, (n, tot) in enumerate(zip(rules.n, [rules.n.sum()] * len(rules))):
        ax.annotate(f"{n:,}", (n, y), textcoords="offset points", xytext=(4, -3))
    ax.set_xscale("log")
    ax.set_xlim(1, rules.n.max() * 4)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor=NC.gray, label="one candidate"),
                       Patch(facecolor=NC.plum, label="contested")],
              loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2)
    stylia.label(ax, xlabel="fills (log)", ylabel="",
                 title=f"{rules.n.sum():,} fills decided", abc=abc)


def plot_identifiers(ax, ids: pd.DataFrame, abc: str) -> None:
    """Coverage of each join key. This is the panel that settles what to join on."""
    species = list(ids.species)
    x = np.arange(len(species))
    for i, c in enumerate(ID_COLS):
        off = (i - (len(ID_COLS) - 1) / 2) * 0.2
        ax.bar(x + off, ids[c], 0.19, color=ID_COLOR[c], label=c)
    ax.set_xticks(x)
    ax.set_xticklabels([SHORT[s] for s in species])
    ax.set_ylim(0, 105)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2)
    stylia.label(ax, xlabel="", ylabel="% of proteome with a value",
                 title="What you can join on", abc=abc)


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(P.SPECIES), choices=list(P.SPECIES))
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    order = tuple(s for s in LABELS if s in args.species)
    tiers = tier_table(order)
    ids = id_table(order)
    rules = rule_table()

    fig, axs = stylia.create_figure(1, 3, width_ratios=[3, 2, 2])
    plot_tiers(axs.next(), tiers, abc="A")
    plot_rules(axs.next(), rules, abc="B")
    plot_identifiers(axs.next(), ids, abc="C")
    out = OUT_DIR / "gene_names.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  NAMING COVERAGE by tier (% of proteome)")
    say(f"     {'species':<14} {'n':>6} {'anchor':>8} {'+exact':>8} {'+uniref90':>10} "
        f"{'= named':>8} {'unnamed':>8}")
    for _, r in tiers.iterrows():
        say(f"     {r.species:<14} {int(r.n):>6} {r['anchor']:>7.1f}% {r.species_exact:>7.1f}% "
            f"{r.species_uniref90:>9.1f}% {r.named:>7.1f}% {r['none']:>7.1f}%")

    say("\n  HOW FILLS WERE DECIDED")
    for _, r in rules.sort_values('n', ascending=False).iterrows():
        say(f"     {r.rule_group:<30} {r.n:>5}   "
            f"{'contested' if r.contested else 'one candidate'}")

    say("\n  IDENTIFIER COVERAGE (% of proteome)")
    say(f"     {'species':<14} " + " ".join(f"{c:>11}" for c in ID_COLS))
    for _, r in ids.iterrows():
        say(f"     {r.species:<14} " + " ".join(f"{r[c]:>10.1f}%" for c in ID_COLS))


if __name__ == "__main__":
    main()
