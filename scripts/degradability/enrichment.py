"""Fisher-exact enrichment of the degradability ranking: COG categories, and the consortium panel.

Asks one question repeatedly: **is this group of proteins over-represented at the top of the
ranking?** Same 2x2 test throughout, so COG categories and the GraDi proteins of interest are
scored on exactly the same footing.

                          in top N     not in top N
    in the group              a              b
    not in the group          c              d

**"The top" is the top 10% of each proteome**, not a fixed count -- 573 Kp / 440 Ec / 289 Sa. The
three proteomes differ by 2x in size, so one absolute N would define "the top" as 1.7% of Kp against
3.5% of Sa and the odds ratios would not be comparable between species. `--top-pct` changes the
fraction; `--top-n` forces an absolute count instead.

`scipy.stats.fisher_exact`, two-sided, then **Benjamini-Hochberg FDR within each
species x activator x test-family block**. The COG categories and the consortium panel are corrected
separately: they are different questions asked of the same ranking, and pooling them would make a
COG category's q-value depend on how many panel families happened to be tested beside it.

Outputs:
    output/results/degradability/enrichment_fisher.tsv   every test, with a, b, c, d
    output/plots/degradability/cog_fisher.png            COG categories, both activators
    output/plots/degradability/interest_panel.png        where the consortium panel actually ranks

Three things to know before reading the numbers
-----------------------------------------------
**One letter per protein.** `cog_category` is the first letter of a possibly multi-letter COG
assignment; 9-16% of classified proteins carry more than one (see `docs/function.md`). `--multi`
re-runs counting a protein in every category it carries, which is the more generous reading; the
default is the single letter because that is the column the rest of the pipeline uses.

**Unclassified is tested, not dropped.** Proteins with no COG letter are their own group. Removing
them would silently change the background for every other test, and their share of the top of the
list is itself a result.

**The odds ratio is reported raw and plotted corrected.** A zero cell makes the ratio 0 or infinite,
which cannot be drawn; the plot uses a Haldane-Anscombe (+0.5) correction while the TSV keeps the
uncorrected value and the raw counts, so nothing is hidden.

Run with the `gradi` env:
    python scripts/degradability/enrichment.py
    python scripts/degradability/enrichment.py --top-pct 5 --multi
    python scripts/degradability/enrichment.py --top-n 100
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import stylia
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import degradability as D  # noqa: E402
from src import function as F  # noqa: E402
from src import interest as I  # noqa: E402
from src import localization as LOC  # noqa: E402
from src import proteomes as P  # noqa: E402

PLOT_DIR = REPO_ROOT / "output" / "plots" / "degradability"
RESULT_DIR = REPO_ROOT / "output" / "results" / "degradability"

# Format: slide | Style: ersilia -- change with stylia.set_format() / stylia.set_style()
stylia.set_format("slide")
stylia.set_style("ersilia")

NC = stylia.NamedColors()
LABELS = {"kpneumoniae": "K. pneumoniae", "ecoli": "E. coli", "saureus": "S. aureus"}
SPECIES_COLOR = {"kpneumoniae": NC.plum, "ecoli": NC.orange, "saureus": NC.mint}

FDR = 0.05
UNCLASSIFIED = "unclassified"

# COGclassifier ships the authoritative letter -> (group, colour, description) table; read it rather
# than hard-coding 26 descriptions that could drift from the tool that produced the calls.
COG_RESOURCE = "cogclassifier/resources/cog_func_category.tsv"


def cog_category_names() -> dict[str, str]:
    """`{letter: description}` from COGclassifier's own bundled resource."""
    import cogclassifier

    path = Path(cogclassifier.__file__).parent / "resources" / "cog_func_category.tsv"
    tab = pd.read_csv(path, sep="\t", header=None,
                      names=["letter", "group", "color", "description"])
    return dict(zip(tab.letter, tab.description))


def frame(species: str) -> pd.DataFrame:
    """Predictions joined to COG letters, compartment, topology, and a gene symbol for the panel."""
    pred = D.load(species)
    cog = F.load_cog(species)[["uniprot_ac", "cog_category", "cog_category_all"]]
    ident = P.load(species)[["uniprot_ac", "gene_name"]]
    egg = F.load_eggnog(species)[["uniprot_ac", "preferred_name"]]
    loc = LOC.load(species)  # compartment + cytoplasmic_fraction + has_signal_peptide
    topo = LOC.load_topology(species)[["uniprot_ac", "n_tm_helix", "n_tm_strand"]]
    df = (pred.merge(cog, on="uniprot_ac", how="left")
              .merge(ident, on="uniprot_ac", how="left")
              .merge(egg, on="uniprot_ac", how="left")
              .merge(loc, on="uniprot_ac", how="left")
              .merge(topo, on="uniprot_ac", how="left"))
    for col in ("cog_category", "cog_category_all", "gene_name", "preferred_name"):
        df[col] = df[col].fillna("").astype(str).str.strip()
    # gene_name first, eggNOG's preferred_name as fallback. Measured: this recovers one extra panel
    # gene on Kp (lpxC) and none at all on Sa, so it is a small improvement, not a fix -- Sa's
    # missing panel members are mostly real absences (no outer membrane) plus a naming gap.
    df["symbol"] = np.where(df.gene_name != "", df.gene_name, df.preferred_name)
    return I.annotate(df, gene_col="symbol")


def fisher_block(df: pd.DataFrame, activator: str, n: int, multi: bool) -> pd.DataFrame:
    """Every test for one species x activator, BH-corrected together."""
    prob = df[f"{activator}_prob"]
    in_top = prob.rank(ascending=False, method="first") <= n

    def members(col_values: pd.Series, key: str) -> pd.Series:
        return col_values == key

    tests = []

    # --- COG categories ---
    if multi:
        letters = sorted({ch for s in df.cog_category_all for ch in s} | {""})
        member_of = {
            lt: (df.cog_category_all.str.contains(lt, regex=False) if lt
                 else df.cog_category_all == "")
            for lt in letters
        }
    else:
        letters = sorted(set(df.cog_category))
        member_of = {lt: members(df.cog_category, lt) for lt in letters}

    names = cog_category_names()
    for lt in letters:
        tests.append(("cog", lt if lt else UNCLASSIFIED,
                      names.get(lt, UNCLASSIFIED if not lt else lt), member_of[lt]))

    # --- localization: the six compartments, then the TMbed features that bear on Clp reach ---
    # DeepLocPro always calls, so the compartments partition the proteome exactly; a compartment
    # absent from a species (cell_wall_surface in a Gram-negative, periplasm/outer_membrane in a
    # Gram-positive) is a structural zero and is skipped by the `a + b == 0` guard below.
    for cls in LOC.LOC_CLASSES:
        tests.append(("loc", LOC.LOC_CLASS_ABBREV[cls], cls, df.localization == cls))
    # Not compartments but the mechanics behind them: an export signal says the protein leaves the
    # cytoplasm, TM helices say it is buried in the membrane, and a beta-barrel is outer-membrane
    # by construction. These are the properties ClpP reachability actually turns on.
    tests.append(("loc", "signal peptide", "has_signal_peptide",
                  df.has_signal_peptide.fillna(False).astype(bool)))
    tests.append(("loc", "TM helix", "n_tm_helix > 0", df.n_tm_helix.fillna(0) > 0))
    tests.append(("loc", "beta-barrel", f"n_tm_strand >= {LOC.BETA_BARREL_MIN_STRANDS}",
                  df.n_tm_strand.fillna(0) >= LOC.BETA_BARREL_MIN_STRANDS))
    tests.append(("loc", "wholly cytoplasmic", "cytoplasmic_fraction == 1",
                  df.cytoplasmic_fraction.fillna(0) >= 1.0))

    # --- the consortium panel, family by family, and as a whole ---
    for fam in I.PANEL:
        tests.append(("panel", fam, fam, df.interest_family == fam))
    tests.append(("panel", "ALL PANEL", "all consortium targets", df.is_interest))

    rows = []
    for kind, key, label, is_member in tests:
        a = int((is_member & in_top).sum())
        b = int((is_member & ~in_top).sum())
        c = int((~is_member & in_top).sum())
        d = int((~is_member & ~in_top).sum())
        if a + b == 0:
            continue  # the group is absent from this proteome -- nothing to test
        odds, p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
        # Haldane-Anscombe, for plotting only; `odds_ratio` above stays uncorrected.
        adj = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append({
            "kind": kind, "key": key, "label": label,
            "a_group_top": a, "b_group_rest": b, "c_other_top": c, "d_other_rest": d,
            "group_n": a + b, "group_top_pct": 100 * a / (a + b),
            "background_top_pct": 100 * (a + c) / (a + b + c + d),
            "odds_ratio": odds, "log2_or_adj": float(np.log2(adj)), "p_value": p,
        })
    out = pd.DataFrame(rows)
    # BH **within each test family**, not across the block. The COG categories and the consortium
    # panel are separate questions asked of the same ranking; pooling them would make a COG
    # category's q-value depend on how many panel families happened to be tested alongside it.
    out["q_value"] = np.nan
    for kind_name, idx in out.groupby("kind").groups.items():
        out.loc[idx, "q_value"] = multipletests(out.loc[idx, "p_value"], method="fdr_bh")[1]
    out["significant"] = out.q_value < FDR
    return out


# ---------------------------------------------------------------- figures


def plot_fisher(ax, res: pd.DataFrame, activator: str, order: list[str], abc: str,
                ylabels: bool, kind: str, title: str,
                row_label=None) -> None:
    """log2 odds ratio per tested group, one point per species.

    Significance is encoded as **marker shape** (`o` significant, `x` not) rather than as fill:
    a hollow marker is nearly invisible once 26 rows are squeezed into one panel, and the
    significant/not distinction is the whole point of the figure.
    """
    from matplotlib.lines import Line2D

    sub = res[(res.activator == activator) & (res.kind == kind)]
    y = {k: i for i, k in enumerate(order)}
    for sp, color in SPECIES_COLOR.items():
        s = sub[(sub.species == sp) & sub.key.isin(y)]
        yy = np.array([y[k] for k in s.key])
        sig = s.significant.to_numpy()
        ax.scatter(s.log2_or_adj[sig], yy[sig], color=color, marker="o", zorder=3)
        ax.scatter(s.log2_or_adj[~sig], yy[~sig], color=color, marker="x", zorder=2)
    ax.axvline(0, color=NC.gray, linestyle="--")
    ax.set_yticks(range(len(order)))
    # Labels on the left panel only -- both panels share the row order, so repeating 26 long
    # labels would halve the space available for the data.
    if ylabels:
        ax.set_yticklabels([row_label(k) for k in order] if row_label else order)
    else:
        ax.set_yticklabels([])
    ax.set_ylim(-0.8, len(order) - 0.2)
    handles = [Line2D([], [], color=c, marker="o", linestyle="none", label=LABELS[s])
               for s, c in SPECIES_COLOR.items()]
    handles += [
        Line2D([], [], color=NC.black, marker="o", linestyle="none", label=f"FDR<{FDR}"),
        Line2D([], [], color=NC.black, marker="x", linestyle="none", label="not significant"),
    ]
    ax.legend(handles=handles, loc="lower right", ncol=2)
    stylia.label(ax, xlabel="log2 odds ratio", ylabel="",
                 title=f"{activator} · {title}", abc=abc)


def plot_interest(ax, data: dict, species: str, activators: list[str], abc: str) -> None:
    """Where each consortium target sits in the ranking, as a percentile, by family."""
    df = data[species]
    fams = [f for f in I.PANEL if (df.interest_family == f).any()]
    markers = {activators[0]: "o", activators[1]: "^"}
    for act in activators:
        pct = df[f"{act}_prob"].rank(pct=True) * 100
        for i, fam in enumerate(fams):
            m = df.interest_family == fam
            jitter = 0.16 if act == activators[1] else -0.16
            ax.scatter(pct[m], np.full(m.sum(), i + jitter),
                       marker=markers[act], color=NC.plum if act == activators[0] else NC.orange,
                       alpha=0.85, label=act if i == 0 else None, zorder=3)
    # Name the few that actually rank high. Without this the reader can see that *something* in
    # the Sec row reaches the 96th percentile but not that it is secA -- and which protein it is
    # the whole question.
    lead = activators[0]
    pct0 = df[f"{lead}_prob"].rank(pct=True) * 100
    marked = df.assign(pct=pct0)[df.is_interest].nlargest(3, "pct")
    # The top few are usually in the SAME family row at near-identical percentiles (secA 96,
    # secB 95, yajC 92 on Kp), so a fixed offset stacks them into an unreadable blur. Step the
    # vertical offset per label within a row instead.
    used: dict[str, int] = {}
    for r in marked.itertuples():
        k = used.get(r.interest_family, 0)
        used[r.interest_family] = k + 1
        ax.annotate(r.symbol, (r.pct, fams.index(r.interest_family) - 0.16),
                    textcoords="offset points", xytext=(4, 5 + 11 * k), color=NC.plum)

    ax.axvline(50, color=NC.gray, linestyle="--")
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels([f if f != "v5 proposal (cytosolic)" else "v5 (cytosolic)" for f in fams])
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.8, len(fams) - 0.2)
    if abc == "A":
        ax.legend(loc="upper left")
    stylia.label(ax, xlabel="Ranking percentile", ylabel="", title=LABELS[species], abc=abc)


# ---------------------------------------------------------------- main


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--species", nargs="+", default=list(D.SPECIES), choices=list(D.SPECIES))
    ap.add_argument("--top-pct", type=float, default=10.0,
                    help="'the top' as a percentage of each proteome (default: 10)")
    ap.add_argument("--top-n", type=int,
                    help="override --top-pct with an absolute list length, the same for every "
                         "species")
    ap.add_argument("--multi", action="store_true",
                    help="count a protein in every COG letter it carries, not just the first")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.mkdir(parents=True, exist_ok=True)

    order = [s for s in LABELS if s in args.species]
    data = {s: frame(s) for s in order}
    activators = list(D.ACTIVATORS)

    # "The top" is a FRACTION of each proteome by default, not a fixed count: the three differ by
    # 2x in size (5,728 / 4,403 / 2,889), so one absolute N would mean the top 1.7% of Kp against
    # the top 3.5% of Sa and the odds ratios would not be comparable across species.
    top_n = {
        sp: (args.top_n if args.top_n else int(round(len(data[sp]) * args.top_pct / 100)))
        for sp in order
    }
    say("  THE TOP" + ("  (absolute)" if args.top_n else f"  (top {args.top_pct:g}%)"))
    for sp in order:
        say(f"     {sp:<14} {top_n[sp]:>5} of {len(data[sp]):>5} proteins")

    blocks = []
    for sp in order:
        for act in activators:
            b = fisher_block(data[sp], act, top_n[sp], args.multi)
            b.insert(0, "species", sp)
            b.insert(1, "activator", act)
            b.insert(2, "top_n", top_n[sp])
            blocks.append(b)
    res = pd.concat(blocks, ignore_index=True)
    out_tsv = RESULT_DIR / "enrichment_fisher.tsv"
    res.to_csv(out_tsv, sep="\t", index=False)
    say(f"  -> {out_tsv.relative_to(REPO_ROOT)}   ({len(res)} tests, "
        f"{int(res.significant.sum())} significant at FDR {FDR})")

    # Panel coverage first: a missing gene is not an absent target.
    say("\n  CONSORTIUM PANEL - coverage (gene symbol match)")
    for sp in order:
        cov = I.coverage(data[sp], gene_col="symbol")
        say(f"     {sp:<14} {int(cov.n_found.sum())}/{int(cov.n_panel.sum())} panel genes found")
        for _, r in cov[cov.n_found < cov.n_panel].iterrows():
            say(f"        {r.family:<26} {r.n_found}/{r.n_panel}   missing: {r.missing}")

    # The COG figure, ordered by mean effect so the reader's eye finds the extremes.
    cog = res[res.kind == "cog"]
    mean_or = cog.groupby("key").log2_or_adj.mean().sort_values()
    keys = [k for k in mean_or.index]
    # Height is set explicitly here, against the usual rule: 26 category rows in a default-height
    # slide panel leaves the tick labels overlapping and unreadable.
    fig, axs = stylia.create_figure(1, len(activators), height=0.62)
    names = cog_category_names()
    for i, (abc, act) in enumerate(zip("AB", activators)):
        plot_fisher(axs.next(), res, act, keys, abc, ylabels=(i == 0), kind="cog",
                    title="COG enrichment",
                    row_label=lambda k: k if k == UNCLASSIFIED else f"{k}  {names.get(k, '')[:30]}")
    out = PLOT_DIR / "cog_fisher.png"
    stylia.save_figure(str(out))
    say(f"\n  -> {out.relative_to(REPO_ROOT)}")

    # The localization figure. Only ~10 rows, so the default height is right.
    loc = res[res.kind == "loc"]
    loc_order = [k for k in loc.groupby("key").log2_or_adj.mean().sort_values().index]
    fig, axs = stylia.create_figure(1, len(activators))
    for i, (abc, act) in enumerate(zip("AB", activators)):
        plot_fisher(axs.next(), res, act, loc_order, abc, ylabels=(i == 0), kind="loc",
                    title="localization enrichment")
    out = PLOT_DIR / "loc_fisher.png"
    stylia.save_figure(str(out))
    say(f"\n  -> {out.relative_to(REPO_ROOT)}")
    for act in activators:
        say(f"     {act}:")
        for _, r in loc[(loc.activator == act)].sort_values("log2_or_adj").iterrows():
            flag = "*" if r.significant else " "
            say(f"       {flag} {r.species:<12} {r.key:<19} OR {r.odds_ratio:>7.2f}  "
                f"q {r.q_value:>8.2e}   {r.a_group_top:>4}/{r.group_n:<5} in top "
                f"{top_n[r.species]}   ({r.label})")

    names = cog_category_names()
    for act in activators:
        sig = cog[(cog.activator == act) & cog.significant].sort_values("log2_or_adj")
        say(f"     {act}: {len(sig)} significant category-species tests")
        for _, r in sig.iterrows():
            arrow = "up  " if r.log2_or_adj > 0 else "down"
            say(f"        {arrow} {r.species:<12} {r.key:<13} OR {r.odds_ratio:>6.2f}  "
                f"q {r.q_value:.2e}   {r.a_group_top:>3}/{r.group_n:<5} in top {top_n[r.species]}"
                f"   {names.get(r.key, UNCLASSIFIED)[:44]}")

    # The panel figure.
    fig, axs = stylia.create_figure(1, len(order))
    for abc, sp in zip("ABC", order):
        plot_interest(axs.next(), data, sp, activators, abc)
    out = PLOT_DIR / "interest_panel.png"
    stylia.save_figure(str(out))
    say(f"\n  -> {out.relative_to(REPO_ROOT)}")

    say("\n  CONSORTIUM PANEL - rank and enrichment")
    for sp in order:
        df = data[sp]
        for act in activators:
            pct = df[f"{act}_prob"].rank(pct=True) * 100
            it = pct[df.is_interest]
            r = res[(res.species == sp) & (res.activator == act)
                    & (res.key == "ALL PANEL")].iloc[0]
            say(f"     {sp:<14} {act:<7} n={int(r.group_n):<3} median percentile "
                f"{it.median():>5.1f}   {r.a_group_top}/{int(r.group_n)} in top {top_n[sp]} "
                f"(expected {r.background_top_pct * r.group_n / 100:.1f})   "
                f"OR {r.odds_ratio:.2f}  q {r.q_value:.3f}"
                f"{'  SIGNIFICANT' if r.significant else ''}")
    # Which individual targets actually rank high -- the question a shortlist asks.
    say("\n  Highest-ranked panel members")
    for sp in order:
        df = data[sp]
        for act in activators:
            df = df.assign(pct=df[f"{act}_prob"].rank(pct=True) * 100)
            t = df[df.is_interest].nlargest(5, "pct")
            say(f"     {sp:<14} {act:<7} " + ", ".join(
                f"{r.symbol} {r.pct:.0f}" for r in t.itertuples()))


if __name__ == "__main__":
    main()
