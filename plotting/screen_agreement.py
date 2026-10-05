"""Three E. coli essentiality screens: what they measured, and how far they agree.

    screen_agreement.png   1  the continuous readout, split by whether the OTHER screen agreed
                           2  pairwise agreement between the three screens
                           3  how many screens call each essential protein essential

**E. coli, not K. pneumoniae, and that is forced.** All three Kp TraDIS screens key onto three
DIFFERENT identifier spaces -- RefSeq `WP_`, EMBL `CCN`, locus tag `KPNRH_` -- and **0 of ~5,000
join to the anchor proteome**, so a per-protein comparison is impossible without the DIAMOND bridge
that `strain_homologs.py` builds and does not persist. The three E. coli screens key on UniProt
accessions and join **100%**. E. coli is also where the measurements actually are: Kp's entire
experimental tier is proxy (see `consensus_evidence.png`).

The three, from `data/processed/essentiality/training_sets/`:

    Keio      arrayed single-gene knockout, K-12   4,190 proteins, 286 essential (base 0.068)
    Goodall   TraDIS, BW25113                      4,056 proteins, 354 essential (base 0.087)
    Choe      Tn-seq, BW25113                      4,272 proteins, 440 essential (base 0.103)

**THE HEADLINE: experimental screens disagree a lot.** Of the 495 proteins called essential by ANY
of the three on the 3,963 all three measured, **only 214 (43%) are called by all three**; 194 are
called by exactly one. Pairwise Jaccard runs **0.443 to 0.727**. "Measured" is not one thing, which
is why `<axis>_evidence` counts sources and concordance rather than treating a single screen as
truth.

**The disagreement is not noise, and the continuous readout proves it.** Goodall ships
`insertion_index` -- the density of transposon insertions, so LOW means essential. Split by Keio's
independent binary call, the medians are **0.0031 (Keio essential) against 0.1248 (not)**, a 40x
separation. The disputed proteins sit exactly where they should: the 90 that only Goodall calls
essential have a median index of **0.0083**, nearly as insertion-free as the agreed essentials,
which says Keio most likely MISSED them rather than Goodall over-calling. The 9 that only Keio calls
sit at **0.0503**, intermediate.

**Why Keio and Goodall agree best (kappa 0.828) while Choe agrees least with both.** Not a quality
ranking -- a transposon screen cannot see a gene whose disruption is merely costly the same way an
arrayed knockout can, and Choe calls the most genes essential (440 against 286). Read the pairwise
panel as "these assays measure overlapping but different things".

Everything is read from the stage's own training sets through `src.essentiality.load_training_set()`.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/` -- this file is
one level deep. Copying `parents[2]` resolves REPO_ROOT to the repo's PARENT and the ImportError
names `src`, not the path, so it reads as a broken conda env.

Colours come from `plotting/palette.py` (stylia's **npg** palette), never from
`stylia.NamedColors()`. `stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units".
`stylia.create_figure(width=, height=)` takes FRACTIONS of the format size, not inches.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/screen_agreement.py
"""

from __future__ import annotations

import argparse
import itertools
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import stylia  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import essentiality as E  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"

# Format: slide | Style: article. "article" not "ersilia" because the ersilia style paints
# every text element and spine plum (#50285A); article gives black. The COLOURS come from
# `plotting/palette.py` -- stylia's npg palette.
stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
SPECIES = "ecoli"

#: short name -> training-set id. All three key on uniprot_ac and join 100% to the anchor.
SCREENS = {
    "Keio": "essential_ecoli_k12_knockout",
    "Goodall": "essential_ecoli_bw25113_tradis_goodall",
    "Choe": "essential_ecoli_bw25113_tnseq_choe",
}
ASSAY = {"Keio": "arrayed knockout", "Goodall": "TraDIS", "Choe": "Tn-seq"}

#: The screen carrying a continuous readout, and its column. LOW = essential.
CONTINUOUS, CONT_COL, CONT_PARTNER = "Goodall", "insertion_index", "Keio"

#: A floor so zeros survive a log axis. 221 proteins sit at exactly 0 insertions -- the most
#: essential signal there is -- and dropping them would delete the strongest evidence in the panel.
FLOOR = 1e-3

AGREE_CLASSES = [
    ("both", True, True, PAL.NPG[4]),
    (f"{CONTINUOUS} only", True, False, PAL.NPG[2]),
    (f"{CONT_PARTNER} only", False, True, PAL.NPG[1]),
    ("neither", False, False, PAL.MUTED),
]


def load_wide() -> pd.DataFrame:
    """One row per anchor protein, one column per screen, NaN where that screen did not measure it.

    A missing value is "this screen had no call", never a measured non-essential, so every
    comparison below drops pairwise rather than filling."""
    w = P.load(SPECIES)[["uniprot_ac"]].copy()
    for short, name in SCREENS.items():
        t = E.load_training_set(name)[["key", "label"]].rename(
            columns={"key": "uniprot_ac", "label": short})
        w = w.merge(t, on="uniprot_ac", how="left")
    cont = E.load_training_set(SCREENS[CONTINUOUS])[["key", CONT_COL]].rename(
        columns={"key": "uniprot_ac"})
    return w.merge(cont, on="uniprot_ac", how="left")


def plot_continuous(ax, w: pd.DataFrame) -> pd.DataFrame:
    """The TraDIS insertion index, split by the four Keio-vs-Goodall agreement classes.

    This is the panel that shows the disagreement is structured rather than random: a protein only
    one screen calls essential still has an insertion index near the agreed essentials, not near
    the agreed non-essentials. Log axis because the classes separate 40-fold and a linear axis
    presses every essential against zero."""
    d = w.dropna(subset=[CONTINUOUS, CONT_PARTNER, CONT_COL]).copy()
    d["_c"] = pd.to_numeric(d[CONT_COL], errors="coerce").clip(lower=FLOOR)
    rows = []
    for i, (label, a, b, color) in enumerate(AGREE_CLASSES):
        v = d.loc[d[CONTINUOUS].astype(bool).eq(a) & d[CONT_PARTNER].astype(bool).eq(b), "_c"]
        rows.append({"klass": label, "n": len(v), "median": float(v.median())})
        if v.empty:
            continue
        parts = ax.violinplot([np.log10(v)], positions=[i], widths=0.78, showextrema=False)
        for body in parts["bodies"]:
            body.set_alpha(1.0)
            body.set_facecolor(color)
            body.set_edgecolor("white")
            body.set_linewidth(0.9)
        ax.scatter([i], [np.log10(v.median())], s=16, color=PAL.INK, zorder=5)

    ax.set_xticks(range(len(AGREE_CLASSES)))
    ax.set_xticklabels([f"{label}\nn={r['n']:,}" for (label, *_r), r in zip(AGREE_CLASSES, rows)],
                       fontsize=SS * 0.78)
    ticks = [0.001, 0.01, 0.1, 0.5]
    ax.set_yticks(np.log10(ticks))
    ax.set_yticklabels([str(t) for t in ticks], fontsize=SS)
    ax.set_xlim(-0.7, len(AGREE_CLASSES) - 0.3)
    stylia.label(ax, xlabel="", ylabel=f"{CONTINUOUS} insertion index (low = essential)",
                 title="Disagreement is structured, not noise")
    return pd.DataFrame(rows)


def _pair_stats(w: pd.DataFrame, a: str, b: str) -> dict:
    """Jaccard and Cohen's kappa on the proteins BOTH screens measured.

    Jaccard because the classes are wildly unbalanced and plain agreement is ~0.95 whatever
    happens; kappa beside it because Jaccard ignores the agreed negatives entirely and kappa is the
    standard chance-corrected number a reviewer will look for."""
    s = w.dropna(subset=[a, b])
    A, B = s[a].astype(bool), s[b].astype(bool)
    tp, fa, fb, tn = int((A & B).sum()), int((A & ~B).sum()), int((~A & B).sum()), int((~A & ~B).sum())
    n = len(s)
    po = (tp + tn) / n
    pe = ((tp + fa) * (tp + fb) + (fb + tn) * (fa + tn)) / n**2
    return {"pair": f"{a} ~ {b}", "n": n, "both": tp, f"only_{a}": fa, f"only_{b}": fb,
            "jaccard": tp / max(tp + fa + fb, 1), "kappa": (po - pe) / (1 - pe), "agreement": po}


def plot_pairwise(ax, w: pd.DataFrame) -> list[dict]:
    """Jaccard and kappa for each pair. Two bars, because they answer different questions."""
    pairs = list(itertools.combinations(SCREENS, 2))
    stats = [_pair_stats(w, a, b) for a, b in pairs]
    y = np.arange(len(pairs))[::-1]
    height = 0.36
    ax.barh(y + height / 2, [s["jaccard"] for s in stats], height=height, color=PAL.PRIMARY,
            label="Jaccard (essentials only)")
    ax.barh(y - height / 2, [s["kappa"] for s in stats], height=height, color=PAL.TERTIARY,
            label="Cohen's kappa")
    for yi, s in zip(y, stats):
        ax.text(s["jaccard"] + 0.015, yi + height / 2, f"{s['jaccard']:.2f}", va="center",
                fontsize=SS * 0.78, color=PAL.INK)
        ax.text(s["kappa"] + 0.015, yi - height / 2, f"{s['kappa']:.2f}", va="center",
                fontsize=SS * 0.78, color=PAL.INK)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{a} ~ {b}" for a, b in pairs], fontsize=SS * 0.85)
    ax.set_xlim(0, 1)
    ax.legend(fontsize=SS * 0.78, frameon=False, loc="lower right", handletextpad=0.5)
    stylia.label(ax, xlabel="Agreement", ylabel="", title="No two screens agree closely")
    return stats


def plot_unanimity(ax, w: pd.DataFrame) -> pd.Series:
    """Of the proteins any screen calls essential, how many screens agree.

    Restricted to the proteins ALL THREE measured, or the counts would mix "no screen called it"
    with "no screen looked" -- the same distinction `<axis>_evidence` exists to keep."""
    both = w.dropna(subset=list(SCREENS))
    n_calls = both[list(SCREENS)].sum(axis=1).astype(int)
    counts = n_calls[n_calls > 0].value_counts().sort_index()
    colors = [PAL.NPG[0], PAL.NPG[2], PAL.NPG[4]]
    x = np.arange(1, 4)
    vals = [int(counts.get(i, 0)) for i in x]
    ax.bar(x, vals, width=0.6, color=colors)
    total = sum(vals)
    for xi, v in zip(x, vals):
        ax.text(xi, v + total * 0.015, f"{v:,}\n{100 * v / total:.0f}%", ha="center", va="bottom",
                fontsize=SS * 0.85, color=PAL.INK)
    ax.set_xticks(x)
    ax.set_xticklabels(["1 screen", "2 screens", "all 3"], fontsize=SS)
    ax.set_ylim(0, max(vals) * 1.35)
    ax.set_xlim(0.4, 3.6)
    stylia.label(ax, xlabel="", ylabel="Proteins called essential",
                 title=f"Only {100 * vals[-1] / total:.0f}% are unanimous")
    return counts


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    say("\n  experimental screen agreement -- E. coli")
    w = load_wide()
    say(f"  {len(w):,} anchor proteins; measured by "
        + ", ".join(f"{k} {int(w[k].notna().sum()):,}" for k in SCREENS))

    fig, axs = stylia.create_figure(1, 3, width_ratios=[1.25, 1.1, 0.85], width=1.0, height=0.40)
    cont = plot_continuous(axs.next(), w)
    stats = plot_pairwise(axs.next(), w)
    counts = plot_unanimity(axs.next(), w)
    out = OUT_DIR / "screen_agreement.png"
    stylia.save_figure(str(out))
    say(f"  -> {out.relative_to(REPO_ROOT)}")

    say("\n  THE THREE SCREENS")
    reg = E.registry()
    for short, name in SCREENS.items():
        r = reg[reg["dataset_id"] == name]
        if len(r):
            r = r.iloc[0]
            say(f"    {short:<9} {ASSAY[short]:<18} {int(r['n_labels']):>5,} proteins  "
                f"{int(r['n_positives']):>4} essential  base {float(r['base_rate']):.4f}")

    say(f"\n  {CONTINUOUS.upper()} INSERTION INDEX BY AGREEMENT CLASS (low = essential)")
    for _, r in cont.iterrows():
        say(f"    {r['klass']:<16} n={int(r['n']):>5,}   median {r['median']:.4f}")
    say("    The disputed proteins sit between the agreed classes, not with the non-essentials:")
    say("    a protein only one screen calls essential is still nearly insertion-free.")

    say("\n  PAIRWISE AGREEMENT")
    for s in stats:
        only = {k: v for k, v in s.items() if k.startswith("only_")}
        say(f"    {s['pair']:<20} n={s['n']:>5,}  both={s['both']:>4}  "
            + "  ".join(f"{k.replace('only_', 'only ')}={v}" for k, v in only.items())
            + f"   Jaccard {s['jaccard']:.3f}   kappa {s['kappa']:.3f}   raw agreement {s['agreement']:.4f}")
    say("    Raw agreement is ~0.94-0.98 for all three pairs and is NOT the number to quote:")
    say("    essentials are <11% of a proteome, so agreeing on the negatives is nearly free.")

    total = int(counts.sum())
    say(f"\n  UNANIMITY  (of the {total:,} proteins any screen called essential, among the "
        f"{len(w.dropna(subset=list(SCREENS))):,} all three measured)")
    for i in (1, 2, 3):
        v = int(counts.get(i, 0))
        say(f"    called by {i} screen(s): {v:>4}  ({100 * v / total:.1f}%)")

    say("\n  CAVEATS")
    say("    - E. coli, not K. pneumoniae, and not by choice: the three Kp screens key onto three")
    say("      different id spaces and 0 of ~5,000 join to the anchor proteome.")
    say("    - A blank is 'this screen had no call', never a measured non-essential. Every")
    say("      comparison drops pairwise; the unanimity panel uses only the 3-way measured set.")
    say("    - Kappa differences are NOT a quality ranking. A transposon screen cannot see a gene")
    say("      whose loss is merely costly the way an arrayed knockout can, and Choe calls the")
    say("      most genes essential (440 vs Keio's 286).")
    say("    - This is why the evidence ladder counts SOURCES and CONCORDANCE rather than")
    say("      treating any one screen as truth. See consensus_evidence.png.")


if __name__ == "__main__":
    main()
