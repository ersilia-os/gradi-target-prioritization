"""E. coli x ASKA: are the proteins that pass the degradability gate really higher-scoring than chance?

The gate is `ASKA AND cytosolic AND (degradability_consensus >= CUT OR (a measured S. aureus ortholog
is a degrader AND consensus >= the proteome median))` (see `ecoli_aska_analysis.py`). This draws the
three scores of the proteins that pass it against the same number of proteins picked AT RANDOM from
the pool the gate was applied to (ASKA and cytosolic) that did NOT pass, so no protein sits in both
groups. A second random draw from the rest of the whole proteome is printed, not drawn.

Writes   output/plots/presentation/ecoli_aska_degradability.png
Prints   medians, the share above each activator's base-rate cut, the AUROC of "passes the gate",
         and a 1,000-draw check that the gap is not one lucky sample.

**Read it for what it is.** The gate is built from these scores, so the passing group is higher by
construction: this shows HOW MUCH higher and where the ortholog route reaches below the cut, not
that the rule works. E. coli scores are extrapolated from S. aureus labels (ranking hypotheses).

`stylia.label(..., xlabel=None)` writes the placeholder "X-axis / Units"; pass "" for no label.
`create_figure(width=, height=)` takes fractions of the format size, not inches. Run plot scripts
ONE AT A TIME (stylia deletes the matplotlib font cache at import).

Run with the `gradi` env:  python plotting/ecoli_aska_degradability.py [-q] [--seed 0]
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import numpy as np  # noqa: E402
import stylia  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

from plotting import palette as PAL  # noqa: E402
from src import degradability as D  # noqa: E402

stylia.set_format("slide")
stylia.set_style("article")
SS = stylia.SLIDE_FONTSIZE_SMALL
OUT = REPO_ROOT / "output" / "plots" / "presentation" / "ecoli_aska_degradability.png"
N_DRAWS = 1000
PANELS = [("degradability_consensus", "Degradability consensus", None),
          ("adep4_prob", "ADEP4 probability", D.BASE_RATE_THRESHOLD["adep4"]),
          ("onc212_prob", "ONC212 probability", D.BASE_RATE_THRESHOLD["onc212"])]


def load():
    spec = importlib.util.spec_from_file_location("ecoli_aska_analysis", REPO_ROOT / "plotting/ecoli_aska_analysis.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.load_ecoli()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    say = (lambda *_: None) if a.quiet else print

    d = load()
    pool = d[d.in_aska & d.cytosolic]
    hit = pool[pool.degradable]
    n = len(hit)
    rng = np.random.default_rng(a.seed)
    # The random draw EXCLUDES the proteins that passed the gate: it is a sample of what was not
    # selected, so a protein can never sit in both groups.
    rest = pool[~pool.degradable]
    rest_all = d[~d.uniprot_ac.isin(hit.uniprot_ac)]
    rand = rest.iloc[rng.choice(len(rest), n, replace=False)]
    rand_all = rest_all.iloc[rng.choice(len(rest_all), n, replace=False)]
    assert not set(rand.uniprot_ac) & set(hit.uniprot_ac) and not set(rand_all.uniprot_ac) & set(hit.uniprot_ac)
    say(f"pool (ASKA and cytosolic) {len(pool):,} | passes the gate {n:,} | not selected {len(rest):,} | CUT {d.attrs['cons_deg_cut']:.4f} | "
        f"ortholog floor {d.attrs['orth_floor']:.4f} | seed {a.seed}")

    fig, axs = stylia.create_figure(1, 3, width_ratios=[1, 1, 1], width=1.0, height=0.36)
    for col, title, cut in PANELS:
        ax = axs.next()  # stylia hands out axes with .next(); zip() over `axs` skips its labelling
        bins = np.linspace(0, 1, 26)
        ax.hist(rand[col], bins=bins, color=PAL.MUTED, alpha=0.75, density=True, label=f"Random {n:,} of the {len(rest):,} not selected")
        ax.hist(hit[col], bins=bins, color=PAL.PRIMARY, alpha=0.75, density=True, label=f"Pass the gate ({n:,})")
        ax.axvline(hit[col].median(), color=PAL.PRIMARY, lw=1.6)
        ax.axvline(rand[col].median(), color=PAL.MUTED, lw=1.6, ls="--")
        if cut is not None:
            ax.axvline(cut, color=PAL.INK, lw=1.0, ls=":")
        stylia.label(ax, xlabel=title, ylabel="Density" if col == PANELS[0][0] else "", title="")
        ax.set_xlim(0, 1)
    first = fig.axes[0]
    first.legend(fontsize=SS, frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(OUT))
    say(f"wrote {OUT.relative_to(REPO_ROOT)}\n")

    y = np.r_[np.ones(n), np.zeros(len(rest))]
    say(f"{'score':26s}{'median pass':>12s}{'median random':>14s}{'median proteome-random':>24s}{'AUROC pass vs rest':>20s}{'>= base-rate cut: pass / random':>34s}")
    for col, title, cut in PANELS:
        auc = roc_auc_score(y, np.r_[hit[col], rest[col]])
        share = f"{100 * (hit[col] >= cut).mean():.0f}% / {100 * (rand[col] >= cut).mean():.0f}%" if cut else "n/a (a percentile)"
        say(f"{title:26s}{hit[col].median():12.3f}{rand[col].median():14.3f}{rand_all[col].median():24.3f}{auc:20.3f}{share:>34s}")
    meds = np.array([rest.iloc[rng.choice(len(rest), n, replace=False)]["degradability_consensus"].median() for _ in range(N_DRAWS)])
    say(f"\nconsensus median over {N_DRAWS} random draws of {n}: {meds.mean():.3f} (range {meds.min():.3f}-{meds.max():.3f}); "
        f"observed {hit['degradability_consensus'].median():.3f}; draws at or above it: {int((meds >= hit['degradability_consensus'].median()).sum())}")
    by = hit.groupby("deg_route")
    say("\nby route (consensus median / ADEP4 median / ONC212 median / n):")
    for r, g in by:
        say(f"  {r:16s}{g.degradability_consensus.median():8.3f}{g.adep4_prob.median():8.3f}{g.onc212_prob.median():8.3f}{len(g):6d}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
