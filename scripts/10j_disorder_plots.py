"""2x3 disorder slide for the degradability axis (docs §3.1) — E. coli K-12.

House style: stylia "slide" format, NPG palette, single organism per slide.

  1  disorder anatomy        where along the chain disorder actually sits, per terminus
  2  engageability           % of proteome with a terminal run >= L, against 5 / 20 / 37 aa
  3  which terminus          N vs C initiation regions, with the named Gr-ADI targets
  4  architecture            the five disorder shapes, by compartment (i.e. by reachable handle)
  5  internal degrons        disorder that is neither N- nor C-terminal — the FtsZ route
  6  ligandable + degradable the bimodal sweet spot: folded core AND engageable tail

The mechanical thresholds in panels 1-3 are measurements, not conventions: ~5 aa for closed-channel
ClpXP recognition (Saunders 2020), ~20 aa for open-channel engagement (Saunders 2020; Fei 2020), and
~37 aa to reach the ClpP active site from the pore (Kenniston 2005). They are what turns "this
protein is disordered" into "this protein can be threaded".

Reads <prefix>_disorder.csv (10d) + localization (09g) + ligandability (06g).
Output: output/plots/10j_disorder_<prefix>.png. Run with `gradi`.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# stylia does an unconditional `shutil.rmtree(mpl.get_cachedir())` at import time. The usual
# `makedirs(get_cachedir())` guard is not enough: matplotlib silently falls back to a temporary
# cache dir when the default is unwritable, so the path we create and the path stylia deletes can
# differ, and the import dies with FileNotFoundError. Pin the cache to a stable repo-local path
# before matplotlib is imported, then re-create it after stylia has wiped it.
_MPL_CACHE = Path(__file__).resolve().parents[1] / "tmp" / "mplcache"
_MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE))

import matplotlib  # noqa: E402
matplotlib.use("Agg")
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import stylia  # noqa: E402
os.makedirs(matplotlib.get_cachedir(), exist_ok=True)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import degradability as D  # noqa: E402
from src import localization as LOC  # noqa: E402

NPG = stylia.CategoricalPalette("npg").colors
SS = stylia.SLIDE_FONTSIZE_SMALL

ORG_COLOR = {"kp": "#E64B35", "ec": "#4DBBD5"}
NTERM_C, CTERM_C = "#3C5488", "#E64B35"
HIT, MISS, GREY, BG = "#00A087", "#C9C9C7", "#555555", "#D8D8D6"
ARCH_COLOR = {"ordered": "#C9C9C7", "internal_linker": "#8491B4", "terminal_idr": "#00A087",
              "multi_idr": "#3C5488", "mostly_disordered": "#E64B35"}
ARCH_ORDER = ["ordered", "internal_linker", "terminal_idr", "multi_idr", "mostly_disordered"]

# Named proteins, so the slide can be read against known answers rather than only as a distribution.
# GroEL is the worked negative control: CLIPPERs depleted it only ~40% because both termini face the
# barrel interior and it is among the most abundant proteins in the cell.
CONTROLS = {"GroEL": ("groL", "groEL"), "AcpP": ("acpP",), "DnaK": ("dnaK",),
            "GyrA": ("gyrA",), "GyrB": ("gyrB",), "FtsZ": ("ftsZ",), "RpoS": ("rpoS",)}
NEGATIVE = {"GroEL"}
PROFILE_LEN = 60
ENGAGE_AA = D.INIT_THRESHOLDS[1]


def _save(out: Path) -> None:
    """stylia.save_figure with real padding.

    `stylia.save_figure(pad=...)` accepts a pad argument but ignores it — both branches call a bare
    `plt.tight_layout()` — so the 2x3 grid ends up with the matplotlib default 1.08 inter-panel pad,
    which is a lot of white at slide scale. Same dpi/bbox/facecolor as stylia so the output is
    otherwise identical.
    """
    import matplotlib.pyplot as plt
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(str(out), dpi=600, transparent=False, bbox_inches="tight")


def _note(ax, msg):
    ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", color="#999",
            fontsize=SS)
    ax.set_xticks([]); ax.set_yticks([])


def _read(org, prefix, name, cols):
    f = D.results_dir(org) / f"{prefix}_{name}.csv"
    if not f.exists():
        return None
    have = pd.read_csv(f, nrows=0).columns
    use = [c for c in cols if c in have]
    return pd.read_csv(f, usecols=use) if "uniprot_accession" in use else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(D.ORGANISMS), default="ecoli")
    args = ap.parse_args()
    org = args.organism
    _, prefix = D.ORGANISMS[org]
    orgname = D.ORG_DISPLAY[org]
    col = ORG_COLOR[prefix]

    d = pd.read_csv(D.results_dir(org) / f"{prefix}_disorder.csv")
    d["gene"] = d["uniprot_accession"].map(D.load_genes(org))
    for name, cols in (("localization", ["uniprot_accession", "localization",
                                         "has_signal_peptide"]),):
        t = _read(org, prefix, name, cols)
        if t is not None:
            d = d.merge(t, on="uniprot_accession", how="left")
    n = len(d)

    by_gene = {}
    sym = d["gene"].astype(str).str.lower()
    for label, aliases in CONTROLS.items():
        hit = d[sym.isin([a.lower() for a in aliases])]
        if not hit.empty:
            by_gene[label] = hit.iloc[0]

    stylia.set_format("slide")
    fig, axs = stylia.create_figure(2, 3, width=1.0, height=0.5625)

    # ---- panel 1: where along the chain disorder sits ----
    # Per position from each terminus, the % of proteins whose residue there is below pLDDT 70.
    # The interior baseline makes the terminal enrichment readable as a multiple, not just a shape.
    ax = axs.next()
    n_hits = np.zeros(PROFILE_LEN); c_hits = np.zeros(PROFILE_LEN); n_tot = 0
    interior_dis = 0.0; interior_len = 0
    for acc in d["uniprot_accession"]:
        p = D.plddt_series(org, acc)
        if len(p) < 2 * PROFILE_LEN + 20:
            continue
        n_tot += 1
        arr = np.asarray(p)
        n_hits += (arr[:PROFILE_LEN] < D.PLDDT_DISORDER_CUT)
        c_hits += (arr[-PROFILE_LEN:][::-1] < D.PLDDT_DISORDER_CUT)
        core = arr[PROFILE_LEN:-PROFILE_LEN]
        interior_dis += float((core < D.PLDDT_DISORDER_CUT).sum()); interior_len += len(core)
    xs = np.arange(1, PROFILE_LEN + 1)
    base = 100 * interior_dis / max(interior_len, 1)
    ax.plot(xs, 100 * n_hits / n_tot, color=NTERM_C, lw=1.6, label="N-terminus")
    ax.plot(xs, 100 * c_hits / n_tot, color=CTERM_C, lw=1.6, label="C-terminus")
    ax.axhline(base, color=GREY, lw=0.9, ls="--")
    ax.text(2, base + 1.2, f"interior baseline {base:.1f}%", fontsize=SS, color=GREY,
            va="bottom", ha="left")
    ax.set_xlim(1, PROFILE_LEN)
    ax.text(0.97, 0.95, f"residue 1–5 is {100*n_hits[:5].mean()/n_tot/base:.1f}× "
                        f"the interior rate", transform=ax.transAxes, ha="right", va="top",
            fontsize=SS, color=GREY)
    stylia.label(ax, xlabel="residues from terminus", ylabel="% of proteins disordered (pLDDT < 70)",
                 title=f"Where disorder sits — {orgname}")
    ax.legend(fontsize=SS, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.86))

    # ---- panel 2: how many proteins clear the mechanical thresholds ----
    ax = axs.next()
    Ls = np.arange(0, 61)
    for t, lab, c, ls in (("n", "N-terminus", NTERM_C, "-"), ("c", "C-terminus", CTERM_C, "-")):
        v = d[f"{t}term_init_len_70"].dropna().to_numpy()
        ax.plot(Ls, [(v >= L).mean() * 100 for L in Ls], ls, color=c, lw=1.6, label=lab)
    either = ((d["nterm_init_len_70"] >= ENGAGE_AA) | (d["cterm_init_len_70"] >= ENGAGE_AA))
    ax.plot(Ls, [((d["nterm_init_len_70"] >= L) | (d["cterm_init_len_70"] >= L)).mean() * 100
                 for L in Ls], color="#1A1A1A", lw=1.2, ls=":", label="either terminus")
    # Only the numbers go on the lines; the mechanism they correspond to goes in the axis label,
    # so three labels cannot collide with each other or with the legend.
    for thr in D.INIT_THRESHOLDS:
        ax.axvline(thr, color=GREY, lw=0.7, ls=":")
        ax.text(thr + 0.8, 98, f"{thr}", fontsize=SS, color=GREY, va="top")
    ax.text(0.97, 0.52, f"{int(either.sum()):,} of {n:,} ({either.mean():.1%}) clear "
                        f"{ENGAGE_AA} aa\nat either terminus",
            transform=ax.transAxes, ha="right", va="top", fontsize=SS, color=HIT,
            fontweight="bold")
    ax.set_ylim(0, 100); ax.set_xlim(0, 60)
    stylia.label(ax, xlabel="initiation-region length L (aa)\n"
                            "5 = recognition · 20 = engagement · 37 = reach to ClpP active site",
                 ylabel="% of proteome with a run ≥ L",
                 title=f"How many can be threaded? — {orgname}")
    ax.legend(fontsize=SS, frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.80))

    # ---- panel 3: which terminus ----
    ax = axs.next()
    s = d.dropna(subset=["nterm_init_len_70", "cterm_init_len_70"])
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("dens", ["#FFFFFF", col])
    ax.hexbin(s["nterm_init_len_70"].clip(upper=60), s["cterm_init_len_70"].clip(upper=60),
              gridsize=30, bins="log", cmap=cmap, mincnt=1, extent=(0, 60, 0, 60))
    ax.axvline(ENGAGE_AA, color=GREY, lw=0.8, ls="--")
    ax.axhline(ENGAGE_AA, color=GREY, lw=0.8, ls="--")
    for g, row in by_gene.items():
        x = min(row.get("nterm_init_len_70") or 0, 60)
        y = min(row.get("cterm_init_len_70") or 0, 60)
        neg = g in NEGATIVE
        ax.scatter([x], [y], s=34, marker="X" if neg else "*",
                   color="#B2182B" if neg else "#1A1A1A", zorder=5, linewidths=0)
        right = x > 45
        ax.annotate(g + (" (neg)" if neg else ""), (x, y), fontsize=SS,
                    color="#B2182B" if neg else "#333",
                    xytext=(-4 if right else 3, 3), textcoords="offset points",
                    ha="right" if right else "left")
    vc = d["engageable_terminus"].value_counts()
    ax.text(0.96, 0.94, f"N only {int(vc.get('n', 0))}   C only {int(vc.get('c', 0))}   "
                        f"both {int(vc.get('both', 0))}",
            transform=ax.transAxes, ha="right", va="top", fontsize=SS, color=GREY)
    stylia.label(ax, xlabel="N-terminal initiation region (aa)",
                 ylabel="C-terminal initiation region (aa)",
                 title=f"Which terminus can be pulled? — {orgname}")

    # ---- panel 4: the shapes of disorder, by compartment ----
    ax = axs.next()
    if "localization" not in d.columns:
        _note(ax, "no localization table")
    else:
        ct = pd.crosstab(d["localization"], d["disorder_architecture"])
        ct = ct.reindex(columns=[c for c in ARCH_ORDER if c in ct.columns], fill_value=0)
        ct = ct.loc[ct.sum(axis=1).sort_values(ascending=False).index]
        pct = ct.div(ct.sum(axis=1), axis=0) * 100
        xi = np.arange(len(pct)); bottom = np.zeros(len(pct))
        for a in pct.columns:
            ax.bar(xi, pct[a].to_numpy(), bottom=bottom, color=ARCH_COLOR[a], width=0.72,
                   label=a.replace("_", " "))
            bottom += pct[a].to_numpy()
        ax.set_xticks(xi)
        ax.set_xticklabels([f"{LOC.LOC_ABBREV.get(i, i)}\n({int(v)})"
                            for i, v in zip(pct.index, ct.sum(axis=1))], fontsize=SS)
        ax.set_ylim(0, 100)
        # The envelope compartments look highly "terminal IDR", but AlphaFold models the PRECURSOR:
        # the disordered N-terminus of a periplasmic or OM protein is overwhelmingly its uncleaved
        # Sec signal peptide, which the mature protein does not have. Measured here rather than
        # assumed. It is not simply an artefact — the spec's pre-export window argument makes a Sec
        # precursor genuinely engageable — but it is a different population and must be named.
        sp_note = ""
        if "has_signal_peptide" in d.columns:
            t = d[d["disorder_architecture"] == "terminal_idr"]
            env = t[t["localization"].isin(["periplasm", "outer_membrane", "extracellular",
                                            "cell_wall_surface"])]
            if len(env):
                frac = (env["has_signal_peptide"] == True).mean()  # noqa: E712
                sp_note = f"envelope terminal-IDR is {frac:.0%} Sec signal peptide"
        stylia.label(ax, xlabel=sp_note, ylabel="% of compartment",
                     title=f"Shapes of disorder, by compartment — {orgname}")
        ax.xaxis.label.set_fontsize(SS)
        ax.xaxis.label.set_color("#777777")
        ax.legend(fontsize=SS, frameon=False, ncol=2, loc="upper center",
                  bbox_to_anchor=(0.5, -0.20), columnspacing=0.9, handlelength=1.1)

    # ---- panel 5: is there any way in at all? ----
    # The previous version of this panel plotted the distribution of internal run lengths, which
    # states the claim without answering it. What matters is the accounting: how much does allowing
    # an internal loop actually add over terminal-only threading? Internal runs are held to the same
    # 20 aa engagement bar as termini, so the comparison is like-for-like.
    ax = axs.next()
    ok = d[d["af_available"] == True]                                        # noqa: E712
    term = ok["max_init_len_70"].fillna(0) >= ENGAGE_AA
    inter = ok["internal_disorder_max_run"].fillna(0) >= ENGAGE_AA
    cats = [("terminal tail only", int((term & ~inter).sum()), HIT),
            ("both tail and loop", int((term & inter).sum()), "#3C5488"),
            ("internal loop only", int((~term & inter).sum()), "#8491B4"),
            ("no route in", int((~term & ~inter).sum()), MISS)]
    ys = np.arange(len(cats))[::-1]
    bars = ax.barh(ys, [c[1] for c in cats], color=[c[2] for c in cats], height=0.68)
    ax.bar_label(bars, labels=[f"{c[1]:,}  ({c[1]/len(ok):.1%})" for c in cats],
                 padding=3, fontsize=SS)
    ax.set_yticks(ys)
    ax.set_yticklabels([c[0] for c in cats], fontsize=SS)
    ax.set_xlim(0, max(c[1] for c in cats) * 1.30)
    ftsz = by_gene.get("FtsZ")
    if ftsz is not None and pd.notna(ftsz.get("internal_disorder_max_run")):
        ax.annotate(f"FtsZ sits here — its {int(ftsz['internal_disorder_start'])}–"
                    f"{int(ftsz['internal_disorder_end'])} loop\ncarries the known ClpXP site 349–358",
                    xy=(cats[2][1], ys[2]), xytext=(0.42, 0.52), textcoords="axes fraction",
                    fontsize=SS, color="#1A1A1A", va="center",
                    arrowprops=dict(arrowstyle="-", color="#999", lw=0.7))
    stylia.label(ax, xlabel=f"proteins (of {len(ok):,} with an AlphaFold model)", ylabel="",
                 title=f"Is there any way in? — {orgname}")

    # ---- panel 6: what the named targets actually look like ----
    # After five distributions, end on the molecules. Each bar is one protein drawn to scale, shaded
    # where AlphaFold pLDDT < 70, so a threadable terminus is something you can see rather than infer.
    ax = axs.next()
    picks = [g for g in ("AcpP", "GroEL", "GyrB", "FtsZ", "DnaK", "GyrA") if g in by_gene]
    for y, g in enumerate(picks):
        r = by_gene[g]
        p = D.plddt_series(org, r["uniprot_accession"])
        if not p:
            continue
        Lp = len(p)
        ax.broken_barh([(0, Lp)], (y - 0.32, 0.64), facecolors="#EAEAE8", edgecolors="none")
        runs, i0 = [], None
        for k, v in enumerate(p):
            if v < D.PLDDT_DISORDER_CUT:
                i0 = k if i0 is None else i0
            elif i0 is not None:
                runs.append((i0, k - 1)); i0 = None
        if i0 is not None:
            runs.append((i0, Lp - 1))
        for a, b in runs:
            ln = b - a + 1
            terminal = (a == 0) or (b == Lp - 1)
            if terminal and ln >= ENGAGE_AA:
                c = HIT                      # long enough to thread
            elif terminal:
                c = "#BFE3D8"                # disordered terminus, but too short
            elif ln >= ENGAGE_AA:
                c = "#8491B4"                # internal candidate (the FtsZ route)
            else:
                c = "#D2D2D0"
            ax.broken_barh([(a, ln)], (y - 0.32, 0.64), facecolors=c, edgecolors="none")
        ax.text(Lp + 12, y, f"{Lp} aa", va="center", ha="left", fontsize=SS, color="#777777")
    ax.set_yticks(range(len(picks)))
    ax.set_yticklabels([g + (" (neg)" if g in NEGATIVE else "") for g in picks],
                       fontsize=SS)
    for tick, g in zip(ax.get_yticklabels(), picks):
        if g in NEGATIVE:
            tick.set_color("#B2182B")
    ax.set_ylim(-0.7, len(picks) - 0.3)
    ax.set_xlim(0, max(int(by_gene[g]["modeled_len"]) for g in picks) * 1.16)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=HIT, label=f"terminal, ≥{ENGAGE_AA} aa"),
                       Patch(color="#BFE3D8", label="terminal, too short"),
                       Patch(color="#8491B4", label=f"internal, ≥{ENGAGE_AA} aa"),
                       Patch(color="#EAEAE8", label="ordered")],
              fontsize=SS, frameon=False, ncol=2, loc="lower right", handlelength=1.1,
              columnspacing=0.9)
    stylia.label(ax, xlabel="residue", ylabel="",
                 title=f"What a threadable terminus looks like — {orgname}")

    out = D.REPO_ROOT / "output" / "plots" / f"10j_disorder_{prefix}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    _save(out)
    print(f"[{org}] wrote {out.relative_to(D.REPO_ROOT)}", flush=True)


if __name__ == "__main__":
    main()
