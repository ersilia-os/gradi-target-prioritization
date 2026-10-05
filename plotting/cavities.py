"""The top-ranked K. pneumoniae pockets, as structures.

    cavities.png   three AlphaFold models, each with its best-scoring P2Rank pocket

Three proteins, side by side: **irp2 · acrB · gyrB** (owner, 2026-10-05). All three sit at the top
of the `pockets_consensus` ranking and all three are `pockets_evidence` 3 -- a drug-like ligand has
actually been seen in that protein's own PDB structures. `gyrB` is in `src/interest.py` and carries
14 of them.

**THE THREE PANELS ARE ON A COMMON SCALE** -- the worker renders twice, once to learn what camera
distance each molecule needs and once with the largest of them for every panel, and the montage
crops all three to a shared HEIGHT with per-panel widths and matching `width_ratios`. Without both, each frame is filled independently and an
805-aa protein comes out larger than a 1,048-aa one, which contradicts the titles.

**READ THE LENGTH IN EACH TITLE.** `pockets_consensus` correlates with protein length at rho 0.66,
so a high rank partly reflects size: `irp2` is 2,035 aa, the 100th percentile of the proteome. The
titles carry the length so a reader can see that rather than take the ranking at face value.

**The cartoon carries pLDDT (AlphaFold palette) and the pocket is a transparent magenta surface**,
with the camera oriented on the pocket. **pLDDT is model confidence, NOT pocket confidence** -- an
orange loop means the fold there is uncertain and says nothing about whether the cavity is real.

Crosses a process boundary into **`gradi-pymol`** (`from pymol import cmd` exists only there), so
the rendering lives in `scripts/pockets/workers/pymol_render.py` and is invoked with `conda run`.
This file selects, tiles and captions; it never imports pymol.

NOTE `REPO_ROOT` is `parents[1]`, not the `parents[2]` used under `scripts/<task>/`.

Run with the `gradi` env, ONE PLOT SCRIPT AT A TIME (stylia rmtree's the matplotlib cache dir):
    python plotting/cavities.py
    python plotting/cavities.py --refresh     # re-render the PNGs instead of reusing the cache
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

os.makedirs(matplotlib.get_cachedir(), exist_ok=True)
import stylia  # noqa: E402

from plotting import palette as PAL  # noqa: E402

from src import pockets as PK  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output" / "plots" / "presentation"
RENDER_DIR = REPO_ROOT / "data" / "processed" / "pockets" / "scratch" / "renders"
CIF_DIR = REPO_ROOT / "data" / "source" / "alphafold"
WORKER = REPO_ROOT / "scripts" / "pockets" / "workers" / "pymol_render.py"
CONDA = Path.home() / "miniconda3" / "bin" / "conda"
PYMOL_ENV = "gradi-pymol"

stylia.set_format("slide")
stylia.set_style("article")

SS = stylia.SLIDE_FONTSIZE_SMALL
SPECIES = "kpneumoniae"
#: The three shown, in this order (owner, 2026-10-05). They are drawn from the top of the
#: `pockets_consensus` ranking among named proteins, and each earns its place differently:
#: `gyrB` is the classic antibacterial target and carries 14 measured PDB ligands, `acrB` is the
#: efflux pump, and `irp2` is the highest-consensus protein in the proteome. Keeping `irp2` is
#: also honest about the confound -- at 2,035 aa it is the 100th percentile of length, and the
#: consensus is rho 0.66 with length.
GENES = ("irp2", "acrB", "gyrB")
N_PANELS = len(GENES)


def select() -> pd.DataFrame:
    """The named proteins in GENES, in that order.

    A fixed list rather than `nlargest`, so the figure does not silently change its subject when
    the pockets table is rebuilt. Every one is still checked against the ranking below and its
    percentile printed, so a drift is visible rather than hidden.
    """
    pk = PK.load(SPECIES)
    pr = P.load(SPECIES)[["uniprot_ac", "gene_name", "protein_name", "sequence"]].copy()
    d = pk.merge(pr, on="uniprot_ac", how="left")
    d["length"] = d["sequence"].str.len()

    rows = []
    for gene in GENES:
        hit = d[d["gene_name"].astype(str).str.strip() == gene]
        if hit.empty:
            sys.exit(f"FATAL no protein named {gene!r} in {SPECIES} -- the gene name moved")
        rows.append(hit.nlargest(1, "pockets_consensus").iloc[0])
    top = pd.DataFrame(rows).reset_index(drop=True)

    # Percentile WITHIN the proteome, not within the three -- the same convention the shortlist
    # heatmap uses, so a reader who has seen that panel reads this the same way.
    top["cons_pct"] = [100 * float((d["pockets_consensus"] < v).mean())
                       for v in top["pockets_consensus"]]
    top["len_pct"] = [100 * float((d["length"] < v).mean()) for v in top["length"]]
    return top


def best_pocket(acc: str) -> list[int]:
    """Residues of this protein's best-scoring ADMITTED P2Rank pocket.

    P2Rank rather than fpocket because it is the column the consensus leans on most (rho 0.83
    against fpocket's 0.75), and `admitted` because the stage already decided which pockets clear
    its pLDDT and size floors -- re-deriving that here would let the figure disagree with the
    table it illustrates.
    """
    p = PK.load_pockets(SPECIES, admitted_only=True)
    g = p[(p["uniprot_ac"] == acc) & (p["tool"] == "p2rank")]
    if not len(g):
        return []
    row = g.sort_values("score", ascending=False).iloc[0]
    return [int(x) for x in str(row["residues"]).split(";") if str(x).strip().isdigit()]


def render(top: pd.DataFrame, refresh: bool, say) -> dict:
    """Render each structure under gradi-pymol. Cached: PNGs are expensive and deterministic."""
    RENDER_DIR.mkdir(parents=True, exist_ok=True)
    jobs, have = [], {}
    for r in top.itertuples():
        png = RENDER_DIR / f"{r.uniprot_ac}.png"
        if png.exists() and not refresh:
            have[r.uniprot_ac] = png
            continue
        cif = CIF_DIR / SPECIES / f"AF-{r.uniprot_ac}-F1-model_v6.cif"
        if not cif.exists():
            say(f"    WARN no AlphaFold model for {r.uniprot_ac} ({r.gene_name}) -- skipped")
            continue
        jobs.append({"acc": r.uniprot_ac, "label": str(r.gene_name), "cif": str(cif),
                     "pocket": best_pocket(r.uniprot_ac)})

    if jobs:
        say(f"    rendering {len(jobs)} structures under {PYMOL_ENV} ...")
        job_path = RENDER_DIR / "jobs.json"
        job_path.write_text(json.dumps({"jobs": jobs}))
        cmd = [str(CONDA), "run", "-n", PYMOL_ENV, "pymol", "-cq", str(WORKER),
               "--", str(job_path), str(RENDER_DIR)]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            sys.exit(f"FATAL pymol worker failed ({res.returncode}):\n{res.stderr[-2000:]}")
        say(f"    {res.stdout.strip()}")
        for j in jobs:
            png = RENDER_DIR / f"{j['acc']}.png"
            if png.exists():
                have[j["acc"]] = png
    return have


def shared_crop(images: list, pad: int = 14) -> list:
    """Crop to a shared HEIGHT but each panel's own WIDTH, keeping one scale.

    Three constraints, and the third is what makes this fiddly:

    1. Scale must be common, or an 805-aa protein looks bigger than a 1,048-aa one. The worker
       already renders every molecule at one camera distance; the montage must not undo that.
    2. Panels must not carry dead space, or the structures read as far apart.
    3. `imshow` STRETCHES an image to fill its axes. So cropping each panel to its own width is
       only safe if each axes is made proportionally wide -- otherwise a narrow image is stretched
       back to full width and the common scale is lost again, silently.

    So: one height for all (rows align, vertical scale shared), each image's own width, and
    `main()` passes the widths as `width_ratios` so pixels-per-angstrom is constant across panels.
    """
    boxes = []
    for img in images:
        a = img[..., :3] if img.ndim == 3 else img
        mask = (a < 0.98).any(axis=2) if a.ndim == 3 else (a < 0.98)
        rows, cols = np.where(mask.any(axis=1))[0], np.where(mask.any(axis=0))[0]
        boxes.append((rows[0], rows[-1], cols[0], cols[-1]) if len(rows) and len(cols) else None)
    if not any(boxes):
        return images

    h = max(b[1] - b[0] for b in boxes if b) + 2 * pad
    out = []
    for img, b in zip(images, boxes):
        if b is None:
            out.append(img)
            continue
        cy = (b[0] + b[1]) // 2
        r0 = min(max(cy - h // 2, 0), max(img.shape[0] - h, 0))
        c0 = max(b[2] - pad, 0)
        c1 = min(b[3] + pad + 1, img.shape[1])
        out.append(img[r0:r0 + h, c0:c1])
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--refresh", action="store_true", help="re-render instead of reusing the cache")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    say = (lambda m: None) if args.quiet else (lambda m: print(m, flush=True))
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    top = select()
    say(f"  {len(top)} K. pneumoniae proteins, fixed set {GENES}:")
    for r in top.itertuples():
        say(f"    {str(r.gene_name):<8} consensus {r.pockets_consensus:.3f} "
            f"(pct {r.cons_pct:>4.1f})   {int(r.length):>5} aa (pct {r.len_pct:>4.1f})   "
            f"evidence {int(r.pockets_evidence)}   PDB ligands {int(r.n_ligands_pdb)}")

    have = render(top, args.refresh, say)
    top = top[top["uniprot_ac"].isin(have)]
    if not len(top):
        sys.exit("FATAL nothing rendered -- cannot draw the figure")

    imgs = shared_crop([mpimg.imread(str(have[r.uniprot_ac])) for r in top.itertuples()])
    # width_ratios = the cropped widths, so each axes is exactly as wide as its own image and
    # imshow does not stretch. This is what keeps one scale while removing the padding.
    fig, axs = stylia.create_figure(1, len(imgs), width=0.92, height=0.42,
                                    width_ratios=[im.shape[1] for im in imgs])
    for r, img in zip(top.itertuples(), imgs):
        ax = axs.next()
        ax.imshow(img)
        ax.set_axis_off()
        # The caption carries LENGTH beside consensus on purpose: the ranking is rho 0.66 with
        # length, and a reader must be able to see that for themselves.
        # ONE LINE, and no in-panel caption. A caption below the axes collides with the next
        # row's title; one inside the axes lands on the molecule, because each cropped render has
        # its own aspect ratio. The consensus value is dropped on purpose -- the panels are
        # ordered by it, so the ranking already says it, and the owner's standing rule is that a
        # figure carries minimal text.
        ax.set_title(f"{r.gene_name}  ·  {int(r.length)} aa  ·  "
                     f"{int(r.n_ligands_pdb)} PDB ligand{'s' if r.n_ligands_pdb != 1 else ''}",
                     fontsize=SS, color=PAL.INK, pad=6)
    for _ in range(N_PANELS - len(top)):
        axs.next().set_axis_off()

    # stylia leaves a generous gutter, which on three image panels is dead white space.
    fig.subplots_adjust(wspace=0.02)
    out = OUT_DIR / "cavities.png"
    stylia.save_figure(str(out))
    say(f"\n  -> {out.relative_to(REPO_ROOT)}")
    say("  cartoon coloured by pLDDT (AlphaFold palette, blue = confident); magenta transparent"
        " surface = the best-scoring admitted P2Rank pocket, oriented toward the camera.")
    say(f"  CAVEAT: pockets_consensus is rho ~0.66 with protein length, so a high rank partly"
        f" reflects size -- lengths are in the titles. All {len(top)} are evidence level 3"
        " (a drug-like ligand measured in their own PDB structures).")


if __name__ == "__main__":
    main()
