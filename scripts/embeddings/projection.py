"""Project the stage-01 ESM-C embeddings to 2D, one map per species.

Reduces each species' 1152-dim ESM-C 600M matrix to two coordinates per protein, so the proteome can
be looked at and so later stages can paint their own scores onto a fixed frame.

The recipe
----------
openTSNE multiscale (perplexities 50/500), cosine, on PCA-50 of the z-scored embeddings, with
`dof=0.8`. This is v1's choice, inherited rather than re-derived: it won a 33-config sweep over
input representation, distance metric, t-SNE structure, and method (vs UMAP and PaCMAP). The
verdict is tabulated in `docs/embeddings.md` Part 2; the original sweep code was deleted, so that
table and `legacy/docs/projection_exploration_log.md` are the only record. Do not re-run the sweep.

`dof` is the load-bearing knob, not perplexity. 0.5-0.6 gives compact clusters flung too far apart;
1.0 keeps them close but smears them; 0.8 is the balance. `--method umap|pacmap` runs a comparison
into `evidence/` and never touches the deliverable.

Scope
-----
The three bacteria. Human has no stage-01 embeddings, so there is nothing to project.

Coordinates are relative and per-species -- three independent embeddings, no shared frame. See
`src/projections.py` before plotting or comparing them.

The control
-----------
`sklearn.manifold.trustworthiness` (k=10, cosine) against the PCA-50 input, with a 0.90 floor. The
floor is calibrated on measured baselines: random 2D scores 0.4996, PCA-2 alone 0.7821, this recipe
0.9752. A t-SNE that diverged or was handed the wrong metric still returns a plausible-looking array
of finite numbers, so a shape check is not enough.

Output
------
    data/processed/embeddings/projection_<species>.tsv
        uniprot_ac      str     UniProt accession, the canonical key
        tsne_x          float   first  t-SNE coordinate (relative, per-species)
        tsne_y          float   second t-SNE coordinate (relative, per-species)

    evidence/projection_manifest.tsv       params, KL, trustworthiness, timings per species
    evidence/projection_<method>_<sp>.tsv  only with --method umap|pacmap

No colouring is stored -- join whatever you want to show on `uniprot_ac`.

Run with the `gradi` env:

    ~/miniconda3/envs/gradi/bin/python scripts/embeddings/projection.py
    ~/miniconda3/envs/gradi/bin/python scripts/embeddings/projection.py --species saureus
    ~/miniconda3/envs/gradi/bin/python scripts/embeddings/projection.py --method umap
"""

from __future__ import annotations

import argparse
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src import embeddings as E  # noqa: E402
from src import projections as PROJ  # noqa: E402
from src import matrices as M  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "processed" / "embeddings"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
DEFAULT_SPECIES = PROJ.SPECIES

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


# ---------------------------------------------------------------- the projection


def prepare(mat: np.ndarray, n_components: int, seed: int) -> np.ndarray:
    """z-score, then PCA. The shared input for the projection AND for the control.

    PCA-50 was measured indistinguishable from the full 1152-d map, so it is the fast default.
    Everything downstream sees this same matrix -- v1 fed t-SNE the PCA but its colour-clustering
    the full 1152-d space, which was both slower and inconsistent.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    z = StandardScaler().fit_transform(mat)
    k = min(n_components, *mat.shape)
    return PCA(n_components=k, random_state=seed).fit_transform(z)


def project_opentsne(
    pca: np.ndarray, seed: int, dof: float, perplexities: list[int]
) -> tuple[np.ndarray, float]:
    """openTSNE multiscale, cosine -- the chosen configuration. Returns `(xy, kl_divergence)`.

    Two optimisation phases. The `exaggeration=12` in phase one is standard *early* exaggeration and
    is not the same knob as the final-phase exaggeration, which stays off: above ~1.5 with heavy
    tails it produces filament artifacts.
    """
    from openTSNE import TSNEEmbedding, affinity, initialization

    aff = affinity.Multiscale(
        pca,
        perplexities=perplexities,
        metric=PROJ.METRIC,
        n_jobs=-1,
        random_state=seed,
    )
    init = initialization.pca(pca, random_state=seed)
    emb = TSNEEmbedding(init, aff, dof=dof, n_jobs=-1, random_state=seed)
    emb = emb.optimize(
        n_iter=PROJ.EARLY_ITER,
        exaggeration=PROJ.EARLY_EXAGGERATION,
        momentum=PROJ.EARLY_MOMENTUM,
    )
    emb = emb.optimize(n_iter=PROJ.LATE_ITER, momentum=PROJ.LATE_MOMENTUM)
    return np.asarray(emb), float(emb.kl_divergence)


def project_umap(pca: np.ndarray, seed: int) -> tuple[np.ndarray, float]:
    """UMAP, for comparison only. Rejected in v1: continents cohere but families blend."""
    import umap

    xy = umap.UMAP(
        n_components=2, n_neighbors=15, min_dist=0.1, metric=PROJ.METRIC, random_state=seed
    ).fit_transform(pca)
    return np.asarray(xy), float("nan")


def project_pacmap(pca: np.ndarray, seed: int) -> tuple[np.ndarray, float]:
    """PaCMAP, for comparison only. Rejected in v1: splits into 2 masses with filament tails."""
    import pacmap

    xy = pacmap.PaCMAP(n_components=2, n_neighbors=None, random_state=seed).fit_transform(pca)
    return np.asarray(xy), float("nan")


def trustworthiness(pca: np.ndarray, xy: np.ndarray) -> float:
    """Fraction of each point's k nearest neighbours that survive the reduction. 0.5 = random."""
    from sklearn.manifold import trustworthiness as tw

    return float(tw(pca, xy, n_neighbors=PROJ.TRUSTWORTHINESS_K, metric=PROJ.METRIC))


# ---------------------------------------------------------------- per species


def run_species(
    species: str,
    method: str,
    n_components: int,
    perplexities: list[int],
    dof: float,
    seed: int,
    limit: int | None,
    refresh: bool,
) -> dict:
    """Project one species and write its table. Returns the manifest row."""
    smoke = limit is not None
    if method == PROJ.METHOD:
        name = f"projection_{species}.tsv"
        out = SCRATCH_DIR / f"smoke_{name}" if smoke else OUT_DIR / name
        cols = list(PROJ.COORD_COLS)
    else:
        # A comparison map never lands beside the deliverable, and its columns are named after the
        # method so it cannot be mistaken for one in a merged frame.
        name = f"projection_{method}_{species}.tsv"
        out = (SCRATCH_DIR if smoke else EVIDENCE_DIR) / (f"smoke_{name}" if smoke else name)
        cols = [f"{method}_x", f"{method}_y"]

    n_expected = E.metadata(species)["n"]

    if out.exists() and not refresh:
        cached = pd.read_csv(out, sep="\t")
        say(f"  {species:<14} cached  {len(cached):>6} rows  ({out.name})")
        return dict(
            species=species, n=len(cached), n_expected=n_expected, method=method,
            pca_components=n_components, perplexities=",".join(map(str, perplexities)),
            metric=PROJ.METRIC, dof=dof, seed=seed, kl_divergence=float("nan"),
            trustworthiness_k10=float("nan"), seconds=0.0, cached=True,
            path=str(out.relative_to(REPO_ROOT)),
            run_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        )

    t0 = time.time()
    accs, mat = E.load(species)
    if smoke:
        accs, mat = accs[:limit], mat[:limit]
    say(f"  {species:<14} {mat.shape[0]:>6} x {mat.shape[1]}  ->  PCA-{n_components} ...")

    pca = prepare(mat, n_components, seed)

    if method == PROJ.METHOD:
        xy, kl = project_opentsne(pca, seed, dof, perplexities)
    elif method == "umap":
        xy, kl = project_umap(pca, seed)
    else:
        xy, kl = project_pacmap(pca, seed)

    tw = trustworthiness(pca, xy)
    elapsed = time.time() - t0

    frame = pd.DataFrame({"uniprot_ac": accs, cols[0]: xy[:, 0], cols[1]: xy[:, 1]})
    frame = M.reindex(frame, species)   # canonical row order -- see src/matrices.py
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp.tsv")
    frame.to_csv(tmp, sep="\t", index=False)
    tmp.rename(out)  # atomic: a killed run never leaves a half-written map behind

    flag = "" if tw >= PROJ.MIN_TRUSTWORTHINESS else "   <- BELOW FLOOR"
    say(
        f"  {'':<14} {len(frame):>6} rows  trustworthiness {tw:.4f}{flag}  "
        f"KL {kl:.3f}  {elapsed:.0f}s  -> {out.name}"
    )

    return dict(
        species=species, n=len(frame), n_expected=len(accs), method=method,
        pca_components=pca.shape[1], perplexities=",".join(map(str, perplexities)),
        metric=PROJ.METRIC, dof=dof, seed=seed, kl_divergence=round(kl, 4),
        trustworthiness_k10=round(tw, 4), seconds=round(elapsed, 1), cached=False,
        path=str(out.relative_to(REPO_ROOT)),
        run_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )


def merge_manifest(rows: list[dict], man_path: Path) -> pd.DataFrame:
    """Carry a cached species' previous manifest row forward instead of overwriting it.

    A cached row has no measurement to report -- nothing was fitted, so there is no trustworthiness,
    no KL and no runtime. Writing it verbatim would replace the real numbers from the run that
    actually produced the table with blanks, so a plain re-run would silently destroy the run record.
    Only a species that was genuinely recomputed overwrites its row.
    """
    if not man_path.exists():
        return pd.DataFrame(rows)
    prior = pd.read_csv(man_path, sep="\t").set_index("species")
    out = []
    for r in rows:
        if r["cached"] and r["species"] in prior.index:
            # `species` first, so carrying a row forward cannot reorder the manifest's columns.
            out.append({"species": r["species"], **prior.loc[r["species"]].to_dict()})
        else:
            out.append(r)
    # And the freshly-built row's key order is the canonical one.
    return pd.DataFrame(out)[list(rows[0])]


def verify(species: str, method: str, limit: int | None) -> list[str]:
    """Re-read what was written and check it against the embeddings. Returns failure strings."""
    if method != PROJ.METHOD or limit is not None:
        return []  # only the deliverable is guarded; comparison maps are exploratory
    fails: list[str] = []
    accs, _ = E.load(species)
    df = PROJ.load(species)

    if len(df) != len(accs):
        fails.append(f"{species}: {len(df)} rows, expected {len(accs)}")
    if not np.isfinite(df[list(PROJ.COORD_COLS)].to_numpy()).all():
        fails.append(f"{species}: non-finite coordinates (t-SNE diverged)")
    if list(df["uniprot_ac"]) != list(accs):
        missing = set(accs) - set(df["uniprot_ac"])
        extra = set(df["uniprot_ac"]) - set(accs)
        fails.append(
            f"{species}: accessions do not match the embeddings row-for-row "
            f"({len(missing)} missing, {len(extra)} unexpected)"
        )
    return fails


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--species", nargs="+", default=list(DEFAULT_SPECIES),
                    choices=list(PROJ.SPECIES),
                    help="species to project (default: all three bacteria)")
    ap.add_argument("--method", default=PROJ.METHOD,
                    choices=[PROJ.METHOD, *PROJ.COMPARISON_METHODS],
                    help="openTSNE is canonical; umap/pacmap are comparison-only "
                         "and write to evidence/")
    ap.add_argument("--pca", type=int, default=PROJ.PCA_COMPONENTS,
                    help="PCA components before projecting (default: 50, measured ~= full 1152-d)")
    ap.add_argument("--perplexities", type=int, nargs="+", default=list(PROJ.PERPLEXITIES),
                    help="multiscale perplexities (default: 50 500 = local + global)")
    ap.add_argument("--dof", type=float, default=PROJ.DOF,
                    help="t-SNE tail weight, the compact <-> detached knob (default: 0.8)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--limit", type=int, help="only the first N proteins per species (smoke test)")
    ap.add_argument("--refresh", action="store_true", help="recompute even if the table exists")
    ap.add_argument("--dry-run", action="store_true", help="print the plan and write nothing")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    rule("=")
    say("STAGE 01 - 2D projections of the ESM-C embeddings")
    rule("=")
    say("  in       : data/processed/embeddings/embeddings_<species>.npz")
    say(f"  out      : {OUT_DIR.relative_to(REPO_ROOT)}/projection_<species>.tsv  (+ evidence/ + scratch/)")
    say(f"  method   : {args.method}"
        + ("  (canonical)" if args.method == PROJ.METHOD else "  (COMPARISON ONLY -> evidence/)"))
    say(f"  recipe   : PCA-{args.pca} of z-scored, multiscale "
        f"{args.perplexities}, {PROJ.METRIC}, dof={args.dof}")
    say(f"  species  : {', '.join(args.species)}")
    say(f"  control  : trustworthiness k={PROJ.TRUSTWORTHINESS_K} "
        f"floor {PROJ.MIN_TRUSTWORTHINESS}")
    say(f"  seed     : {args.seed}   (structure is reproducible; floats are not, n_jobs=-1)")
    if args.limit:
        say(f"  limit    : {args.limit} proteins per species (SMOKE TEST -> scratch/smoke_*)")
    say()

    if args.dry_run:
        say("  --dry-run: nothing projected, nothing written.")
        for sp in args.species:
            n = E.metadata(sp)["n"]
            target = (f"projection_{sp}.tsv" if args.method == PROJ.METHOD
                      else f"evidence/projection_{args.method}_{sp}.tsv")
            say(f"    {sp:<14} {n:>6} proteins  ->  {target}")
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    rule()
    say("PROJECT")
    rule()
    rows = [
        run_species(sp, args.method, args.pca, args.perplexities, args.dof,
                    args.seed, args.limit, args.refresh)
        for sp in args.species
    ]
    say()

    rule()
    say("VERIFY")
    rule()
    fails: list[str] = []
    for sp in args.species:
        fails.extend(verify(sp, args.method, args.limit))
    if args.method != PROJ.METHOD:
        say("  skipped: comparison maps are exploratory, not guarded.")
    elif args.limit:
        say("  skipped: smoke test.")
    elif fails:
        for f in fails:
            say(f"  FAIL  {f}")
    else:
        say(f"  {len(args.species)} species: row counts, finite coordinates and accession "
            "order all match the embeddings.")
    say()

    # A smoke test and a comparison run each get their own manifest, so neither can overwrite the
    # deliverable's -- the `smoke_`-prefix rule the other stages use.
    stem = "projection_manifest.tsv" if args.method == PROJ.METHOD \
        else f"projection_{args.method}_manifest.tsv"
    man_path = (SCRATCH_DIR if args.limit else EVIDENCE_DIR) / (f"smoke_{stem}" if args.limit else stem)
    man = merge_manifest(rows, man_path)
    man.to_csv(man_path, sep="\t", index=False)
    # The summary reports what this run measured; the manifest keeps the run that measured it.
    measured = {r["species"]: r for r in man.to_dict("records")}

    rule()
    say("SUMMARY")
    rule()
    say(f"  {'species':<14} {'n':>7} {'expected':>9} {'trustworth':>11} {'KL':>8} {'sec':>7}")
    for r in rows:
        # For a cached species these come from the run that actually fitted the map.
        m = measured.get(r["species"], r)
        tw, kl, secs = m["trustworthiness_k10"], m["kl_divergence"], m["seconds"]
        flag = ""
        if r["n"] != r["n_expected"]:
            flag = "   <- MISMATCH"
        elif tw == tw and tw < PROJ.MIN_TRUSTWORTHINESS:  # tw == tw is False for NaN
            flag = "   <- BELOW FLOOR"
        elif r["cached"]:
            flag = "   (cached)"
        tw_s = f"{tw:.4f}" if tw == tw else "n/a"
        kl_s = f"{kl:.3f}" if kl == kl else "-"
        say(f"  {r['species']:<14} {r['n']:>7} {r['n_expected']:>9} "
            f"{tw_s:>11} {kl_s:>8} {secs:>7.1f}{flag}")
    total = sum(measured.get(r["species"], r)["seconds"] for r in rows)
    say(f"\n  total {total / 60:.1f} min   manifest -> "
        f"{man_path.relative_to(REPO_ROOT)}")
    say()

    rule()
    say("OUTPUTS")
    rule()
    for r in rows:
        p = REPO_ROOT / r["path"]
        say(f"  {r['path']:<58} {p.stat().st_size / 1e3:>8.1f} kB")
    say("  no colour column: join COG / localization / degradability on `uniprot_ac`.")
    say("  load via src/projections.py -- load, load_all, coords_for, manifest.")
    say()

    # Guards last, so everything above is printed before we exit.
    if fails:
        sys.exit(
            "FAILED: the projection does not match the embeddings it came from.\n  "
            + "\n  ".join(fails)
            + "\n  Every downstream landscape joins on these accessions; do not trust this run."
        )
    # Checked against `measured`, so a cached re-run still enforces the floor on the stored score
    # rather than skipping the check because this run computed nothing.
    low = [
        m for m in (measured.get(r["species"], r) for r in rows)
        if m["trustworthiness_k10"] == m["trustworthiness_k10"]
        and m["trustworthiness_k10"] < PROJ.MIN_TRUSTWORTHINESS
    ]
    if low and args.method == PROJ.METHOD:
        sys.exit(
            "FAILED trustworthiness floor:\n  "
            + "\n  ".join(
                f"{r['species']} scored {r['trustworthiness_k10']:.4f}, "
                f"below the {PROJ.MIN_TRUSTWORTHINESS} floor"
                for r in low
            )
            + "\n  For scale: random 2D coordinates score ~0.50 and PCA-2 alone ~0.78, while this"
              "\n  recipe measured 0.975. A score this low means the neighbourhood structure did not"
              "\n  survive the reduction -- wrong metric, collapsed PCA, or a diverged optimisation."
        )
    rule("=")
    say("stage 01 projections complete.")
    rule("=")


if __name__ == "__main__":
    main()
