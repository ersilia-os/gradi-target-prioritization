"""What does a raw ESM-C cosine to an S. aureus protein actually mean, biologically?

`nn_similarity` is the cosine to the nearest S. aureus training protein. This calibrates it against
orthology, which is independent of the embedding, and checks the ortholog table itself on a short
list of proteins that must have an S. aureus ortholog.

Reads   data/processed/embeddings/embeddings_<species>.npz     ESM-C, via src.embeddings
        data/processed/orthology/orthologs.tsv                  OrthoFinder + reciprocal best hit
Writes  output/results/degradability/cosine_calibration.tsv     one row per species x similarity bin
        output/results/degradability/conserved_gene_check.tsv   one row per species x gene

**Measured result (2026-11-05): raw cosine barely separates relatives from strangers.** Random
Ec/Kp-Sa pairs have median 0.775/0.761 and a 95th percentile of 0.90/0.89; a true ortholog pair has
median 0.906/0.894; and a protein's nearest S. aureus neighbour is a true ortholog only 1.7-6.9% of
the time below 0.97. That is why `degradability_evidence` (which cut at 0.90) is no longer shipped.

Run with the `gradi` env:  python scripts/degradability/cosine_calibration.py [-q]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import embeddings as E  # noqa: E402
from src import orthology as O  # noqa: E402
from src import proteomes as P  # noqa: E402

OUT_DIR = REPO_ROOT / "output/results/degradability"
SPECIES = ("ecoli", "kpneumoniae")
BINS = [(0, .70), (.70, .80), (.80, .85), (.85, .90), (.90, .93), (.93, .95), (.95, .97), (.97, 1.01)]
# Proteins that every bacterium needs, so each MUST have an S. aureus ortholog. parC/parE are
# topoisomerase IV, close paralogs of gyrA/gyrB, where a clean one-to-one call is not expected.
CONSERVED = ["rpoA", "rpoB", "rpoC", "gyrA", "gyrB", "parC", "parE", "ftsZ", "secA", "secY", "dnaK",
             "groEL", "dnaA", "dnaN", "recA", "tufA", "fusA", "infB", "rplB", "rpsC", "atpD", "ffh",
             "ftsY", "map", "pyrG", "tsf", "lepA", "ileS", "leuS", "gltX"]
N_RANDOM = 200_000


def say(msg: str = "", quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def unit(x: np.ndarray) -> np.ndarray:
    return x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-9)


def embed(sp: str) -> tuple[list[str], np.ndarray]:
    acc, x = E.vectors_for(sp, P.load(sp)["uniprot_ac"].tolist())
    return list(acc), unit(x)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    q = ap.parse_args().quiet

    say("cosine_calibration.py  raw ESM-C cosine vs orthology, S. aureus as the target", q)
    sa_acc, sa = embed("saureus")
    sa_idx = {a: i for i, a in enumerate(sa_acc)}
    orth = O.load_orthologs()
    orth = orth[orth.target_species == "saureus"]
    rng = np.random.default_rng(0)
    bins, genes = [], []

    for sp in SPECIES:
        acc, x = embed(sp)
        idx = {a: i for i, a in enumerate(acc)}
        sim = x @ sa.T
        nn, nn_i = sim.max(1), sim.argmax(1)
        rnd = sim.ravel()[rng.integers(0, sim.size, N_RANDOM)]
        o = orth[(orth.query_species == sp) & orth.query_ac.isin(idx) & orth.target_ac.isin(sa_idx)]
        pair = np.array([sim[idx[a], sa_idx[b]] for a, b in zip(o.query_ac, o.target_ac)])
        say(f"\n  {sp}: random-pair cosine median {np.median(rnd):.3f}  p95 {np.percentile(rnd, 95):.3f}"
            f"  p99 {np.percentile(rnd, 99):.3f}", q)
        say(f"    ortholog pairs {len(o):,} ({int(o.is_rbh.sum()):,} RBH): median {np.median(pair):.3f}  "
            f"p10 {np.percentile(pair, 10):.3f}  p90 {np.percentile(pair, 90):.3f}", q)

        pairs = set(zip(o.query_ac, o.target_ac))
        has_orth = set(o.query_ac)
        nn_is_orth = np.array([(acc[i], sa_acc[nn_i[i]]) in pairs for i in range(len(acc))])
        has = np.array([a in has_orth for a in acc])
        say(f"    {'nn_similarity bin':<20}{'n':>6}{'has an ortholog':>17}{'nearest is ortholog':>21}", q)
        for lo, hi in BINS:
            m = (nn >= lo) & (nn < hi)
            if not m.sum():
                continue
            say(f"    [{lo:.2f},{hi:.2f}){'':<9}{int(m.sum()):>6}{has[m].mean():>16.1%}{nn_is_orth[m].mean():>20.1%}", q)
            bins.append({"species": sp, "bin": f"[{lo:.2f},{hi:.2f})", "n": int(m.sum()),
                         "share_with_ortholog": round(float(has[m].mean()), 4),
                         "share_nearest_is_ortholog": round(float(nn_is_orth[m].mean()), 4),
                         "random_pair_median": round(float(np.median(rnd)), 4),
                         "random_pair_p95": round(float(np.percentile(rnd, 95)), 4),
                         "ortholog_pair_median": round(float(np.median(pair)), 4)})

        # does the ortholog table recover proteins that must have an S. aureus ortholog?
        prot = P.load(sp).set_index("uniprot_ac")["gene_name"].astype(str).str.lower()
        miss = []
        for g in CONSERVED:
            acs = prot.index[prot == g.lower()]
            hit = o[o.query_ac.isin(acs)] if len(acs) else o.iloc[0:0]
            status = ("no_gene_name_in_proteome" if not len(acs)
                      else "no_ortholog" if hit.empty else "rbh" if hit.is_rbh.any() else "orthofinder_only")
            genes.append({"species": sp, "gene": g, "status": status})
            if status in ("no_ortholog", "no_gene_name_in_proteome"):
                miss.append(f"{g}({status})")
        n_ok = sum(1 for r in genes if r["species"] == sp and r["status"] in ("rbh", "orthofinder_only"))
        say(f"    conserved genes with an S. aureus ortholog: {n_ok}/{len(CONSERVED)}   missing: {miss or 'none'}", q)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(bins).to_csv(OUT_DIR / "cosine_calibration.tsv", sep="\t", index=False)
    pd.DataFrame(genes).to_csv(OUT_DIR / "conserved_gene_check.tsv", sep="\t", index=False)
    say("\n  wrote output/results/degradability/cosine_calibration.tsv and conserved_gene_check.tsv", q)
    return 0


if __name__ == "__main__":
    sys.exit(main())
