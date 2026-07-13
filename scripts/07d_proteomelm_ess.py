"""ProteomeLM-Ess essentiality prediction (docs §4.3a, PRIMARY predictor).

ProteomeLM (Bitbol-Lab, Apache-2.0) is a proteome-scale transformer that contextualises a whole
proteome at once from per-protein ESM-C 600M embeddings. We reuse the ESM-C embeddings already
computed by 01a as its input, run the backbone to get contextualised per-protein embeddings, and
then score essentiality.

The published `-Ess` head is NOT released (the repo README is a TODO), so we train our own linear
head — exactly the pattern ProteomeLM ships for its PPI task. The training labels are the curated
E. coli essential set (`ec_transfer_essential` from 07c, i.e. the EcoGene/Keio 299-gene set), read
off the E. coli ProteomeLM embeddings. We report a cross-validated AUROC on E. coli, then freeze the
head and apply it to the requested organism. Because the E. coli labels are curated (not OGEE-derived),
this side-steps the OGEE-parroting concern the spec raises for ProteomeLM.

Steps (all cached / resumable):
  1. load ESM-C 600M embeddings (data/processed/<org>/embeddings/<prefix>_esmc600m_embeddings.npz)
  2. ProteomeLM forward over the whole proteome -> contextualised embeddings, cached to .pt
  3. train (or load) the logistic essentiality head on E. coli embeddings + labels
  4. predict -> output/results/<org>/<prefix>_ess_proteomelm.csv
     columns: uniprot_accession, proteomelm_ess_score, proteomelm_ess_essential, proteomelm_ess_status

Run with the `gradi` conda env interpreter (CPU fp32 by default; the forward is ~seconds–minutes).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import essentiality as E  # noqa: E402

HEAD_DIR = E.REPO_ROOT / "data" / "processed" / "other" / "proteomelm"
ESSENTIAL_PROB_CUTOFF = 0.5  # head probability at/above -> called essential


def _esmc_npz(org: str) -> Path:
    _, prefix = E.ORGANISMS[org]
    return E.REPO_ROOT / "data" / "processed" / org / "embeddings" / f"{prefix}_esmc600m_embeddings.npz"


def get_ctx_embeddings(org: str, size: str, device: str):
    """ProteomeLM contextualised embeddings for a whole proteome; cached to a .pt file.

    Returns (accessions: np.ndarray[str], ctx: np.ndarray[float32] (N, dim)).
    """
    import torch

    _, prefix = E.ORGANISMS[org]
    npz_path = _esmc_npz(org)
    if not npz_path.exists():
        raise SystemExit(f"Missing ESM-C embeddings {npz_path}. Run scripts/01a_esmc_embeddings.py --organism {org}.")
    z = np.load(npz_path, allow_pickle=True)
    accs = z["accessions"].astype(str)
    esm = z["embeddings"].astype(np.float32)

    cache = E.essentiality_processed_dir(org, "proteomelm") / f"{prefix}_proteomelm_{size}_ctx.npz"
    if cache.exists():
        c = np.load(cache, allow_pickle=True)
        print(f"[{org}] using cached ProteomeLM ctx {cache.relative_to(E.REPO_ROOT)} {c['ctx'].shape}", flush=True)
        return c["accessions"].astype(str), c["ctx"].astype(np.float32)

    from proteomelm import ProteomeLMForMaskedLM

    print(f"[{org}] loading ProteomeLM-{size} ...", flush=True)
    t = time.time()
    model = ProteomeLMForMaskedLM.from_pretrained(f"Bitbol-Lab/ProteomeLM-{size}").to(device).float().eval()
    emb = torch.tensor(esm, device=device).unsqueeze(0)  # (1, N, 1152)
    with torch.no_grad():
        out = model(inputs_embeds=emb, group_embeds=None,
                    output_hidden_states=True, output_attentions=False)
    ctx = out.hidden_states[-1].squeeze(0).float().cpu().numpy()
    print(f"[{org}] ProteomeLM forward: {ctx.shape} in {time.time()-t:.1f}s", flush=True)
    np.savez_compressed(cache, accessions=accs, ctx=ctx)
    return accs, ctx


def _ecoli_labels() -> pd.DataFrame:
    """E. coli essential labels (from 07c). Columns: uniprot_accession, essential (bool)."""
    p = E.results_dir("ecoli") / "ec_ess_ecoli.csv"
    if not p.exists():
        raise SystemExit("Missing output/results/ecoli/ec_ess_ecoli.csv. Run scripts/07c_ecoli_transfer.py --organism ecoli.")
    df = pd.read_csv(p)[["uniprot_accession", "ec_transfer_essential"]]
    return df.rename(columns={"ec_transfer_essential": "essential"})


def train_or_load_head(size: str, device: str):
    """Train (cached) the logistic essentiality head on E. coli ctx embeddings + labels.

    Returns (pipeline, cv_auroc). The pipeline standardises features then logistic-regresses.
    """
    import joblib
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    HEAD_DIR.mkdir(parents=True, exist_ok=True)
    head_path = HEAD_DIR / f"ess_head_{size}.joblib"
    if head_path.exists():
        d = joblib.load(head_path)
        print(f"using cached head {head_path.relative_to(E.REPO_ROOT)} (cv AUROC={d['cv_auroc']:.3f})", flush=True)
        return d["pipe"], d["cv_auroc"]

    accs, ctx = get_ctx_embeddings("ecoli", size, device)
    lab = _ecoli_labels().set_index("uniprot_accession")["essential"].reindex(accs).fillna(False).to_numpy()
    print(f"training head on E. coli: {len(accs)} proteins, {int(lab.sum())} essential", flush=True)

    pipe = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0),
    )
    # honest cross-validated AUROC before fitting on everything
    proba = cross_val_predict(pipe, ctx, lab, cv=5, method="predict_proba")[:, 1]
    cv_auroc = roc_auc_score(lab, proba)
    pipe.fit(ctx, lab)
    joblib.dump({"pipe": pipe, "cv_auroc": cv_auroc, "size": size}, head_path)
    print(f"head trained: 5-fold CV AUROC = {cv_auroc:.3f}; cached to {head_path.relative_to(E.REPO_ROOT)}", flush=True)
    return pipe, cv_auroc


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(E.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--size", choices=["XS", "S", "M", "L"], default="M")
    ap.add_argument("--device", choices=["cpu", "mps"], default="cpu")
    args = ap.parse_args()
    org, size = args.organism, args.size
    _, prefix = E.ORGANISMS[org]

    try:
        pipe, cv_auroc = train_or_load_head(size, args.device)
        accs, ctx = get_ctx_embeddings(org, size, args.device)
        score = pipe.predict_proba(ctx)[:, 1]
        status = f"ok (head cv_auroc={cv_auroc:.3f})"
    except Exception as exc:  # noqa: BLE001 — never let the predictor crash the pipeline
        print(f"[{org}] ProteomeLM failed ({exc}); emitting deferred placeholder", flush=True)
        accs = np.array(E.load_accessions(org))
        score = np.full(len(accs), np.nan)
        cv_auroc = float("nan")
        status = f"deferred: {type(exc).__name__}"

    df = pd.DataFrame({
        "uniprot_accession": accs,
        "proteomelm_ess_score": np.round(score, 4),
        "proteomelm_ess_essential": score >= ESSENTIAL_PROB_CUTOFF if np.isfinite(score).all() else pd.NA,
        "proteomelm_ess_status": status,
    })
    out = E.results_dir(org) / f"{prefix}_ess_proteomelm.csv"
    df.to_csv(out, index=False)
    n_ess = int((df["proteomelm_ess_score"] >= ESSENTIAL_PROB_CUTOFF).sum()) if np.isfinite(score).all() else 0
    print(f"[{org}] wrote {out.relative_to(E.REPO_ROOT)} ({len(df)} proteins; "
          f"{n_ess} predicted essential; {status})", flush=True)


if __name__ == "__main__":
    main()
