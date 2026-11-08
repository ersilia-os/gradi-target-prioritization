"""Write webapp/data/roc.json: the ROC curves shown at the top of the Methods dialog.

Scores are the shipped fused degradability probabilities (out-of-fold for every *S. aureus*
protein in a screen), labels are the measured ADEP4 / ONC212 calls. The curve is a plotted
summary of the same numbers the cross-validation reports (fusion_eval.tsv): the script
exits non-zero if the AUROC it computes disagrees with the stated one by more than 0.02.
Env: gradi.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve

REPO_ROOT = Path(__file__).resolve().parents[2]
D = REPO_ROOT / "data" / "processed" / "degradability"
OUT = REPO_ROOT / "webapp" / "data" / "roc.json"
EXPECTED = {"adep4": 0.882, "onc212": 0.769}   # fused, cluster-grouped CV (CLAUDE.md, degradability)


def main():
    labels = pd.read_csv(D / "evidence" / "labels_saureus.tsv", sep="\t")
    scores = pd.read_csv(D / "degradability_saureus.tsv", sep="\t")
    d = labels.merge(scores, on="uniprot_ac", how="inner")
    out = {}
    for act in ("adep4", "onc212"):
        m = d[[act, f"{act}_prob"]].dropna()
        y, p = m[act].astype(int).to_numpy(), m[f"{act}_prob"].to_numpy()
        auc = roc_auc_score(y, p)
        fpr, tpr, _ = roc_curve(y, p)
        grid = np.linspace(0, 1, 81)
        out[act] = {"auc": round(float(auc), 3), "n": int(len(y)), "n_pos": int(y.sum()),
                    "fpr": [round(float(x), 4) for x in grid],
                    "tpr": [round(float(np.interp(x, fpr, tpr)), 4) for x in grid]}
        print(f"{act}: n={len(y)} pos={int(y.sum())} AUROC={auc:.4f} (stated {EXPECTED[act]})")
        if abs(auc - EXPECTED[act]) > 0.02:
            sys.exit(f"{act}: AUROC {auc:.3f} disagrees with the stated {EXPECTED[act]}")
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    print(f"wrote {OUT.relative_to(REPO_ROOT)} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
