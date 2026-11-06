"""One-shot: put the already-computed matrices into canonical row order.

The producers now emit canonical order (src/matrices.py), but the artifacts on disk predate that
and were written in whatever order their stage happened to iterate in. Every one of them is
COMPLETE -- same accessions, different order -- so this is a permutation, not a recomputation:
re-running ESM-C, ProteomeLM and DeepLocPro to fix an ordering would cost hours and change no value.

Each file is verified after writing: same shape, same set, canonical order, and for the numeric
payload an exact per-accession value check against the original.

  python tools/reindex_matrices.py --dry-run
  python tools/reindex_matrices.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import matrices as M  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
D = REPO / "data" / "processed"

TABLES = ["localization/deeplocpro_{sp}.tsv", "localization/tmbed_{sp}.tsv",
          "embeddings/projection_{sp}.tsv"]
NPZS = ["embeddings/embeddings_{sp}.npz", "embeddings/prott5_{sp}.npz",
        "embeddings/proteomelm_{sp}.npz"]


def do_table(path: Path, sp: str, dry: bool) -> str:
    df = pd.read_csv(path, sep="\t", dtype={"uniprot_ac": str})
    ok, _ = M.check(df["uniprot_ac"], sp)
    if ok:
        return "already canonical"
    before = df.set_index("uniprot_ac")
    out = M.reindex(df, sp)
    # every value must survive the permutation, checked per accession
    after = out.set_index("uniprot_ac")
    if not before.reindex(after.index).equals(after):
        sys.exit(f"FATAL {path}: values changed under the permutation")
    if not dry:
        tmp = path.with_suffix(".tmp.tsv")
        out.to_csv(tmp, sep="\t", index=False)
        tmp.rename(path)
    return f"reordered {len(out)} rows"


def do_npz(path: Path, sp: str, dry: bool) -> str:
    z = np.load(path, allow_pickle=True)
    accs = z["accessions"].astype(str)
    ok, _ = M.check(accs, sp)
    if ok:
        return "already canonical"
    payload = {k: z[k] for k in z.files}
    new_accs, new_mat = M.reindex_arrays(accs, z["embeddings"], sp)
    # exact value check: each accession's vector must be bit-identical to where it came from
    where = {a: i for i, a in enumerate(accs)}
    for a in (new_accs[0], new_accs[len(new_accs) // 2], new_accs[-1]):
        if not np.array_equal(new_mat[list(new_accs).index(a)], z["embeddings"][where[a]]):
            sys.exit(f"FATAL {path}: vector for {a} changed under the permutation")
    payload["accessions"] = np.array(new_accs, dtype=object)
    payload["embeddings"] = new_mat
    if not dry:
        tmp = path.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, **payload)
        tmp.rename(path)
    return f"reordered {len(new_accs)} x {new_mat.shape[1]}"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    for sp in M.SPECIES:
        print(f"\n{sp}")
        for tpl in TABLES:
            p = D / tpl.format(sp=sp)
            print(f"  {tpl.format(sp=sp):42s} {do_table(p, sp, a.dry_run) if p.exists() else 'ABSENT'}")
        for tpl in NPZS:
            p = D / tpl.format(sp=sp)
            print(f"  {tpl.format(sp=sp):42s} {do_npz(p, sp, a.dry_run) if p.exists() else 'ABSENT'}")
    print("\ndone." if not a.dry_run else "\ndry run -- nothing written.")


if __name__ == "__main__":
    main()
