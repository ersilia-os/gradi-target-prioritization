"""Pockets on every AlphaFold model: fpocket 4.0 and P2Rank 2.5.1, side by side, never merged.

    data/processed/pockets/evidence/pockets_<species>.tsv      LONG: one row per pocket per tool
    data/processed/pockets/scratch/{pdb,fpocket,p2rank}/<species>/   caches, safe to delete

Long table columns:

    uniprot_ac · tool · pocket · score · raw_score · volume · n_residues · mean_plddt ·
    admitted · residues · model

* `score` is the number the deliverable reads: fpocket's **druggability score** (Schmidtke &
  Barril 2010, 0-1) and P2Rank's **calibrated probability** (0-1, `-c alphafold`).
* `raw_score` is each tool's own ranking score (fpocket `Score`, P2Rank `score`), kept so the
  ranking can be reproduced, but not used.
* `volume` is fpocket's (P2Rank does not report one).
* `model` is the predictor of the structure the pocket was found on: `alphafold_db_v6`, or
  `esmfold_v1` for the 32 proteins AlphaFold DB does not cover (see `esmfold.py`).
* `mean_plddt` is the mean AlphaFold pLDDT over the pocket's lining residues.
* `admitted` = `mean_plddt >= 70`. **This is the ONLY place model confidence enters the axis.**
  v1 used pLDDT twice -- inside a pocket consensus AND again as a whole-protein disorder penalty
  (`legacy/scripts/06e_pockets.py:264`, `06g:129`). Neither tool reads it: P2Rank's `alphafold`
  config drops B-factor as a feature and fpocket ignores it, so a pocket in a low-confidence
  region is scored as if it were real and must be filtered after the fact. 70 is AlphaFold's own
  "confident" band, not a tuned number.

A protein with a model and no pocket contributes no rows here; `merge.py` turns that into a
measured 0. A protein with no model is absent here and carries `evidence = none` in the
deliverable -- the distinction between "looked, found nothing" and "could not look" is kept.

Tools live in `gradi-pockets` (osx-64, Rosetta): `FPOCKET_BIN`, `P2RANK_DIR`, `POCKETS_JAVA_HOME`
override. Run THIS script with the `gradi` env; it reaches the tools across a process boundary.
  python scripts/pockets/predict.py
  python scripts/pockets/predict.py --species saureus --workers 8 --limit 50     # smoke
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from structures import model_ca, model_file  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
TASK_DIR = REPO_ROOT / "data" / "processed" / "pockets"
EVIDENCE_DIR = TASK_DIR / "evidence"
SCRATCH_DIR = TASK_DIR / "scratch"

FPOCKET_BIN = os.environ.get(
    "FPOCKET_BIN", str(Path.home() / "miniconda3/envs/gradi-pockets/bin/fpocket"))
P2RANK_DIR = Path(os.environ.get("P2RANK_DIR", str(REPO_ROOT / "tmp/tools/p2rank_2.5.1")))
JAVA_HOME = os.environ.get(
    "POCKETS_JAVA_HOME", str(Path.home() / "miniconda3/envs/gradi-pockets/lib/jvm"))
FPOCKET_VERSION = "4.0"
P2RANK_VERSION = "2.5.1"

PLDDT_ADMIT = 70.0

COLS = ["uniprot_ac", "tool", "pocket", "score", "raw_score", "volume", "n_residues",
        "mean_plddt", "admitted", "residues", "model"]

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


# --------------------------------------------------------------------------- cif -> pdb
def cif_to_pdb(cif: Path, out: Path) -> None:
    """fpocket and P2Rank both want PDB. biotite keeps residue numbering (1..L for AlphaFold)."""
    import biotite.structure.io.pdb as pdb
    import biotite.structure.io.pdbx as pdbx

    arr = pdbx.get_structure(pdbx.CIFFile.read(str(cif)), model=1, extra_fields=["b_factor"])
    f = pdb.PDBFile()
    pdb.set_structure(f, arr)
    out.parent.mkdir(parents=True, exist_ok=True)
    f.write(str(out))


# --------------------------------------------------------------------------- fpocket
def _parse_info(path: Path) -> dict[int, dict]:
    pockets: dict[int, dict] = {}
    cur = None
    for line in path.read_text().splitlines():
        s = line.strip()
        if s.startswith("Pocket") and s.endswith(":"):
            cur = int(s.split()[1])
            pockets[cur] = {}
        elif cur is not None and ":" in s:
            k, _, v = s.partition(":")
            try:
                pockets[cur][k.strip().lower()] = float(v.strip())
            except ValueError:
                pass
    return pockets


def _residues(atm: Path) -> list[int]:
    nums = set()
    for line in atm.read_text().splitlines():
        if line.startswith("ATOM"):
            nums.add(int(line[22:26]))
    return sorted(nums)


def fpocket_one(acc: str, pdb_path: str, out_json: str) -> str:
    """Run fpocket on one structure, keep a compact json, delete the bulky output dir."""
    pdb_path, out_json = Path(pdb_path), Path(out_json)
    out_dir = pdb_path.with_name(pdb_path.stem + "_out")
    shutil.rmtree(out_dir, ignore_errors=True)
    try:
        subprocess.run([FPOCKET_BIN, "-f", str(pdb_path)], check=True,
                       capture_output=True, timeout=900)
        info = out_dir / f"{pdb_path.stem}_info.txt"
        if not info.exists():
            raise RuntimeError("fpocket wrote no _info.txt")
        rec = []
        for n, props in _parse_info(info).items():
            atm = out_dir / "pockets" / f"pocket{n}_atm.pdb"
            rec.append({"pocket": n,
                        "score": props.get("druggability score"),
                        "raw_score": props.get("score"),
                        "volume": props.get("volume"),
                        "residues": _residues(atm) if atm.exists() else []})
        out_json.write_text(json.dumps({"uniprot_ac": acc, "pockets": rec}))
        return "ok"
    except Exception as e:  # noqa: BLE001 -- recorded, not swallowed: the json says why
        out_json.write_text(json.dumps({"uniprot_ac": acc, "error": str(e)[:300]}))
        return "error"
    finally:
        shutil.rmtree(out_dir, ignore_errors=True)


# --------------------------------------------------------------------------- P2Rank
def p2rank_batch(pdbs: list[Path], out_dir: Path, threads: int) -> None:
    """One JVM over a dataset list: model loading dominates P2Rank's cost on small proteins."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = out_dir / "dataset.ds"
    ds.write_text("\n".join(str(p.resolve()) for p in pdbs) + "\n")
    cmd = [str(P2RANK_DIR / "prank"), "predict", str(ds), "-c", "alphafold",
           "-o", str(out_dir), "-threads", str(threads), "-visualizations", "0"]
    log = out_dir / "p2rank.log"
    with open(log, "w") as fh:
        subprocess.run(cmd, check=True, env=dict(os.environ, JAVA_HOME=JAVA_HOME),
                       stdout=fh, stderr=subprocess.STDOUT)


def parse_p2rank(csv_path: Path) -> list[dict]:
    d = pd.read_csv(csv_path, skipinitialspace=True)
    d.columns = [c.strip() for c in d.columns]
    out = []
    for _, r in d.iterrows():
        res = sorted({int(t.split("_")[-1]) for t in str(r["residue_ids"]).split()
                      if t.split("_")[-1].lstrip("-").isdigit()})
        out.append({"pocket": int(r["rank"]), "score": float(r["probability"]),
                    "raw_score": float(r["score"]), "volume": None, "residues": res})
    return out


# --------------------------------------------------------------------------- per species
def run(species: str, workers: int, threads: int, limit: int | None) -> pd.DataFrame:
    af = pd.read_csv(EVIDENCE_DIR / f"alphafold_{species}.tsv", sep="\t")
    models = af[af["af_status"] == "model"]
    accs = models["uniprot_ac"].tolist()
    if limit:
        accs = accs[:limit]
    say(f"  {len(accs):,} proteins with a validated model "
        + ", ".join(f"({k} {v:,})" for k, v in
                    models["model_source"].value_counts().items()))

    pdb_dir = SCRATCH_DIR / "pdb" / species
    fp_dir = SCRATCH_DIR / "fpocket" / species
    p2_dir = SCRATCH_DIR / "p2rank" / species
    for d in (pdb_dir, fp_dir, p2_dir):
        d.mkdir(parents=True, exist_ok=True)

    # One structure file per protein -- the AFDB cif, or the ESMFold PDB -- cached under <acc>.
    stems = [(a, model_file(species, a)) for a in accs]

    # 1. structure -> pdb
    t0 = time.time()
    todo = [(st, f) for st, f in stems if not (pdb_dir / f"{st}.pdb").exists()]
    for i, (st, f) in enumerate(todo, 1):
        if f.suffix == ".cif":
            cif_to_pdb(f, pdb_dir / f"{st}.pdb")
        else:
            shutil.copyfile(f, pdb_dir / f"{st}.pdb")
        if i % 1000 == 0:
            say(f"    converted {i:,}/{len(todo):,}")
    say(f"  [1/3] -> pdb: {len(todo):,} converted, {len(stems) - len(todo):,} cached "
        f"({time.time() - t0:.0f}s)")

    # 2. fpocket, parallel, resumable per unit
    t0 = time.time()
    todo = [st for st, _ in stems if not (fp_dir / f"{st}.json").exists()]
    status = {"ok": 0, "error": 0}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(fpocket_one, st, str(pdb_dir / f"{st}.pdb"), str(fp_dir / f"{st}.json"))
                for st in todo]
        for i, fut in enumerate(as_completed(futs), 1):
            status[fut.result()] += 1
            if i % 500 == 0:
                say(f"    fpocket {i:,}/{len(todo):,}  ({time.time() - t0:.0f}s)")
    say(f"  [2/3] fpocket {FPOCKET_VERSION}: {len(todo):,} run ({status['ok']:,} ok, "
        f"{status['error']:,} error), {len(stems) - len(todo):,} cached ({time.time() - t0:.0f}s)")

    # 3. P2Rank, one batch over whatever has no predictions yet
    t0 = time.time()
    todo = [st for st, _ in stems if not (p2_dir / f"{st}.pdb_predictions.csv").exists()]
    if todo:
        p2rank_batch([pdb_dir / f"{st}.pdb" for st in todo], p2_dir, threads)
    say(f"  [3/3] P2Rank {P2RANK_VERSION} (-c alphafold): {len(todo):,} run, "
        f"{len(stems) - len(todo):,} cached ({time.time() - t0:.0f}s)")

    # 4. long table, with pocket pLDDT from the model's own CA B-factors (0-100 on both predictors)
    rows, errors, missing_p2 = [], [], []
    for a in accs:
        _, plddt, source = model_ca(species, a)
        fp = json.loads((fp_dir / f"{a}.json").read_text())
        if "error" in fp:
            errors.append(a)
        per_tool = {"fpocket": fp.get("pockets", [])}
        p2csv = p2_dir / f"{a}.pdb_predictions.csv"
        if p2csv.exists():
            per_tool["p2rank"] = parse_p2rank(p2csv)
        else:
            missing_p2.append(a)
        for tool, pockets in per_tool.items():
            for p in pockets:
                res = [r for r in p["residues"] if 1 <= r <= len(plddt)]
                mp = float(np.mean(plddt[[r - 1 for r in res]])) if res else np.nan
                rows.append({"uniprot_ac": a, "tool": tool, "pocket": p["pocket"],
                             "score": p["score"], "raw_score": p["raw_score"],
                             "volume": p["volume"], "n_residues": len(res),
                             "mean_plddt": round(mp, 2),
                             "admitted": bool(mp >= PLDDT_ADMIT) if res else False,
                             "residues": ";".join(map(str, res)), "model": source})
    if errors:
        say(f"  WARNING fpocket failed on {len(errors):,} structures (e.g. {errors[:3]})")
    if missing_p2:
        sys.exit(f"FATAL {species}: P2Rank produced no predictions file for "
                 f"{len(missing_p2):,} structures (e.g. {missing_p2[:3]}). An absent file is not "
                 "'no pocket' -- P2Rank writes an empty table for that. Read p2rank.log.")
    return pd.DataFrame(rows, columns=COLS)


def summarise(species: str, long: pd.DataFrame, n_model: int) -> list[dict]:
    out = []
    for tool in ("fpocket", "p2rank"):
        t = long[long["tool"] == tool]
        adm = t[t["admitted"]]
        out.append({
            "species": species, "tool": tool, "proteins": n_model,
            "pockets": len(t), "admitted": len(adm),
            "pct_admitted": 100 * len(adm) / max(len(t), 1),
            "proteins_any": t["uniprot_ac"].nunique(),
            "proteins_admitted": adm["uniprot_ac"].nunique(),
            "pct_proteins_admitted": 100 * adm["uniprot_ac"].nunique() / n_model,
        })
    return out


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", choices=SPECIES, default=list(SPECIES))
    ap.add_argument("--workers", type=int, default=8, help="parallel fpocket processes")
    ap.add_argument("--threads", type=int, default=8, help="P2Rank threads")
    ap.add_argument("--limit", type=int, default=None,
                    help="first N proteins; writes scratch/smoke_pockets_<species>.tsv only")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet

    say("=" * 92)
    say(f"pockets/predict.py -- fpocket {FPOCKET_VERSION} + P2Rank {P2RANK_VERSION} on AlphaFold "
        f"models; pockets admitted at mean pLDDT >= {PLDDT_ADMIT:.0f}")
    say("  in : data/processed/pockets/evidence/alphafold_<species>.tsv + data/source/alphafold/")
    say("  out: data/processed/pockets/evidence/pockets_<species>.tsv (long, one row per pocket)")
    say("=" * 92)

    summary = []
    for sp in args.species:
        say(f"\n[{sp}]")
        long = run(sp, args.workers, args.threads, args.limit)
        if args.limit:
            out = SCRATCH_DIR / f"smoke_pocket_list_{sp}.tsv"
        else:
            out = EVIDENCE_DIR / f"pocket_list_{sp}.tsv"
        long.to_csv(out, sep="\t", index=False)
        af = pd.read_csv(EVIDENCE_DIR / f"alphafold_{sp}.tsv", sep="\t")
        n_model = int((af["af_status"] == "model").sum()) if not args.limit else args.limit
        summary += summarise(sp, long, n_model)
        say(f"  wrote {out.relative_to(REPO_ROOT)}  ({len(long):,} pockets)")

    say("\n" + "-" * 92)
    say("SUMMARY -- pockets per tool, and how many survive the pLDDT admission")
    say(pd.DataFrame(summary).to_string(index=False, float_format=lambda x: f"{x:.1f}"))


if __name__ == "__main__":
    main()
