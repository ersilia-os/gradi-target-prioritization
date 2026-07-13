"""Structure-based pocket detection on AlphaFold models (docs §2.3b).

Runs BOTH fpocket (Voronoi/alpha-sphere druggability) and P2Rank (ML ligand-site predictor,
AlphaFold-tuned `-c alphafold` config) on every AlphaFold model of the proteome, then summarises
the best pocket per protein. AlphaFold pLDDT (stored in the model B-factor) is used to
**down-weight pockets that sit in low-confidence / disordered regions** — a pocket lining
made of low-pLDDT residues is likely an artefact, so the consensus score multiplies each
tool's score by the mean pLDDT fraction of its top pocket's residues (couples §2.3b with §2.4).

Pipeline (all cached + resumable, keyed by UniProt accession):
  1. cif -> pdb  (biotite; pLDDT preserved in B-factor)  -> data/processed/<org>/pockets/pdb/
  2. fpocket per structure   -> parsed json in       data/processed/<org>/pockets/fpocket/
  3. P2Rank in batch (one JVM over a dataset list) -> data/processed/<org>/pockets/p2rank/
  4. merge -> output/results/<org>/<prefix>_pockets.csv

External tools (osx-64 `gradi-pockets` env + P2Rank tarball); override paths via env vars
FPOCKET_BIN, P2RANK_DIR, POCKETS_JAVA_HOME. Run the SCRIPT itself with the `gradi` env.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src import ligandability as L  # noqa: E402

import biotite.structure.io.pdb as pdb  # noqa: E402
import biotite.structure.io.pdbx as pdbx  # noqa: E402

FPOCKET_BIN = os.environ.get(
    "FPOCKET_BIN", str(Path.home() / "miniconda3/envs/gradi-pockets/bin/fpocket")
)
P2RANK_DIR = Path(
    os.environ.get("P2RANK_DIR", str(L.REPO_ROOT / "tmp/tools/p2rank_2.5.1"))
)
JAVA_HOME = os.environ.get(
    "POCKETS_JAVA_HOME", str(Path.home() / "miniconda3/envs/gradi-pockets/lib/jvm")
)


# --------------------------------------------------------------------------- cif -> pdb
def cif_to_pdb(cif: Path, out_pdb: Path) -> bool:
    """Convert an AlphaFold cif to pdb, preserving pLDDT in the B-factor. False on failure."""
    try:
        f = pdbx.CIFFile.read(str(cif))
        arr = pdbx.get_structure(f, model=1, extra_fields=["b_factor"])
        out = pdb.PDBFile()
        pdb.set_structure(out, arr)
        out.write(str(out_pdb))
        return True
    except Exception as e:  # noqa: BLE001
        print(f"  [convert-fail] {cif.name}: {e}", flush=True)
        return False


def residue_plddt(pdb_path: Path) -> dict[int, float]:
    """Mean B-factor (== pLDDT) per residue number, parsed straight from the pdb."""
    acc: dict[int, list[float]] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    resnum = int(line[22:26])
                    b = float(line[60:66])
                except ValueError:
                    continue
                acc.setdefault(resnum, []).append(b)
    return {r: float(np.mean(v)) for r, v in acc.items() if v}


# --------------------------------------------------------------------------- fpocket
def parse_fpocket_info(info_path: Path) -> dict[int, dict]:
    """pocket_number -> {drug_score, volume, hydrophobicity} from <name>_info.txt."""
    pockets: dict[int, dict] = {}
    cur: int | None = None
    with open(info_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("Pocket"):
                try:
                    cur = int(s.split()[1])
                except (IndexError, ValueError):
                    cur = None
                if cur is not None:
                    pockets[cur] = {}
            elif cur is not None and ":" in s:
                key, _, val = s.partition(":")
                key = key.strip().lower()
                val = val.strip()
                try:
                    num = float(val)
                except ValueError:
                    continue
                if "druggability score" in key:
                    pockets[cur]["drug_score"] = num
                elif key.startswith("volume"):
                    pockets[cur]["volume"] = num
                elif "hydrophobicity score" in key:
                    pockets[cur]["hydrophobicity"] = num
    return pockets


def pocket_residue_nums(atm_pdb: Path) -> list[int]:
    nums: set[int] = set()
    with open(atm_pdb) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                try:
                    nums.add(int(line[22:26]))
                except ValueError:
                    pass
    return sorted(nums)


def run_fpocket_one(acc: str, pdb_path: Path, out_json: Path) -> None:
    """Run fpocket, summarise, cache a compact json, and delete the bulky _out dir."""
    out_dir = pdb_path.with_name(pdb_path.stem + "_out")
    if out_dir.exists():
        shutil.rmtree(out_dir, ignore_errors=True)
    try:
        subprocess.run(
            [FPOCKET_BIN, "-f", str(pdb_path)],
            check=True,
            capture_output=True,
            timeout=600,
        )
    except Exception as e:  # noqa: BLE001
        out_json.write_text(json.dumps({"accession": acc, "error": str(e)[:200], "pockets": []}))
        shutil.rmtree(out_dir, ignore_errors=True)
        return

    info = out_dir / f"{pdb_path.stem}_info.txt"
    rec = {"accession": acc, "pockets": []}
    if info.exists():
        plddt = residue_plddt(pdb_path)
        pmap = parse_fpocket_info(info)
        pock_dir = out_dir / "pockets"
        for pnum, props in pmap.items():
            atm = pock_dir / f"pocket{pnum}_atm.pdb"
            resnums = pocket_residue_nums(atm) if atm.exists() else []
            pl = [plddt[r] for r in resnums if r in plddt]
            rec["pockets"].append(
                {
                    "n": pnum,
                    "drug_score": props.get("drug_score", 0.0),
                    "volume": props.get("volume", 0.0),
                    "hydrophobicity": props.get("hydrophobicity", 0.0),
                    "n_residues": len(resnums),
                    "mean_plddt": round(float(np.mean(pl)), 2) if pl else None,
                }
            )
    out_json.write_text(json.dumps(rec))
    shutil.rmtree(out_dir, ignore_errors=True)


def _fpocket_worker(args):
    acc, pdb_path, out_json = args
    try:
        run_fpocket_one(acc, Path(pdb_path), Path(out_json))
        return acc, True
    except Exception as e:  # noqa: BLE001
        return acc, f"{e}"


# --------------------------------------------------------------------------- P2Rank
def run_p2rank_batch(pdb_paths: list[Path], out_dir: Path, threads: int) -> Path:
    """Run P2Rank once over a dataset list (amortises JVM/model warmup). Returns out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    ds = out_dir / "dataset.ds"
    ds.write_text("\n".join(str(p.resolve()) for p in pdb_paths) + "\n")
    env = dict(os.environ, JAVA_HOME=JAVA_HOME)
    cmd = [
        str(P2RANK_DIR / "prank"),
        "predict",
        str(ds),  # dataset list is a positional argument
        "-c", "alphafold",
        "-o", str(out_dir),
        "-threads", str(threads),
        "-visualizations", "0",
    ]
    print(f"  P2Rank: {len(pdb_paths)} structures, {threads} threads ...", flush=True)
    subprocess.run(cmd, check=True, env=env)
    return out_dir


def parse_p2rank_predictions(csv_path: Path, plddt: dict[int, float]) -> dict:
    """Top-pocket summary from a P2Rank *_predictions.csv (header has padded names)."""
    import csv as _csv

    rec = {"pockets": []}
    if not csv_path.exists():
        return rec
    with open(csv_path) as fh:
        reader = _csv.reader(fh)
        header = [h.strip() for h in next(reader, [])]
        idx = {h: i for i, h in enumerate(header)}
        for row in reader:
            if not row or len(row) < len(header):
                continue
            try:
                score = float(row[idx["score"]].strip())
                prob = float(row[idx["probability"]].strip())
            except (ValueError, KeyError):
                continue
            resids = row[idx["residue_ids"]].strip() if "residue_ids" in idx else ""
            nums = []
            for tok in resids.split():
                # tokens look like "A_123"
                part = tok.split("_")[-1]
                try:
                    nums.append(int(part))
                except ValueError:
                    pass
            pl = [plddt[n] for n in nums if n in plddt]
            rec["pockets"].append(
                {
                    "score": score,
                    "probability": prob,
                    "n_residues": len(nums),
                    "mean_plddt": round(float(np.mean(pl)), 2) if pl else None,
                }
            )
    return rec


# --------------------------------------------------------------------------- summarise
def summarise(acc: str, fp: dict, p2: dict, af_plddt: float | None) -> dict:
    row = {"uniprot_accession": acc}
    fpock = fp.get("pockets", []) if fp else []
    row["fpocket_n_pockets"] = len(fpock)
    if fpock:
        best = max(fpock, key=lambda p: p.get("drug_score", 0.0))
        row["fpocket_max_drug_score"] = best.get("drug_score", 0.0)
        row["fpocket_best_volume"] = best.get("volume", 0.0)
        row["fpocket_best_hydrophobicity"] = best.get("hydrophobicity", 0.0)
        row["fpocket_best_pocket_plddt"] = best.get("mean_plddt")
    else:
        row.update(
            fpocket_max_drug_score=0.0,
            fpocket_best_volume=0.0,
            fpocket_best_hydrophobicity=0.0,
            fpocket_best_pocket_plddt=None,
        )

    p2pock = p2.get("pockets", []) if p2 else []
    row["p2rank_n_pockets"] = len(p2pock)
    if p2pock:
        best = max(p2pock, key=lambda p: p.get("probability", 0.0))
        row["p2rank_top_score"] = best.get("score", 0.0)
        row["p2rank_top_prob"] = best.get("probability", 0.0)
        row["p2rank_top_pocket_plddt"] = best.get("mean_plddt")
    else:
        row.update(p2rank_top_score=0.0, p2rank_top_prob=0.0, p2rank_top_pocket_plddt=None)

    # pLDDT-weighted consensus: each tool's best-pocket score scaled by its pocket pLDDT
    # fraction; average the two (both in [0,1]).
    fp_pl = (row["fpocket_best_pocket_plddt"] or 0.0) / 100.0
    p2_pl = (row["p2rank_top_pocket_plddt"] or 0.0) / 100.0
    fp_w = float(row["fpocket_max_drug_score"]) * fp_pl
    p2_w = float(row["p2rank_top_prob"]) * p2_pl
    row["pocket_consensus_score"] = round((fp_w + p2_w) / 2.0, 4)
    row["af_mean_plddt"] = af_plddt
    return row


# --------------------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--organism", choices=list(L.ORGANISMS), default="kpneumoniae")
    ap.add_argument("--limit", type=int, default=None, help="first N proteins (debug)")
    ap.add_argument("--workers", type=int, default=6, help="fpocket parallel processes")
    ap.add_argument("--p2rank-threads", type=int, default=6)
    ap.add_argument("--skip-p2rank", action="store_true")
    ap.add_argument("--skip-fpocket", action="store_true")
    args = ap.parse_args()

    org = args.organism
    _, prefix = L.ORGANISMS[org]
    accs = L.load_accessions(org)
    if args.limit:
        accs = accs[: args.limit]

    pdb_dir = L.processed_dir(org, "pockets", "pdb")
    fp_dir = L.processed_dir(org, "pockets", "fpocket")
    p2_dir = L.processed_dir(org, "pockets", "p2rank")
    p2_run_dir = L.processed_dir(org, "pockets", "p2rank_run")

    # ---- 1. cif -> pdb (only those with an AF model + not already converted)
    pdb_paths: dict[str, Path] = {}
    to_convert = []
    for acc in accs:
        cif = L.af_cif_path(org, acc)
        if not cif.exists():
            continue
        pp = pdb_dir / f"{acc}.pdb"
        pdb_paths[acc] = pp
        if not pp.exists():
            to_convert.append((acc, cif, pp))
    print(f"[{org}] {len(pdb_paths)} AF models; {len(to_convert)} to convert", flush=True)
    for i, (acc, cif, pp) in enumerate(to_convert, 1):
        cif_to_pdb(cif, pp)
        if i % 500 == 0:
            print(f"  converted {i}/{len(to_convert)}", flush=True)
    pdb_paths = {a: p for a, p in pdb_paths.items() if p.exists()}

    # ---- 2. fpocket (parallel, resumable via per-acc json)
    if not args.skip_fpocket:
        todo = [
            (a, str(pdb_paths[a]), str(fp_dir / f"{a}.json"))
            for a in pdb_paths
            if not (fp_dir / f"{a}.json").exists()
        ]
        print(f"[{org}] fpocket: {len(todo)} to run ({args.workers} workers)", flush=True)
        done = 0
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = [ex.submit(_fpocket_worker, t) for t in todo]
            for fut in as_completed(futs):
                done += 1
                if done % 250 == 0:
                    print(f"  fpocket {done}/{len(todo)}", flush=True)

    # ---- 3. P2Rank (single batch over the still-missing structures)
    if not args.skip_p2rank:
        missing = [pdb_paths[a] for a in pdb_paths if not (p2_dir / f"{a}.json").exists()]
        if missing:
            run_p2rank_batch(missing, p2_run_dir, args.p2rank_threads)
            pred_root = p2_run_dir
            for a in list(pdb_paths):
                if (p2_dir / f"{a}.json").exists():
                    continue
                # P2Rank writes <name>_predictions.csv (name = pdb filename)
                csv_path = _find_predictions(pred_root, f"{a}.pdb")
                plddt = residue_plddt(pdb_paths[a]) if csv_path else {}
                rec = parse_p2rank_predictions(csv_path, plddt) if csv_path else {"pockets": []}
                rec["accession"] = a
                (p2_dir / f"{a}.json").write_text(json.dumps(rec))
            print(f"[{org}] P2Rank parsed -> {p2_dir}", flush=True)

    # ---- 4. merge
    af_plddt = _load_af_plddt(org, prefix)
    rows = []
    for acc in pdb_paths:
        fp = _read_json(fp_dir / f"{acc}.json")
        p2 = _read_json(p2_dir / f"{acc}.json")
        rows.append(summarise(acc, fp, p2, af_plddt.get(acc)))
    import pandas as pd

    out = L.results_dir(org) / f"{prefix}_pockets.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"[{org}] wrote {out} ({len(rows)} proteins)", flush=True)


def _find_predictions(root: Path, pdb_name: str) -> Path | None:
    for cand in root.rglob(f"{pdb_name}_predictions.csv"):
        return cand
    return None


def _read_json(p: Path) -> dict:
    try:
        return json.loads(p.read_text())
    except Exception:  # noqa: BLE001
        return {}


def _load_af_plddt(org: str, prefix: str) -> dict[str, float]:
    import pandas as pd

    p = L.results_dir(org) / f"{prefix}_alphafold_structure.csv"
    if not p.exists():
        return {}
    df = pd.read_csv(p, usecols=["uniprot_accession", "af_mean_plddt"])
    return dict(zip(df["uniprot_accession"], df["af_mean_plddt"]))


if __name__ == "__main__":
    main()
