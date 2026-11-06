"""Can a de-novo function predictor annotate the 2,996 proteins stage 02 cannot?

Stage 02 annotates GO in two tiers -- `curated` (UniProt) then `eggnog` -- and leaves **2,996 of
13,020 proteins (23%) with no GO term at all** (Kp 1,506 · Ec 496 · Sa 994). The eggNOG tier bought
only **322 proteins for 50.6 GB**, and `docs/function.md` records why: ~73% of the unannotated DO
get an orthologous group, but **only 2-18% of those groups carry any GO**. Re-verified for the
current gap: ~70% have an eggNOG OG and just 18 / 26 / 8 of those groups carry GO.

**Every route tried so far asks "what is this protein similar to?"** -- eggNOG, InterPro2GO and the
rejected ESM-C k-NN all transfer annotation from something else, so all three inherit the same
failure. This stage tests predictors that assign function **de novo from sequence**, needing no
lookup source to carry GO.

`--tool sprofgo` (default) runs **SPROF-GO** (Yuan et al., *Brief Bioinform* 2023): ProtT5-XL-U50
embeddings + self-attention pooling, then homology-based label diffusion. `--tool deepgometa` is
reserved for the escalation (DeepGOMeta, prokaryote-trained) and reuses this entire harness.

THIS SCRIPT DECIDES NOTHING BY ITSELF -- it writes to `evidence/ + scratch/` only. Nothing touches
`goslim_<species>.tsv` until the numbers below justify a third tier.

THE CONTROL IS K. PNEUMONIAE, NOT E. COLI -- THIS INVERTS STAGE 02's CONVENTION
------------------------------------------------------------------------------
These tools train on Swiss-Prot-derived data. Measured reviewed status of our proteomes:
**Ec 4,403/4,403 (100%) · Sa 815 (28.2%) · Kp 7 (0.1%)**. Stage 02 controls on E. coli *because* it
is best annotated; here that is backwards -- E. coli proteins are the likeliest to sit inside the
training set, and accuracy measured there would be inflated exactly as Geptop's self-scoring was.
**Kp is the clean control**: its curated GO is electronic (InterPro2GO-derived), independent of an
experimental training signal, and it is the anchor organism. E. coli is still reported, labelled
LIKELY CIRCULAR, so the inflation is visible rather than hidden.

THE THREE NUMBERS
-----------------
  (a) accuracy on the Kp control, head-to-head with **eggNOG scored on the same proteins** -- stage
      02 only ever scored eggNOG on E. coli (MF 80.8 / BP 82.4 / CC 81.1), so the like-for-like Kp
      number is computed here. That comparison, not a published Fmax, is the decision.
  (b) yield: how many of the 2,996 gain a slim term. eggNOG's answer was 322.
  (c) **the subgroup that has defeated everything**: accuracy and yield split by whether the protein
      has an eggNOG OG that carries no GO (~70% of the gap). A tool that fails there has inherited
      the same limitation and the download buys nothing, whatever its headline.

RUNNING IT
----------
SPROF-GO lives in `gradi-sprofgo` and is driven across a process boundary (`GRADI_SPROFGO_BIN`) --
the mandatory env split; it pins transformers 4.17, which would wreck `gradi`. Three patches to the
vendored copy, each recorded in `data/source/sprofgo/SOURCE.md` and each validated against
the authors' own bundled demo, which our output reproduces **byte-identically**:
its hard-coded `ProtTrans_path`, CUDA-only device selection (now reaches MPS), and a float64
constraint matrix MPS cannot hold.

**It caps at 5,000 sequences per call**, so input is batched. Measured ~1.6 s/protein on MPS
(1.8x faster than CPU, and verified to give the same terms -- max score delta 0.001, which is the
output format's own rounding). Expect ~5-6 h for 13,020.

    python scripts/function/deepgo.py --predict          # the long part, batched, resumable
    python scripts/function/deepgo.py --analyse          # the three numbers, seconds
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
from src import function as F  # noqa: E402
from src import proteomes as P  # noqa: E402

RAW = REPO_ROOT / "data" / "source"
SPROFGO_DIR = RAW / "sprofgo" / "SPROF-GO"
PROTT5_DIR = RAW / "sprofgo" / "prot_t5_xl_uniref50"
OUT_DIR = REPO_ROOT / "data" / "processed" / "function"
EVIDENCE_DIR = OUT_DIR / "evidence"
SCRATCH_DIR = OUT_DIR / "scratch"
RUN_DIR = SCRATCH_DIR / "deepgo"
DEFAULT_SPROFGO_BIN = Path.home() / "miniconda3" / "envs" / "gradi-sprofgo" / "bin" / "python"

SPECIES = ("kpneumoniae", "ecoli", "saureus")
BATCH = 4000              # the tool refuses >5000 per call
# T5 self-attention is O(L^2). Measured: a 3,163-aa protein runs fine on MPS, but S. aureus Q2FYJ6
# (ebh, 9,535 aa) asks for a 23.85 GiB buffer and kills the run. Exactly ONE protein of 13,020
# exceeds 4,000 aa, and it is already curated -- so it is EXCLUDED and reported, not truncated:
# a truncated prediction would be a silent fabrication for a protein we do not even need.
MAX_LEN = 4000
ASPECTS = ("MF", "BP", "CC")
# Stage 02's measured eggNOG control on E. coli -- the bar, and the reason we recompute it on Kp.
EGGNOG_ECOLI_CONTROL = {"MF": 0.808, "BP": 0.824, "CC": 0.811}
EGGNOG_YIELD = 322

VERBOSE = True


def say(msg: str = "") -> None:
    if VERBOSE:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    say(char * width)


def stage02():
    """Import function/goslim.py for its goatools machinery -- do NOT reimplement the slim map."""
    spec = importlib.util.spec_from_file_location(
        "s02", REPO_ROOT / "scripts" / "function/goslim.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    for name in ("build_slim_mapper", "split_go"):
        if not hasattr(m, name):
            sys.exit(f"FATAL function/goslim.py has no {name} -- it was refactored")
    return m


def sprofgo_bin() -> str:
    p = Path(os.environ.get("GRADI_SPROFGO_BIN", DEFAULT_SPROFGO_BIN))
    if not p.exists():
        sys.exit(f"FATAL no sprofgo interpreter at {p}. See install.sh. Do NOT install it into "
                 "`gradi`: it pins transformers 4.17.")
    if not (PROTT5_DIR / "pytorch_model.bin").exists():
        sys.exit(f"FATAL ProtT5 weights missing at {PROTT5_DIR}. SPROF-GO needs the FULL "
                 "prot_t5_xl_uniref50 (11.3 GB); our half encoder-only build is NOT the same "
                 "weights and must not be substituted.")
    return str(p)


def predict(species: tuple[str, ...], use_gpu: bool, refresh: bool) -> None:
    """Run SPROF-GO over every proteome, in batches of BATCH, skipping batches already done."""
    binary = sprofgo_bin()
    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)
    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    prot = pd.concat([P.load(sp)[["uniprot_ac", "sequence"]].assign(species=sp)
                      for sp in species], ignore_index=True)
    too_long = prot[prot.sequence.str.len() > MAX_LEN]
    if len(too_long):
        say(f"  EXCLUDED {len(too_long)} protein(s) longer than {MAX_LEN:,} aa -- T5 attention is "
            f"O(L^2) and these blow past addressable memory:")
        for _, r in too_long.iterrows():
            say(f"    {r.uniprot_ac}  {len(r.sequence):,} aa  ({r.species})")
        too_long[["uniprot_ac", "species"]].assign(
            length=too_long.sequence.str.len(), reason=f"longer than {MAX_LEN} aa").to_csv(
            EVIDENCE_DIR / "deepgo_excluded_too_long.tsv", sep="\t", index=False)
        prot = prot[prot.sequence.str.len() <= MAX_LEN].reset_index(drop=True)
    say(f"  {len(prot):,} proteins over {len(species)} species, batches of {BATCH:,}")

    n_batches = (len(prot) + BATCH - 1) // BATCH
    for b in range(n_batches):
        chunk = prot.iloc[b * BATCH:(b + 1) * BATCH]
        outdir = RUN_DIR / f"run_{b:02d}"
        done = outdir / f"batch_{b:02d}_all_preds.txt"
        manifest = outdir / f"batch_{b:02d}.accessions.txt"
        want = "\n".join(chunk.uniprot_ac)
        # Verify the cached batch holds the SAME proteins. Excluding a protein shifts every batch
        # after it, so an index-only check would happily reuse the wrong output.
        if done.exists() and manifest.exists() and manifest.read_text() == want and not refresh:
            say(f"    batch {b + 1}/{n_batches}: cached ({len(chunk):,} proteins, membership verified)")
            continue
        if done.exists() and not refresh:
            say(f"    batch {b + 1}/{n_batches}: cached output exists but membership CHANGED -- rerunning")
        outdir.mkdir(parents=True, exist_ok=True)
        manifest.write_text(want)
        fa = outdir / f"batch_{b:02d}.fa"
        with open(fa, "w") as fh:
            for a, s in zip(chunk.uniprot_ac, chunk.sequence):
                fh.write(f">{a}\n{s}\n")
        say(f"    batch {b + 1}/{n_batches}: {len(chunk):,} proteins -> running")
        cmd = [binary, "./script/predict.py", "--fasta", str(fa.resolve()),
               "--outpath", str(outdir.resolve()) + "/"]
        if use_gpu:
            cmd.append("--gpu")
        env = dict(os.environ, GRADI_PROTT5_PATH=str(PROTT5_DIR.resolve()))
        r = subprocess.run(cmd, cwd=SPROFGO_DIR, env=env, capture_output=True, text=True)
        if r.returncode != 0 or not done.exists():
            sys.exit(f"FATAL SPROF-GO failed on batch {b} (exit {r.returncode}):\n"
                     f"{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
        say(f"      done -> {done.relative_to(REPO_ROOT)}")


def parse_all_preds(path: Path, floor: float) -> dict[str, dict[str, float]]:
    """Parse SPROF-GO's `*_all_preds.txt` into {protein: {GO: score}}, sparsified at `floor`.

    Format: a header giving the per-aspect GO vocabulary (MF 790, BP 4,766, CC 667 terms), then four
    lines per protein -- id, then `MF:`/scores, `BP:`/scores, `CC:`/scores. Sparsifying on read
    matters: 13,020 x 6,223 dense floats is ~81M numbers, almost all zero.
    """
    lines = path.read_text().splitlines()
    vocab: dict[str, list[str]] = {}
    i = 0
    while i < len(lines):
        if lines[i].rstrip(":") in ASPECTS and i + 1 < len(lines) and lines[i + 1].startswith("GO:"):
            vocab[lines[i].rstrip(":")] = [x.strip() for x in lines[i + 1].split(";")]
            i += 3            # aspect, ids, names
        else:
            i += 1
        if len(vocab) == len(ASPECTS):
            break
    if len(vocab) != len(ASPECTS):
        sys.exit(f"FATAL could not read the GO vocabulary header from {path}")

    out: dict[str, dict[str, float]] = {}
    cur, aspect = None, None
    for line in lines[i:]:
        line = line.rstrip("\n")
        if not line:
            continue
        stripped = line.rstrip(":")
        if stripped in ASPECTS and len(line) == len(stripped) + 1:
            aspect = stripped
        elif line.startswith("GO:"):
            continue                      # a vocabulary line, already consumed
        elif ";" in line and cur is not None and aspect is not None:
            vals = line.split(";")
            terms = vocab[aspect]
            if len(vals) != len(terms):
                sys.exit(f"FATAL {path}: {len(vals)} scores for {len(terms)} {aspect} terms")
            for go, v in zip(terms, vals):
                v = float(v)
                if v >= floor:
                    out[cur][go] = v
        else:
            cur = line.strip()
            out.setdefault(cur, {})
    return out


# ---------------------------------------------------------------- the three numbers

def _slim(map_terms, go_ids):
    """GO ids -> per-aspect slim terms, via stage 02's own goatools mapper."""
    return map_terms(sorted(set(go_ids)))


def _score(pairs, aspects) -> dict:
    """Stage 02's control metric, unchanged: `any slim term shared` and mean Jaccard.

    Counted ONLY where the truth has a term in that aspect -- identical to `control_eggnog`, so
    these numbers sit directly beside eggNOG's and the rejected ESM-C tier's.
    """
    st = {a: {"n": 0, "exact": 0, "jac": 0.0} for a in aspects}
    for truth, pred in pairs:
        for a in aspects:
            if not truth[a]:
                continue
            st[a]["n"] += 1
            st[a]["exact"] += int(bool(set(pred[a]) & set(truth[a])))
            st[a]["jac"] += _jaccard(pred[a], truth[a])
    return st


def _jaccard(a, b) -> float:
    sa, sb = set(a), set(b)
    return len(sa & sb) / len(sa | sb) if (sa or sb) else 0.0


def _show(label: str, st: dict, aspects) -> None:
    say(f"    {label}")
    say(f"      {'aspect':<6} {'n':>7} {'any slim term shared':>22} {'mean Jaccard':>14}")
    for a in aspects:
        s_ = st[a]
        if s_["n"]:
            say(f"      {a:<6} {s_['n']:>7,} {100 * s_['exact'] / s_['n']:>21.1f}% "
                f"{s_['jac'] / s_['n']:>14.3f}")


def analyse(preds: dict[str, dict[str, float]]) -> None:
    """Accuracy on the clean control, yield on the gap, and the subgroup that decides it."""
    s02 = stage02()
    aspects = list(s02.ASPECTS.values())
    # build_slim_mapper returns (map_terms, meta, dag, slim_ids) -- its docstring says two, so
    # take the first element by index rather than unpacking a shape that may change again.
    _mapper = s02.build_slim_mapper(RAW / "go" / "go-basic.obo",
                                    RAW / "go" / "goslim_prokaryote.obo")
    map_terms = _mapper[0] if isinstance(_mapper, tuple) else _mapper
    split_go = s02.split_go

    rows = []
    for sp in SPECIES:
        gs = F.load_goslim(sp)[["uniprot_ac", "goslim_source"]]
        ann = P.load_annotation(sp)[["uniprot_ac", "go_id"]]
        eg = pd.read_csv(EVIDENCE_DIR / f"eggnog_{sp}.tsv", sep="\t", dtype=str,
                         keep_default_na=False)[["uniprot_ac", "eggnog_ogs", "gos"]]
        d = gs.merge(ann, on="uniprot_ac", how="left").merge(eg, on="uniprot_ac", how="left")
        d["species"] = sp
        rows.append(d)
    tab = pd.concat(rows, ignore_index=True).fillna("")

    # ---- (a) accuracy, on the CLEAN control -------------------------------------------------
    rule()
    say("(a) ACCURACY -- K. pneumoniae is the control; E. coli is reported but LIKELY CIRCULAR")
    rule()
    say("    Both tools are Swiss-Prot-trained. Reviewed fraction: Kp 0.1%, Sa 28.2%, Ec 100%.")
    say("    Scoring on E. coli would largely be scoring on the training set.")
    say()
    # **The metric MUST be applied at a threshold.** SPROF-GO emits a score for all 6,223 GO
    # terms, so at a low cut it predicts ~342 terms per protein and "any slim term shared" is
    # trivially 100% -- measured, and meaningless. eggNOG emits a specific set, so it is scored
    # once. Sweeping the cut is the only way to compare them fairly; Jaccard is the metric that
    # penalises over-prediction and is therefore the one to read.
    for sp, note in (("kpneumoniae", "THE CONTROL (clean)"), ("ecoli", "LIKELY CIRCULAR")):
        sub = tab[(tab.species == sp) & (tab.go_id != "")]
        sub = sub[sub.uniprot_ac.isin(preds)]
        if sub.empty:
            say(f"  {sp}: no overlap between curated GO and predictions"); continue
        say(f"  {sp} -- {note}   n={len(sub):,} proteins with curated GO and a prediction")
        truths = {ac: _slim(map_terms, split_go(g)) for ac, g in zip(sub.uniprot_ac, sub.go_id)}

        say(f"      {'cut':>5} {'terms/prot':>11} " + " ".join(
            f"{a + ' shared':>12} {a + ' Jacc':>10}" for a in aspects))
        for cut in (0.01, 0.1, 0.3, 0.5, 0.7, 0.9):
            pairs, nterm = [], []
            for ac in sub.uniprot_ac:
                keep = [g for g, v in preds[ac].items() if v >= cut]
                pred = _slim(map_terms, keep)
                nterm.append(sum(len(v) for v in pred.values()))
                pairs.append((truths[ac], pred))
            st = _score(pairs, aspects)
            cells = []
            for a in aspects:
                s_ = st[a]
                cells.append(f"{100 * s_['exact'] / s_['n']:>11.1f}% {s_['jac'] / s_['n']:>10.3f}"
                             if s_["n"] else f"{'-':>12} {'-':>10}")
            say(f"      {cut:>5.2f} {np.mean(nterm):>11.1f} " + " ".join(cells))

        egg = [(truths[ac], _slim(map_terms, split_go(g)))
               for ac, g in zip(sub.uniprot_ac, sub.gos) if g]
        if egg:
            st = _score(egg, aspects)
            nterm = np.mean([sum(len(v) for v in p.values()) for _t, p in egg])
            cells = []
            for a in aspects:
                s_ = st[a]
                cells.append(f"{100 * s_['exact'] / s_['n']:>11.1f}% {s_['jac'] / s_['n']:>10.3f}"
                             if s_["n"] else f"{'-':>12} {'-':>10}")
            say(f"      {'eggNOG':>5} {nterm:>11.1f} " + " ".join(cells)
                + f"   (n={len(egg):,}, recomputed on THESE proteins)")
        say()
    say("    reference: eggNOG on E. coli MF 80.8 / BP 82.4 / CC 81.1;")
    say("               the REJECTED ESM-C k-NN, reweighted: MF 68 / BP 53 / CC 66")

    # ---- (b) yield, with a threshold sweep --------------------------------------------------
    rule()
    say("(b) YIELD -- of the 2,996 with no GO at all, how many gain a slim term?")
    rule()
    dark = tab[(tab.goslim_source == "") & (tab.uniprot_ac.isin(preds))]
    say(f"    {len(dark):,} of 2,996 unannotated proteins have a prediction")
    say(f"      {'cut':>6} {'gain a slim term':>18} {'vs eggNOG 322':>16}")
    sweep = []
    for cut in (0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7):
        n = sum(1 for ac in dark.uniprot_ac
                if any(_slim(map_terms, [g for g, v in preds[ac].items() if v >= cut]).values()))
        sweep.append({"cut": cut, "n_gained": n})
        say(f"      {cut:>6.2f} {n:>18,} {n - EGGNOG_YIELD:>+16,}")
    pd.DataFrame(sweep).to_csv(EVIDENCE_DIR / "deepgo_yield_sweep.tsv", sep="\t", index=False)

    # ---- (c) the subgroup that has defeated everything --------------------------------------
    rule()
    say("(c) THE DECIDING SUBGROUP -- proteins whose eggNOG group carries NO GO")
    rule()
    say("    ~70% of the gap. eggNOG, InterPro2GO and the ESM-C k-NN all failed here; a tool that")
    say("    fails too has inherited the same limitation, whatever its headline number.")
    say()
    has_og = dark[(dark.eggnog_ogs != "") & (dark.gos == "")]
    no_og = dark[dark.eggnog_ogs == ""]
    say(f"      {'subgroup':<34} {'n':>7} {'gain a slim term @0.1':>23}")
    for label, grp in (("has eggNOG OG, OG carries no GO", has_og),
                       ("no eggNOG OG at all", no_og)):
        if not len(grp):
            continue
        n = sum(1 for ac in grp.uniprot_ac
                if any(_slim(map_terms, [g for g, v in preds[ac].items() if v >= 0.1]).values()))
        say(f"      {label:<34} {len(grp):>7,} {n:>22,} ({100 * n / len(grp):.1f}%)")

    EVIDENCE_DIR.mkdir(parents=True, exist_ok=True)

    SCRATCH_DIR.mkdir(parents=True, exist_ok=True)
    say()
    say(f"  wrote {(RUN_DIR / 'yield_sweep.tsv').relative_to(REPO_ROOT)}")


def load_predictions(floor: float) -> dict[str, dict[str, float]]:
    """Every batch's predictions, merged. Fails loudly if a batch is missing."""
    runs = sorted(RUN_DIR.glob("run_*/batch_*_all_preds.txt"))
    if not runs:
        sys.exit(f"FATAL no predictions under {RUN_DIR}. Run with --predict first.")
    preds: dict[str, dict[str, float]] = {}
    for r in runs:
        part = parse_all_preds(r, floor)
        preds.update(part)
        say(f"    {r.parent.name}: {len(part):,} proteins")
    return preds


def main() -> None:
    global VERBOSE
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tool", default="sprofgo", choices=["sprofgo", "deepgometa"])
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--predict", action="store_true", help="run the predictor (the long part)")
    ap.add_argument("--analyse", action="store_true", help="the three numbers, from cached output")
    ap.add_argument("--floor", type=float, default=0.01,
                    help="sparsify scores below this when parsing (does not set the decision cut)")
    ap.add_argument("--no-gpu", action="store_true", help="force CPU (MPS is ~1.8x faster)")
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    VERBOSE = not args.quiet
    if args.tool == "deepgometa":
        sys.exit("deepgometa is the ESCALATION path and is not wired yet -- run sprofgo first and "
                 "let its numbers decide (see the module docstring).")

    rule("=")
    say("STAGE 02 - de novo GO prediction for the proteins homology cannot reach")
    rule("=")
    say(f"  tool         {args.tool} (SPROF-GO: ProtT5-XL-U50 + attention + label diffusion)")
    say(f"  target       2,996 proteins with NO GO after curated+eggnog (Kp 1,506 · Ec 496 · Sa 994)")
    say(f"  control      K. pneumoniae -- NOT E. coli, which is 100% reviewed and likely circular")
    say(f"  bar          eggNOG on E. coli: MF {EGGNOG_ECOLI_CONTROL['MF']:.3f} "
        f"BP {EGGNOG_ECOLI_CONTROL['BP']:.3f} CC {EGGNOG_ECOLI_CONTROL['CC']:.3f}; "
        f"yield {EGGNOG_YIELD}")
    say("  writes       evidence/ + scratch/ only -- this script promotes nothing to a tier")
    rule("=")
    if args.dry_run:
        say("dry run: nothing computed.")
        return

    if args.predict:
        rule()
        say("PREDICT - batched, resumable; ~1.6 s/protein on MPS")
        rule()
        predict(tuple(args.species), use_gpu=not args.no_gpu, refresh=args.refresh)
        say("  prediction complete.")

    if args.analyse:
        rule()
        say("LOAD")
        rule()
        preds = load_predictions(args.floor)
        say(f"  {len(preds):,} proteins with predictions (scores >= {args.floor})")
        analyse(preds)

    if not (args.predict or args.analyse):
        say("nothing to do: pass --predict and/or --analyse")


if __name__ == "__main__":
    main()
