"""ProteomeLM-Ess: the paper's own essentiality head, applied to our three anchors.

    data/processed/essentiality/proteomelm_ess_<species>.tsv
        uniprot_ac
        proteomelm_ess            p(essential), 0-1, NEVER null -- the head scores every protein
        proteomelm_ess_rank       1 = most essential, within this proteome only
        proteomelm_ess_evidence   held_out | in_training | unseen_species   <- READ THIS FIRST

This is the fourth evidence source on the axis, beside `geptop_` (prediction), `deg_` (measured)
and `ogee_` (prediction). It was deferred for months because the weights were unreleased; Cyril
Malbranke released them on 2026-10-01 (`Bitbol-Lab/ProteomeLM-ess`, **Apache-2.0** -- note that is
a far easier licence than TabPFN's non-commercial weights, so this column can ship).

WHAT THE NUMBER MEANS DEPENDS ON THE SPECIES, which is why `proteomelm_ess_evidence` is not
decoration. From the authors' own `data/genomes.tsv` (staged at `data/source/proteomelm/`):

    ecoli        held_out        taxid 83333, 290 E / 3,969 NE -- their Fig. 5B held-out genome,
                                 so their reported AUROC 0.952 is honest ON OUR EXACT ANCHOR.
    saureus      in_training     taxid 93061, 2,889 proteins (our exact proteome), 395 E / 2,484 NE,
                                 in their cross-validation set. The column is closer to RECALL than
                                 prediction here -- do not quote it as out-of-sample performance.
    kpneumoniae  unseen_species  no K. pneumoniae anywhere in their 89 genomes (the only Klebsiella
                                 is K. michiganensis M5al, and it is positives-only). This is a
                                 genuine out-of-distribution prediction -- and it is the anchor.

So `proteomelm_ess` is comparable WITHIN a species and not across one, the same rule `ogee_ess`
carries for the same reason.

WHY THIS DOES NOT REUSE `proteomelm_<species>.npz`
--------------------------------------------------
It is the same backbone and the same layer 8, and our matrices were built with the tool's own code
path -- but they are NOT what this head was fitted on, and feeding them would be the silent failure
CLAUDE.md's *External models* section exists to prevent. Three differences, all measured by reading
the authors' `proteomelm/essentiality.py`:

    z-scoring     ours z-scores genome-wide; the head takes RAW `hidden_states[8]`
    truncation    ours WINDOWS sequences above 4,096 aa; the head TRUNCATES to the first 4,096
    dtype         the head runs the backbone in bfloat16

(The one difference that turned out NOT to exist: their `group_embeds=x` and our `self` mode are
the same input, because `modeling_proteomelm.py` does `if group_embeds is None: group_embeds =
inputs_embeds.clone()`. Checked rather than assumed, since a mismatch there would have been
invisible.)

So the worker recomputes ESM-C and the backbone pass itself, end to end, from our FASTA.

ENV: the head ships in the authors' **git** build, which pulls torch 2.14 against `gradi`'s 2.12 --
the collision that breaks stage 01's ESM-C. It therefore lives in **`gradi-plm-ess`** and is reached
across a process boundary (`scripts/essentiality/workers/proteomelm_ess.py`, `GRADI_PLM_ESS_BIN`
overrides). Never activate that env to run this stage.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src import matrices as M  # noqa: E402
from src import proteomes as P  # noqa: E402

SPECIES = ("kpneumoniae", "ecoli", "saureus")
OUT_DIR = REPO_ROOT / "data" / "processed" / "essentiality"
EVIDENCE = OUT_DIR / "evidence"
SCRATCH = OUT_DIR / "scratch"
WORKER = REPO_ROOT / "scripts" / "essentiality" / "workers" / "proteomelm_ess.py"
DEFAULT_BIN = Path.home() / "miniconda3" / "envs" / "gradi-plm-ess" / "bin"
HEAD = "Bitbol-Lab/ProteomeLM-ess"

# From the authors' own data/genomes.tsv -- staged at data/source/proteomelm/ess_genomes.tsv.
# DECLARED, never inferred: getting this wrong silently turns a recall number into a prediction.
ROLE = {"ecoli": "held_out", "saureus": "in_training", "kpneumoniae": "unseen_species"}

# The biological control is IMPORTED from summary.py rather than restated, so the two cannot drift
# -- the same reason `ligands/effort.py` spec-loads its predicates from `chembl.py`.
_spec = importlib.util.spec_from_file_location(
    "_ess_summary", REPO_ROOT / "scripts" / "essentiality" / "summary.py")
_summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_summary)
RIBOSOME: re.Pattern = _summary.RIBOSOME
DISPENSABLE: set = _summary.DISPENSABLE

# Calibrated on the first full run, NOT set to 1.0. See the run log in docs/essentiality.md: a
# continuous score cannot put every ribosomal protein in the top decile, and demanding it would
# fail a working model -- the mistake `summary.py` already records for the screens.
MIN_RIBOSOME_AUROC = 0.80


def say(msg: str = "", quiet: bool = False) -> None:
    if not quiet:
        print(msg, flush=True)


def rule(char: str = "-", width: int = 92) -> None:
    print(char * width, flush=True)


def worker_bin() -> Path:
    """The `gradi-plm-ess` interpreter, or exit telling the user how to make one."""
    bin_dir = Path(os.environ.get("GRADI_PLM_ESS_BIN", DEFAULT_BIN))
    python = bin_dir / "python"
    if not python.exists():
        sys.exit(
            f"FAILED: no interpreter at {python}.\n"
            "  ProteomeLM-Ess needs the authors' git build, which pulls torch 2.14 and must NOT go\n"
            "  into `gradi`. Create it with:\n"
            "    conda create -y -n gradi-plm-ess python=3.11\n"
            '    ~/miniconda3/envs/gradi-plm-ess/bin/pip install "proteomelm @ '
            'git+https://github.com/Bitbol-Lab/ProteomeLM" httpx\n'
            "  or point GRADI_PLM_ESS_BIN at an existing bin/ directory.")
    return python


def write_fasta(df: pd.DataFrame, path: Path) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for acc, seq in zip(df["uniprot_ac"], df["sequence"]):
            fh.write(f">{acc}\n{seq}\n")
    return len(df)


def score_fasta(fasta: Path, out: Path, manifest: Path, device: str, quiet: bool) -> pd.DataFrame:
    """Run the worker across the process boundary and read its TSV back."""
    cmd = [str(worker_bin()), str(WORKER), "--fasta", str(fasta), "--out", str(out),
           "--manifest", str(manifest), "--head", HEAD, "--device", device]
    say(f"    -> {' '.join(cmd[:2])} ... --device {device}", quiet)
    proc = subprocess.run(cmd, capture_output=quiet)
    if proc.returncode != 0:
        if quiet and proc.stderr:
            sys.stderr.write(proc.stderr.decode(errors="replace"))
        sys.exit(f"FAILED: the gradi-plm-ess worker exited {proc.returncode} for {fasta.name}.")
    return pd.read_csv(out, sep="\t")


def control(df: pd.DataFrame, species: str) -> dict:
    """Is the column the right way up? Ribosome vs textbook dispensables, on the shipped score.

    An inverted essentiality column is well-formed, confident and exactly wrong, and nothing about
    its distribution says so -- this is the same polarity test `summary.py` applies to every screen,
    reusing that file's own panels.
    """
    names = P.load(species)[["uniprot_ac", "gene_name"]]
    d = df.merge(names, on="uniprot_ac", how="left")
    sym = d["gene_name"].fillna("").astype(str)

    ribo = d.loc[sym.str.match(RIBOSOME), "proteomelm_ess"]
    disp = d.loc[sym.str.lower().isin(DISPENSABLE), "proteomelm_ess"]
    top = d["proteomelm_ess"] >= d["proteomelm_ess"].quantile(0.90)

    auroc = ""
    if len(ribo) and len(disp):
        # Mann-Whitney U / (n1*n2): P(a ribosomal protein outranks a dispensable one).
        both = pd.concat([ribo, disp])
        ranks = both.rank()
        u = ranks[: len(ribo)].sum() - len(ribo) * (len(ribo) + 1) / 2
        auroc = round(float(u / (len(ribo) * len(disp))), 4)

    return {
        "species": species,
        "role": ROLE[species],
        "n": len(d),
        "n_ribosomal": len(ribo),
        "ribosome_median": round(float(ribo.median()), 4) if len(ribo) else "",
        "ribosome_top_decile": round(float(top[sym.str.match(RIBOSOME)].mean()), 3) if len(ribo) else "",
        "n_dispensable": len(disp),
        "dispensable_median": round(float(disp.median()), 4) if len(disp) else "",
        "dispensable_top_decile": round(float(top[sym.str.lower().isin(DISPENSABLE)].mean()), 3) if len(disp) else "",
        "ribosome_vs_dispensable_auroc": auroc,
    }


# Measured screens that key DIRECTLY onto an anchor proteome (`features_from == <species>`), so
# they can be scored with no join at all. The Kp screens are not here because they key onto their
# own strain proteomes -- `--validate-strains` scores those strains themselves instead, which is
# the honest test and needs no identifier mapping either.
ANCHOR_SCREENS = {
    "ecoli": ["essential_ecoli_k12_knockout", "essential_ecoli_bw25113_tradis_goodall",
              "essential_ecoli_bw25113_tnseq_choe", "core_essential_gammaproteobacteria"],
}

# Kp screens: (training set, the strain proteome its keys live in). Scoring the STRAIN rather than
# transferring its labels to HS11286 is deliberate -- CLAUDE.md records that exact-sequence label
# transfer recovers only 13.7-71.6% of rows, so it would price the column with a decimated join.
STRAIN_SCREENS = [
    ("essential_kpneumoniae_ecl8_tradis", "kpneumoniae__ecl8__GCA_000315385.1"),
    ("essential_kpneumoniae_atcc43816_tradis", "kpneumoniae__kppr1__GCF_000742755.1"),
    ("essential_kpneumoniae_rh201207_tradis", "rh201207_bruchmann"),
]


def _auc(y: pd.Series, score: pd.Series) -> tuple[float, float]:
    from sklearn.metrics import average_precision_score, roc_auc_score
    return (round(float(roc_auc_score(y, score)), 4),
            round(float(average_precision_score(y, score)), 4))


def validate_anchor(species: str, scores: pd.DataFrame) -> list[dict]:
    """Score the shipped column against measured screens that key onto this anchor directly."""
    from src import essentiality as Es
    rows = []
    for name in ANCHOR_SCREENS.get(species, []):
        lab = Es.load_training_set(name)[["key", "label"]].dropna()
        m = scores.merge(lab, left_on="uniprot_ac", right_on="key", how="inner")
        if m["label"].nunique() < 2:
            continue
        auroc, aupr = _auc(m["label"], m["proteomelm_ess"])
        rows.append({"species": species, "role": ROLE[species], "scored": species,
                     "screen": name, "n": len(m), "n_pos": int(m["label"].sum()),
                     "base_rate": round(float(m["label"].mean()), 4),
                     "auroc": auroc, "aupr": aupr, "join": "direct (uniprot_ac)"})

    # DEG, where it exists on this anchor: measured, already keyed on uniprot_ac. Null means
    # unmeasured and is dropped -- never read as non-essential.
    try:
        deg = Es.load_deg(species)[["uniprot_ac", "deg_ess"]].dropna()
    except Exception:
        deg = pd.DataFrame()
    if len(deg):
        m = scores.merge(deg, on="uniprot_ac", how="inner")
        m = m[m["deg_ess"].isin([0.0, 1.0])]          # 0.5 is two screens disagreeing, not a label
        if len(m) and m["deg_ess"].nunique() == 2:
            auroc, aupr = _auc(m["deg_ess"], m["proteomelm_ess"])
            rows.append({"species": species, "role": ROLE[species], "scored": species,
                         "screen": "deg_ess (measured)", "n": len(m),
                         "n_pos": int(m["deg_ess"].sum()),
                         "base_rate": round(float(m["deg_ess"].mean()), 4),
                         "auroc": auroc, "aupr": aupr, "join": "direct (uniprot_ac)"})
    return rows


def validate_strains(device: str, quiet: bool) -> list[dict]:
    """Score each K. pneumoniae SCREEN STRAIN itself and evaluate on its own measured labels.

    This is the number that prices the Kp column: K. pneumoniae is absent from the head's training
    data entirely, and these TraDIS calls are labels the authors never saw. No identifier mapping
    is involved -- the strain's own proteome is scored and joined to its own screen on its own ids.
    """
    from src import essentiality as Es
    from src import strains as S
    rows = []
    for name, label in STRAIN_SCREENS:
        say(f"  strain {label}", quiet)
        frame = S.load_frame(label)
        idcol = next((c for c in ("protein_id", "id", "uniprot_ac", "key") if c in frame.columns),
                     frame.columns[0])
        seqcol = next((c for c in frame.columns if "seq" in c.lower()), None)
        if seqcol is None:
            say(f"    no sequence column in {label} -- skipped", quiet)
            continue
        frame = frame[[idcol, seqcol]].rename(columns={idcol: "uniprot_ac", seqcol: "sequence"})
        frame = frame[frame["sequence"].astype(str).str.strip().ne("")].drop_duplicates("uniprot_ac")

        fasta = SCRATCH / f"proteomelm_ess_strain_{label}.faa"
        raw = SCRATCH / f"proteomelm_ess_strain_{label}.tsv"
        man = SCRATCH / f"proteomelm_ess_strain_{label}.json"
        if not raw.exists():
            write_fasta(frame, fasta)
            say(f"    {len(frame):,} sequences -> scoring", quiet)
            score_fasta(fasta, raw, man, device, quiet)
        sc = pd.read_csv(raw, sep="\t").rename(
            columns={"protein_id": "uniprot_ac", "p_essential": "proteomelm_ess"})
        sc["uniprot_ac"] = sc["uniprot_ac"].astype(str)

        lab = Es.load_training_set(name)[["key", "label"]].dropna()
        lab["key"] = lab["key"].astype(str)
        m = sc.merge(lab, left_on="uniprot_ac", right_on="key", how="inner")
        if m.empty or m["label"].nunique() < 2:
            say(f"    joined {len(m)} rows -- cannot evaluate", quiet)
            continue
        auroc, aupr = _auc(m["label"], m["proteomelm_ess"])
        rows.append({"species": "kpneumoniae", "role": "unseen_species", "scored": label,
                     "screen": name, "n": len(m), "n_pos": int(m["label"].sum()),
                     "base_rate": round(float(m["label"].mean()), 4),
                     "auroc": auroc, "aupr": aupr,
                     "join": f"strain's own ids ({len(m)}/{len(lab)} of the screen)"})
        say(f"    AUROC {auroc}  AUPR {aupr}  on {len(m):,} joined rows", quiet)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", nargs="+", default=list(SPECIES), choices=list(SPECIES))
    ap.add_argument("--device", default="cpu", choices=["cpu", "mps", "cuda"])
    ap.add_argument("--validate", action="store_true",
                    help="evaluate the shipped column against measured screens on the anchors")
    ap.add_argument("--validate-strains", action="store_true",
                    help="also score the three K. pneumoniae screen strains (slow, and the only "
                         "honest out-of-distribution number this axis can get for Kp)")
    ap.add_argument("--refresh", action="store_true", help="re-score even if the TSV exists")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("-q", "--quiet", action="store_true")
    args = ap.parse_args()
    q = args.quiet

    rule("=")
    say("scripts/essentiality/proteomelm_ess.py -- ProteomeLM-Ess, the paper's own head", q)
    say(f"  head      {HEAD}  (Apache-2.0)", q)
    say(f"  env       gradi-plm-ess, across a process boundary ({WORKER.name})", q)
    say(f"  in        data/processed/proteomes/proteome_<species>.tsv  (sequences)", q)
    say(f"  out       {OUT_DIR.relative_to(REPO_ROOT)}/proteomelm_ess_<species>.tsv", q)
    say(f"  species   {', '.join(args.species)}   device {args.device}", q)
    rule("=")

    if args.dry_run:
        for sp in args.species:
            n = len(P.load(sp))
            say(f"  would score {sp:13s} {n:,} proteins  (role: {ROLE[sp]})", q)
        say("\ndry run -- nothing written.", q)
        return

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    SCRATCH.mkdir(parents=True, exist_ok=True)

    controls, manifests = [], []
    for sp in args.species:
        rule()
        say(f"  {sp}  ({ROLE[sp]})", q)
        final = OUT_DIR / f"proteomelm_ess_{sp}.tsv"
        raw = SCRATCH / f"proteomelm_ess_raw_{sp}.tsv"
        manifest_json = SCRATCH / f"proteomelm_ess_{sp}.json"
        fasta = SCRATCH / f"proteomelm_ess_{sp}.faa"

        prot = P.load(sp)[["uniprot_ac", "sequence"]]
        prot = prot[prot["sequence"].astype(str).str.strip().ne("")]
        if len(prot) != len(P.load(sp)):
            sys.exit(f"FAILED: {sp} has proteins with no sequence; the head needs all of them.")

        if raw.exists() and not args.refresh:
            say(f"    cached: {raw.relative_to(REPO_ROOT)}", q)
            scores = pd.read_csv(raw, sep="\t")
        else:
            n = write_fasta(prot, fasta)
            say(f"    {n:,} sequences -> {fasta.relative_to(REPO_ROOT)}", q)
            started = time.time()
            scores = score_fasta(fasta, raw, manifest_json, args.device, q)
            say(f"    scored in {time.time() - started:.0f}s", q)

        df = (scores.rename(columns={"protein_id": "uniprot_ac",
                                     "p_essential": "proteomelm_ess",
                                     "rank": "proteomelm_ess_rank"})
              [["uniprot_ac", "proteomelm_ess", "proteomelm_ess_rank"]])
        df["uniprot_ac"] = df["uniprot_ac"].astype(str)
        df["proteomelm_ess"] = df["proteomelm_ess"].round(6)
        df["proteomelm_ess_evidence"] = ROLE[sp]

        # COMPLETE and CANONICAL, like every matrix this project ships. `reindex` raises rather
        # than inventing a row, which is what we want: the head scores every protein, so a missing
        # one means the worker dropped it and that must not pass silently.
        df = M.reindex(df, sp)
        df.to_csv(final, sep="\t", index=False)
        say(f"    wrote {final.relative_to(REPO_ROOT)}  ({len(df):,} rows, canonical)", q)

        c = control(df, sp)
        controls.append(c)
        say(f"    control: ribosome median {c['ribosome_median']} (n={c['n_ribosomal']})  vs  "
            f"dispensable {c['dispensable_median']} (n={c['n_dispensable']})  "
            f"AUROC {c['ribosome_vs_dispensable_auroc']}", q)

        if manifest_json.exists():
            m = json.loads(manifest_json.read_text())
            m["species"] = sp
            m["role"] = ROLE[sp]
            manifests.append(m)

    # ---- validation against MEASURED labels. The control above only says the column is the right
    # way up; this says whether it is any good, and it is the only place the Kp number can come
    # from -- K. pneumoniae has no measured essentiality on HS11286 at all.
    val: list[dict] = []
    if args.validate or args.validate_strains:
        rule()
        say("  VALIDATION against measured screens", q)
        scored = {sp: pd.read_csv(OUT_DIR / f"proteomelm_ess_{sp}.tsv", sep="\t")
                  for sp in args.species}
        for sp in args.species:
            val += validate_anchor(sp, scored[sp])
        if args.validate_strains:
            val += validate_strains(args.device, q)
        if val:
            vdf = pd.DataFrame(val)
            vdf.to_csv(EVIDENCE / "proteomelm_ess_validation.tsv", sep="\t", index=False)
            say(vdf.to_string(index=False), q)
            say("\n  Read `role` with every row: `in_training` is recall, not performance.", q)

    cdf = pd.DataFrame(controls)
    cdf.to_csv(EVIDENCE / "proteomelm_ess_control.tsv", sep="\t", index=False)
    if manifests:
        pd.DataFrame(manifests).to_csv(EVIDENCE / "proteomelm_ess_manifest.tsv", sep="\t", index=False)

    rule("=")
    say("  CONTROL -- ribosome must outrank the textbook dispensables", q)
    say(cdf.to_string(index=False), q)
    rule("=")

    failed = [c for c in controls
              if c["ribosome_vs_dispensable_auroc"] != ""
              and c["ribosome_vs_dispensable_auroc"] < MIN_RIBOSOME_AUROC]
    if failed:
        sys.exit(f"FAILED: {len(failed)} species below the {MIN_RIBOSOME_AUROC} polarity floor "
                 f"({', '.join(c['species'] for c in failed)}). The column may be INVERTED -- "
                 "class 0 is essential in the head's output.")
    say("  polarity OK on every species.", q)


if __name__ == "__main__":
    main()
