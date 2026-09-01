# `legacy/` — the v1 pipeline, frozen

This directory holds the **first** version of the GraDi target-prioritization pipeline, archived in
September 2026 when the work restarted at the repository root.

**Start with [`HISTORY.md`](HISTORY.md).** It is the retrospective: what was built and what was not,
the methodological decisions and why, ~70 documented traps, which data sources were obtained and how,
and — importantly — which artifacts here look like data but are not.

## Rules

- **Frozen.** Do not extend anything in here. New work goes at the repository root.
- **Do not trust three things** without reading `HISTORY.md` §7 first: the Kp legacy degradability TSV
  (four verified defects), the E. coli degradability values in the webapp (a deterministic MD5 mock),
  and `data/raw/legacy/clp_substrates/` (45- and 35-row hand-curated substitutes, not the papers'
  tables).
- The webapp under `app/` is still **deployed** from here by `.github/workflows/pages.yml`. It serves
  v1 numbers, including the two untrustworthy degradability columns above.

## Layout

```
HISTORY.md   the retrospective — read this first
scripts/     the numbered v1 pipeline, 00a .. 10o (81 files)
src/         the four shared modules (ligandability, essentiality, localization, degradability)
docs/        the five axis specs plus run logs, reports and references (24 files)
app/         the target-selector webapp (still live)
assets/      the five axis diagrams
data -> ../data      symlink
output -> ../output  symlink
```

## Running something in here

The scripts still work. They resolve paths with `Path(__file__).resolve().parents[1]`, which after the
move points at `legacy/` — so `legacy/data` and `legacy/output` are **symlinks to the real trees one
level up**. That is what keeps these scripts runnable without editing 81 files, and it is why the
archive is a faithful snapshot of what actually ran.

Without those symlinks a script would resolve its paths under `legacy/`, find no inputs, and **exit 0
having written nothing** — a silent failure rather than an error. If a legacy script ever appears to
do nothing, check the symlinks first.

Run them with the `gradi` conda env, from inside this directory:

```bash
cd legacy
python scripts/09h_localization_plots.py --organism ecoli
```

Four scripts need a different env — `03a`/`03c` need `gradi-ortho`, `06e` needs `gradi-pockets`,
`06n` needs `gradi-pymol`, and `09c`/`09d` need `gradi-loc`. See `HISTORY.md` §8 for why each exists;
the `gradi-loc` split in particular is mandatory, not cosmetic.
