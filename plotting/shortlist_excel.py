"""The E. coli x ASKA shortlist as a colour-coded Excel workbook: the function-balanced list of 200 and its random control.

Reads the tables `plotting/shortlist_ecoli.py` writes (nothing is re-derived here, so the workbook cannot
disagree with them):

    shortlist_ecoli_aska_balanced200.tsv         Sheet "Shortlist"      proteins with slots proportional to the cytosolic ASKA
                                                                        proteome by COG function AND length band
    shortlist_ecoli_aska_balanced_control50.tsv  Sheet "Random control" random cytosolic ASKA proteins below the median
                                                                        degradability, not on the list, same function mix
                                                                        and essential share per function
                                                                        Sheet "Function and length mix"    proteome vs list vs control, per COG letter
                                                                        Sheet "Method"          how it was built, caveats, column guide

NO composite score. The order within each sheet is `degradability_consensus`, descending, a stated tiebreak.
Selectivity, ligandability, pockets and studiedness are decorators, shown and never used to select or rank.

Writes   output/results/shortlist/shortlist_ecoli_aska.xlsx

Run with the `gradi` env:  python plotting/shortlist_excel.py [-q]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import pandas as pd  # noqa: E402
from openpyxl import Workbook  # noqa: E402
from openpyxl.formatting.rule import CellIsRule, ColorScaleRule  # noqa: E402
from openpyxl.styles import Alignment, Border, Font, PatternFill, Side  # noqa: E402
from openpyxl.utils import get_column_letter  # noqa: E402

from plotting import palette as PAL  # noqa: E402

OUT = REPO_ROOT / "output" / "results" / "shortlist"
SPECIES = "ecoli"
SOURCE_LABEL = {"essential_rule": "Meets essential rule", "degradability_only": "Degradability only",
                "random_control_below_median": "Random control"}
ROUTE_LABEL = {"measured_and_consensus": "Measured + consensus", "measured_only": "Measured only",
               "consensus_only": "Consensus only", "neither": "Neither"}


def hexa(rgb) -> str:
    """palette.py roles are RGB float tuples; openpyxl wants RRGGBB."""
    return "".join(f"{round(c * 255):02X}" for c in rgb[:3])


def describe_function(cog_categories: pd.Series) -> pd.Series:
    """COG category names, `; `-joined, from the `;`-joined letters in function_<species>.tsv."""
    cog = pd.read_csv(REPO_ROOT / "data/source/cdd/cog_func_category.tsv", sep="\t", header=None,
                      names=["letter", "group", "colour", "name"]).set_index("letter")["name"]
    return cog_categories.fillna("").map(lambda s: "; ".join(cog.get(x, x) for x in str(s).split(";") if x))


def load_sheet(name: str) -> pd.DataFrame:
    from src import function

    d = pd.read_csv(OUT / name, sep="\t")
    fn = function.load(SPECIES)[["uniprot_ac", "cog_categories"]]
    d = d.merge(fn, on="uniprot_ac", how="left", validate="1:1")
    d["function"] = describe_function(d["cog_categories"])
    d["essential_any"] = d["deg_essential_any"] == 1
    d["essential_all"] = d["deg_essential_all"] == 1
    d["human_ortholog"] = d["has_human_ortholog"].astype(bool)
    d["recorded_drug"] = d["n_drugs"].fillna(0) > 0
    d["list_source"] = d["list_source"].map(SOURCE_LABEL)
    if "essential_route" in d:
        d["essential_route"] = d["essential_route"].map(ROUTE_LABEL)
        d["below_gate"] = ~d["passes_degradability_gate"].astype(bool)
    d["cog_letter"] = d["cog"].replace("none", "-")
    d = d.sort_values("degradability_consensus", ascending=False).reset_index(drop=True)
    d.insert(0, "rank", range(1, len(d) + 1))
    return d


COLUMNS = [  # (column, header, width, kind)
    ("rank", "#", 5, "int"),
    ("uniprot_ac", "UniProt", 11, "text"),
    ("gene_name", "Gene", 10, "text"),
    ("protein_name", "Protein", 42, "text"),
    ("function", "Function (COG category)", 40, "text"),
    ("cog_letter", "COG letter (- = unclassified)", 9, "text"),
    ("length_band", "Length band (aa)", 10, "text"),
    ("deg_pct_in_len", "Degradability percentile within its length decile", 14, "dec"),
    ("list_source", "Why on the list", 20, "text"),
    ("essential_route", "Essential: how (measured / consensus >= 0.7)", 18, "text"),
    ("top50", "In the absolute top 50", 9, "bool"),
    ("below_gate", "Below degradability gate", 10, "bool"),
    ("degradability_consensus", "Degradability consensus", 14, "scale_deg"),
    ("adep4_prob", "ADEP4 prob", 10, "dec"),
    ("onc212_prob", "ONC212 prob", 10, "dec"),
    ("deg_route", "Degradability route", 14, "text"),
    ("ortholog_support", "S. aureus ortholog", 14, "text"),
    ("essential_any", "Measured essential (any screen)", 12, "bool"),
    ("essential_all", "Measured essential (all screens)", 12, "bool"),
    ("deg_n_datasets", "DEG screens covering it", 10, "int"),
    ("essentiality_consensus", "Essentiality consensus", 14, "scale_ess"),
    ("essentiality_evidence", "Ess. evidence", 10, "ev"),
    ("human_ortholog", "Human ortholog", 10, "bool"),
    ("bacterial_panel_orthologs", "Conserved in bacterial panel (fraction)", 12, "dec"),
    ("ligands_consensus", "Ligands consensus", 12, "scale_lig"),
    ("ligands_evidence", "Lig. evidence", 10, "ev"),
    ("n_ligands_bacterial", "Potent ligands (bacterial)", 12, "int"),
    ("ligand_evidence", "Ligand evidence (potent ChEMBL, PDB ligand or drug)", 14, "bool"),
    ("recorded_drug", "Recorded drug", 10, "bool"),
    ("n_drugs", "Drugs (n)", 8, "int"),
    ("n_ligands_pdb", "PDB ligands", 9, "int"),
    ("n_pdb_structures", "PDB structures", 9, "int"),
    ("pockets_evidence", "Pockets evidence", 10, "ev"),
    ("p2rank_score", "P2Rank score (length-driven)", 12, "dec"),
    ("length", "Length (aa)", 9, "int"),
    ("ribosomal", "Ribosomal", 9, "bool"),
    ("on_panel", "On consortium panel", 10, "bool"),
    ("n_papers_uniprot_prokaryotic", "Papers", 8, "int"),
    ("aska_id", "ASKA clone", 12, "text"),
    ("aska_confidence", "ASKA mapping confidence (3 best)", 10, "int"),
]


def write_sheet(ws, s: pd.DataFrame, hdr: PatternFill, border: Border, columns=COLUMNS) -> None:
    cols = [c for c in columns if c[0] in s.columns]
    for j, (_, head, w, _) in enumerate(cols, 1):
        c = ws.cell(row=1, column=j, value=head)
        c.font, c.fill = Font(bold=True, color="FFFFFF"), hdr
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[1].height = 58
    for i, row in enumerate(s.itertuples(index=False), 2):
        r = row._asdict()
        for j, (col, _, _, kind) in enumerate(cols, 1):
            v = r[col]
            if kind == "bool":
                v = "yes" if bool(v) else ""
            elif v is None or (not isinstance(v, str) and pd.isna(v)):
                v = None
            elif kind in ("int", "ev"):
                v = int(v)
            elif kind == "dec" or kind.startswith("scale"):
                v = round(float(v), 4)
            c = ws.cell(row=i, column=j, value=v)
            c.border = border
            if kind != "text":
                c.alignment = Alignment(horizontal="center")
    last = len(s) + 1
    idx = {c: get_column_letter(j) for j, (c, *_) in enumerate(cols, 1)}

    def scale(col, hi):
        ws.conditional_formatting.add(
            f"{idx[col]}2:{idx[col]}{last}",
            ColorScaleRule(start_type="num", start_value=0, start_color="FFFFFF",
                           end_type="num", end_value=1, end_color=hi))

    scale("essentiality_consensus", hexa(PAL.SECONDARY))
    scale("degradability_consensus", hexa(PAL.PRIMARY))
    scale("ligands_consensus", hexa(PAL.TERTIARY))
    scale("adep4_prob", "9DB1D9")
    scale("onc212_prob", "9DB1D9")
    for col in ("essentiality_evidence", "ligands_evidence", "pockets_evidence"):
        for val, colour in ((1, "F4B6A9"), (2, "FDE9A9"), (3, "A9DDB8")):
            ws.conditional_formatting.add(
                f"{idx[col]}2:{idx[col]}{last}",
                CellIsRule(operator="equal", formula=[str(val)], fill=PatternFill("solid", bgColor=colour)))
    for col, colour in (("essential_any", "A9DDB8"), ("essential_all", "A9DDB8"), ("recorded_drug", "A9DDB8"),
                        ("human_ortholog", "F4B6A9"), ("top50", "A9DDB8"), ("below_gate", "F4B6A9"), ("ligand_evidence", "A9DDB8")):
        if col not in idx:
            continue
        ws.conditional_formatting.add(
            f"{idx[col]}2:{idx[col]}{last}",
            CellIsRule(operator="equal", formula=['"yes"'], fill=PatternFill("solid", bgColor=colour)))
    ws.freeze_panes = "E2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{last}"


def write_mix(ws, mix: pd.DataFrame, lmix: pd.DataFrame, hdr: PatternFill, border: Border) -> None:
    heads = ["COG letter (- = unclassified)", "Function", "Cytosolic ASKA proteome (n)", "Proteome %", "List (n)", "List %",
             "List minus proteome (points)", "Control (n)"]
    widths = [14, 38, 16, 11, 9, 9, 16, 11]
    for j, (h, w) in enumerate(zip(heads, widths), 1):
        c = ws.cell(row=1, column=j, value=h)
        c.font, c.fill = Font(bold=True, color="FFFFFF"), hdr
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[1].height = 45
    for i, r in enumerate(mix.itertuples(index=False), 2):
        vals = [r.letter, r.name, int(r.proteome_n), round(r.proteome_pct, 1), int(r.list_n), round(r.list_pct, 1),
                round(r.list_pct - r.proteome_pct, 1), int(r.control_n)]
        for j, v in enumerate(vals, 1):
            c = ws.cell(row=i, column=j, value=v)
            c.border = border
            if j > 2:
                c.alignment = Alignment(horizontal="center")
    r0 = len(mix) + 4
    for j, h in enumerate(["Length band (aa)", "", "Cytosolic ASKA proteome (n)", "Proteome %", "List (n)", "List %", "List minus proteome (points)", "Control (n)",
                           "Target % (mean of 3 distributions)", "List minus target (points)"], 1):
        c = ws.cell(row=r0, column=j, value=h)
        c.font, c.fill = Font(bold=True, color="FFFFFF"), hdr
        c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
    for i, r in enumerate(lmix.itertuples(index=False), r0 + 1):
        vals = [r.letter, "", int(r.proteome_n), round(r.proteome_pct, 1), int(r.list_n), round(r.list_pct, 1),
                round(r.list_pct - r.proteome_pct, 1), int(r.control_n), round(r.target_pct, 1), round(r.list_pct - r.target_pct, 1)]
        for j, v in enumerate(vals, 1):
            c = ws.cell(row=i, column=j, value=v)
            c.border = border
            if j > 2:
                c.alignment = Alignment(horizontal="center")
    for j, w in ((9, 20), (10, 16)):
        ws.column_dimensions[get_column_letter(j)].width = w
    ws.freeze_panes = "A2"


def write(sheets: dict[str, pd.DataFrame], mix: pd.DataFrame, lmix: pd.DataFrame, params: dict, path: Path) -> None:
    hdr = PatternFill("solid", fgColor=hexa(PAL.PRIMARY))
    thin = Side(style="thin", color="D9D9D9")
    border = Border(top=thin, bottom=thin, left=thin, right=thin)
    wb = Workbook()
    first = True
    for name, s in sheets.items():
        ws = wb.active if first else wb.create_sheet()
        first = False
        ws.title = name
        write_sheet(ws, s, hdr, border)
    write_mix(wb.create_sheet("Function and length mix"), mix, lmix, hdr, border)

    p = params
    m = wb.create_sheet("Method")
    m.column_dimensions["A"].width = 40
    m.column_dimensions["B"].width = 120
    rows = [
        ("ORGANISM", "E. coli K-12 MG1655, restricted to proteins with a clone in the ASKA (-) library and predicted cytosolic"),
        ("", ""),
        ("HOW THE LIST WAS BUILT", "Not a score: stated steps. Slots follow the cytosolic ASKA proteome on two margins, COG function and length band."),
        ("1. Quotas", f"{p['n']} slots, matched as closely as possible on COG function (first letter) to the {p['pool']:,} cytosolic ASKA proteins AND on "
                      f"length band ({', '.join(p['bands'])} aa) to a TARGET mix: the equal-weight average of three distributions ({p['target_srcs']}). No single "
                      "distribution is right: degradation references (S. aureus hits, essential proteins) run shorter and ligandable proteins longer; this sits "
                      "in between. Each step adds a protein to the function and length band furthest below their targets; exact function x length cells are not "
                      "targeted (too small for 200 proteins)."),
        ("2. Eligibility", f"A protein is eligible if its degradability consensus is in the top {p['elig_pct']:.0f}% of its own length decile ({p['eligible']:,} proteins). "
                           "The consensus is strongly length-driven, so a raw cut would only admit short proteins; this asks whether a protein is unusually "
                           f"degradable for its size. The {p['elig_pct']:.0f}% is a CHOSEN value. The raw degradability gate is NOT applied: {p['below_gate']} of the "
                           "200 fall below it and are flagged."),
        ("3. Order inside a cell", f"(a) meets the essential rule: measured essential in any DEG/Keio screen OR essentiality consensus >= {p['or_cut']} (a CHOSEN cut; the "
                                   f"derived base-rate cut is {p['ess_cut']:.4f}); (b) has ligand evidence; (c) degradability percentile within its length decile. "
                                   "Ribosomal proteins never get the ligand preference: a recorded drug on the ribosome belongs to the whole complex."),
        ("Result", f"{p['n']} proteins; {p['ess_rule']} meet the essential rule, {p['ess_any']} are measured essential in any screen, {p['ess_strict']} under the strict rule "
                   f"(measured AND consensus >= {p['ess_cut']:.4f}); {p['ligand']} have ligand evidence; {p['top50']} of the absolute top 50 by raw degradability are "
                   f"in it; {p['ribo']} are ribosomal; median length {p['list_len']:.0f} aa against {p['pool_len']:.0f} for the cytosolic proteome."),
        ("Balance", "Sheet 'Function and length mix' compares the list with the cytosolic proteome and the length targets. Largest function deviation: "
                    f"{p['max_dev']:.1f} points; length bands: {p['band_txt']}."),
        ("", ""),
        ("RANDOM CONTROL", f"{p['n_control']} random cytosolic ASKA proteins, not on the list, with raw degradability consensus BELOW THE PROTEOME MEDIAN ({p['floor']:.4f}), drawn from "
                           f"{p['control_pool']:,}. It follows the list's function mix and measured-essential share ({p['control_ess']} of {p['n_control']}; the pool's own share is "
                           f"lower). Seed {p['seed']}. Highest raw consensus {p['control_raw_max']:.2f}, median {p['control_raw_med']:.2f}. LIMITATION: no protein up to 100 aa and only one of 286 at "
                           f"100-150 aa scores below the median, so the short bands cannot be filled. Length bands {p['ctrl_band_txt']}; median length {p['control_len']:.0f} aa against "
                           f"{p['list_len']:.0f} for the list. Compare with the list within a length band, and note that the list prefers ligand evidence by design, so ligand "
                           "evidence is not matched."),
        ("ORDER", "Each sheet is ordered by degradability consensus, descending. A tiebreak, not a score."),
        ("", ""),
        ("COLUMN GUIDE", ""),
        ("COG letter", "First COG category letter of the protein (9-16% carry several); - means unclassified."),
        ("Degradable gate", f"Degradability consensus >= {p['cut']:.4f} (derived: the quantile where it passes the same share as either activator above its "
                            f"base-rate cut), OR a measured S. aureus ortholog is a degrader AND the consensus >= {p['floor']:.4f} (the proteome median)."),
        ("Degradability consensus", "Mean within-species percentile of ADEP4 and ONC212 probabilities, trained on S. aureus screens and EXTRAPOLATED to "
                                    "E. coli: ranking hypotheses, not measurements. The two activators are close to one opinion (rho ~0.97). "
                                    "Strongly length-driven (small proteins score high)."),
        ("Essential: how", "Measured + consensus / Measured only / Consensus only / Neither (top-ups). Consensus cut 0.7, chosen. Essentiality is a FLAG here, "
                           "not a gate; 'Neither' proteins were selected for degradability alone."),
        ("Measured essential (any / all)", "Any = at least one DEG screen calls it essential; all = every screen covering it does. Rich-medium screens."),
        ("Essentiality consensus", "Mean within-species percentile of Geptop, ProteomeLM-Ess and the screen-trained models. On E. coli it is built from the same "
                                   "screens as the measured calls, so agreement is not independent confirmation. Compare within a species only."),
        ("Human ortholog", "Union of OrthoFinder and reciprocal best hit. A decorator, never a gate: yes means a selectivity liability."),
        ("Ligand evidence", "Yes = a potent ChEMBL ligand (this protein or a bacterial homolog), a drug-like ligand on its own PDB structure, or a recorded drug. "
                            "Used only to order the top-ups inside a function, never as a gate or a score; P2Rank is not used because it mostly ranks protein length."),
        ("Ligands / pockets", "Decorators only. Ligands consensus is ChEMBL potent-ligand evidence (0 with evidence 1 = nobody screened it). Recorded drug = "
                              "a drug in any ChEMBL mechanism category against a >=95% identical target (target-level, mostly ribosome/PBP/gyrase). "
                              "P2Rank score is largely a ranking by protein length."),
        ("Evidence (1/2/3)", "Per-axis provenance ladder: red 1 = prediction only / could not look, amber 2, green 3 = best corroborated. Not comparable across axes."),
        ("Degradability route", "consensus_only / both / ortholog_only: which part of the degradable rule the protein passes. S. aureus ortholog: "
                                "degrader / not a degrader / unmeasured / no ortholog (shown, never a veto)."),
    ]
    for i, (a, b) in enumerate(rows, 1):
        ca, cb = m.cell(row=i, column=1, value=a), m.cell(row=i, column=2, value=b)
        if a.isupper() or a == "ORGANISM":
            ca.font = Font(bold=True, color="FFFFFF")
            ca.fill = cb.fill = hdr
            cb.font = Font(bold=True, color="FFFFFF")
        cb.alignment = Alignment(wrap_text=True, vertical="top")
        ca.alignment = Alignment(wrap_text=True, vertical="top")
    wb.save(path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("-q", "--quiet", action="store_true")
    ap.add_argument("--n", type=int, default=200)
    ap.add_argument("--control", type=int, default=50)
    a = ap.parse_args()
    say = (lambda *_: None) if a.quiet else print

    say("shortlist_excel.py  E. coli x ASKA: function-balanced list, random control")
    lst = load_sheet(f"shortlist_ecoli_aska_balanced{a.n}.tsv")
    ctrl = load_sheet(f"shortlist_ecoli_aska_balanced_control{a.control}.tsv")
    assert lst.uniprot_ac.is_unique and not set(lst.uniprot_ac) & set(ctrl.uniprot_ac), "list and control overlap"

    # The cuts and the proteome mix come from the same loader the shortlist uses, so the workbook states what ran.
    from plotting import shortlist_ecoli as S

    d = S.load()
    pool = d[d.in_aska & d.cytosolic]
    cog_names = pd.read_csv(REPO_ROOT / "data/source/cdd/cog_func_category.tsv", sep="\t", header=None,
                            names=["letter", "group", "colour", "name"]).set_index("letter")["name"]
    letters = pool.cog.value_counts().index
    mix = pd.DataFrame({"letter": ["-" if c == "none" else c for c in letters],
                        "name": ["Unclassified" if c == "none" else cog_names.get(c, c) for c in letters],
                        "proteome_n": pool.cog.value_counts().reindex(letters).values,
                        "list_n": lst.cog.value_counts().reindex(letters).fillna(0).values,
                        "control_n": ctrl.cog.value_counts().reindex(letters).fillna(0).values})
    mix["proteome_pct"] = mix.proteome_n / mix.proteome_n.sum() * 100
    mix["list_pct"] = mix.list_n / mix.list_n.sum() * 100
    for df in (mix,):
        df["dev"] = (df.list_pct - df.proteome_pct)
    pool_b = pd.cut(pool.length, S.LENGTH_EDGES, labels=S.LENGTH_LABELS)
    lmix = pd.DataFrame({"letter": S.LENGTH_LABELS, "proteome_n": pool_b.value_counts().reindex(S.LENGTH_LABELS).values,
                         "list_n": lst.length_band.value_counts().reindex(S.LENGTH_LABELS).fillna(0).values,
                         "control_n": ctrl.length_band.value_counts().reindex(S.LENGTH_LABELS).fillna(0).values})
    lmix["proteome_pct"] = lmix.proteome_n / lmix.proteome_n.sum() * 100
    lmix["list_pct"] = lmix.list_n / lmix.list_n.sum() * 100
    targets = pd.read_csv(OUT / "shortlist_ecoli_aska_length_targets.tsv", sep="\t", index_col=0).reindex(S.LENGTH_LABELS)
    lmix["target_pct"] = targets.target.values * 100
    assert mix.dev.abs().max() < 2 and (lmix.list_pct - lmix.target_pct).abs().max() < 2, "list is more than 2 points from its target mixes"

    params = {"cut": d.attrs["cons_deg_cut"], "floor": d.attrs["orth_floor"], "ess_cut": d.attrs["cons_cut"], "or_cut": 0.7,
              "n": len(lst), "pool": len(pool), "seed": 0, "bands": S.LENGTH_LABELS, "elig_pct": 30, "eligible": int((pool.length > 0).sum() * 0 + (
                  pool.groupby(pd.qcut(pool.length, 10, labels=False)).degradability_consensus.rank(pct=True) >= 0.7).sum()),
              "below_gate": int(lst.below_gate.sum()), "ess_rule": int((lst.list_source == "Meets essential rule").sum()),
              "ess_any": int(lst.essential_any.sum()), "ess_strict": int(lst.essential_strict.sum()), "ligand": int(lst.ligand_evidence.sum()),
              "top50": int(lst.top50.sum()), "ribo": int(lst.ribosomal.sum()), "max_dev": float(mix.dev.abs().max()),
              "band_txt": ", ".join(f"{r.letter} {r.list_pct:.0f}% (target {r.target_pct:.0f}%, proteome {r.proteome_pct:.0f}%)" for r in lmix.itertuples()),
              "target_srcs": ", ".join(f"{k.replace('_', ' ')}" for k in targets.columns if k != "target"),
              "ctrl_band_txt": ", ".join(f"{r.letter} {int(r.control_n)}" for r in lmix.itertuples()),
              "n_control": len(ctrl), "control_ess": int(ctrl.essential_any.sum()),
              "control_pool": int(((pool.degradability_consensus < d.attrs["orth_floor"]) & ~pool.uniprot_ac.isin(lst.uniprot_ac)).sum()),
              "control_len": float(ctrl.length.median()),  "control_raw_med": float(ctrl.degradability_consensus.median()), "control_raw_max": float(ctrl.degradability_consensus.max()),
              "list_len": float(lst.length.median()), "pool_len": float(pool.length.median())}

    path = OUT / "shortlist_ecoli_aska.xlsx"
    write({f"Shortlist ({len(lst)})": lst, f"Random control ({len(ctrl)})": ctrl}, mix, lmix, params, path)
    say(f"\n  wrote {path.relative_to(REPO_ROOT)}")
    say(f"  list {len(lst)}: {lst.list_source.value_counts().to_dict()} | measured essential {params['ess_any']} | strict {params['ess_strict']} | "
        f"top 50 inside {params['top50']} | below gate {params['below_gate']} | ribosomal {params['ribo']} | human ortholog {int(lst.human_ortholog.sum())} | "
        f"ligand evidence {params['ligand']} | median length {params['list_len']:.0f}")
    say(f"  control {len(ctrl)}: measured essential {params['control_ess']} | length bands {params['ctrl_band_txt']} | median length {params['control_len']:.0f} | "
        f"raw consensus median {params['control_raw_med']:.2f}, max {params['control_raw_max']:.2f}")
    say(f"  largest function deviation {params['max_dev']:.1f} points; length bands {params['band_txt']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
