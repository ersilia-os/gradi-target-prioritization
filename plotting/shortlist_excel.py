"""The E. coli x ASKA shortlist as a colour-coded Excel workbook: the list of 200, its random control, the overflow.

Reads the tables `plotting/shortlist_ecoli.py` writes (nothing is re-derived here, so the workbook cannot
disagree with them):

    shortlist_ecoli_aska_top200.tsv           Sheet "Shortlist"       the absolute top 50 by degradability among the
                                                                      proteins passing ASKA, cytosolic and degradable,
                                                                      topped up with the essential-gate survivors
    shortlist_ecoli_aska_control50.tsv        Sheet "Random control"  random cytosolic ASKA proteins below the median
                                                                      degradability, not on the list
    shortlist_ecoli_aska_top200_overflow.tsv  Sheet "Left out"        essential-gate proteins that did not fit
    shortlist_ecoli_aska_gates.tsv            Sheet "Method"          funnel and leave-one-out

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
SOURCE_LABEL = {"top_degradability": "Top 50 degradability", "essential_gate": "Essential gate",
                "both": "Both", "random_control_below_median": "Random control"}


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
    if "list_source" in d:
        d["list_source"] = d["list_source"].map(SOURCE_LABEL)
    d = d.sort_values("degradability_consensus", ascending=False).reset_index(drop=True)
    d.insert(0, "rank", range(1, len(d) + 1))
    return d


COLUMNS = [  # (column, header, width, kind)
    ("rank", "#", 5, "int"),
    ("uniprot_ac", "UniProt", 11, "text"),
    ("gene_name", "Gene", 10, "text"),
    ("protein_name", "Protein", 42, "text"),
    ("function", "Function (COG category)", 40, "text"),
    ("list_source", "Why on the list", 20, "text"),
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
                        ("human_ortholog", "F4B6A9")):
        ws.conditional_formatting.add(
            f"{idx[col]}2:{idx[col]}{last}",
            CellIsRule(operator="equal", formula=['"yes"'], fill=PatternFill("solid", bgColor=colour)))
    ws.freeze_panes = "E2"
    ws.auto_filter.ref = f"A1:{get_column_letter(len(cols))}{last}"


def write(sheets: dict[str, pd.DataFrame], gates: pd.DataFrame, params: dict, path: Path) -> None:
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

    m = wb.create_sheet("Method")
    m.column_dimensions["A"].width = 40
    m.column_dimensions["B"].width = 120
    rows = [
        ("ORGANISM", "E. coli K-12 MG1655, restricted to proteins with a clone in the ASKA (-) library"),
        ("", ""),
        ("GATES (in order)", "proteins remaining"),
        *[(f"{r.gate}: {r.rule}", f"{r.n_after:,}  (removes {r.n_removed:,}; survivors with this gate dropped: {int(r.n_without_gate):,})")
          for r in gates.itertuples()],
        ("", ""),
        ("THE LIST OF 200", "Not a score. The absolute top 50 by degradability consensus among the proteins that pass the first three gates "
                            "(ASKA, cytosolic, degradable), topped up with the proteins that also pass the essential gate, in degradability "
                            "order, to 200. 'Why on the list' says which part: Top 50 degradability / Essential gate / Both. "
                            f"{params['n_left_out']} essential-gate proteins with the lowest degradability did not fit (sheet 'Left out')."),
        ("RANDOM CONTROL", f"{params['n_control']} random cytosolic ASKA proteins below the proteome-median degradability consensus "
                           f"({params['floor']:.4f}) that are not on the list; {params['control_ess']} measured essential, "
                           f"{params['n_control'] - params['control_ess']} not (the list's essential share, 80%; the pool's own is 12%). Seed 0. "
                           "It is much longer than the list (median 382 vs 156 aa), so a comparison mixes degradability with protein size."),
        ("ORDER", "Each sheet is ordered by degradability consensus, descending. A tiebreak, not a score."),
        ("", ""),
        ("COLUMN GUIDE", ""),
        ("Cytosolic", "DeepLocPro class cytoplasm OR TMbed cytoplasmic fraction = 1 (two independent predictors; DeepLocPro mislabels many ribosomal proteins)."),
        ("Degradable", f"Degradability consensus >= {params['cut']:.4f} (derived: the quantile where it passes the same share as either activator above its "
                       f"base-rate cut), OR a measured S. aureus ortholog is a degrader AND the consensus >= {params['floor']:.4f} (the proteome median)."),
        ("Degradability consensus", "Mean within-species percentile of ADEP4 and ONC212 probabilities, trained on S. aureus screens and EXTRAPOLATED to "
                                    "E. coli: ranking hypotheses, not measurements. The two activators are close to one opinion (rho ~0.97). "
                                    "Strongly length-driven (small proteins score high)."),
        ("Essential gate", f"Measured essential in any DEG/Keio screen AND essentiality consensus >= {params['ess_cut']:.4f} (the quantile at 1 minus the "
                           "measured base rate, top 15.8%). On E. coli the consensus is built from the same screens as the measured calls, so this reads "
                           "'the screens agree and nothing contradicts them', not independent confirmation. Only ~55% of measured essentials clear the cut."),
        ("Measured essential (any / all)", "Any = at least one DEG screen calls it essential; all = every screen covering it does. Rich-medium screens."),
        ("Essentiality consensus", "Mean within-species percentile of Geptop, ProteomeLM-Ess and the screen-trained models. Compare within a species only."),
        ("Human ortholog", "Union of OrthoFinder and reciprocal best hit. A decorator, never a gate: yes means a selectivity liability."),
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
    a = ap.parse_args()
    say = (lambda *_: None) if a.quiet else print

    say("shortlist_excel.py  E. coli x ASKA: list of 200, random control, left out")
    lst = load_sheet("shortlist_ecoli_aska_top200.tsv")
    ctrl = load_sheet("shortlist_ecoli_aska_control50.tsv")
    over = load_sheet("shortlist_ecoli_aska_top200_overflow.tsv")
    gates = pd.read_csv(OUT / "shortlist_ecoli_aska_gates.tsv", sep="\t")
    assert lst.uniprot_ac.is_unique and not set(lst.uniprot_ac) & set(ctrl.uniprot_ac), "list and control overlap"

    # The cuts are re-derived by the same loader the shortlist uses, so the Method sheet states what ran.
    import importlib.util

    spec = importlib.util.spec_from_file_location("ecoli_aska_analysis", REPO_ROOT / "plotting/ecoli_aska_analysis.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    d = mod.load_ecoli()
    params = {"cut": d.attrs["cons_deg_cut"], "floor": d.attrs["orth_floor"], "ess_cut": d.attrs["cons_cut"],
              "n_left_out": len(over), "n_control": len(ctrl), "control_ess": int(ctrl.essential_any.sum())}

    path = OUT / "shortlist_ecoli_aska.xlsx"
    write({f"Shortlist ({len(lst)})": lst, f"Random control ({len(ctrl)})": ctrl, f"Left out ({len(over)})": over},
          gates, params, path)
    say(f"\n  wrote {path.relative_to(REPO_ROOT)}")
    say(f"  list {len(lst)}: {lst.list_source.value_counts().to_dict()} | measured essential {int(lst.essential_any.sum())} | "
        f"human ortholog {int(lst.human_ortholog.sum())} | recorded drug {int(lst.recorded_drug.sum())}")
    say(f"  control {len(ctrl)}: measured essential {int(ctrl.essential_any.sum())} | median consensus {ctrl.degradability_consensus.median():.3f}")
    say(f"  left out {len(over)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
