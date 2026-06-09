"""2x2 descriptor figure for the project's three reference proteomes.

Everything here is derivable from the *raw* proteome download (sequence + FASTA-header metadata),
before any annotation step — the proteome QC/overview companion to scripts/00_fetch_proteomes.py.
Panels (across E. coli K-12 / K. pneumoniae HS11286 / human reviewed-canonical):

  top-left     reviewed vs unreviewed   (>sp| / >tr| counts; Kp ~entirely unreviewed)
  top-right    protein-existence level  (normalised PE1..PE5; Kp ~99% predicted/inferred)
  bottom-left  gene-symbol coverage     (% with a real gene symbol; Kp ~18%)
  bottom-right protein length           (violin, log10 y; bacteria compact vs human long-tailed)

Panels TL/TR/BR come from the local FASTAs; gene-symbol coverage needs `gene_primary` (the FASTA
`GN=` field conflates real symbols with `KPHS_*` locus tags), fetched live from UniProt.

Styling via stylia (ersilia-os/stylia). One output file (script-number plot-naming rule):
  output/plots/00b_proteome_descriptors.png
Run with the `gradi` conda env interpreter.
"""

from __future__ import annotations

import io
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np
import pandas as pd
import requests
import stylia
from tenacity import retry, stop_after_attempt, wait_exponential

REPO_ROOT = Path(__file__).resolve().parents[1]
PLOT_PATH = REPO_ROOT / "output" / "plots" / "00b_proteome_descriptors.png"

# (display name, strain label, UniProt proteome ID, local FASTA) — order as the user listed them
PROTEOMES = [
    (
        "E. coli K-12",
        "",
        "UP000000625",
        REPO_ROOT / "data/raw/ecoli/proteome/UP000000625_EcoliK12.fasta",
    ),
    (
        "K. pneumoniae",
        "HS11286",
        "UP000007841",
        REPO_ROOT / "data/raw/kpneumoniae/proteome/UP000007841_HS11286.fasta",
    ),
    (
        "Human",
        "",
        "UP000005640",
        REPO_ROOT / "data/raw/human/proteome/UP000005640_Human.fasta",
    ),
]
PE_LEVELS = ["1", "2", "3", "4", "5"]
PE_LABELS = {
    "1": "PE1 protein",
    "2": "PE2 transcript",
    "3": "PE3 homology",
    "4": "PE4 predicted",
    "5": "PE5 uncertain",
}


def iter_fasta(path: Path):
    header, seq = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line.rstrip("\n"), []
            else:
                seq.append(line.strip())
    if header is not None:
        yield header, "".join(seq)


@retry(stop=stop_after_attempt(5), wait=wait_exponential(multiplier=2, min=2, max=60))
def gene_symbol_fraction(proteome_id: str) -> float:
    """% of entries with a real primary gene symbol (live UniProt; FASTA GN= conflates locus tags)."""
    url = (
        f"https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=tsv"
        f"&query=proteome:{proteome_id}&fields=accession,gene_primary"
    )
    df = pd.read_csv(io.StringIO(requests.get(url, timeout=300).text), sep="\t")
    col = df.columns[-1]
    nonempty = df[col].fillna("").astype(str).str.strip().ne("").sum()
    return 100 * nonempty / len(df)


def style_violin(parts, colors):
    for body, c in zip(parts["bodies"], colors):
        body.set_facecolor(c)
        body.set_edgecolor("none")
        body.set_alpha(0.85)
    for key in ("cbars", "cmins", "cmaxes", "cmedians"):
        if key in parts:
            parts[key].set_color("#444444")
            parts[key].set_linewidth(1.0)


def main() -> None:
    labels, reviewed, unreviewed, lengths, pe_counts, gene_pct = [], [], [], [], [], []
    for name, strain, pid, fasta in PROTEOMES:
        rev = unrev = 0
        L = []
        pe = {lvl: 0 for lvl in PE_LEVELS}
        for header, seq in iter_fasta(fasta):
            if header.startswith(">sp|"):
                rev += 1
            elif header.startswith(">tr|"):
                unrev += 1
            if seq:
                L.append(len(seq))
            m = re.search(r"PE=(\d)", header)
            if m and m.group(1) in pe:
                pe[m.group(1)] += 1
        labels.append(f"{name}\n{strain}" if strain else name)
        reviewed.append(rev)
        unreviewed.append(unrev)
        lengths.append(np.array(L))
        pe_counts.append(pe)
        gene_pct.append(gene_symbol_fraction(pid))
        print(
            f"{name:14s}: {rev + unrev} proteins | reviewed {rev}, unreviewed {unrev} | "
            f"median len {int(np.median(L))} aa | PE1 {pe['1']} | gene-symbol {gene_pct[-1]:.1f}%"
        )

    org_colors = stylia.CategoricalPalette().get(len(PROTEOMES))
    rev_colors = stylia.CategoricalPalette().get(2)
    pe_colors = stylia.CategoricalPalette().get(len(PE_LEVELS))
    x = list(range(len(PROTEOMES)))

    fig, axs = stylia.create_figure(
        2, 2, width=1.9, height=1.0
    )  # flatter, slide-friendly

    # ---- TL: reviewed vs unreviewed (absolute counts) ----
    ax = axs.next()
    totals = [r + u for r, u in zip(reviewed, unreviewed)]
    thr = max(totals) * 0.03
    b1 = ax.bar(x, reviewed, color=rev_colors[0], label="Reviewed (Swiss-Prot)")
    b2 = ax.bar(
        x, unreviewed, bottom=reviewed, color=rev_colors[1], label="Unreviewed (TrEMBL)"
    )
    ax.bar_label(
        b1,
        labels=[f"{v:,}" if v >= thr else "" for v in reviewed],
        label_type="center",
        color="white",
    )
    ax.bar_label(
        b2,
        labels=[f"{v:,}" if v >= thr else "" for v in unreviewed],
        label_type="center",
        color="white",
    )
    ax.bar_label(b2, labels=[f"{t:,}" for t in totals], padding=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("")
    ax.set_ylabel("Number of proteins")
    ax.set_ylim(0, max(totals) * 1.15)
    ax.set_title("Reviewed vs unreviewed")
    ax.legend(loc="upper left")

    # ---- TR: protein-existence level (normalised stacked) ----
    ax = axs.next()
    pe_totals = [sum(pe.values()) for pe in pe_counts]
    bottoms = [0.0] * len(PROTEOMES)
    for lvl, color in zip(PE_LEVELS, pe_colors):
        fracs = [
            100 * pe_counts[i][lvl] / pe_totals[i] if pe_totals[i] else 0 for i in x
        ]
        ax.bar(x, fracs, bottom=bottoms, color=color, label=PE_LABELS[lvl])
        bottoms = [b + f for b, f in zip(bottoms, fracs)]
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("")
    ax.set_ylabel("Share of proteome (%)")
    ax.set_ylim(0, 100)
    ax.set_title("Protein-existence evidence")
    ax.legend(fontsize="small", loc="lower left")

    # ---- BL: gene-symbol coverage ----
    ax = axs.next()
    bars = ax.bar(x, gene_pct, color=org_colors)
    ax.bar_label(bars, labels=[f"{p:.0f}%" for p in gene_pct], padding=3)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("")
    ax.set_ylabel("Entries with a gene symbol (%)")
    ax.set_ylim(0, 110)
    ax.set_title("Gene-symbol coverage")

    # ---- BR: protein length (log scale, no numeric annotation) ----
    ax = axs.next()
    vp = ax.violinplot([np.log10(v) for v in lengths], positions=x, showmedians=True)
    style_violin(vp, org_colors)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("")
    ax.set_ylabel("Protein length (aa, log scale)")
    yt = [1, 2, 3, 4]
    ax.set_yticks(yt)
    ax.set_yticklabels([f"$10^{t}$" for t in yt])
    ax.set_title("Protein length")

    PLOT_PATH.parent.mkdir(parents=True, exist_ok=True)
    stylia.save_figure(str(PLOT_PATH))
    print(f"\nWrote {PLOT_PATH.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
