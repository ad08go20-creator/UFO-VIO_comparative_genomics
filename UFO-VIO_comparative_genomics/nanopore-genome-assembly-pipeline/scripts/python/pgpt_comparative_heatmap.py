"""
Comparative PGPT (PLaBAse) gene-content heatmap across multiple genomes.

Input: a tab-separated raw gene-count table (rows = PGPT category names,
columns = genome names, values = gene counts), as produced by PGPg_finder's
`gene_counts` output for each genome in the comparison set, merged into one
table.

PGPg_finder's own `heatmap.py` normalizes each genome's counts against the
grand total of *every* genome submitted in that batch run — this makes
outputs from separate batch runs (e.g., one genome run alone vs. a larger
multi-genome run) numerically incompatible. If your genomes were annotated
in separate PGPg_finder runs, recompute raw counts into a single merged
table first (do not merge the tool's own "normalized" outputs across runs).

Usage:
    python pgpt_comparative_heatmap.py raw_gene_counts.tsv output_prefix
"""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams['font.family'] = 'Arial'


def clean_label(x: str) -> str:
    return x.replace('|', ' / ').replace('_', ' ').title()


def build_heatmap(counts_path: str, output_prefix: str, genome_order: list[str] | None = None,
                   highlight_genome: str | None = None, annotate_counts: bool = False) -> None:
    df = pd.read_csv(counts_path, sep='\t', index_col=0)
    if genome_order:
        df = df[genome_order]
    df['total'] = df.sum(axis=1)
    df = df.sort_values('total', ascending=False).drop(columns='total')

    row_labels = [clean_label(c) for c in df.index]
    n_rows, n_cols = df.shape
    raw = df.values
    logvals = np.log10(raw + 1)

    cmap = LinearSegmentedColormap.from_list("pgpt_green", ["#F7FCF5", "#A1D99B", "#41AB5D", "#00441B"])

    fig, ax = plt.subplots(figsize=(9.5, 0.25 * n_rows + 1))
    im = ax.imshow(logvals, cmap=cmap, aspect='auto')

    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticklabels(row_labels, fontsize=7)

    if highlight_genome and highlight_genome in df.columns:
        idx = list(df.columns).index(highlight_genome)
        ax.get_xticklabels()[idx].set_fontstyle('italic')
        ax.get_xticklabels()[idx].set_fontweight('bold')

    if annotate_counts:
        vmax = logvals.max()
        for i in range(n_rows):
            for j in range(n_cols):
                val = logvals[i, j]
                color = 'white' if val > vmax * 0.6 else 'black'
                ax.text(j, i, str(int(raw[i, j])), ha='center', va='center', fontsize=5.5, color=color)

    # thin white gridlines between cells + thin black border frame
    # (kept consistent across all comparative heatmaps in this project: ANI, KEGG, PGPT)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.6)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=2, labelsize=7)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color('black')

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('log10(gene count + 1)', fontsize=9)
    cbar.ax.tick_params(labelsize=7.5)

    plt.tight_layout()
    fig.savefig(f'{output_prefix}.png', dpi=400, facecolor='white', bbox_inches='tight')
    fig.savefig(f'{output_prefix}.svg', facecolor='white', bbox_inches='tight')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    build_heatmap(sys.argv[1], sys.argv[2])
