"""
Comparative KEGG module completeness heatmap across multiple genomes.

Input: a tab-separated completeness matrix (rows = KEGG module IDs or
curated module names, columns = genome names, values = completeness
fraction 0-1), as produced by Anvi'o's `anvi-estimate-metabolism`
(`*-module_stepwise_completeness-MATRIX.txt` / `*-module_pathwise_completeness-MATRIX.txt`).

If your matrix is indexed by raw KEGG module IDs (M00001, M00002, ...),
map them to readable names first using Anvi'o's own `module_names.txt`
(shipped alongside the completeness matrix under
ANVIO/METABOLISM/STEPWISE/), e.g.:

    names = dict(line.strip().split('\t', 1) for line in open('module_names.txt'))
    df.index = [names.get(m, m) for m in df.index]

Usage:
    python kegg_comparative_heatmap.py completeness_matrix.tsv output_prefix
"""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams['font.family'] = 'Arial'


def build_heatmap(matrix_path: str, output_prefix: str, genome_order: list[str] | None = None,
                   highlight_genome: str | None = None) -> None:
    df = pd.read_csv(matrix_path, sep='\t', index_col=0)
    if genome_order:
        df = df[genome_order]
    row_labels = list(df.index)
    n_rows, n_cols = df.shape
    vals = df.values.astype(float)

    cmap = LinearSegmentedColormap.from_list("kegg_blue", ["#F7FBFF", "#6BAED6", "#2171B5", "#08306B"])

    fig, ax = plt.subplots(figsize=(9, 0.5 * n_rows + 1))
    im = ax.imshow(vals, cmap=cmap, vmin=0, vmax=1, aspect='auto')

    ax.set_xticks(range(n_cols))
    ax.set_yticks(range(n_rows))
    ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=9)
    ax.set_yticklabels(row_labels, fontsize=9)

    if highlight_genome and highlight_genome in df.columns:
        idx = list(df.columns).index(highlight_genome)
        ax.get_xticklabels()[idx].set_fontstyle('italic')
        ax.get_xticklabels()[idx].set_fontweight('bold')

    # thin white gridlines between cells + thin black border frame
    # (kept consistent across all comparative heatmaps in this project: ANI, KEGG, PGPT)
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which='minor', color='white', linewidth=0.8)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=2, labelsize=9)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.8)
        spine.set_color('black')

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label('Completeness', fontsize=10)
    cbar.ax.tick_params(labelsize=8)

    plt.tight_layout()
    fig.savefig(f'{output_prefix}.png', dpi=400, facecolor='white', bbox_inches='tight')
    fig.savefig(f'{output_prefix}.svg', facecolor='white', bbox_inches='tight')


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(1)
    build_heatmap(sys.argv[1], sys.argv[2])
