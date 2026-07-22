#!/usr/bin/env python
# Copyright (C) <2015> EMBL-European Bioinformatics Institute
#
# Modified version of roary_plots.py (Marco Galardini) producing
# publication-ready SVG output with a cleaner visual style.
#
# Usage:
#   python roary_plots_custom.py <tree.newick> <gene_presence_absence.csv>

__author__ = "Marco Galardini"
__version__ = "0.9.0"


def get_options():
    """Parse command-line arguments."""
    import argparse

    description = "Create enhanced plots from Roary outputs."
    parser = argparse.ArgumentParser(description=description,
                                      prog="roary_plots_custom.py")

    parser.add_argument("tree", action="store",
                         help="Newick tree file (e.g., accessory_binary_genes.fa.newick)")
    parser.add_argument("spreadsheet", action="store",
                         help="Roary gene presence/absence spreadsheet (e.g., gene_presence_absence.csv)")

    parser.add_argument("--labels", action="store_true",
                         default=False,
                         help="Add node labels to the tree (up to 10 chars)")
    parser.add_argument("--format",
                         choices=("png", "tiff", "pdf", "svg"),
                         default="svg",
                         help="Output format [default: svg]")
    parser.add_argument("-N", "--skipped-columns", action="store",
                         type=int,
                         default=14,
                         help="First N columns of Roary's output to exclude [default: 14]")

    parser.add_argument("--version", action="version",
                         version="%(prog)s " + __version__)

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    from Bio import Phylo
    from matplotlib.colors import ListedColormap

    # --- Global style configuration ---
    sns.set_theme(style="ticks")
    plt.rc("text", color="black")
    plt.rc("axes", labelcolor="black", edgecolor="black")
    plt.rc("xtick", color="black")
    plt.rc("ytick", color="black")
    plt.rc("figure", dpi=300)

    # --- Load data ---
    print("-> Loading tree and Roary spreadsheet...")
    t = Phylo.read(options.tree, "newick")
    roary = pd.read_csv(options.spreadsheet, low_memory=False)
    roary.set_index("Gene", inplace=True)
    roary.drop(list(roary.columns[:options.skipped_columns - 1]), axis=1, inplace=True)
    roary.replace(".{2,100}", 1, regex=True, inplace=True)
    roary.replace(np.nan, 0, regex=True, inplace=True)
    print("-> Data loaded successfully.")

    # --- Pangenome category thresholds (computed once, reused by all plots) ---
    num_genomes = roary.shape[1]
    core_threshold = num_genomes * 0.99
    soft_core_threshold = num_genomes * 0.95
    shell_threshold = num_genomes * 0.15

    gene_sums = roary.sum(axis=1)

    # --- Plot 1: pangenome gene-frequency histogram ---
    print("-> Generating pangenome frequency plot...")

    freq_df = pd.DataFrame({"frequency": gene_sums})

    def assign_category_hist(freq):
        if freq >= core_threshold:
            return "Core genes"
        elif freq >= soft_core_threshold:
            return "Soft-core genes"
        elif freq >= shell_threshold:
            return "Shell genes"
        else:
            return "Cloud genes"

    freq_df["Category"] = freq_df["frequency"].apply(assign_category_hist)

    comp_palette = sns.color_palette("Pastel1", 4)
    color_map = {
        "Cloud genes": comp_palette[0],
        "Shell genes": comp_palette[1],
        "Soft-core genes": comp_palette[2],
        "Core genes": comp_palette[3],
    }

    plt.figure(figsize=(10, 6))
    sns.histplot(data=freq_df, x="frequency", hue="Category",
                 multiple="stack",
                 bins=num_genomes,
                 palette=color_map,
                 edgecolor="white",
                 linewidth=0.5,
                 hue_order=["Cloud genes", "Shell genes", "Soft-core genes", "Core genes"])

    plt.xlim(-0.5, num_genomes + 0.5)

    plt.title("Pangenome Gene Frequency", fontsize=16, fontweight="bold")
    plt.xlabel("Number of Genomes", fontsize=12, fontweight="bold")
    plt.ylabel("Number of Genes", fontsize=12, fontweight="bold")
    sns.despine()
    plt.tight_layout()
    plt.savefig(f"pangenome_frequency.{options.format}")
    plt.clf()

    # --- Plot 2: phylogenetic tree + presence/absence matrix ---
    print("-> Generating tree and presence/absence matrix plot...")
    roary_sorted = roary.loc[roary.sum(axis=1).sort_values(ascending=False).index]
    roary_sorted = roary_sorted[[x.name for x in t.get_terminals()]]

    fig = plt.figure(figsize=(17, 10))
    gs = fig.add_gridspec(1, 40, wspace=0)

    ax_tree = fig.add_subplot(gs[0, 0:10])
    with plt.rc_context({"lines.linewidth": 0.5}):
        Phylo.draw(t, axes=ax_tree, show_confidence=False, label_func=lambda x: None,
                   xticks=([],), yticks=([],), ylabel=("",), xlabel=("",),
                   axis=("off",), title=(f"Phylogenetic Tree\n({roary.shape[1]} strains)",),
                   do_show=False)

    ax_matrix = fig.add_subplot(gs[0, 10:40])
    # White for absence (0), pastel green for presence (1)
    pastel_cmap = ListedColormap(["#e9f7f1ff", "#51be8eff"])
    ax_matrix.matshow(roary_sorted.T, cmap=pastel_cmap, aspect="auto", interpolation="none")
    ax_matrix.set_yticks([])
    ax_matrix.set_xticks([])
    ax_matrix.axis("off")
    ax_matrix.set_title(f"Gene Presence/Absence\n({roary.shape[0]} gene clusters)")

    plt.savefig(f"pangenome_matrix.{options.format}")
    plt.clf()

    # --- Donut chart categories ---
    core_genes_count = (gene_sums >= core_threshold).sum()
    soft_core_genes_count = ((gene_sums >= soft_core_threshold) & (gene_sums < core_threshold)).sum()
    shell_genes_count = ((gene_sums >= shell_threshold) & (gene_sums < soft_core_threshold)).sum()
    cloud_genes_count = (gene_sums < shell_threshold).sum()

    pangenome_counts = {
        "Core genes": core_genes_count,
        "Soft-core genes": soft_core_genes_count,
        "Shell genes": shell_genes_count,
        "Cloud genes": cloud_genes_count,
    }
    pangenome_df = pd.DataFrame(list(pangenome_counts.items()), columns=["Category", "Count"])
    total = pangenome_df["Count"].sum()

    # --- Plot 3: pangenome composition donut chart ---
    print("-> Generating pangenome donut chart...")
    plt.figure(figsize=(10, 8))
    colors = sns.color_palette("Pastel1", 4)
    plt.pie(pangenome_df["Count"], labels=pangenome_df["Category"],
            autopct=lambda p: "{:.0f}".format(p * total / 100),
            colors=colors, wedgeprops={"edgecolor": "white", "linewidth": 1.5},
            pctdistance=0.85, startangle=90,
            textprops={"color": "black"})

    my_circle = plt.Circle((0, 0), 0.7, color="white")
    p = plt.gcf()
    p.gca().add_artist(my_circle)

    plt.title("Pangenome Composition", fontsize=16, fontweight="bold")
    plt.savefig(f"pangenome_pie.{options.format}")
    plt.clf()

    print("\nAll plots generated successfully!")