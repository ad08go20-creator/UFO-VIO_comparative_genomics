#!/usr/bin/env Rscript
# UpSet plot of genes shared between a genome of interest and a set of
# reference genomes, based on Roary's gene_presence_absence.csv.
#
# Required packages: tidyverse, UpSetR, svglite

library(tidyverse)
library(UpSetR)
library(svglite)

# --- Configuration ---
genome_cols <- c(
    "A56AF",
    "BASUSDA_45",
    "CV1",
    "CV8",
    "Cv1",
    "DSM26508",
    "SAM215",
    "chromo_utfo"
)
genome_of_interest <- "UFO_VIO"
tab20_colors <- c(
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
    "#8C564B", "#E377C2", "#7F7F7F", "#9467BD"
)

# --- Load and prepare data ---
df <- read_csv("gene_presence_absence.csv")

df_binary <- df %>%
    select(all_of(genome_cols)) %>%
    mutate(across(everything(), ~ as.integer(!is.na(.) & . != ""))) %>%
    rename(!!genome_of_interest := chromo_utfo)

df_filtered <- df_binary %>%
    filter(!!sym(genome_of_interest) == 1)

if (nrow(df_filtered) == 0) {
    stop(paste0("No genes found for '", genome_of_interest, "'."), call. = FALSE)
}

# --- Build the plot ---
p <- upset(
    as.data.frame(df_filtered),
    nsets = ncol(df_filtered),
    order.by = "freq",
    decreasing = TRUE,
    main.bar.color = tab20_colors[8],
    sets.bar.color = tab20_colors,
    matrix.color = tab20_colors[8],
    mainbar.y.label = paste("Genes in intersection (with", genome_of_interest, ")"),
    sets.x.label = paste("Genes shared with", genome_of_interest),
    text.scale = 1.2
)

# --- Save to SVG ---
svglite::svglite("pangenome_upset.svg", width = 12, height = 7)
print(p)
dev.off()

message("Plot saved as pangenome_upset.svg")