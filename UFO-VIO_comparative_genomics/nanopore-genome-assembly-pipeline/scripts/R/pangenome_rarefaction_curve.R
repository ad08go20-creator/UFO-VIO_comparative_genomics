#!/usr/bin/env Rscript
# Rarefaction curve of conserved vs. total pangenome genes from Roary output.
# Inputs: number_of_conserved_genes.Rtab, number_of_genes_in_pan_genome.Rtab

library(ggplot2)

conserved_file <- "number_of_conserved_genes.Rtab"
total_file <- "number_of_genes_in_pan_genome.Rtab"

conserved <- colMeans(read.table(conserved_file))
total <- colMeans(read.table(total_file))

data <- data.frame(
    genes = c(conserved, total),
    genomes = c(seq_along(conserved), seq_along(total)),
    type = c(
        rep("Conserved genes", length(conserved)),
        rep("Total genes", length(total))
    )
)

p <- ggplot(data = data, aes(x = genomes, y = genes, color = type)) +
    geom_line(size = 1.2) +
    theme_classic() +
    xlab("Number of genomes") +
    ylab("Number of genes") +
    theme_bw(base_size = 16) +
    theme(legend.position = "right")

ggsave(filename = "pangenome_rarefaction.png", plot = p, width = 10, height = 10)
ggsave(filename = "pangenome_rarefaction.svg", plot = p, width = 10, height = 10)
