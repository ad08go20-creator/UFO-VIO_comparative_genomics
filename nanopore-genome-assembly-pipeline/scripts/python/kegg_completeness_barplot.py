#!/usr/bin/env python
"""Bar plot of KEGG pathway completeness from a KEGGaNOG *_pathways.tsv file.

TODO: extend to export SVG output and expose a configurable color palette.
"""
import kegganog as kgn
import pandas as pd

INPUT_FILE = "SAMPLE_pathways.tsv"
OUTPUT_PNG = "KEGG_barplot.png"

df = pd.read_csv(INPUT_FILE, sep="\t")

kgnbar = kgn.barplot(
    df,
    figsize=(8, 10),
    sort_order="descending",
    yticks_fontsize=8,
    cmap="Blues",
)

fig = kgnbar.plotfig()
fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")

print(f"Saved '{OUTPUT_PNG}'")