#!/usr/bin/env python
"""Horizontal bar chart of gene counts per COG functional category.

Input: the .emapper.cog.tsv file produced by eggnogCOGextractor.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
INPUT_FILE = "MM_wd3sxapu.emapper.cog.tsv"
OUTPUT_PLOT_PNG = "cog_bar_chart.png"
OUTPUT_PLOT_SVG = "cog_bar_chart.svg"

# --- Load data ---
print(f"Loading data from '{INPUT_FILE}'...")
try:
    df = pd.read_csv(INPUT_FILE, sep="\t")
except FileNotFoundError:
    print(f"ERROR: file not found '{INPUT_FILE}'.")
    exit()
except pd.errors.EmptyDataError:
    print(f"ERROR: file '{INPUT_FILE}' is empty.")
    exit()

df = df[df["COGCount"] > 0]
df = df.sort_values("COGCount", ascending=False)

print("Generating bar chart...")
fig, ax = plt.subplots(figsize=(12, 8))

sns.barplot(data=df,
            x="COGCount",
            y="COGDescription",
            hue="COGDescription",
            palette="tab20",
            orient="h",
            legend=False,
            ax=ax)

# Exact count at the end of each bar
for index, value in enumerate(df["COGCount"]):
    ax.text(value, index, f" {value}", va="center", ha="left")

# COG category letter on the right-hand side
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax2 = ax.twinx()
ax2.spines["right"].set_visible(False)
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(ax.get_yticks())
ax2.set_yticklabels(df["COGID"], fontsize=10, weight="bold")

ax.set_title("Number of Genes per COG Functional Category", fontsize=16, weight="bold")
ax.set_xlabel("Number of Genes", fontsize=12)
ax.set_ylabel("COG Category", fontsize=12)

plt.tight_layout()
plt.savefig(OUTPUT_PLOT_PNG, dpi=300, bbox_inches="tight")
plt.savefig(OUTPUT_PLOT_SVG, bbox_inches="tight")

print(f"Saved '{OUTPUT_PLOT_PNG}' and '{OUTPUT_PLOT_SVG}'")