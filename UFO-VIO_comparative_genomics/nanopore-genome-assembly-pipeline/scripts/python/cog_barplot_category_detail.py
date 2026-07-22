#!/usr/bin/env python
"""Bar chart of the top-scoring gene functions within specific COG categories.

Input: the raw .emapper.annotations file produced by eggNOG-mapper.
Edit COLOR_MAP below to select which COG categories to plot.
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- Configuration ---
EGGNOG_ANNOTATIONS_FILE = "MM_wd3sxapu.emapper.annotations"

# Color per category of interest
COLOR_MAP = {
    "Q": "#2eabb8ff",  # Secondary metabolites biosynthesis, transport and catabolism
    "V": "#d1d197ff",  # Defense mechanisms
}
CATEGORIES_TO_PLOT = list(COLOR_MAP.keys())


def generate_detailed_plot(df, cog_letter, plot_color):
    """Filter by COG category and plot the top annotated functions."""
    print(f"\n--- Processing category: {cog_letter} ---")

    category_df = df[df["Cleaned_COG"] == cog_letter].copy()

    if category_df.empty:
        print(f"No genes found for category '{cog_letter}'.")
        return

    plot_data = category_df.loc[category_df.groupby("Description")["score"].idxmax()]
    plot_data = plot_data.sort_values("score", ascending=False).head(25)

    print(f"Found {len(plot_data)} unique functions to plot.")

    plt.figure(figsize=(12, 10))
    sns.barplot(data=plot_data, x="score", y="Description", color=plot_color, orient="h")

    plt.title(f'Detailed Functions for COG Category: "{cog_letter}"', fontsize=16, weight="bold")
    plt.xlabel("Annotation Score (Bit-score)", fontsize=12)
    plt.ylabel("Function Description", fontsize=12)
    plt.tight_layout()

    output_filename_svg = f"detailed_plot_category_{cog_letter}.svg"
    plt.savefig(output_filename_svg, bbox_inches="tight")

    print(f"Plot saved as '{output_filename_svg}'")
    plt.close()


print(f"Loading data from '{EGGNOG_ANNOTATIONS_FILE}'...")
try:
    header_row = 0
    with open(EGGNOG_ANNOTATIONS_FILE, "r") as f:
        for i, line in enumerate(f):
            if line.startswith("#query"):
                header_row = i
                break
    main_df = pd.read_csv(EGGNOG_ANNOTATIONS_FILE, sep="\t", header=header_row)
except FileNotFoundError:
    print(f"ERROR: file not found '{EGGNOG_ANNOTATIONS_FILE}'.")
    exit()

print("Cleaning COG categories...")
main_df.dropna(subset=["COG_category"], inplace=True)
main_df["Cleaned_COG"] = main_df["COG_category"].apply(lambda x: str(x)[0] if str(x) != "-" else None)

for category, color in COLOR_MAP.items():
    generate_detailed_plot(main_df, category, color)

print("\nDone.")