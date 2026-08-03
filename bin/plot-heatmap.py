#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Plot a figure with specified width and height.')
parser.add_argument('--matrix', default='pocp-matrix.tsv', help='TSV file holding the pairwise POCP matrix')
parser.add_argument('--width', type=float, required=True, help='Width of the figure')
parser.add_argument('--height', type=float, required=True, help='Height of the figure')

# Parse the command-line arguments
args = parser.parse_args()

# Read the TSV file
data = pd.read_csv(args.matrix, sep="\t", index_col=0)

# Convert values to numeric, pairs that were not compared (one-vs-all mode) become NaN
# and are drawn as empty cells
data = data.apply(pd.to_numeric, errors="coerce")

# Check if conversion is successful
if data.isnull().values.all():
    raise ValueError("Unable to convert any value to numeric!")

# Create a heatmap using seaborn
sns.set_theme(font_scale=1.0)
plt.figure(figsize=(args.width, args.height))

heatmap = sns.heatmap(data, cmap="viridis", annot=True, fmt=".1f", linewidths=.5, square=True, cbar_kws={"shrink": 0.6})

# Move the x-axis labels to the top
plt.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

# Decide whether to rotate top x-axis labels
renderer = plt.gcf().canvas.get_renderer()
label_widths = [label.get_window_extent(renderer=renderer).width for label in heatmap.get_xticklabels()]
total_label_width = sum(label_widths)
fig_width_px = plt.gcf().bbox.width

if total_label_width * 2 > fig_width_px:
    plt.xticks(rotation=90)  # Rotate to vertical if overlapping
else:
    plt.xticks(rotation=0)

# Rotate y-axis labels for better readability
plt.yticks(rotation=0)

# Set axis labels and plot title
plt.xlabel("")
plt.ylabel("")
plt.title("Pairwise percentage of conserved proteins (POCP)", fontweight="bold")

# Save the figure as an SVG and PDF file
heatmap.get_figure().savefig("pocp-heatmap.svg", format="svg", bbox_inches="tight")
heatmap.get_figure().savefig("pocp-heatmap.pdf", format="pdf", bbox_inches="tight")

# Optional: Show the plot interactively
# plt.show()
