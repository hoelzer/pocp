import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description='Plot a figure with specified width and height.')
parser.add_argument('--width', type=float, required=True, help='Width of the figure')
parser.add_argument('--height', type=float, required=True, help='Height of the figure')

# Parse the command-line arguments
args = parser.parse_args()

# Read the TSV file
data = pd.read_csv("pocp-matrix.tsv", sep="\t", index_col=0)

# Convert values to numeric
data = data.apply(pd.to_numeric, errors="coerce")

# Check if conversion is successful
if data.isnull().values.any():
    raise ValueError("Unable to convert all values to numeric!")

# Use upper or lower triangle of the matrix
upper_triangle = np.triu(data)
lower_triangle = np.tril(data)

# Combine upper and lower triangles
combined_data = upper_triangle + lower_triangle - np.diag(np.diag(data))

# Create a new DataFrame with the original labels
combined_df = pd.DataFrame(combined_data, index=data.index, columns=data.columns)

# Create a heatmap using seaborn
sns.set(font_scale=1.0)
plt.figure(figsize=(args.width, args.height))

heatmap = sns.heatmap(combined_df, cmap="viridis", annot=True, fmt=".1f", linewidths=.5, square=True, cbar_kws={"shrink": 0.6})

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
