import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Read the TSV file
#data = pd.read_csv("pocp-diff.tsv", sep="\t", index_col=0)
data = pd.read_csv("pocp-diff-old-vs-new-blast.tsv", sep="\t", index_col=0)
#data = pd.read_csv("brucella-pocp-diff.tsv", sep="\t", index_col=0)
#data = pd.read_csv("enterococcus-pocp-diff.tsv", sep="\t", index_col=0)

# Convert values to numeric
data = data.apply(pd.to_numeric, errors="coerce")

# Check if conversion is successful
#if data.isnull().values.any():
#    raise ValueError("Unable to convert all values to numeric!")

# Use upper or lower triangle of the matrix
upper_triangle = np.triu(data)
lower_triangle = np.tril(data)

# Combine upper and lower triangles
combined_data = upper_triangle + lower_triangle - np.diag(np.diag(data))

# Create a new DataFrame with the original labels
combined_df = pd.DataFrame(combined_data, index=data.index, columns=data.columns)

# Create a heatmap using seaborn
sns.set(font_scale=1.0)
plt.figure(figsize=(16, 8))

heatmap = sns.heatmap(combined_df, cmap="viridis", annot=True, fmt=".1f", linewidths=.5, square=True, cbar_kws={"shrink": 0.6})

# Move the x-axis labels to the top
plt.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)

# Set axis labels and plot title
plt.xlabel("")
plt.ylabel("")
plt.title("")

# Rotate y-axis labels for better readability
plt.yticks(rotation=0)

# Save the figure as an SVG file
#heatmap.get_figure().savefig("pocp-diff.svg", format="svg", bbox_inches="tight")
heatmap.get_figure().savefig("pocp-diff-old-vs-new-blast.svg", format="svg", bbox_inches="tight")
#heatmap.get_figure().savefig("brucella-pocp-diff.svg", format="svg", bbox_inches="tight")
#heatmap.get_figure().savefig("enterococcus-pocp-diff.svg", format="svg", bbox_inches="tight")

# Show the plot
plt.show()
