import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
from mpl_toolkits.mplot3d import Axes3D

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update(
    {
        "font.size": 12,
        "font.family": "serif",
        "text.usetex": True,
        "axes.labelsize": 14,
        "axes.titlesize": 15,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "lines.linewidth": 1.9,
        "lines.markersize": 0,
    }
)

colors = ["#377eb8", "#4daf4a", "#e41a1c", "#984ea3", "#ff7f00"]
markers = ["o", "s", "D", "^", "v"]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

csv_files = sorted(glob.glob("degree_*.csv"))

for idx, csv_file in enumerate(csv_files):
    degree = int(csv_file.split("_")[1].split(".")[0])
    data = pd.read_csv(csv_file, skipinitialspace=True)

    methods = data["Method"].unique()

    for i, method in enumerate(methods):
        subset = data[data["Method"] == method]
        xs = subset["Control Points"]
        ys = np.full_like(xs, degree)
        zs = subset["Time (s)"]

        ax.plot(
            xs,
            ys,
            zs,
            marker=markers[i % len(markers)],
            linestyle="-",
            color=colors[i % len(colors)],
            label=f"{method}" if idx == 0 else "",
        )

ax.set_xlabel("Number of Control Points")
ax.set_ylabel("Degree")
ax.set_zlabel("Time (s)")
ax.set_title("Benchmark Results")

# Handle legend uniquely to avoid duplicate entries
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), title="Method", loc="best")

plt.tight_layout()
plt.savefig("benchmark.png", format="png", dpi=300)
