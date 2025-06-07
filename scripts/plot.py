import matplotlib.pyplot as plt
import pandas as pd
import glob

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'text.usetex': True,
    'axes.labelsize': 14,
    'axes.titlesize': 15,
    'legend.fontsize': 12,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 1.8,
    'lines.markersize': 3
})

colors = ['#377eb8', '#4daf4a', '#e41a1c', '#984ea3', '#ff7f00']
markers = ['o', 's', 'D', '^', 'v']

csv_files = sorted(glob.glob("degree_*.csv"))

for csv_file in csv_files:
    degree = csv_file.split('_')[1].split('.')[0]
    data = pd.read_csv(csv_file, skipinitialspace=True)

    plt.figure(figsize=(8, 6))

    methods = data['Method'].unique()

    for i, method in enumerate(methods):
        subset = data[data['Method'] == method]
        plt.plot(subset['Control Points'], subset['Time (s)'], 
                 marker=markers[i % len(markers)], linestyle='-', color=colors[i % len(colors)], label=method)

    plt.title(f'Benchmark Results for Degree {degree}')
    plt.xlabel('Number of Control Points')
    plt.ylabel('Time (s)')
    plt.legend(title='Method', loc='best', frameon=True)
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f'benchmark_degree_{degree}.png', format='png', dpi=300)
    plt.close()
