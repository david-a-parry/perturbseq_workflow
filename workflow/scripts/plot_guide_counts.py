import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
matplotlib.use("Agg")

df = pd.DataFrame()
for csv_input in snakemake.input:
    tmp_df = pd.read_csv(csv_input)
    tmp_df['filename'] = csv_input
    df = df.append(tmp_df)
df['sample_name'] = df.filename.apply(
    lambda x: "-".join(os.path.basename(x).split('-')[:-1]))
count_df = df.groupby(
    ['sample_name', 'Target', 'Gene'])['Coverage'].sum().to_frame()
count_df.reset_index(inplace=True)
g = sns.FacetGrid(count_df, row='sample_name')
g.map(sns.distplot, 'Coverage')
plt.savefig("results/plots/read_counts.pdf")
