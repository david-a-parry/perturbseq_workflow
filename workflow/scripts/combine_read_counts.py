import os
import pandas as pd

df = pd.DataFrame()
for csv_input in snakemake.input:
    tmp_df = pd.read_csv(csv_input)
    tmp_df['filename'] = csv_input
    df = df.append(tmp_df)
df['sample_name'] = df.filename.apply(
    lambda x: "-".join(os.path.basename(x).split('-')[:-1]))
df.to_csv("results/read_counts_per_alignment.csv.gz", index=False)
count_df = df.groupby(
    ['sample_name', 'Target', 'Gene'])['Coverage'].sum().to_frame()
count_df.reset_index(inplace=True)
count_df.to_csv("results/read_counts_per_sample.csv.gz", index=False)
