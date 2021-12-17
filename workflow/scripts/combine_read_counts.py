import os
import pandas as pd

# read in per-bam read counts
df = pd.DataFrame()
for csv_input in snakemake.input:
    tmp_df = pd.read_csv(csv_input)
    tmp_df['filename'] = csv_input
    df = df.append(tmp_df)
df['sample_name'] = df.filename.apply(
    lambda x: "-".join(os.path.basename(x).split('-')[:-1]))
df.to_csv(snakemake.output[0], index=False)
count_df = df.groupby(
    ['sample_name', 'Target', 'Gene'])['Coverage'].sum().to_frame()
count_df.reset_index(inplace=True)
samp2mapped = count_df.groupby(['sample_name'])['Coverage'].sum().to_dict()
count_df['fraction_reads'] = count_df.apply(
    lambda x: x.Coverage/samp2mapped[x.sample_name],
    axis=1)
count_df.to_csv(snakemake.output[1], index=False)
pivot = df.pivot_table(index=["Target", 'Gene'],
                       columns='sample_name')['Coverage']
pivot.to_csv(snakemake.output[2], sep='\t')
# get sample information, T0 timepoint and test sample IDs
sample_df = pd.read_csv(snakemake.config['samples'], sep='\t')
# write out read counts for each genotype
for gt in sample_df.genotype.unique():
    sample_ids = list(sample_df[sample_df.genotype == gt].sample_name.values)
    pivot[sample_ids].to_csv("results/read_counts/{}_read_counts.tsv".format(gt),
                             sep='\t')
