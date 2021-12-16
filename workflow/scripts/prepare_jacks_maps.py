import os
import pandas as pd


df = pd.DataFrame()
sample_df = pd.read_csv(snakemake.config['samples'], sep='\t')
with open(snakemake.output[0], 'wt') as fh:
    fh.write("\t".join(["Replicate", "Sample", "Control"]) + "\n")
    for gt in sample_df.genotype.unique():
        gt_df = sample_df[sample_df.genotype == gt]
        t0 = gt_df.timepoint.min()
        t0_smp = "{}-T{}".format(gt, t0)
        for tp in gt_df.timepoint.unique():
            for smpl in gt_df[gt_df.timepoint == tp].sample_name.unique():
                fh.write("\t".join(
                    [smpl, "{}-T{}".format(gt, tp), t0_smp]) + "\n")
# write sgrna map
read_counts = pd.read_csv(snakemake.input[0], sep='\t')
read_counts[['Target', 'Gene']].to_csv(snakemake.output[1],
                                       sep='\t',
                                       index=False)
