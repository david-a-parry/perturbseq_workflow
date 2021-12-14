import pandas as pd


df = pd.read_excel(snakemake.input[0])
df['reference_name'] = df.GUIDE_ID.apply(
    lambda x: x.replace(':', '_').replace('-', '_', 1).replace(
        '-', 'minus', 1).replace('+', 'plus', 1)) + "_" + df['TARGET EXON']
with open(snakemake.output[0], 'wt') as fh:
    df.apply(lambda x: fh.write(">{}\n{}\n".format(x.reference_name,
                                                   x.SEQUENCE)),
             axis=1)
