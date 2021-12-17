import pandas as pd


df = pd.read_excel(snakemake.input[0])
exon_col = 'TARGET EXON' if 'TARGET EXON' in df.columns else 'EXON'
df['reference_name'] = df.GUIDE_ID.apply(
    lambda x: x.replace(':', '_').replace('-', '_', 1).replace(
        '-', 'minus', 1).replace('+', 'plus', 1)) + "_" + df[exon_col]
with open(snakemake.output[0], 'wt') as fh:
    df.apply(lambda x: fh.write(">{}\n{}\n".format(x.reference_name,
                                                   x.SEQUENCE)),
             axis=1)
