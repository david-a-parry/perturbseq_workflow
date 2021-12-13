import pysam
import os
import glob
import logging
import pandas as pd
from collections import defaultdict


logger = logging.getLogger("Read Counts")
logger.setLevel(logging.INFO)
formatter = logging.Formatter(
    '[%(asctime)s] %(name)s - %(levelname)s - %(message)s')
ch = logging.StreamHandler()
ch.setLevel(logger.level)
ch.setFormatter(formatter)
logger.addHandler(ch)


def target_from_ref(seqname):
    '''Return the gene, strand and exon from guide fasta sequence name'''
    try:
        chrom, start, end, gene, strand, _, exon = seqname.split('_')
    except ValueError:
        try:
            chrom, exon, end, gene, strand, _ = seqname.split('_')
        except ValueError:
            raise RuntimeError(
                "Could not parse sequence name: {}".format(seqname))
    strand = strand.replace('minus', '-').replace('plus', '+')
    return (gene, strand, exon)


def get_guide_coverage(bam, ref, minimum_aligned_length=0,
                       minimum_fraction_aligned=0.0):
    '''
        Determine the number of reads aligned to sequences in FASTA reference
        file in BAM alignment file. Returns a pandas dataframe of Target names
        and read counts.

        Args:
            bam:    Alignment file in BAM/SAM format

            ref:    FASTA reference file (same as used for creating BAM/SAM
                    alignment)

            minimum_aligned_length:
                    Minimum number of bases that must have been aligned in
                    order to count a record. Default=0

            minimum_fraction_aligned:
                    Minimum fraction of bases that must have been aligned in
                    order to count a record. Default=0.0

    '''
    bamfile = pysam.AlignmentFile(bam)
    fasta = pysam.FastaFile(ref)
    total_reads = 0
    counts = defaultdict(int)
    mapping = {'mapped': [0], 'unmapped': [0], 'filtered': [0], 'passed': [0]}
    for read in bamfile.fetch(until_eof=True):
        total_reads += 1
        if read.is_unmapped:
            mapping['unmapped'][0] += 1
        else:
            mapping['mapped'][0] += 1
            if minimum_aligned_length:
                if read.query_alignment_length < minimum_aligned_length:
                    mapping['filtered'][0] += 1
                    continue
            if minimum_fraction_aligned:
                if (read.query_alignment_length/read.query_length <
                        minimum_fraction_aligned):
                    mapping['filtered'][0] += 1
                    continue
            counts[read.reference_name] += 1
            mapping['passed'][0] += 1
    logger.info("{:,}/{:,} mapped reads ({:g})".format(
        mapping['mapped'][0], total_reads,
        mapping['mapped'][0]/total_reads))
    count_rows = {'Target': [], 'Coverage': []}
    for seq in fasta.references:
        count_rows['Target'].append(seq)
        if seq in counts:
            count_rows['Coverage'].append(counts[seq])
        else:
            count_rows['Coverage'].append(0)
    df = pd.DataFrame.from_dict(count_rows)
    df['Gene'] = df.Target.apply(lambda x: target_from_ref(x)[0])
    mapping_df = pd.DataFrame.from_dict(mapping)
    return df, mapping_df


count_df, map_df = get_guide_coverage(snakemake.input[0], snakemake.input[1])
count_df.to_csv(snakemake.output['read_counts'], index=False)
map_df.to_csv(snakemake.output['map_counts'], index=False)
logger.info("Read counts written to {}".format(
    snakemake.output['read_counts']))
