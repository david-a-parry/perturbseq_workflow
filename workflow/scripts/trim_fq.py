import os
import glob
import logging
import gzip
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import defaultdict, OrderedDict


logging.basicConfig(filename=snakemake.log[0],
                    level=logging.INFO,
                    format='%(asctime)s %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    )


def _get_open_function(f):
    if f.endswith('.gz'):
        return gzip.open
    else:
        return open


def trim_fastq(fastq,
               output,
               five_prime_trim,
               seq_length,
               progress_interval=5e6):
    '''
        Remove five_prime_trim bases from start of each sequence
        in fastq and retain seq_length bases from sequence.

        Args:
            fastq:  filename for FASTQ input file

            output: output filename for trimmed FASTQ

            five_prime_trim:
                    number of nucleotides to trim from beginning of each
                    read.

            seq_length:
                    length of sequence to keep. Each sequence will be
                    sliced as follows:
                        seq[five_prime_trim:five_prime_trim+seq_length]
    '''
    fastq_func = _get_open_function(fastq)
    out_funct = _get_open_function(output)
    with fastq_func(fastq, 'rt') as infile, out_funct(output, 'wt') as outfile:
        n = 0
        written = False
        head = None
        for line in infile:
            n += 1
            if n % 4 == 1:
                head = line
            if n % 4 == 2:
                seq = line.rstrip()[five_prime_trim:five_prime_trim +
                                    seq_length]
                if seq:  #prevent writing empty sequences
                    outfile.write(head)
                    outfile.write(seq + '\n')
                    written = True
            else:
                if written:
                    if n % 4 == 0:
                        qual = line.rstrip()[five_prime_trim:five_prime_trim +
                                             seq_length]
                        outfile.write(qual + '\n')
                    else:
                        outfile.write(line)
                if n % 4 == 0:  #reset at end of fastq record
                    written = False
                    head = None
                    if progress_interval is not None:
                        records = int(n / 4)
                        if records % progress_interval == 0:
                            logger.info(
                                "{:,} records processed".format(records))
        records = int(n / 4)
        logger.info("Finished processing {:,} records".format(records))



trim_fastq(snakemake.input[0],
           snakemake.output[0],
           five_prime_trim=snakemake.config['trim'],
           seq_length=snakemake.config['max_length'])
