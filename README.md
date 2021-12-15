# Perturbseq Workflow

[Snakemake](snakemake.github.io) workflow for generating read count data from
CRISPR screens.

## Setup

### Snakemake

You will need to initialise a conda environment with Snakemake and pandas
available to run this workflow, ideally using [mamba](https://github.com/mamba-org/mamba):

    mamba create --name snakemake  'snakemake-minimal>=5.24.1' 'pandas>=1.1'

The necessary environements for the different stages of this workflw  will then
be created once you run the workflow (see below).

### Your Data

Edit the config files `config/config.yaml`, `config/samples.tsv` and
`config/units.tsv` according to your needs.

#### config/config.yaml

* samples:
path to your samples table. Does not need to be edited if you intend to edit
config/samples.tsv file instead.

* units:
path to your units table. Does not need to be edited if you intend to edit
config/units.tsv file instead.

* guides:
URL (omitting the http:// prefix) to excel spreadsheet of guide sequences.
Currently those in the format of the Toronto knockout libraries are supported
(e.g. [TKOv3](media.addgene.org/cms/filer_public/71/a8/71a81179-7a62-4d75-9b53-236e6f6b7d4d/tkov3_guide_sequence.xlsx)
and [mTKO](media.addgene.org/cms/filer_public/1a/c4/1ac4f468-fc05-4c49-9d36-b61ec18ed759/mtko_library.xlsx)
).

* trim:
Optional parameter to specify the number of nucleotides to trim from the start
of each read in your FASTQ files.

* max_length:
Optional parameter to specify the length of sequence to keep from your FASTQ
files (i.e. if you need to trim from the 3' ends of your sequences).


#### config/samples.tsv

Edit the example file to specify the individual samples and timepoints in your
project.

#### config/units.tsv

Edit the example file to specify the individual run files in your project. A
sample may have multiple FASTQ files associated with it distinguished by the
'unit_name' parameter (which may be any arbitrary field as long as it does not
contain a hyphen).

## Running

Invoke snakemake with a number of cores that suits your hardware. The 
`--use-conda` flag is required as cluster execution is not yet supported.

`snakemake --use-conda --cores 4`

## Author

Written by David A. Parry at the University of Edinburgh.
