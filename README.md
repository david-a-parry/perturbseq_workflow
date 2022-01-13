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
`config/units.tsv` according to your needs. See schemas in `workflow/schemas`.

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

* control_genotype

control genotype for synthetic lethal/suppressor gene identification. If provided, a regression analysis of other genotypes will be performed comparing BAGEL and JACKS scores against this genotype. Assumes matching timepoints for each genotype.

* bagel_zip:

description: URL of BAGEL github repository.

* trim:

Optional parameter to specify the number of nucleotides to trim from the start
of each read in your FASTQ files.

* max_length:

Optional parameter to specify the length of sequence to keep from your FASTQ
files (i.e. if you need to trim from the 3' ends of your sequences).

* species:

Used to automatically determine which core essential and negative control gene sets to use. Either human or mouse are supported. Ignored if bagel_ess_genes and bagel_neg_genes are specified. Default = human.

* bagel_ess_genes:

Path to text file containing core essential genes relevant to your study. File must contain gene symbols as the first (or only) whitespace separated column.

* bagel_neg_genes:

Path to text file containing negative control non-essential genes relevant to your study. File must contain gene symbols as the first (or only) whitespace separated column.


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
`--use-conda` flag will create all the necessary conda environments for the
workflow.

    $ snakemake --use-conda --cores 4

To run on an SGE cluster, example configurations are provided in `cluster-qsub`
and `cluster_config.yaml`. The workflow can be run on an SGE cluster with a
command like the following:

    $snakemake --profile cluster-qsub --cluster-config cluster_config.yaml --use-conda --cores 8

## Author

Written by David A. Parry at the University of Edinburgh.
