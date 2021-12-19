def all_read_counts(wildcards):
    return ["results/read_counts/per_bam/{}-{}_counts.csv.gz".format(sample,
                                                                     unit)
            for (sample, unit) in units.index]


rule count_reads:
    input:
        "alignments/{sample_name}-{unit_name}.bam",
        guides_fasta,
        guides_fasta + '.fai'
    output:
        read_counts = "results/read_counts/per_bam/{sample_name}-{unit_name}_counts.csv.gz",
        map_counts = "results/read_counts/per_bam/{sample_name}-{unit_name}_mapped.csv.gz"
    log:
        "logs/read_counts/{sample_name}_{unit_name}.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/read_counts_from_bam.py"

rule combine_read_counts:
    input:
        all_read_counts
    output:
        "results/read_counts/read_counts_per_alignment.csv.gz",
        "results/read_counts/read_counts_per_sample.csv.gz",
        "results/read_counts/all_read_counts.tsv",
        expand("results/read_counts/{genotype}_read_counts.tsv", genotype=genotypes)
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/combine_read_counts.py"

rule plot_guide_counts:
    input:
        "results/read_counts/read_counts_per_sample.csv.gz",
    output:
        "results/plots/lorenz_curves_t0.pdf",
        "results/plots/guide_counts_t0.pdf",
        "results/plots/guide_counts_all.pdf",
        "results/read_counts/norm_counts.tsv.gz",
        "results/read_counts/fold_change.tsv"
    conda:
        "../envs/stats.yaml"
    notebook:
        "../notebooks/guide_plots.py.ipynb"

