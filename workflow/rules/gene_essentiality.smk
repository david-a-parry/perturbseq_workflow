def get_t0_sample(wildcards):
    smpls = samples[(samples.genotype == wildcards.genotype) &
                    (samples.timepoint == t0)].sample_name.values
    if len(smpls) != 1:
        raise ValueError("Expected exactly 1 T0 sample for genotype " +
                         "{}, got {}".format(wildcards.genotype,
                                             ", ".join(smpls) if smpls else 0))
    return smpls[0]


def get_bagel_samples(wildcards):
    return ','.join(
        samples[(samples.genotype == wildcards.genotype) &
                (samples.timepoint == int(wildcards.timepoint))
        ].sample_name.unique())


def get_bagel_ref_genes():
    ess = config.get("bagel_ess_genes")
    neg = config.get("bagel_neg_genes")
    if ess is None or neg is None:
        species = config.get("species")
        if species is None or species.lower() == 'human':
            ess = ess if ess is not None else "resources/bagel/CEGv2.txt"
            neg = neg if neg is not None else "resources/bagel/NEGv1.txt"
        elif species.lower() == 'mouse':
            ess = ess if ess is not None else "resources/bagel/CEG_mouse.txt"
            neg = neg if neg is not None else "resources/bagel/NEG_mouse.txt"
        else:
            raise ValueError("Unsupported species: {}".format(species))
        return ess, neg


bagel_ess_genes, bagel_neg_genes = get_bagel_ref_genes()


rule download_bagel:
    input:
        HTTP.remote(config["bagel_zip"])
    output:
        "resources/bagel/BAGEL.py"
    log:
        "logs/download_bagel.log"
    run:
        shell("unzip -o {config[bagel_zip]} -d resources 2>&1 >{log}")
        shell("mkdir -p resources/bagel 2>> {log}")
        shell("mv resources/bagel-*/* resources/bagel 2>> {log}")
        shell("rmdir resources/bagel-* 2>> {log}")

rule run_bagel_fc:
    input:
        "results/read_counts/{genotype}_read_counts.tsv",
        "resources/bagel/BAGEL.py"
    output:
        "results/bagel/{genotype}.foldchange"
    log:
        "logs/bagel-fc-{genotype}.log"
    conda:
        "../envs/stats.yaml"
    params:
        control_sample = get_t0_sample,
        prefix = "results/bagel/{genotype}"
    shell:
        "python resources/bagel/BAGEL.py fc -i {input[0]} -o {params.prefix} "
        "-c {params.control_sample} 2> {log}"


rule run_bagel_bf:
    input:
        "results/bagel/{genotype}.foldchange"
    output:
        "results/bagel/{genotype}-{timepoint}.bf"
    log:
        "logs/bagel-bf-{genotype}-{timepoint}.log"
    conda:
        "../envs/stats.yaml"
    params:
        sample_cols = get_bagel_samples,
        neg = bagel_neg_genes,
        ess = bagel_ess_genes
    shell:
        "python resources/bagel/BAGEL.py bf -i {input[0]} -o {output} "
        "-n {params.neg} -e {params.ess} -c {params.sample_cols} 2>&1 >{log}"

rule run_bagel_pr:
    input:
        "results/bagel/{genotype}-{timepoint}.bf"
    output:
        "results/bagel/{genotype}-{timepoint}.pr"
    log:
        "logs/bagel-pr-{genotype}-{timepoint}.log"
    conda:
        "../envs/stats.yaml"
    params:
        neg = bagel_neg_genes,
        ess = bagel_ess_genes
    shell:
        "python resources/bagel/BAGEL.py pr -i {input} -o {output} "
        "-n {params.neg} -e {params.ess} 2> {log}"

rule plot_bagel_results:
    input:
        expand("results/bagel/{gt}-{time}.pr", gt=genotypes, time=timepoints)
    output:
        "results/plots/bagel_pr_curve.pdf",
        "results/plots/bagel_ess_neg_dist.pdf"
    conda:
        "../envs/stats.yaml"
    log:
        "logs/plot_bagel_results.log"
    params:
        neg = bagel_neg_genes,
        ess = bagel_ess_genes
    notebook:
        "../notebooks/plot_bagel_results.py.ipynb"

rule prepare_jacks_maps:
    input:
        "results/read_counts/norm_counts.tsv.gz"
    output:
        replicate_map = "results/jacks_analysis/replicate_map.txt",
        sgrna_map = "results/jacks_analysis/sgrna_map.txt"
    log:
        "logs/prepare_jacks_maps.log"
    conda:
        "../envs/stats.yaml"
    script:
        "../scripts/prepare_jacks_maps.py"

rule run_jacks:
    input:
        read_counts = "results/read_counts/all_read_counts.tsv",
        replicate_map = "results/jacks_analysis/replicate_map.txt",
        sgrna_map = "results/jacks_analysis/sgrna_map.txt",
        bagel = "resources/bagel/BAGEL.py"  # run after bagel install in case we need NEG genes
    output:
        "results/jacks_analysis/jacks_gene_JACKS_results.txt"
    conda:
        "../envs/jacks.yaml"
    log:
        "logs/run_jacks.log"
    params:
        neg = bagel_neg_genes
    shell:
        "JACKS.py {input.read_counts} "
        "{input.replicate_map} {input.sgrna_map} "
        "--ctrl_sample_hdr=Control --gene_hdr=Gene --sgrna_hdr Target "
        "--outprefix=results/jacks_analysis/jacks "
        "--ctrl_genes={params.neg} 2> {log}"

rule plot_jacks_results:
    input:
        "results/jacks_analysis/jacks_gene_JACKS_results.txt",
    output:
        "results/plots/jacks_pr_curve.pdf",
        "results/plots/jacks_ess_neg_dist.pdf"
    conda:
        "../envs/stats.yaml"
    log:
        "logs/plot_jacks_results.log"
    params:
        neg = bagel_neg_genes,
        ess = bagel_ess_genes
    notebook:
        "../notebooks/plot_jacks_results.py.ipynb"

rule regression_analysis:
    input:
        "results/jacks_analysis/jacks_gene_JACKS_results.txt",
        expand("results/bagel/{gt}-{time}.pr", gt=genotypes, time=timepoints)
    output:
        "results/regression_analyses/regression_hits.xlsx"
    conda:
        "../envs/stats.yaml"
    params:
        neg = bagel_neg_genes,
        ess = bagel_ess_genes
    notebook:
        "../notebooks/regression_analysis.py.ipynb"
