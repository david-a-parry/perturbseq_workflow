local_guides_xlsx = os.path.join("resources",
                                 os.path.basename(config["guides"]))


def get_original_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    return units.loc[(wildcards.sample_name, wildcards.unit_name), ["fq1"]].dropna()


def get_trimmed_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    fq = get_original_fastq(wildcards)
    return os.path.join("results",
                        "trimmed_fastq",
                        "{}-{}-trimmed.fastq.gz".format(wildcards.sample_name,
                                                        wildcards.unit_name))


def get_fastq(wildcards):
    """Get fastq files of given sample-unit."""
    if config.get('trim', 0) == 0 and config.get('max_length', 20) >= 20:
        return get_original_fastq(wildcards)
    return get_trimmed_fastq(wildcards)


rule download_guides:
    input:
        HTTP.remote(config["guides"], keep_local=True)
    output:
        local_guides_xlsx
    log:
        "logs/download_guides.log"
    run:
        shell("mv {input} {local_guides_xlsx} 2> {log}")

rule build_reference:
    input:
        local_guides_xlsx
    output:
        guides_fasta
    conda:
        "../envs/stats.yaml"
    log:
        "logs/build_reference.log"
    resources:
        mem_gb=4
    script:
        "../scripts/ref_from_xlsx.py"

rule bt_index_reference:
    input:
        guides_fasta
    params:
        prefix = bt_idx,
    output:
        ebwt1 = bt_idx + '.1.ebwt',
        ebwt2 = bt_idx + '.2.ebwt',
        ebwt3 = bt_idx + '.3.ebwt',
        ebwt4 = bt_idx + '.4.ebwt',
        ebwtrev1 = bt_idx + '.rev.1.ebwt',
        ebwtrev2 = bt_idx + '.rev.2.ebwt'
    log:
        "logs/bt_index_reference.log"
    conda:
        "../envs/bowtie.yaml"
    resources:
        mem_gb=16
    shell:
        "bowtie-build --threads {threads} -f {input} {params.prefix} 2>&1 >{log}"

rule faidx_index_reference:
    input:
        guides_fasta
    output:
        guides_fasta + '.fai'
    log:
        "logs/faidx_index_reference.log"
    conda:
        "../envs/bowtie.yaml"
    shell:
        "samtools faidx {input} 2>> {log}"

rule trim_fastq:
    input:
        get_original_fastq
    output:
        temp("results/trimmed_fastq/{sample_name}-{unit_name}-trimmed.fastq.gz")
    log:
        "logs/trim_fq/{sample_name}_{unit_name}.log"
    conda:
        "../envs/stats.yaml"
    resources:
        mem_gb=4
    script:
        "../scripts/trim_fq.py"

rule bowtie_map:
    input:
        fq = get_fastq,
        ebwt1 = bt_idx + '.1.ebwt',
        ebwt2 = bt_idx + '.2.ebwt',
        ebwt3 = bt_idx + '.3.ebwt',
        ebwt4 = bt_idx + '.4.ebwt',
        ebwtrev1 = bt_idx + '.rev.1.ebwt',
        ebwtrev2 = bt_idx + '.rev.2.ebwt'
    params:
        ref_idx = bt_idx
    output:
        bam=temp("results/alignments/{sample_name}-{unit_name}.bam"),
        bai=temp("results/alignments/{sample_name}-{unit_name}.bam.bai")
    log:
        "logs/bowtie_map/{sample_name}_{unit_name}.log"
    conda:
        "../envs/bowtie.yaml"
    resources:
        tmpdir="tmp"
    threads: 8
    resources:
        mem_gb=16,
        runtime_min=120
    shell:
        "(bowtie -m 1 -v 2 -p {threads} {params.ref_idx} {input.fq} -S | "
        "samtools sort -O BAM "
        "-T tmp/{wildcards.sample_name}_{wildcards.unit_name} - )"
        "> {output.bam} 2> {log} && samtools index {output.bam} 2>> {log}"

