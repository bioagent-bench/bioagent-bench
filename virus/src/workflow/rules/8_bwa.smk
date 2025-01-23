rule bwa_index:
    input:
        scaffolds = rules.length_filter.output.scaffolds_length_filtered,
    output:
        idx = multiext(rules.length_filter.output.scaffolds_length_filtered, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/{fraction}/{sample}/bwa_index/{sample}.log",
    params:
        algorithm="bwtsw",
    threads:
        config["BWA"]["threads"]
    wrapper:
        "v1.7.1/bio/bwa/index"

rule bwa_mem:
    input:
        reads = get_sample_fastqs,
        idx = multiext(rules.length_filter.output.scaffolds_length_filtered, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = "results/{fraction}/{sample}/bwa/{sample}.bam",
    log:
        "logs/{fraction}/{sample}/bwa_mem/{sample}.log",
    params:
        extra=r"-R '@RG\tID:{sample}\tSM:{sample}'",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
    threads:
        config["BWA"]["threads"]
    wrapper:
        "v1.10.0/bio/bwa/mem"

rule samtools_stats:
    input:
        bam = rules.bwa_mem.output.bam,
    output:
        samtools_stats = "results/{fraction}/{sample}/bwa/{sample}.samtools_stats.txt",
    params:
        extra="",  # Optional: extra arguments.
    log:
        "logs/{fraction}/{sample}/samtools_stats/{sample}.log",
    wrapper:
        "v1.7.1/bio/samtools/stats"
