rule assembly:
    input:
        clean_paired_1 = rules.fastp.output.paired_1,
        clean_paired_2 = rules.fastp.output.paired_2,
        clean_unpaired_1 = rules.fastp.output.unpaired_1,
        clean_unpaired_2 = rules.fastp.output.unpaired_2
    output:
        scaffolds="results/{fraction}/{sample}/assembly/scaffolds.fasta"
    log:
        "logs/{fraction}/{sample}/assembly/{sample}.assembly.log"
    threads:
        config["SPADES"]["threads"]
    params:
        k = config["SPADES"]["k"],
        assembly_dir = "results/{fraction}/{sample}/assembly"
    conda:
        "../envs/assembly.yaml"
    shell:
        "spades.py "
        "-t {threads} "
        "--meta "
        "-k {params.k} "
        "--pe1-1 {input.clean_paired_1} "
        "--pe1-2 {input.clean_paired_2} "
        "--pe1-s {input.clean_unpaired_1} "
        "--pe1-s {input.clean_unpaired_2} "
        "-o {params.assembly_dir} "
        "&>{log}"

rule scaffolds_header_fix:
    input:
        scaffolds = rules.assembly.output.scaffolds
    output:
        scaffolds_header_fixed = "results/{fraction}/{sample}/scaffolds/{sample}_scaffolds.fasta"
    log:
        "logs/{fraction}/{sample}/hearder_fix/{sample}.header_fix.log"
    params:
        prefix = "{sample}_"
    shell:
        "perl -p -e 's/^>/>{params.prefix}/g' {input} > {output} "
        "2> {log}"
