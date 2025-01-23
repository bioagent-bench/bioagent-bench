LENGTH =str(config['SEQTK']['length'])

rule length_filter:
    input:
        scaffolds_header_fixed = rules.scaffolds_header_fix.output.scaffolds_header_fixed
    output:
        scaffolds_length_filtered = "results/{fraction}/{sample}/scaffolds/{sample}_scaffolds_gt" + LENGTH + ".fasta"
    params:
        length = LENGTH
    conda:
        "../envs/seqtk.yaml"
    shell:
        "seqtk seq -L {params.length} {input.scaffolds_header_fixed} > {output.scaffolds_length_filtered}"

rule scaffolds_concatenate:
    input:
        get_scaffolds_by_fraction
    output:
        scaffolds_concatenated = "results/{fraction}/concatenated_scaffolds/{fraction}_scaffolds_gt"+ LENGTH + ".fasta"
    shell:
        "cat {input} > {output}"
