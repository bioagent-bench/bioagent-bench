LENGTH = str(config['SEQTK']['length'])

rule cat_contigs:
    input:
        contigs = "results/{fraction}/concatenated_scaffolds/{fraction}_scaffolds_gt"+ LENGTH + ".fasta"
    output:
        cat_contig_classification = "results/{fraction}/cat/{fraction}_scaffolds_gt" + LENGTH + ".contig2classification.txt"
    threads:
        config['CAT']['threads']
    params:
        database = config['CAT']['database'],
        taxonomy = config['CAT']['taxonomy'],
        out_prefix = "results/{fraction}/cat/{fraction}_scaffolds_gt" + LENGTH,
        path_to_diamond = config['CAT']['path_to_diamond']
    conda:
        "../envs/cat.yaml"
    shell:
        "CAT contigs "
        "-n {threads} "
        "-c {input} "
        "-d {params.database} "
        "-t {params.taxonomy} "
        "-o {params.out_prefix} "
        "--path_to_diamond {params.path_to_diamond}"
        
rule cat_add_names:
    input: 
        cat_contig_classification = rules.cat_contigs.output.cat_contig_classification,
    output:
         cat_add_official_names = "results/{fraction}/cat/{fraction}_scaffolds_gt" + LENGTH + "_cat_taxonomy_official.txt"
    params:
        taxonomy = config['CAT']['taxonomy'],
    conda:
        "../envs/cat.yaml"
    shell:
        "CAT add_names "
        "-i {input} "
        "-t {params.taxonomy} "
        "-o {output} "
        "--only_official"
        
rule cat_summarise:
    input: 
        cat_add_official_names = rules.cat_add_names.output.cat_add_official_names,
        contigs = "results/{fraction}/concatenated_scaffolds/{fraction}_scaffolds_gt"+ LENGTH + ".fasta"
    output:
         cat_summarise_output = "results/{fraction}/cat/{fraction}_scaffolds_gt" + LENGTH + "_cat_summarise.txt"
    conda:
        "../envs/cat.yaml"
    shell:
        "CAT summarise "
        "-i {input.cat_add_official_names} "
        "-c {input.contigs} "
        "-o {output}"
