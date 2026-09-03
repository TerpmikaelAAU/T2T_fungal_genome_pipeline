configfile: "config/config.yaml"
rule dorado_polish:
    input:
        a = rules.dorado_align.output.a,
        b = rules.contig_count.output.a,
    output:
        a = "data/dorado_polish/{input}.fasta"
    threads:
        16
    resources:
        mem_mb=resources["porechop_api"]["mem_mb"],
        runtime=resources["porechop_api"]["time"],
        partition="gpu",
        gpus=1,
    shell:
        """
        "{config[dorado]}" polish --batchsize 8 --device cuda:all {input.a} {input.b} > {output.a}
        
        """

