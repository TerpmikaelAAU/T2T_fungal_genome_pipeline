configfile: "config/config.yaml"


rule dorado:
    input:
        a = rules.rasusa.output.a,
    output:
        a = "data/dorado/{input}.fasta"
    threads:
        16
    resources:
        mem_mb=resources["dorado"]["mem_mb"],
        runtime=resources["dorado"]["runtime"],
        gres=GPU_GRES,
    shell:
        """
        "{config[dorado]}" correct --index-size 2G --device cuda:all {input.a} > {output.a}
        
        """
