configfile: "config/config.yaml"

rule rasusa:
    input:
        a = rules.chopper_dorado.output.a,
    output:
        a = "data/rasusa/Coverage/{input}.fastq"
    threads:
        12  
    resources:
        mem_mb=resources["rasusa"]["mem_mb"],
        runtime=resources["rasusa"]["time"]
    conda:
        "rasusa.yml"
    shell:
        """
        rasusa reads --coverage 100 --genome-size 40mb --seed 1 -o {output.a} {input.a}
        """