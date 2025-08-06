configfile: "config/config.yaml"

rule chopper_dorado:
    input:
        a = rules.porechop_abi.output.a
    output:
        a = "data/chopper/L10kbQ10/{input}.fastq"
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["time"]
    conda:
        "chopper.yml"
    shell:
        """
        chopper -q 10 -l 10000 --threads $(nproc) < {input.a} > {output.a}
        """
