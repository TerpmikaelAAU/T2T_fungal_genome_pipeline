configfile: "config/config.yaml"


rule chopper_flye:
    input:
        a = rules.porechop_abi.output.a
    output:
        a = "data/chopper/Flye/{input}.fastq"
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["time"]
    conda:
        "chopper.yml"
    shell:
        """
        chopper -q 20 -l 20000 --threads $(nproc) < {input.a} > {output.a}
        """

