configfile: "config/config.yaml"

rule chopper_ultralong:
    input:
        a = rules.porechop_abi.output.a
    output:
        a = "data/chopper/ultralong/{input}.fastq"
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["time"]
    conda:
        "chopper.yml"
    shell:
        """
        chopper -q 10 -l 50000 --threads $(nproc) < {input.a} > {output.a}
        """

