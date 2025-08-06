configfile: "config/config.yaml"


rule porechop_abi:
    input:
        a = rules.bam_to_fastq.output.a
    output:
        a = "data/porechopped/{input}.fastq"
    threads:
        75
    resources:
        mem_mb=resources["porechop_api"]["mem_mb"],
        runtime=resources["porechop_api"]["time"]
    conda:
        "porechop_abi.yml"
    shell:
        """
        porechop_abi --ab_initio -i {input} -o {output.a}
        """
