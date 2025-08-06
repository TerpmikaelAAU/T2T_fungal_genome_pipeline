configfile: "config/config.yaml"
rule dorado_align:
    input:
        a = rules.dorado_align.output.a,
    output: 
        a = "data/dorado_align/{input}_indexed.bam
    threads:
        12
    resources:
        mem_mb=resources["porechop_api"]["mem_mb"],
        runtime=resources["porechop_api"]["time"]
        partition="general"
    conda:
        "samtools.yml"
    shell:
        """
        samtools index --threads $(nproc) {input.a} > {output.a}

        """

