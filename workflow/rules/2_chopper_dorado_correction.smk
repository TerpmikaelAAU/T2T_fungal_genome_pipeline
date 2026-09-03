rule chopper_dorado:
    input:
        a = rules.porechop_abi.output.a
    output:
        a = temp("data/chopper/L10kbQ10/{input}.fastq")
    threads:
        12
    resources:
        mem_mb=scaled_mem(0.15, 8000),
        runtime=resources["chopper"]["runtime"]
    conda:
        "../envs/chopper.yml"
    shell:
        """
        chopper -q 10 -l 10000 --threads $(nproc) < {input.a} > {output.a}
        """