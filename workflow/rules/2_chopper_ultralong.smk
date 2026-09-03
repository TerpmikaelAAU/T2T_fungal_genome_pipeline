rule chopper_ultralong:
    input:
        a = rules.porechop_abi.output.a
    output:
        a = temp("data/chopper/ultralong/{input}.fastq")
    threads:
        12
    resources:
        mem_mb=scaled_mem(0.15, 8000),
        runtime=resources["chopper"]["runtime"]
    conda:
        "../envs/chopper.yml"
    shell:
        """
        chopper -q 10 -l 50000 --threads $(nproc) < {input.a} > {output.a}
        """