rule porechop_abi:
    input:
        a = get_raw_fastq
    output:
        a = temp("data/porechopped/{input}.fastq")
    threads:
        75
    resources:
        mem_mb=scaled_mem(1.5, 32000),
        runtime=scaled_time(0.05, 600),
    conda:
        "../envs/porechop_abi.yml"
    shell:
        """
        porechop_abi --ab_initio -i {input.a} -o {output.a}
        """