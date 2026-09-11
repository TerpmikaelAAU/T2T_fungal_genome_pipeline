rule chopper_flye:
    input:
        a = get_trimmed_fastq
    output:
        a = temp("data/chopper/Flye/{input}.fastq")
    threads:
        12
    resources:
        mem_mb=scaled_mem(0.15, 8000),
        runtime=resources["chopper"]["runtime"]
    conda:
        "../envs/chopper.yml"
    shell:
        """
        if [[ "{input.a}" == *.gz ]]; then zcat "{input.a}"; else cat "{input.a}"; fi \
            | chopper -q 20 -l 20000 --threads {threads} > {output.a}
        """