# Separate, fixed Q20/l20000 filter feeding the organelle-recovery path
# (flye -> getorganelle), independent of the main (min_q, min_len) grid in
# config.yaml `filter:` -- that grid is tuned for hifiasm's nuclear assembly,
# not for finding one mitochondrial contig.
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
    container:
        "docker://quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0"
    conda:
        "../envs/chopper.yml"
    shell:
        """
        if [[ "{input.a}" == *.gz ]]; then zcat "{input.a}"; else cat "{input.a}"; fi \
            | chopper -q 20 -l 20000 --threads {threads} > {output.a}
        """