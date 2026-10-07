# Separate, fixed Q20/l20000 filter feeding the organelle-recovery path
# (flye -> getorganelle), independent of the main (min_q, min_len) grid in
# config.yaml `filter:` -- that grid is tuned for hifiasm's nuclear assembly,
# not for finding one mitochondrial contig.
# Reads the prefilter cell's output when Q20/l20000 is stricter than it.
rule chopper_flye:
    input:
        a = lambda w: chopper_source(w, 20, 20000)
    output:
        a = temp("data/chopper/Flye/{input}.fastq")
    # Same sizing as rule chopper (2_chopper.smk).
    threads:
        2
    resources:
        mem_mb=scaled_mem(0, 4000),
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