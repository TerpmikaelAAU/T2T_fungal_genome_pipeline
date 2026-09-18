# ============================================================================
#  Read-filtering grid
# ============================================================================
# Applies BOTH cutoffs in one pass -- one rule instance per (min_q, min_len)
# cell of the grid in config.yaml `filter:`. Downstream (assembly, or dorado
# correct if enabled) runs independently for every cell; contig_count picks
# the winner. Also reused for hifiasm's --ul input (config `ultralong:`),
# which is just another (fixed) point on the same (minq, minlen) plane.
#
# input.a is always plain text now (see 0_3_decompress.smk) -- no more
# per-grid-point zcat of the same multi-GB .gz file.
rule chopper:
    input:
        a = get_trimmed_fastq
    output:
        a = temp("data/chopper/{input}_q{minq}_l{minlen}.fastq")
    threads:
        12
    resources:
        # Floor raised from 8000: a real run's efficiency report showed
        # chopper sitting at 98.9% memory utilization, one bad read away
        # from an OOM.
        mem_mb=scaled_mem(0.15, 16000),
        runtime=resources["chopper"]["runtime"]
    container:
        "docker://quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0"
    conda:
        "../envs/chopper.yml"
    shell:
        """
        chopper -q {wildcards.minq} -l {wildcards.minlen} --threads {threads} \
            < {input.a} > {output.a}
        """
