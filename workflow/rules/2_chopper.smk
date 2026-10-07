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
# per-grid-point zcat of the same multi-GB .gz file. Cells stricter than
# config `prefilter` read the prefilter cell's output, not the full read set
# (see chopper_source in the Snakefile).
rule chopper:
    input:
        a = get_chopper_input
    output:
        a = temp("data/chopper/{input}_q{minq}_l{minlen}.fastq")
    # chopper is I/O bound: ~0.6 cores busy whatever --threads is (5% CPU
    # efficiency on 12 threads). Its real memory is ~1-2 GB; the efficiency
    # report's "99% memory" is page cache from streaming the FASTQ, which
    # fills any limit (it read 99% at 8 GB and again at 16 GB, never OOM).
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
        chopper -q {wildcards.minq} -l {wildcards.minlen} --threads {threads} \
            < {input.a} > {output.a}
        """
