# ============================================================================
#  FASTQ samples: decompress once, drop duplex parents
# ============================================================================
# The filter grid means every (min_q, min_len) chopper job used to zcat the
# same multi-GB .fastq.gz independently -- measured at 12.8% CPU efficiency,
# 98.9% memory on a real run: every chopper instance was starved waiting on
# single-threaded gzip, not doing filtering work. This rule decompresses it
# exactly once; get_raw_fastq() points every grid point at the plain-text
# result instead of the original .gz.
#
# It also removes the simplex parents of duplex reads (dorado duplex output),
# so each duplex molecule enters the assembly once -- see
# workflow/scripts/drop_duplex_parents.sh. A no-op copy for simplex data.
# BAM samples get the same filter (by the dx tag) in bam_to_fastq.
rule prepare_fastq:
    input:
        a = lambda w: SAMPLES[w.input]["path"]
    output:
        a = temp("data/reads/{input}.fastq")
    threads:
        8
    resources:
        mem_mb=scaled_mem(0.05, 4000),
        runtime=scaled_time(0.1, 120),
    params:
        script = os.path.join(workflow.basedir, "workflow/scripts/drop_duplex_parents.sh"),
    log:
        "logs/prepare_fastq/{input}.log"
    shell:
        """
        bash "{params.script}" "{input.a}" "{output.a}" {threads} > {log} 2>&1
        """
