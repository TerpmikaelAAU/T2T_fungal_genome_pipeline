# ============================================================================
#  Decompress once
# ============================================================================
# The filter grid means every (min_q, min_len) chopper job used to zcat the
# same multi-GB .fastq.gz independently -- 12.8% CPU efficiency, starved by
# single-threaded gzip (see Handoff.md). This rule decompresses it exactly
# once; get_raw_fastq() points every grid point at the plain-text result
# instead of the original .gz.
rule decompress_once:
    input:
        a = lambda w: SAMPLES[w.input]["path"]
    output:
        a = temp("data/decompressed/{input}.fastq")
    threads:
        8
    resources:
        mem_mb=scaled_mem(0.05, 4000),
        runtime=scaled_time(0.05, 120),
    shell:
        """
        if command -v pigz >/dev/null 2>&1; then
            pigz -dc -p {threads} "{input.a}" > "{output.a}"
        else
            zcat "{input.a}" > "{output.a}"
        fi
        """
