# Only runs for samples with type: pod5. GPU-only (bio-node10); output BAM
# carries the move table dorado_align/dorado_polish need later, and is
# protected() since basecalling is the most expensive step to redo.
rule dorado_basecall:
    input:
        dorado = dorado_bin,
        a = lambda w: SAMPLES[w.input]["path"]
    output:
        a = protected("data/dorado_basecall/{input}.bam")
    threads:
        16
    resources:
        mem_mb=resources["dorado_basecall"]["mem_mb"],
        runtime=resources["dorado_basecall"]["runtime"],
        gres=GPU_GRES,
    log:
        "logs/dorado_basecall/{input}.log"
    shell:
        """
        "{input.dorado}" basecaller sup --emit-moves --device cuda:all {input.a} > {output.a} 2> {log}
        """