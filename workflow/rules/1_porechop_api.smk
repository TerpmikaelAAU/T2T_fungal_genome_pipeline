# Adapter trimming, opt-out per sample (config `porechop: false`, default on
# -- see trim_adapters()). --ab_initio auto-detects adapters instead of using
# a fixed database. CAUTION: OOM'd at 95 GB on a large test sample, because the
# ab-initio phase feeds the WHOLE file to its internal approx_counter even
# though it only samples 40k reads from it. If re-enabling for a similarly
# large sample, use a two-phase pattern instead of this single-shot one:
# `porechop_abi -go` (guess adapters only) on a small subset of reads, then
# `porechop_abi -cap adapters.txt -ddb` on the full file using that guess --
# plus `-tmp` pointed at node-local scratch (it defaults to `./tmp/`, which
# on this cluster is network storage) and an explicit `--threads`.
rule porechop_abi:
    input:
        a = get_raw_fastq
    output:
        a = temp("data/porechopped/{input}.fastq")
    # Mostly single-threaded (0.45% CPU efficiency on 64 threads without an
    # explicit --threads); only the adapter alignment uses --threads. Memory:
    # on a ~100 GB input, 1.5x (150 GB) was OOM-killed and the 3x retry
    # finished, so start at 2.5x (5x on retry).
    threads:
        16
    resources:
        mem_mb=scaled_mem(2.5, 32000),
        runtime=scaled_time(0.05, 600),
    container:
        "docker://quay.io/biocontainers/porechop_abi:0.5.0--py312h5e9d817_5"
    conda:
        "../envs/porechop_abi.yml"
    shell:
        """
        porechop_abi --ab_initio --threads {threads} -i {input.a} -o {output.a}
        """