# Adapter trimming, opt-out per sample (config `porechop: false`, default on
# -- see trim_adapters()). --ab_initio auto-detects adapters instead of using
# a fixed database. CAUTION: OOM'd at 95 GB on p_infestans_88069, because the
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