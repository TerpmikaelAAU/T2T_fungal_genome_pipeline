# The assembly grid's min_len sweep, applied to dorado-corrected reads.
# Correction discards real per-base quality (see 3_seqtk_fasta_to_fastq.smk),
# so only min_len is swept here -- min_q plays no further part once
# correction is on for a sample (see filter_grid() in Snakefile).
def get_corrected_fastq(wildcards):
    return f"data/seqtk/fasta_to_fastq/{wildcards.input}.fastq"

rule length_filter_corrected:
    input:
        a = get_corrected_fastq
    output:
        a = temp("data/dorado_filtered/{input}_l{minlen}.fastq")
    threads:
        12
    resources:
        mem_mb=scaled_mem(0.15, 16000),
        runtime=resources["chopper"]["runtime"]
    container:
        "docker://quay.io/biocontainers/chopper:0.13.0--h7f49ad2_0"
    conda:
        "../envs/chopper.yml"
    shell:
        """
        chopper -l {wildcards.minlen} --threads {threads} \
            < {input.a} > {output.a}
        """
