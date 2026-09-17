# ============================================================================
#  Statistics and final deliverables
# ============================================================================
# These rules CONSUME the temp() intermediates, so Snakemake keeps each file
# alive until its stats have been collected, then deletes it. That is why the
# read/assembly stats rules list files that are otherwise transient.

def read_stage_files(wildcards):
    """Every read file worth measuring, for this sample's entry point.

    Filtering/correction now fork into a (min_q, min_len) grid (see
    config.yaml `filter:`), so there is no longer a single "the chopped
    file" or "the corrected file" -- those per-grid-cell stats live in
    assembly_stats.tsv instead, alongside each candidate's contig stats."""
    n = wildcards.input
    f = {"raw": get_raw_fastq(wildcards)}
    if trim_adapters(n):
        f["porechopped"] = f"data/porechopped/{n}.fastq"
    return f


rule read_stats:
    input:
        unpack(read_stage_files)
    output:
        a = "results/{input}/read_stats.tsv"
    threads:
        4
    resources:
        mem_mb = scaled_mem(0.1, 8000),
        runtime = 120,
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input} > {output.a}
        """


def assembly_stats_candidates(wildcards):
    grid = filter_grid(wildcards.input)
    return expand("data/hifiasm/{input}_q{minq}_l{minlen}/{input}_q{minq}_l{minlen}.fa",
                  input=wildcards.input,
                  minq=grid["min_q"],
                  minlen=grid["min_len"])

rule assembly_stats:
    input:
        candidates = assembly_stats_candidates,
        best = "data/contig/{input}_lowest_contig_file.fa",
    output:
        a = "results/{input}/assembly_stats.tsv"
    threads:
        4
    resources:
        mem_mb = 16000,
        runtime = 120,
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input.candidates} {input.best} > {output.a}
        """


rule final_genome:
    """The deliverable: polished where possible, best raw assembly otherwise."""
    input:
        a = lambda w: (f"data/dorado_polish/{w.input}.fasta" if has_bam(w.input)
                       else f"data/contig/{w.input}_lowest_contig_file.fa")
    output:
        a = protected("results/{input}/{input}_final.fasta")
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = 30,
    shell:
        """
        cp {input.a} {output.a}
        """


rule summary:
    input:
        reads = "results/{input}/read_stats.tsv",
        asm   = "results/{input}/assembly_stats.tsv",
        final = "results/{input}/{input}_final.fasta",
        busco = "data/busco/{input}/BUSCO",
    output:
        a = "results/{input}/summary.txt"
    params:
        sample   = lambda w: w.input,
        polished = lambda w: has_bam(w.input),
        lineage  = lambda w: busco_lineage(w.input),
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = 30,
    script:
        "../scripts/summarize.py"