# ============================================================================
#  Statistics and final deliverables
# ============================================================================
# These rules CONSUME the temp() intermediates, so Snakemake keeps each file
# alive until its stats have been collected, then deletes it. That is why the
# read/assembly stats rules list files that are otherwise transient.

def read_stage_files(wildcards):
    """Every read file worth measuring, for this sample's entry point."""
    n = wildcards.input
    f = {"porechopped": f"data/porechopped/{n}.fastq",
         "chopped_L10kbQ10": f"data/chopper/L10kbQ10/{n}.fastq",
         "corrected": f"data/dorado/{n}.fasta"}
    if subsampled(n):
        f["subsampled"] = f"data/rasusa/Coverage/{n}.fastq"
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


rule assembly_stats:
    input:
        candidates = expand("data/hifiasm/{{input}}_{length}/{{input}}_{length}.fa",
                            length=config["length"]),
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