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
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input} > {output.a}
        """


def assembly_stats_candidates(wildcards):
    return grid_assemblies(wildcards.input)

rule assembly_stats:
    input:
        candidates = assembly_stats_candidates,
        best = "data/contig/{input}_{selector}_file.fa",
    output:
        a = "results/{input}/{selector}/assembly_stats.tsv"
    threads:
        4
    resources:
        mem_mb = 16000,
        runtime = 120,
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input.candidates} {input.best} > {output.a}
        """


rule final_genome:
    """The deliverable: polished where possible, best raw assembly otherwise --
    one per selector (lowest_contig / highest_busco, see Snakefile SELECTORS)."""
    input:
        a = lambda w: (f"data/dorado_polish/{w.input}_{w.selector}.fasta" if has_bam(w.input)
                       else f"data/contig/{w.input}_{w.selector}_file.fa")
    output:
        a = protected("results/{input}/{selector}/{input}_{selector}_final.fasta")
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = 30,
    shell:
        """
        cp {input.a} {output.a}
        """


def summary_optional_inputs(wildcards):
    """Inputs for summary.txt's optional sections: BUSCO before polishing
    (pod5/bam samples), TELOMERES (telomere: enabled) and CONTAMINATION
    (contamination: enabled)."""
    d = {}
    if has_bam(wildcards.input):
        d["busco_unpolished"] = f"data/busco_unpolished/{wildcards.input}_{wildcards.selector}/BUSCO"
    if TELOMERE_ON:
        d["telomeres"] = f"results/{wildcards.input}/{wildcards.selector}/telomeres.tsv"
        d["motif"] = f"data/telomere/{wildcards.input}/motif.tsv"
    if CONTAMINATION_ON:
        d["tiara"] = f"results/{wildcards.input}/{wildcards.selector}/contamination/tiara.tsv"
    return d

rule summary:
    input:
        unpack(summary_optional_inputs),
        reads = "results/{input}/read_stats.tsv",
        asm   = "results/{input}/{selector}/assembly_stats.tsv",
        final = "results/{input}/{selector}/{input}_{selector}_final.fasta",
        busco = "data/busco/{input}_{selector}/BUSCO",
    output:
        a = "results/{input}/{selector}/summary.txt"
    params:
        sample   = lambda w: w.input,
        selector = lambda w: w.selector,
        polished = lambda w: has_bam(w.input),
        lineage  = lambda w: busco_lineage(w.input),
        tiara_min_len = CONTAM.get("min_len", 3000),
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = 30,
    script:
        "../scripts/summarize.py"