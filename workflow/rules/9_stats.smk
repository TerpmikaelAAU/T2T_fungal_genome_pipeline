# ============================================================================
#  Statistics and final deliverables
# ============================================================================
# These rules CONSUME the temp() intermediates, so Snakemake keeps each file
# alive until its stats have been collected, then deletes it. That is why the
# read/assembly stats rules list files that are otherwise transient.

def read_stage_files(wildcards):
    """The reads before the filtering grid: the input reads, and the
    adapter-trimmed reads when porechop is on for this sample."""
    n = wildcards.input
    f = [get_raw_fastq(wildcards)]
    if trim_adapters(n):
        f.append(f"data/porechopped/{n}.fastq")
    return f


def filtered_read_stats_files(wildcards):
    """One stats row per grid cell, in the same order as grid_assemblies()."""
    n = wildcards.input
    grid = filter_grid(n)
    return [f"data/read_stats/{n}_q{q}_l{l}.tsv"
            for q in grid["min_q"] for l in grid["min_len"]]


def dorado_read_stats_files(wildcards):
    """With dorado_correct on: the reads going into dorado correct and the
    corrected reads before the min_len sweep. Empty otherwise."""
    n = wildcards.input
    if not wants_dorado_correct(n):
        return []
    return [f"data/read_stats/{n}_dorado_{stage}.tsv" for stage in ("in", "out")]


def dorado_stage_reads(wildcards):
    n = wildcards.input
    if wildcards.stage == "in":
        return get_dorado_correct_input(wildcards)
    return f"data/seqtk/fasta_to_fastq/{n}.fastq"


# seqkit stats -j only runs separate FILES in parallel, so a one-file job
# gets one thread (they ran at 3% CPU efficiency on 4); ~2.5 GB peak.
#
# Measures exactly what hifiasm gets for one grid cell (chopper output,
# rasusa-capped, or length-filtered dorado-corrected reads). One small job per
# cell rather than one job over all of them, so each cell's temp() FASTQ can
# be deleted as soon as its own hifiasm and stats are done instead of every
# cell's reads piling up on disk at once.
rule filtered_read_stats:
    input:
        a = get_assembly_input
    output:
        a = "data/read_stats/{input}_q{minq}_l{minlen}.tsv"
    threads:
        1
    resources:
        mem_mb = scaled_mem(0, 4000),
        runtime = 120,
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input.a} > {output.a}
        """


# How much dorado correct keeps: its input (one fixed min_q/min_len cutoff)
# and its output. Separate small jobs, like filtered_read_stats, so neither
# temp() FASTQ is kept on disk waiting for read_stats.
rule dorado_read_stats:
    input:
        a = dorado_stage_reads
    output:
        a = "data/read_stats/{input}_dorado_{stage}.tsv"
    wildcard_constraints:
        stage = r"in|out"
    threads:
        1
    resources:
        mem_mb = scaled_mem(0, 4000),
        runtime = 120,
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input.a} > {output.a}
        """


rule read_stats:
    input:
        reads = read_stage_files,
        dorado = dorado_read_stats_files,
        cells = filtered_read_stats_files,
    output:
        a = "results/{input}/read_stats.tsv"
    threads:
        2
    resources:
        mem_mb = scaled_mem(0, 4000),
        runtime = 120,
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit stats -a -T -j {threads} {input.reads} > {output.a}
        for f in {input.dorado} {input.cells}; do tail -n +2 "$f" >> {output.a}; done
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
        mem_mb = 2000,
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
        mem_mb = 1000,
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
        mem_mb = 1000,
        runtime = 30,
    script:
        "../scripts/summarize.py"