# ============================================================================
#  Per-sample HTML report: results/{sample}/report.html
# ============================================================================
# One page per sample for reading the filter grid without the tools: what
# each read filter kept, hifiasm's k-mer histogram per grid cell, BUSCO per
# grid assembly and for the finals, and which cell each selector picked
# (make_report.py fills report_template.html; the charts draw in the
# browser, so the page is one self-contained file). On unless config
# `report: enabled: false`.
#
# Like 9_stats.smk, report_reads and report_busco_grid CONSUME temp()
# intermediates (the reads going into filtering, the grid's BUSCO folders),
# so on a run that has already deleted them they can't be made without
# rebuilding half the pipeline. final_targets() therefore only asks for the
# report where those inputs are still there or still to come -- see
# wants_report() in the Snakefile.

def report_corrected_reads(wildcards):
    return f"data/seqtk/fasta_to_fastq/{wildcards.input}.fastq"


# Length and mean quality of every read going into the filter grid, so the
# page can redraw what any (min_q, min_len) cutoff keeps. Mean Q is
# computed the way chopper computes it (from the mean error probability).
# About 10 bytes per read gzipped. One thread: fx2tab -j only parallelises
# across files.
rule report_reads:
    input:
        a = get_trimmed_fastq
    output:
        a = "data/report/{input}/reads.tsv.gz"
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = 120,
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit fx2tab -n -i -l -q {input.a} | cut -f 2,3 | gzip -c > {output.a}
        """


# dorado_correct samples: lengths of the corrected reads, before the
# min_len sweep (their quality values are made up, see 3_seqtk_fasta_to_fastq).
rule report_corrected_reads:
    input:
        a = report_corrected_reads
    output:
        a = "data/report/{input}/corrected.tsv.gz"
    threads:
        1
    resources:
        mem_mb = 2000,
        runtime = 120,
    container:
        "docker://quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit fx2tab -n -i -l {input.a} | cut -f 2 | gzip -c > {output.a}
        """


# One line per grid assembly: its busco_grid folder and BUSCO's
# "C:..%[S:..%,D:..%],F:..%,M:..%,n:.." result (n:0 for the placeholder of
# an assembly too small to search), kept after the folders are deleted.
rule report_busco_grid:
    input:
        dirs = busco_grid_candidates,
    output:
        a = "data/report/{input}/busco_grid.tsv"
    threads:
        1
    resources:
        mem_mb = 500,
        runtime = 10,
    shell:
        """
        for d in {input.dirs}; do
            f=$(find "$d" -name 'short_summary*.txt' | head -n 1)
            line=$( [ -n "$f" ] && grep -m 1 -o 'C:[^ ]*n:[0-9]*' "$f" || true)
            printf '%s\\t%s\\n' "$d" "$line"
        done > {output.a}
        """


def report_inputs(wildcards):
    n = wildcards.input
    d = {
        "reads_hist": f"data/report/{n}/reads.tsv.gz",
        "read_stats": f"results/{n}/read_stats.tsv",
        # Every selector's assembly_stats.tsv lists all grid assemblies.
        "asm_stats":  f"results/{n}/{SELECTORS[0]}/assembly_stats.tsv",
        "busco_grid": f"data/report/{n}/busco_grid.tsv",
        # summary.txt comes last for each selector, so its pick is final.
        "summaries":  [f"results/{n}/{sel}/summary.txt" for sel in SELECTORS],
        "busco_final": [f"data/busco/{n}_{sel}/BUSCO" for sel in SELECTORS],
    }
    if wants_dorado_correct(n):
        d["corrected_hist"] = f"data/report/{n}/corrected.tsv.gz"
    if can_polish(n):
        d["busco_unpolished"] = [f"data/busco_unpolished/{n}_{sel}/BUSCO" for sel in SELECTORS]
    return d


def report_cells(wildcards):
    """Per read cell: its grid cutoffs, the reads hifiasm got (as named in
    read_stats.tsv) and the pass-1 hifiasm log holding its k-mer histogram."""
    n = wildcards.input
    grid = filter_grid(n)
    out = []
    for q in grid["min_q"]:
        for l in grid["min_len"]:
            w = SimpleNamespace(input=n, minq=str(q), minlen=str(l))
            out.append({"q": q, "l": l, "reads": get_assembly_input(w),
                        "hifiasm_log": f"logs/hifiasm/{n}_q{q}_l{l}.log"})
    return out


def report_assemblies(wildcards):
    """Per grid assembly, in grid_assemblies() order: pass, cell, FASTA and
    the busco_grid folder it was searched in."""
    n = wildcards.input
    grid = filter_grid(n)
    passes = ["hifiasm", "hifiasm_telo"] if TELO_RERUN else ["hifiasm"]
    return [{"pass": asm, "q": q, "l": l,
             "fasta": f"data/{asm}/{n}_q{q}_l{l}/{n}_q{q}_l{l}.fa",
             "busco_dir": f"data/busco_grid/{asm}/{n}_q{q}_l{l}/BUSCO"}
            for asm in passes for q in grid["min_q"] for l in grid["min_len"]]


def report_select_logs(wildcards):
    n = wildcards.input
    logs = {"lowest_contig": f"logs/select/{n}_lowest_contig.log",
            "highest_busco": f"logs/select/{n}_highest_busco.log"}
    if TELOMERE_ON:
        logs["most_t2t"] = f"logs/telomere/{n}_most_t2t.log"
    return logs


rule report:
    input:
        unpack(report_inputs)
    output:
        a = "results/{input}/report.html"
    params:
        sample         = lambda w: w.input,
        selectors      = SELECTORS,
        grid           = lambda w: filter_grid(w.input),
        cells          = report_cells,
        assemblies     = report_assemblies,
        select_logs    = report_select_logs,
        dorado_correct = lambda w: wants_dorado_correct(w.input),
        dorado_q       = config["dorado_correct"]["min_q"],
        dorado_len     = config["dorado_correct"]["min_len"],
        capped         = lambda w: subsampled(w.input),
        polished       = lambda w: can_polish(w.input),
        lineage        = lambda w: busco_lineage(w.input),
    threads:
        1
    resources:
        mem_mb = 4000,
        runtime = 60,
    script:
        "../scripts/make_report.py"
