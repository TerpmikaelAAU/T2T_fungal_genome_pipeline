# Picks the winning assembly across a sample's (min_q, min_len) grid: the
# fewest-contigs candidate, on the assumption that's the most contiguous
# assembly. Cells smaller than the sample's min_assembly_mb (config.yaml)
# don't count -- a tiny fragment would otherwise win with 1 contig.
# CAVEAT for repeat-rich genomes (e.g. oomycetes): fewest contigs can mean
# COLLAPSED repeats rather than a genuinely better assembly -- cross-check
# the winner's total length in assembly_stats.tsv against the expected
# genome size before trusting this pick.
def contig_count_candidates(wildcards):
    # Both hifiasm passes when the telomere re-run is on (see Snakefile
    # grid_assemblies()); its placeholder cells are empty, so skipped.
    return grid_assemblies(wildcards.input)

rule contig_count:
    input:
        fasta_files = contig_count_candidates
    output:
        a = temp("data/contig/{input}_lowest_contig_file.fa")
    params:
        min_bp = lambda w: min_assembly_bp(w.input),
    threads:
        1
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["runtime"],
    log:
        "logs/select/{input}_lowest_contig.log"
    script:
        "../scripts/pick_lowest_contig.py"
