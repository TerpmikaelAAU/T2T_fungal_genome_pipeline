# ============================================================================
#  Alternative grid winner: highest BUSCO completeness, not fewest contigs
# ============================================================================
# contig_count (6_00) picks the fewest-contigs assembly across the grid,
# which can reward COLLAPSED repeats on a repeat-rich genome (see its own
# comment). This picks a second, independent winner using an orthogonal,
# biology-grounded signal instead: BUSCO completeness. Both winners go
# through the same rest of the pipeline (see Snakefile SELECTORS /
# final_targets()), so results/{sample}/ ends up with two comparable
# deliverables -- lowest_contig/ and highest_busco/.
#
# COST WARNING: this runs BUSCO on every grid cell, not just the eventual
# winner -- for a 3x4=12-cell grid that's 12 BUSCO runs instead of 1.

rule busco_grid:
    input:
        a = "data/hifiasm/{input}_q{minq}_l{minlen}/{input}_q{minq}_l{minlen}.fa"
    output:
        # temp(): only the winning cell's BUSCO result matters afterwards,
        # and that gets re-run anyway by the final `busco` rule (7_BUSCO.smk)
        # on the possibly-polished assembly -- this is purely for selection.
        dir = temp(directory("data/busco_grid/{input}_q{minq}_l{minlen}/BUSCO")),
    params:
        lineage = lambda w: busco_lineage(w.input)
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
    container:
        "docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2"
    conda:
        "../envs/BUSCO.yml"
    shell:
        """
        mkdir -p {output.dir}
        if [ ! -s {input.a} ]; then
            # Empty/placeholder grid cell (see 5_hifiasm.smk's near-empty-input
            # guard) -- BUSCO can't run on nothing, so score it 0%% instead of
            # crashing, matching how contig_count treats these cells (ignored).
            echo "WARNING: {input.a} is empty; skipping BUSCO, writing a placeholder 0%% summary" >&2
            echo "C:0.0%[S:0.0%,D:0.0%],F:0.0%,M:100.0%,n:0" > {output.dir}/short_summary.placeholder.txt
        else
            busco -i {input.a} -o {output.dir} -l {params.lineage} -m geno -f -c {threads} --metaeuk --tar
        fi
        """


def busco_grid_candidates(wildcards):
    grid = filter_grid(wildcards.input)
    return expand("data/busco_grid/{input}_q{minq}_l{minlen}/BUSCO",
                  input=wildcards.input,
                  minq=grid["min_q"],
                  minlen=grid["min_len"])


rule highest_busco:
    input:
        # contig_count_candidates is defined in 6_00_lowest_contig_count.smk,
        # included just before this file -- both expand() over the exact
        # same grid in the same order, so busco_dirs[i] and fasta_files[i]
        # always refer to the same grid cell (see pick_highest_busco.py).
        busco_dirs  = busco_grid_candidates,
        fasta_files = contig_count_candidates,
    output:
        a = temp("data/contig/{input}_highest_busco_file.fa")
    threads:
        1
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["runtime"],
    script:
        "../scripts/pick_highest_busco.py"
