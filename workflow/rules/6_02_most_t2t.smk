# ============================================================================
#  Third grid winner: most telomere-to-telomere contigs
# ============================================================================
# With `telomere: enabled: true`, alongside fewest contigs (6_00) and
# highest BUSCO (6_01): count, in every grid cell of both hifiasm passes,
# the contigs with a telomere at BOTH ends (the sample's motif, see
# 5_1_telomere_motif.smk) and pick the cell with the most. This is the most
# direct measure of what the pipeline is for -- complete chromosomes -- and
# goes through the same polish/BUSCO/stats as the other selectors.

def telomere_grid_searches(wildcards):
    # data/<asm>/<cell>/<cell>.fa -> data/telomere_grid/<asm>/<cell>.tsv,
    # in grid_assemblies() order so pick_most_t2t.py can pair by position.
    return [f"data/telomere_grid/{fa.split('/')[1]}/{fa.split('/')[2]}.tsv"
            for fa in grid_assemblies(wildcards.input)]

rule telomere_grid_search:
    input:
        fa    = "data/{asm}/{input}_q{minq}_l{minlen}/{input}_q{minq}_l{minlen}.fa",
        motif = "data/telomere/{input}/motif.tsv",
    output:
        tsv = temp("data/telomere_grid/{asm}/{input}_q{minq}_l{minlen}.tsv"),
    params:
        window = TELO.get("window", 2000),
        tmp    = "data/telomere_grid/{asm}/{input}_q{minq}_l{minlen}_tmp",
    threads:
        1
    resources:
        mem_mb  = resources["tidk"]["mem_mb"],
        runtime = resources["tidk"]["runtime"],
    container:
        TIDK_CONTAINER
    conda:
        "../envs/tidk.yml"
    shell:
        """
        motif=$(awk -F'\\t' '$1 == "motif" {{ print $2 }}' {input.motif})
        if [ -s {input.fa} ]; then
            rm -rf {params.tmp}
            tidk search --string "$motif" --window {params.window} \
                --output cell --dir {params.tmp} {input.fa}
            mv {params.tmp}/cell_telomeric_repeat_windows.tsv {output.tsv}
            rm -rf {params.tmp}
        else
            : > {output.tsv}  # empty/placeholder cell -- skipped by the picker
        fi
        """


rule most_t2t:
    input:
        # Both in grid_assemblies() order -- paired by position.
        searches    = telomere_grid_searches,
        fasta_files = contig_count_candidates,
    output:
        a = temp("data/contig/{input}_most_t2t_file.fa")
    params:
        min_repeats = TELO.get("min_repeats", 10),
        min_bp      = lambda w: min_assembly_bp(w.input),
    threads:
        1
    resources:
        mem_mb=resources["fga"]["mem_mb"],
        runtime=resources["fga"]["runtime"],
    log:
        "logs/telomere/{input}_most_t2t.log"
    script:
        "../scripts/pick_most_t2t.py"
