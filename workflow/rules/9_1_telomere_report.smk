# ============================================================================
#  Telomeres on each final assembly
# ============================================================================
# With `telomere: enabled: true`, every selector's final (possibly polished)
# assembly is searched for the sample's telomere motif (5_1_telomere_motif.smk):
#   results/{sample}/{selector}/telomeres.tsv         -- per contig: repeats at
#                                                        each end, telomeric?, T2T?
#   results/{sample}/{selector}/telomeres/*.svg       -- tidk plot along contigs
# and summary.txt gets a TELOMERES section from it.

rule telomere_final_search:
    input:
        fa    = "results/{input}/{selector}/{input}_{selector}_final.fasta",
        motif = "data/telomere/{input}/motif.tsv",
    output:
        tsv = "results/{input}/{selector}/telomeres/{input}_{selector}_telomeric_repeat_windows.tsv",
        svg = "results/{input}/{selector}/telomeres/{input}_{selector}.svg",
    params:
        window = TELO.get("window", 2000),
        dir    = "results/{input}/{selector}/telomeres",
        name   = "{input}_{selector}",
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
        mkdir -p {params.dir}
        motif=$(awk -F'\\t' '$1 == "motif" {{ print $2 }}' {input.motif})
        if [ -s {input.fa} ]; then
            tidk search --string "$motif" --window {params.window} \
                --output {params.name} --dir {params.dir} {input.fa}
            tidk plot --tsv {output.tsv} --output {params.dir}/{params.name}
        else
            : > {output.tsv}
            : > {output.svg}
        fi
        """


rule telomere_report:
    input:
        search = "results/{input}/{selector}/telomeres/{input}_{selector}_telomeric_repeat_windows.tsv",
        motif  = "data/telomere/{input}/motif.tsv",
    output:
        a = "results/{input}/{selector}/telomeres.tsv",
    params:
        min_repeats = TELO.get("min_repeats", 10),
        window      = TELO.get("window", 2000),
    threads:
        1
    resources:
        mem_mb  = 2000,
        runtime = 30,
    script:
        "../scripts/telomere_report.py"
