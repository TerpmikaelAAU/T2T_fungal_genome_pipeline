# ============================================================================
#  Telomere motif discovery (tidk) + telomere-aware hifiasm re-run
# ============================================================================
# Only reached with `telomere: enabled: true` (config.yaml). Works for any
# genome -- nothing assumes TTAGGG:
#   1. telomere_explore  -- `tidk explore` on this sample's pass-1 grid cell
#                           with the fewest contigs lists candidate repeat
#                           units near contig ends; the top `top_motifs` are
#                           each checked with `tidk search`.
#   2. telomere_motif    -- picks the motif that makes the most contig ends
#                           telomeric (pick_telomere_motif.py), and decides
#                           whether it explains enough of the telomeres
#                           (`rerun_min_fraction`, `rerun_min_ends`) to
#                           re-run hifiasm with it.
#   3. hifiasm_telo      -- the same grid cell re-assembled with
#                           `hifiasm --telo-m <motif>`, which lets hifiasm
#                           keep telomeres on contig ends. If the motif didn't
#                           qualify this writes an empty placeholder instead,
#                           which every selector skips. The selectors choose
#                           across BOTH passes (Snakefile grid_assemblies()).
#
# COST: a qualifying motif re-runs hifiasm for every grid cell (the pass-1
# overlap caches are deleted to save disk, so it's a full run each).

TIDK_CONTAINER = "docker://quay.io/biocontainers/tidk:0.2.7--h6872113_0"

def telomere_explore_candidates(wildcards):
    # Pass 1 only: the motif decides whether pass 2 runs at all.
    return grid_assemblies(wildcards.input, passes=["hifiasm"])

rule telomere_explore:
    input:
        assemblies = telomere_explore_candidates,
    output:
        candidates = temp(directory("data/telomere/{input}/candidates")),
    params:
        top       = TELO.get("top_motifs", 10),
        min_len   = TELO.get("explore_min_len", 5),
        max_len   = TELO.get("explore_max_len", 12),
        threshold = TELO.get("explore_threshold", 10),
        window    = TELO.get("window", 2000),
    threads:
        2
    resources:
        mem_mb  = resources["tidk"]["mem_mb"],
        runtime = resources["tidk"]["runtime"],
    log:
        "logs/telomere/{input}_explore.log"
    container:
        TIDK_CONTAINER
    conda:
        "../envs/tidk.yml"
    shell:
        """
        mkdir -p {output.candidates}
        : > {log}
        # The pass-1 cell with the fewest contigs (same pick as contig_count).
        best=$(for f in {input.assemblies}; do
            n=$(grep -c "^>" "$f" || true)
            if [ "$n" -gt 0 ]; then echo "$n $f"; fi
        done | sort -n | awk 'NR == 1 {{ print $2 }}')
        printf 'rank\\tmotif\\texplore_count\\n' > {output.candidates}/candidates.tsv
        if [ -z "$best" ]; then
            echo "No non-empty pass-1 assembly to explore" >> {log}
            exit 0
        fi
        echo "Exploring $best" >> {log}
        tidk explore --minimum {params.min_len} --maximum {params.max_len} \
            --threshold {params.threshold} "$best" \
            > {output.candidates}/explore.tsv 2>> {log}
        # explore sorts by count; keep the top N (rank, motif, count).
        awk -F'\\t' -v OFS='\\t' -v top={params.top} \
            'NR > 1 && NR <= top + 1 {{ print NR - 1, $1, $2 }}' \
            {output.candidates}/explore.tsv >> {output.candidates}/candidates.tsv
        tail -n +2 {output.candidates}/candidates.tsv | while IFS=$'\\t' read -r rank motif count; do
            tidk search --string "$motif" --window {params.window} \
                --output "rank$rank" --dir {output.candidates} "$best" 2>> {log}
        done
        """


rule telomere_motif:
    input:
        candidates = "data/telomere/{input}/candidates",
    output:
        motif = "data/telomere/{input}/motif.tsv",
        table = "results/{input}/telomere_motif_candidates.tsv",
    params:
        min_repeats   = TELO.get("min_repeats", 10),
        min_fraction  = TELO.get("rerun_min_fraction", 0.8),
        min_ends      = TELO.get("rerun_min_ends", 4),
        fallback      = TELO.get("fallback_motif", "TTAGGG"),
        rerun_enabled = TELO_RERUN,
    threads:
        1
    resources:
        mem_mb  = 2000,
        runtime = 30,
    script:
        "../scripts/pick_telomere_motif.py"


rule hifiasm_telo:
    input:
        unpack(hifiasm_inputs),
        motif = "data/telomere/{input}/motif.tsv",
    output:
        fa = temp("data/hifiasm_telo/{input}_q{minq}_l{minlen}/{input}_q{minq}_l{minlen}.fa"),
    threads:
        32
    resources:
        mem_mb=scaled_mem(4.0, 64000),
        runtime=scaled_time(0.3, 720),
    params:
        ul_flag   = lambda w, input: f"--ul {input.b}" if wants_ultralong(w.input) else "",
        min_bytes = int(config.get("hifiasm_min_input_mb", 1)) * 1_000_000,
        prefix    = lambda w: f"data/hifiasm_telo/{w.input}_q{w.minq}_l{w.minlen}/prefix",
    log:
        "logs/telomere/{input}_q{minq}_l{minlen}_hifiasm_telo.log"
    container:
        "docker://quay.io/biocontainers/hifiasm:0.25.0--h5ca1c30_0"
    conda:
        "../envs/hifiasm.yml"
    shell:
        """
        mkdir -p $(dirname {params.prefix})
        rerun=$(awk -F'\\t' '$1 == "rerun" {{ print $2 }}' {input.motif})
        motif5=$(awk -F'\\t' '$1 == "motif_5prime" {{ print $2 }}' {input.motif})
        size=$(stat -c%s {input.a} 2>/dev/null || echo 0)
        if [ "$rerun" != "yes" ]; then
            echo "Telomere motif didn't qualify for a re-run (see {input.motif}); empty placeholder" > {log}
            : > {output.fa}
        elif [ "$size" -lt {params.min_bytes} ]; then
            # Same near-empty-input guard as the pass-1 hifiasm rule.
            echo "{input.a} is ${{size}} bytes; empty placeholder" > {log}
            : > {output.fa}
        else
            echo "Re-running hifiasm with --telo-m $motif5" > {log}
            hifiasm --ont -t {threads} --primary {params.ul_flag} --telo-m "$motif5" \
                -o "{params.prefix}" {input.a} 2>> {log}
            awk '/^S/ {{ print ">"$2; print $3 }}' {params.prefix}.p_ctg.gfa > {output.fa}
            # Only the primary contigs are used; drop graphs, beds and the
            # large overlap caches (see the pass-1 rule).
            rm -f {params.prefix}.*
        fi
        """
