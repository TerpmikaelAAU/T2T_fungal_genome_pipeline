# Genome completeness check on the deliverable assembly (polished if the
# sample could be polished -- can_polish() in the Snakefile -- else the winning raw grid-cell assembly) -- run once per
# selector (lowest_contig and highest_busco, see Snakefile SELECTORS), since
# each selector can pick a different assembly and, if polished, a
# differently-polished one too. `lineage` is per-sample (busco_lineage(),
# default fungi_odb12) -- pick the dataset that matches the organism, e.g.
# stramenopiles_odb* for oomycetes like Phytophthora, which are NOT fungi
# despite this pipeline's name.
rule busco:
    input:
        a = lambda w: (f"data/dorado_polish/{w.input}_{w.selector}.fasta" if can_polish(w.input)
             else f"data/contig/{w.input}_{w.selector}_file.fa")
    output:
        # Only the short summary gets pulled into results/summary.txt (see
        # summarize.py); the full tables/logs/per-gene predictions here are
        # not, so once summary.txt exists this whole directory is dead
        # weight -- temp() it like everything else that isn't a deliverable.
        dir = temp(directory("data/busco/{input}_{selector}/BUSCO")),
    params:
        lineage = lambda w: busco_lineage(w.input),
        min_bp  = lambda w: min_assembly_bp(w.input),
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
    log:
        "logs/busco/{input}_{selector}.log"
    container:
        "docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2"
    conda:
       "../envs/BUSCO.yml"
    shell:
        """
        mkdir -p {output.dir}
        bp=$(grep -v '^>' {input.a} | tr -d '\\n' | wc -c || true)
        if [ "$bp" -lt {params.min_bp} ]; then
            # Only reached when no grid cell reached min_assembly_mb (the
            # selectors' fallback) -- BUSCO can crash on a tiny, gene-less
            # assembly, so note it in summary.txt instead.
            echo "WARNING: {input.a} is $bp bp (< {params.min_bp}); skipping BUSCO" > {log}
            printf '# BUSCO skipped: assembly is %s bp, below min_assembly_mb\\n\\tC:0.0%%[S:0.0%%,D:0.0%%],F:0.0%%,M:100.0%%,n:0\\n' "$bp" \
                > {output.dir}/short_summary.placeholder.txt
        else
            busco -i {input.a} -o {output.dir} -l {params.lineage} -m geno -f -c $(nproc) --metaeuk --tar > {log} 2>&1
        fi
        """

# Same BUSCO run on the UNPOLISHED winner, for polished samples only (an
# unpolished sample's `busco` above already is its "before"). With
# both, summary.txt shows BUSCO before and after polishing side by side, to
# see what dorado polish gained (or lost) for each selector.
use rule busco as busco_unpolished with:
    input:
        a = "data/contig/{input}_{selector}_file.fa"
    output:
        dir = temp(directory("data/busco_unpolished/{input}_{selector}/BUSCO")),
    log:
        "logs/busco/{input}_{selector}_unpolished.log"
