# Genome completeness check on the deliverable assembly (polished if a BAM
# was available, else the winning raw grid-cell assembly) -- run once per
# selector (lowest_contig and highest_busco, see Snakefile SELECTORS), since
# each selector can pick a different assembly and, if polished, a
# differently-polished one too. `lineage` is per-sample (busco_lineage(),
# default fungi_odb12) -- pick the dataset that matches the organism, e.g.
# stramenopiles_odb* for oomycetes like Phytophthora, which are NOT fungi
# despite this pipeline's name.
rule busco:
    input:
        a = lambda w: (f"data/dorado_polish/{w.input}_{w.selector}.fasta" if has_bam(w.input)
             else f"data/contig/{w.input}_{w.selector}_file.fa")
    output:
        # Only the short summary gets pulled into results/summary.txt (see
        # summarize.py); the full tables/logs/per-gene predictions here are
        # not, so once summary.txt exists this whole directory is dead
        # weight -- temp() it like everything else that isn't a deliverable.
        dir = temp(directory("data/busco/{input}_{selector}/BUSCO")),
    params:
        lineage = lambda w: busco_lineage(w.input)
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
        #
    container:
        "docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2"
    conda:
       "../envs/BUSCO.yml"
    shell:
        """
        
        busco -i {input.a} -o {output.dir} -l {params.lineage} -m geno -f -c $(nproc) --metaeuk --tar
        
        
        """

# Same BUSCO run on the UNPOLISHED winner, for pod5/bam samples only (fastq
# samples aren't polished, so `busco` above already is their "before"). With
# both, summary.txt shows BUSCO before and after polishing side by side, to
# see what dorado polish gained (or lost) for each selector.
use rule busco as busco_unpolished with:
    input:
        a = "data/contig/{input}_{selector}_file.fa"
    output:
        dir = temp(directory("data/busco_unpolished/{input}_{selector}/BUSCO")),
