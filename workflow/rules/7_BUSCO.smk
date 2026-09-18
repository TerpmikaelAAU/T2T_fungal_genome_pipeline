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