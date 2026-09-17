rule busco:
    input:
        a = lambda w: (f"data/dorado_polish/{w.input}.fasta" if has_bam(w.input)
             else f"data/contig/{w.input}_lowest_contig_file.fa")
    output:
        # Only the short summary gets pulled into results/summary.txt (see
        # summarize.py); the full tables/logs/per-gene predictions here are
        # not, so once summary.txt exists this whole directory is dead
        # weight -- temp() it like everything else that isn't a deliverable.
        dir = temp(directory("data/busco/{input}/BUSCO")),
    params:
        lineage = lambda w: busco_lineage(w.input)
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
        #
    conda:
       "../envs/BUSCO.yml"
    shell:
        """
        
        busco -i {input.a} -o {output.dir} -l {params.lineage} -m geno -f -c $(nproc) --metaeuk --tar
        
        
        """