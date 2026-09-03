rule busco:
    input:
        a = lambda w: (f"data/dorado_polish/{w.input}.fasta" if has_bam(w.input)
             else f"data/contig/{w.input}_lowest_contig_file.fa")
    output:
        dir = directory("data/busco/{input}/BUSCO"),
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