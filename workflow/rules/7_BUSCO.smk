configfile: "config/config.yaml"

rule busco:
    input:
        a = rules.dorado_polish.output.a
        #Use the line below instead if you do not do dorado polishing
        #a = "data/contig/{input}_lowest_contig_file.fa"
    output:
        dir = directory("data/busco/{input}/BUSCO"),
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["time"],
        partition="general"
        #
    conda:
       "BUSCO.yml"
    shell:
        """
        
        busco -i {input.a} -o {output.dir} -l fungi_odb10 -m geno -f -c $(nproc) --metaeuk --tar
        
        
        """