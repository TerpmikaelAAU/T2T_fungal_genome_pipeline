configfile: "config/config.yaml"

rule flye:
    input:
        a = rules.chopper_flye.output.a
    output:
        a = directory("data/flye/{input}/flye_assembly"),
        b = "data/flye/{input}/flye_assembly/assembly.fasta",
        c = "data/flye/{input}/flye_assembly/flye.log",
        d = "data/flye/{input}/flye_assembly/assembly_graph.gfa",
    threads:
        12
    resources:
        mem_mb=resources["flye"]["mem_mb"],
        runtime=resources["flye"]["time"],
        partition="general"
    conda:
        "Flye.yml"
    shell:
        """
       
        flye -t $(nproc) --nano-hq {input.a} -o {output.a} --no-alt-contigs --scaffold --read-error 0.03
        
        """