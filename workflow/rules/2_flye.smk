# Only built when a sample wants organelle recovery (wants_organelle()).
# flye is NOT the main nuclear assembler -- hifiasm is -- this is purely the
# input flye assembly graph (output.d) that getorganelle mines for a
# mitochondrial contig.
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
        runtime=resources["flye"]["runtime"],
    container:
        "docker://quay.io/biocontainers/flye:2.9.6--py312h734f728_1"
    conda:
        "../envs/Flye.yml"
    shell:
        """
       
        flye -t $(nproc) --nano-hq {input.a} -o {output.a} --no-alt-contigs --scaffold --read-error 0.03
        
        """