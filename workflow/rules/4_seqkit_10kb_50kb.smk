configfile: "config/config.yaml"

rule seqkit_10_50:
    input:
        a = rules.seqtk_fasta_to_fastq.output.a
    output:
        a = "data/seqkit/dorado/{input}_{length}.fastq"
        
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["time"]
    conda:
        "seqkit"
    shell:
        """
        seqkit seq -m {wildcards.length} {input.a} > {output.a}
  
        """



























