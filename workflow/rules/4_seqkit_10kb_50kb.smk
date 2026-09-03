rule seqkit_10_50:
    input:
        a = rules.seqtk_fasta_to_fastq.output.a
    output:
        a = temp("data/seqkit/dorado/{input}_{length}.fastq")
        
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["runtime"]
    conda:
        "../envs/seqkit.yml"
    shell:
        """
        seqkit seq -m {wildcards.length} {input.a} > {output.a}
  
        """