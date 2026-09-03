rule seqtk_fasta_to_fastq:
    input:
        a = "data/dorado/{input}.fasta"
    output:
        a = temp("data/seqtk/fasta_to_fastq/{input}.fastq")
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["runtime"]
    conda:
        "../envs/seqtk.yml"
    shell:
        """
        seqtk seq -F '#' {input.a} > {output.a}
        """