rule bam_to_fastq:
    input:
        a = get_bam
    output:
        a = temp("data/samtools/Fastq/{input}.fastq")
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["runtime"]
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools fastq --threads $(nproc) {input.a} > {output.a}  

        """