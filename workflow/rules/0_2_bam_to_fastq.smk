configfile: "config/config.yaml"

rule bam_to_fastq:
    input:
        a = "data/dorado_basecall/{input}.bam"
    output:
        a = "data/samtools/Fastq/{input}.fastq"
    threads:
        12
    resources:
        mem_mb=resources["chopper"]["mem_mb"],
        runtime=resources["chopper"]["time"]
    conda:
        "samtools.yml"
    shell:
        """
        samtools fastq --threads $(nproc) {input.a} > {output.a}  

        """
