configfile: "config/config.yaml"

rule seqtk_fasta_to_fastq:
    input:
        a = rules.dorado.output.a
    output:
        a = "data/seqtk/fasta_to_fastq/{input}.fastq"
    threads:
        10
    resources:
        mem_mb=resources["seqkit"]["mem_mb"],
        runtime=resources["seqkit"]["time"]
    conda:
        "seqtk.yml"
    shell:
        """
        seqtk seq -F '#' {input.a} > {output.a}
        """
























