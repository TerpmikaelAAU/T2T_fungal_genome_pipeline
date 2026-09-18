# dorado correct outputs FASTA (no quality scores); hifiasm --ont wants
# FASTQ. Only reached when dorado_correct is enabled for this sample -- see
# wants_dorado_correct()/get_assembly_input() -- padding every base with a
# fake '#' quality score. Runs once per sample; the min_len sweep happens
# afterward in 3_2_length_filter_corrected.smk.
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
    container:
        "docker://quay.io/biocontainers/seqtk:1.5--h577a1d6_1"
    conda:
        "../envs/seqtk.yml"
    shell:
        """
        seqtk seq -F '#' {input.a} > {output.a}
        """