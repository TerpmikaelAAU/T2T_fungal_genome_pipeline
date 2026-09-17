# dorado correct outputs FASTA (no quality scores); hifiasm --ont wants
# FASTQ. Only reached when config `dorado_correct.enabled: true` -- see
# get_assembly_input() -- padding every base with a fake '#' quality score.
rule seqtk_fasta_to_fastq:
    input:
        a = "data/dorado/{input}_q{minq}_l{minlen}.fasta"
    output:
        a = temp("data/seqtk/fasta_to_fastq/{input}_q{minq}_l{minlen}.fastq")
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