# Reached for both pod5 (our own basecalled BAM) and bam (user-supplied)
# entry points -- either way the BAM needs converting to plain FASTQ before
# it can enter porechop/chopper. get_bam() resolves which BAM applies.
#
# Duplex BAMs (dorado duplex) hold each duplex read (dx:i:1) AND its two
# simplex parents (dx:i:-1); the parents are dropped so each molecule enters
# the assembly once. Reads with dx:i:0 or no dx tag (all simplex data) are
# kept. The log ends with how many reads were kept and how many are duplex.
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
    log:
        "logs/bam_to_fastq/{input}.log"
    container:
        "docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools view -u --threads $(nproc) -e '!([dx]==-1)' {input.a} 2> {log} \
            | samtools fastq --threads $(nproc) - > {output.a} 2>> {log}
        awk 'NR % 4 == 1 {{ n++; if (index($1, ";")) d++ }}
             END {{ print n + 0 " reads kept (duplex parents removed), " d + 0 " of them duplex" }}' \
            {output.a} >> {log}
        """
