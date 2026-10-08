# Only reached for pod5/bam entry points (has_bam() -- polishing needs a
# BAM). Aligns the original basecalled reads back onto a WINNING grid-cell
# assembly, producing the input dorado_polish needs. Runs once per DISTINCT
# winner: {selector} here is the first selector that picked it (see
# 6_05_polish_once.smk), so identical winners are aligned only once.
#
# Two rules, because dorado and samtools can't share one container: the
# dorado binary needs the host's system libraries (libz etc.), which the
# samtools biocontainer doesn't have ("libz.so.1: cannot open shared object
# file"). So dorado_align_raw runs dorado on the host, like every other
# dorado rule, and dorado_align does the samtools work in its container.
rule dorado_align_raw:
    input:
        dorado = dorado_bin,
        a = unique_winner,
        b = get_bam
    output:
        a = temp("data/dorado_align_raw/{input}_{selector}.bam"),
    threads:
        50
    resources:
        mem_mb=scaled_mem(1.0, 64000),
        runtime=resources["dorado_align"]["runtime"],
    log:
        "logs/dorado_align/{input}_{selector}_aligner.log"
    shell:
        """
        "{input.dorado}" aligner --threads {threads} {input.a} {input.b} > {output.a} 2> {log}
        """

# dorado polish refuses duplex data, so duplex reads (dx:i:1) are dropped
# here and the duplex ("stereo") read group is removed from the header. The
# simplex reads left -- including the parents of every duplex read -- still
# cover all the data. A no-op for simplex BAMs.
rule dorado_align:
    input:
        a = rules.dorado_align_raw.output.a,
    output:
        a = temp("data/dorado_align/{input}_{selector}.bam"),
        bai = temp("data/dorado_align/{input}_{selector}.bam.bai"),
    threads:
        50
    resources:
        mem_mb=scaled_mem(1.0, 64000),
        runtime=resources["dorado_align"]["runtime"],
    container:
        "docker://quay.io/biocontainers/samtools:1.24--h9dcdb79_1"
    conda:
        "../envs/samtools.yml"
    log:
        "logs/dorado_align/{input}_{selector}.log"
    shell:
        """
        samtools view -u -e '!([dx]==1)' {input.a} 2> {log} \
            | samtools sort --threads {threads} -o {output.a}.tmp.bam - 2>> {log}
        samtools view --no-PG -H {output.a}.tmp.bam > {output.a}.header.sam
        if grep -q '^@RG.*stereo' {output.a}.header.sam; then
            echo "duplex BAM: removing the duplex read group(s) for dorado polish" >> {log}
            grep -v '^@RG.*stereo' {output.a}.header.sam > {output.a}.simplex.sam
            samtools reheader {output.a}.simplex.sam {output.a}.tmp.bam > {output.a}
            rm {output.a}.tmp.bam {output.a}.simplex.sam
        else
            mv {output.a}.tmp.bam {output.a}
        fi
        rm {output.a}.header.sam
        samtools index {output.a}
        """
