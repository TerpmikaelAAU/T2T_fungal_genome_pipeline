# Only reached for pod5/bam entry points (has_bam() -- polishing needs a
# BAM). Aligns the original basecalled reads back onto a WINNING grid-cell
# assembly, producing the input dorado_polish needs. Runs once per DISTINCT
# winner: {selector} here is the first selector that picked it (see
# 6_05_polish_once.smk), so identical winners are aligned only once.
rule dorado_align:
    input:
        dorado = dorado_bin,
        a = unique_winner,
        b = get_bam
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
    shell:
        """
        "{input.dorado}" aligner {input.a} {input.b} | samtools sort --threads $(nproc) > {output.a}  
        echo "align done!"
        samtools index {output.a}

        """