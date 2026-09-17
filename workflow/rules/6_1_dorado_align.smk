# Only reached for pod5/bam entry points (has_bam() -- polishing needs a
# BAM's move table). Aligns the original basecalled reads back onto the
# WINNING grid-cell assembly (contig_count's output), producing the input
# dorado_polish needs.
rule dorado_align:
    input:
        dorado = dorado_bin,
        a = rules.contig_count.output.a,
        b = get_bam
    output: 
        a = temp("data/dorado_align/{input}.bam")
    threads:
        50
    resources:
        mem_mb=scaled_mem(1.0, 64000),
        runtime=resources["dorado_align"]["runtime"],
    conda:
        "../envs/samtools.yml"
    shell:
        """
        "{input.dorado}" aligner {input.a} {input.b} | samtools sort --threads $(nproc) > {output.a}  
        echo "align done!"
        samtools index {output.a}

        """