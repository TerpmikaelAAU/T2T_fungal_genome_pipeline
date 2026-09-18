# Only reached for pod5/bam entry points (has_bam() -- polishing needs a
# BAM's move table). Aligns the original basecalled reads back onto the
# WINNING grid-cell assembly -- whichever selector (lowest_contig or
# highest_busco, see Snakefile SELECTORS) is being built -- producing the
# input dorado_polish needs.
rule dorado_align:
    input:
        dorado = dorado_bin,
        a = "data/contig/{input}_{selector}_file.fa",
        b = get_bam
    output:
        a = temp("data/dorado_align/{input}_{selector}.bam")
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