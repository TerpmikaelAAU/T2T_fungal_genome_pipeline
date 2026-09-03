rule dorado_align:
    input:
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
        "{config[dorado]}" aligner {input.a} {input.b} | samtools sort --threads $(nproc) > {output.a}  
        echo "align done!"
        samtools index {output.a}

        """