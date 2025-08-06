configfile: "config/config.yaml"
rule dorado_align:
    input:
        a = rules.contig_count.output.a,
        b = rules.dorado_basecall.output.a
    output: 
        a = "data/dorado_align/{input}.bam"
    threads:
        50
    resources:
        mem_mb=resources["dorado_align"]["mem_mb"],
        runtime=resources["dorado_align"]["time"],
        partition="shared"
    conda:
        "samtools.yml"
    shell:
        """
        "{config[dorado]}" aligner {input.a} {input.b} | samtools sort --threads $(nproc) > {output.a}  
        echo "align done!"
        samtools index {output.a}

        """
