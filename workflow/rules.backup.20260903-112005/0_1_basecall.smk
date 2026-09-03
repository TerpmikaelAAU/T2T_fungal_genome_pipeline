configfile: "config/config.yaml"
import os

rule dorado_basecall:
    input:
        a = lambda wildcards: os.path.join(config["pod5_folder"], wildcards.input)
    output:
        a = "data/dorado_basecall/{input}.bam"
    threads:
        16
    resources:
        mem_mb=resources["dorado_basecall"]["mem_mb"],
        runtime=resources["dorado_basecall"]["time"],
        partition="gpu",
        gpus=1,
    shell:
        """
        "{config[dorado]}" basecaller sup --emit-moves --device cuda:all {input.a} > {output.a}
        
        """


