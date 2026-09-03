# Optional coverage cap. Only built when a sample sets `subsample:` in config.
# NOTE: this is now a deliberate coverage choice, NOT a RAM workaround --
# dorado correct handles memory via blocks (see 3_dorado_correct.smk).
rule rasusa:
    input:
        a = rules.chopper_dorado.output.a,
    output:
        a = temp("data/rasusa/Coverage/{input}.fastq")
    params:
        cov = lambda w: SAMPLES[w.input]["subsample"]["coverage"],
        gsize = lambda w: SAMPLES[w.input]["subsample"]["genome_size"],
    threads:
        12
    resources:
        mem_mb=resources["rasusa"]["mem_mb"],
        runtime=resources["rasusa"]["runtime"]
    conda:
        "../envs/rasusa.yml"
    log:
        "logs/rasusa/{input}.log"
    shell:
        """
        rasusa reads --coverage {params.cov} --genome-size {params.gsize} \
            --seed 1 -o {output.a} {input.a} 2> {log}
        """