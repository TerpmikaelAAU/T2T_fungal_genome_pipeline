# Optional coverage cap. Only built when a sample sets `subsample:` in config.
# NOTE: this is now a deliberate coverage choice, NOT a RAM workaround --
# dorado correct handles memory via blocks (see 3_dorado_correct.smk).
rule rasusa:
    input:
        a = "data/chopper/{input}_q{minq}_l{minlen}.fastq",
    output:
        a = temp("data/rasusa/Coverage/{input}_q{minq}_l{minlen}.fastq")
    params:
        cov = lambda w: SAMPLES[w.input]["subsample"]["coverage"],
        gsize = lambda w: SAMPLES[w.input]["subsample"]["genome_size"],
    threads:
        12
    resources:
        mem_mb=resources["rasusa"]["mem_mb"],
        runtime=resources["rasusa"]["runtime"]
    container:
        "docker://quay.io/biocontainers/rasusa:5.1.0--hfa8f182_0"
    conda:
        "../envs/rasusa.yml"
    log:
        "logs/rasusa/{input}_q{minq}_l{minlen}.log"
    shell:
        """
        rasusa reads --coverage {params.cov} --genome-size {params.gsize} \
            --seed 1 -o {output.a} {input.a} 2> {log}
        """