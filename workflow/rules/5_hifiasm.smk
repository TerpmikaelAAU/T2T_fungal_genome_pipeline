configfile: "config/config.yaml"

rule hifiasm:
    input:
        a = rules.seqkit_10_50.output.a,
        b = rules.chopper_ultralong.output.a,
    output:
        a = "data/hifiasm/{input}_{length}/prefix.p_ctg.gfa",
        b = ["data/hifiasm/{input}_{length}/prefix.a_ctg.gfa",
        "data/hifiasm/{input}_{length}/prefix.a_ctg.lowQ.bed",
        "data/hifiasm/{input}_{length}/prefix.a_ctg.noseq.gfa",
        "data/hifiasm/{input}_{length}/prefix.p_ctg.lowQ.bed",
        "data/hifiasm/{input}_{length}/prefix.p_ctg.noseq.gfa",
        "data/hifiasm/{input}_{length}/prefix.p_utg.gfa",
        "data/hifiasm/{input}_{length}/prefix.p_utg.lowQ.bed",
        "data/hifiasm/{input}_{length}/prefix.p_utg.noseq.gfa",
        "data/hifiasm/{input}_{length}/prefix.r_utg.gfa",
        "data/hifiasm/{input}_{length}/prefix.r_utg.lowQ.bed",
        "data/hifiasm/{input}_{length}/prefix.r_utg.noseq.gfa",
        ]

    threads:
        12
    resources:
        mem_mb=resources["hifiasm"]["mem_mb"],
        runtime=resources["hifiasm"]["time"],
        partition="general"
    conda:
        "hifiasm.yml"
        
    shell:
        """
        hifiasm --ont -t $(nproc) --primary --ul {input.b} -o "data/hifiasm/{wildcards.input}_{wildcards.length}/prefix" {input.a}

        """