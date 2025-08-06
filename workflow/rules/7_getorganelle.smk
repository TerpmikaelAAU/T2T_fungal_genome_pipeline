configfile: "config/config.yaml"

rule getorganelle:
    input:
        a = rules.flye.output.d
    output:
        dir = directory("data/getorganelle/{input}/Mitochondria"),
    threads:
        12
    resources:
        mem_mb=resources["hifiasm"]["mem_mb"],
        runtime=resources["busco"]["time"],
        partition="shared"
    conda:
       "getorganelle.yml"
    shell:
        """
        get_organelle_from_assembly.py -F fungus_mt -g {input.a} --config-dir "0.0.1" -o {output.dir} -t $(nproc) --overwrite

        """