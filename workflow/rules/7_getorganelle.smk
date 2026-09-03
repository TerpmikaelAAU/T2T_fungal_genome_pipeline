rule getorganelle:
    input:
        a = rules.flye.output.d
    output:
        dir = directory("data/getorganelle/{input}/Mitochondria"),
    params:
        db = lambda w: organelle_db(w.input)
    threads:
        12
    resources:
        mem_mb=resources["hifiasm"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
    conda:
       "../envs/getorganelle.yml"
    shell:
        """
        get_organelle_from_assembly.py -F {params.db} -g {input.a} --config-dir "0.0.1" -o {output.dir} -t $(nproc) --overwrite

        """