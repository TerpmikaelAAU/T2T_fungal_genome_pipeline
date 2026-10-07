# Mines flye's assembly graph for an organelle (typically mitochondrial)
# contig. Only built when a sample sets `organelle: true` in config.yaml
# (wants_organelle()); `db` is that sample's required `organelle_type`
# (GetOrganelle's -F flag, e.g. fungus_mt).
# --config-dir is GetOrganelle's reference database, downloaded once by
# rule getorganelle_database (7_getorganelle_database.smk).
# NOTE: this output isn't copied into results/{sample}/ -- it's the only
# copy of the organelle assembly, so don't temp() it.
rule getorganelle:
    input:
        a = rules.flye.output.d,
        db = rules.getorganelle_database.output.dir,
    output:
        dir = directory("data/getorganelle/{input}/Mitochondria"),
    params:
        db = lambda w: organelle_type(w.input)
    threads:
        12
    resources:
        mem_mb=resources["hifiasm"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
    container:
        "docker://quay.io/biocontainers/getorganelle:1.7.7.1--pyhdfd78af_0"
    conda:
       "../envs/getorganelle.yml"
    shell:
        """
        get_organelle_from_assembly.py -F {params.db} -g {input.a} --config-dir "{input.db}" -o {output.dir} -t $(nproc) --overwrite

        """