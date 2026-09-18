# Mines flye's assembly graph for an organelle (typically mitochondrial)
# contig. Only built when a sample sets `organelle: true` in config.yaml
# (wants_organelle()); `db` is that sample's required `organelle_type`
# (GetOrganelle's -F flag, e.g. fungus_mt).
# NOTE: --config-dir "0.0.1" is NOT a declared Snakemake input -- the
# "0.0.1" directory must already exist (run rule getorganelle_database
# first) or this fails with a missing-database error the DAG won't predict.
# NOTE: this output isn't copied into results/{sample}/ -- it's the only
# copy of the organelle assembly, so don't temp() it.
rule getorganelle:
    input:
        a = rules.flye.output.d
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
        get_organelle_from_assembly.py -F {params.db} -g {input.a} --config-dir "0.0.1" -o {output.dir} -t $(nproc) --overwrite

        """