# One-time fetch of GetOrganelle's reference database (v0.0.1) into a fixed
# "0.0.1" directory at the repo root. Not wired as a Snakemake dependency of
# rule getorganelle (which reads it via --config-dir without declaring it as
# input) -- run this manually once, before the first getorganelle build.
rule getorganelle_database:
    output:
        dir = directory("0.0.1"),
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["runtime"],
    shell:
        """
        curl -L https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz | tar zx

        """