configfile: "config/config.yaml"

rule getorganelle_database:
    output:
        dir = directory("0.0.1"),
    threads:
        12
    resources:
        mem_mb=resources["busco"]["mem_mb"],
        runtime=resources["busco"]["time"],
        partition="shared"
    shell:
        """
        curl -L https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz | tar zx

        """