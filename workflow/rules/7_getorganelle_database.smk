# Fetches GetOrganelle's reference database (v0.0.1) into a fixed "0.0.1"
# directory at the repo root (SeedDatabase/ + LabelDatabase/, one fasta per
# organelle type). rule getorganelle takes it as an input, so it is
# downloaded automatically the first time an organelle sample needs it and
# reused by every later run.
rule getorganelle_database:
    output:
        dir = directory("0.0.1"),
    threads:
        1
    resources:
        mem_mb=2000,
        runtime=60,
    log:
        "logs/getorganelle_database.log"
    shell:
        """
        curl -fsSL https://github.com/Kinggerm/GetOrganelleDB/releases/download/0.0.1/v0.0.1.tar.gz 2> {log} \
            | tar zx 2>> {log}
        """
