# ============================================================================
#  Contamination screen: Tiara on every final assembly
# ============================================================================
# With `contamination: enabled: true` (the default). Tiara (Karlicki et al.
# 2022) is a small deep-learning classifier that labels each contig as
# eukarya, bacteria, archaea, organelle (mitochondrion/plastid) or unknown
# -- enough to spot bacterial contigs in a fungal culture's assembly. Its
# models ship with the tool, so unlike NCBI FCS-GX there is no ~500 GB
# database, and it runs in minutes.
#
# Report only: nothing is removed. results/{sample}/{selector}/contamination/
# tiara.tsv has every contig's class and probabilities, and summary.txt gets
# a CONTAMINATION section listing any bacterial/archaeal contigs. Contigs
# shorter than `min_len` aren't classified (Tiara's own limit).
#
# Tiara isn't on bioconda (no biocontainer) and pins old dependencies
# (Python <= 3.9, torch 1.7), so get_tiara installs it once into resources/
# inside the official python:3.9-slim image -- like get_geneml/get_dorado.

TIARA_VERSION = "1.0.3"
TIARA_ENV = f"resources/tiara-{TIARA_VERSION}"
TIARA_CONTAINER = "docker://python:3.9-slim"


rule get_tiara:
    output:
        tiara = protected(f"{TIARA_ENV}/bin/tiara"),
    params:
        env = TIARA_ENV,
    threads:
        2
    resources:
        mem_mb  = 8000,
        runtime = 120,
    log:
        "logs/contamination/get_tiara.log"
    container:
        TIARA_CONTAINER
    shell:
        """
        rm -rf {params.env}
        python -m venv {params.env} > {log} 2>&1
        # setuptools<70: Tiara's pinned skorch 0.9 still imports pkg_resources.
        {params.env}/bin/pip install --no-cache-dir \
            "tiara=={TIARA_VERSION}" "setuptools<70" >> {log} 2>&1
        {output.tiara} --help >> {log} 2>&1
        """


rule tiara:
    input:
        tiara = f"{TIARA_ENV}/bin/tiara",
        fa    = "results/{input}/{selector}/{input}_{selector}_final.fasta",
    output:
        tsv = "results/{input}/{selector}/contamination/tiara.tsv",
    params:
        min_len = CONTAM.get("min_len", 3000),
    threads:
        8
    resources:
        mem_mb  = 16000,
        runtime = 240,
    log:
        "logs/contamination/{input}_{selector}_tiara.log"
    container:
        TIARA_CONTAINER
    shell:
        """
        if [ -s {input.fa} ]; then
            {input.tiara} -i {input.fa} -o {output.tsv} -t {threads} \
                -m {params.min_len} --probabilities > {log} 2>&1
        else
            echo "Empty assembly -- nothing to classify" > {log}
            printf 'sequence_id\\tclass_fst_stage\\tclass_snd_stage\\n' > {output.tsv}
        fi
        """
