# ============================================================================
#  Fetch the dorado binary
# ============================================================================
# Dorado is not on conda, so it used to be a manual install with an absolute
# path hardcoded in config.yaml. This rule fetches a pinned version instead,
# so a fresh clone of the pipeline is self-contained.
#
# To use a dorado you installed yourself, set an absolute path in config:
#     dorado:
#       path: /my/own/dorado/bin/dorado
# and this rule will never run.

rule get_dorado:
    output:
        binary = protected(f"resources/dorado-{DORADO_VERSION}-linux-x64/bin/dorado")
    params:
        url = ("https://cdn.oxfordnanoportal.com/software/analysis/"
               f"dorado-{DORADO_VERSION}-linux-x64.tar.gz"),
        outdir = "resources",
        tarball = f"resources/dorado-{DORADO_VERSION}-linux-x64.tar.gz",
    threads:
        2
    resources:
        # ~3.5 GB download, ~8 GB unpacked
        mem_mb = 4000,
        runtime = 120,
        disk_mb = 20000,
    log:
        "logs/get_dorado.log"
    shell:
        """
        mkdir -p {params.outdir}
        echo "Fetching {params.url}" > {log}
        curl -fsSL --retry 3 --retry-delay 10 "{params.url}" -o "{params.tarball}" 2>> {log}
        tar -xzf "{params.tarball}" -C "{params.outdir}" 2>> {log}
        rm -f "{params.tarball}"
        chmod +x "{output.binary}"
        "{output.binary}" --version >> {log} 2>&1
        """