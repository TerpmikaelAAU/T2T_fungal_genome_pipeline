# ============================================================================
#  Annotation (BAGS): geneML genes -> antiSMASH BGCs -> BUSCO on proteins
# ============================================================================
# Only with `bags: enabled: true` (config.yaml). The same steps as the BAGS
# pipeline (github.com/TerpmikaelAAU/BAGS), run on every selector's final
# assembly:
#   1. geneml         -- deep-learning fungal gene prediction (GFF3 + proteins)
#   2. antismash      -- secondary-metabolite gene clusters, using geneML's
#                        genes (--genefinding-tool none)
#   3. busco_proteins -- BUSCO in protein mode: how complete the ANNOTATION
#                        is (the genome-mode BUSCO in 7_BUSCO.smk measures
#                        the assembly)
#   4. bags_summary   -- one row per assembly: genes, BUSCO, BGCs by type;
#                        bags_report collects them into results/annotation_summary.tsv
#
# geneML is only on PyPI (no container), so get_geneml installs it once
# into resources/ inside the official python image -- like get_dorado -- and
# geneml runs from there in that same image. The antiSMASH databases are
# likewise downloaded once into resources/.

GENEML_VERSION = config.get("bags", {}).get("geneml_version", "1.1.0")
GENEML_ENV = f"resources/geneml-{GENEML_VERSION}"
PYTHON_CONTAINER = "docker://python:3.12-slim"
ANTISMASH_CONTAINER = "docker://quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1"


rule get_geneml:
    output:
        geneml = protected(f"{GENEML_ENV}/bin/geneml"),
    params:
        env = GENEML_ENV,
    threads:
        2
    resources:
        mem_mb  = 8000,
        runtime = 120,
    log:
        "logs/bags/get_geneml.log"
    container:
        PYTHON_CONTAINER
    shell:
        """
        rm -rf {params.env}
        python -m venv {params.env} > {log} 2>&1
        {params.env}/bin/pip install --no-cache-dir "geneml=={GENEML_VERSION}" >> {log} 2>&1
        {output.geneml} --version >> {log} 2>&1
        """


rule geneml:
    input:
        geneml = f"{GENEML_ENV}/bin/geneml",
        fa     = "results/{input}/{selector}/{input}_{selector}_final.fasta",
    output:
        gff      = "results/{input}/{selector}/annotation/{input}_{selector}.gff3",
        proteins = "results/{input}/{selector}/annotation/{input}_{selector}.faa",
    threads:
        30
    resources:
        mem_mb  = resources["geneml"]["mem_mb"],
        runtime = resources["geneml"]["runtime"],
    log:
        "logs/bags/{input}_{selector}_geneml.log"
    container:
        PYTHON_CONTAINER
    shell:
        """
        {input.geneml} {input.fa} --output {output.gff} --proteins {output.proteins} \
            --cores {threads} --cpu-only > {log} 2>&1
        """


rule antismash_database:
    output:
        db = protected(directory("resources/antismash_databases")),
    threads:
        4
    resources:
        mem_mb  = 8000,
        runtime = 480,
    log:
        "logs/bags/antismash_database.log"
    container:
        ANTISMASH_CONTAINER
    conda:
        "../envs/antismash.yml"
    shell:
        """
        download-antismash-databases --database-dir {output.db} > {log} 2>&1
        """


rule antismash:
    input:
        fa  = "results/{input}/{selector}/{input}_{selector}_final.fasta",
        gff = "results/{input}/{selector}/annotation/{input}_{selector}.gff3",
        db  = "resources/antismash_databases",
    output:
        dir = directory("results/{input}/{selector}/annotation/antismash"),
    threads:
        8
    resources:
        mem_mb  = resources["antismash"]["mem_mb"],
        runtime = resources["antismash"]["runtime"],
    log:
        "logs/bags/{input}_{selector}_antismash.log"
    container:
        ANTISMASH_CONTAINER
    conda:
        "../envs/antismash.yml"
    shell:
        """
        antismash --cpus {threads} \
            --databases {input.db} \
            --taxon fungi \
            --genefinding-tool none \
            --genefinding-gff3 {input.gff} \
            --cc-mibig --cb-general \
            --allow-long-headers \
            --output-dir {output.dir} \
            {input.fa} > {log} 2>&1
        """


rule busco_proteins:
    input:
        proteins = "results/{input}/{selector}/annotation/{input}_{selector}.faa",
    output:
        summary = "results/{input}/{selector}/annotation/busco_proteins_short_summary.txt",
    params:
        lineage = lambda w: busco_lineage(w.input),
        out     = "data/busco_proteins/{input}_{selector}",
    threads:
        8
    resources:
        mem_mb  = resources["busco"]["mem_mb"],
        runtime = resources["busco"]["runtime"],
    log:
        "logs/bags/{input}_{selector}_busco_proteins.log"
    container:
        "docker://quay.io/biocontainers/busco:6.1.0--pyhdfd78af_2"
    conda:
        "../envs/BUSCO.yml"
    shell:
        """
        rm -rf {params.out}
        busco -i {input.proteins} -o {params.out} -l {params.lineage} \
            -m proteins -c {threads} -f > {log} 2>&1
        cp $(find {params.out} -name 'short_summary*.txt' | head -n 1) {output.summary}
        rm -rf {params.out}
        """


rule bags_summary:
    input:
        gff       = "results/{input}/{selector}/annotation/{input}_{selector}.gff3",
        antismash = "results/{input}/{selector}/annotation/antismash",
        busco     = "results/{input}/{selector}/annotation/busco_proteins_short_summary.txt",
    output:
        a = "results/{input}/{selector}/annotation/annotation_summary.tsv",
    params:
        name = "{input}_{selector}",
    threads:
        1
    resources:
        mem_mb  = 2000,
        runtime = 30,
    script:
        "../scripts/bags_summary.py"


rule bags_report:
    input:
        [f"results/{n}/{sel}/annotation/annotation_summary.tsv"
         for n in SAMPLES for sel in SELECTORS],
    output:
        a = "results/annotation_summary.tsv",
    threads:
        1
    resources:
        mem_mb  = 1000,
        runtime = 10,
    shell:
        """
        head -n 1 {input[0]} > {output.a}
        for f in {input}; do tail -n +2 "$f" >> {output.a}; done
        """
