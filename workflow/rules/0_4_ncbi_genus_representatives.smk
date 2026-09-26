# ============================================================================
#  Reference genomes for the phylogeny: one NCBI assembly per genus
# ============================================================================
# Downloads ONE representative genome for every genus under
# config['phylogeny']['taxon'] (default: Fungi) with NCBI datasets. These are
# the backbone of the tree built in 10_phylogeny.smk, which places each
# sample's two final assemblies (lowest_contig + highest_busco) among them.
#
# None of these rules depend on any sample, so Snakemake schedules them
# straight away, in parallel with basecalling/assembly -- by the time the
# assemblies are done the reference set is normally already on disk.
#
# Three steps:
#   1. ncbi_genome_summary          -- metadata for every assembly in the
#                                      taxon + the NCBI taxonomy (lineage) of
#                                      every organism in it (no sequences yet)
#   2. pick_genus_representatives   -- best assembly per genus (see
#                                      pick_genus_representatives.py for the
#                                      ranking)
#   3. download_genus_representatives -- fetch just those, renamed to
#                                      readable tree tip labels
#
# DISK: all fungal genera is a few thousand genomes -- expect ~100 GB in
# data/phylogeny/reference_genomes/. It is kept (not temp()) so the download
# only ever happens once; delete it by hand to force a fresh download (e.g.
# after changing `phylogeny.taxon`).
#
# Optional: export NCBI_API_KEY=<your key> before running Snakemake -- the
# datasets CLI picks it up and NCBI allows it more requests per second.

rule ncbi_genome_summary:
    output:
        genomes  = "data/phylogeny/ncbi/genome_summary.jsonl",
        taxonomy = "data/phylogeny/ncbi/taxonomy_summary.jsonl",
    params:
        taxon  = lambda w: config.get("phylogeny", {}).get("taxon", "Fungi"),
        taxids = "data/phylogeny/ncbi/taxids.txt",
    threads:
        1
    resources:
        mem_mb  = resources["ncbi_datasets"]["mem_mb"],
        runtime = 120,
    log:
        "logs/phylogeny/ncbi_genome_summary.log"
    container:
        "docker://staphb/ncbi-datasets:18.37.0"
    conda:
        "../envs/ncbi_datasets.yml"
    shell:
        """
        mkdir -p $(dirname {output.genomes})
        # GenBank only: every RefSeq (GCF_) assembly also exists as a GenBank
        # (GCA_) one, so this avoids counting each genome twice.
        # --exclude-atypical drops MAGs/chimeric/low-quality flagged assemblies.
        datasets summary genome taxon "{params.taxon}" \
            --assembly-source genbank \
            --exclude-atypical \
            --as-json-lines > {output.genomes} 2> {log}
        echo "$(wc -l < {output.genomes}) assemblies found under '{params.taxon}'" >> {log}

        # Lineage of every organism, so genus comes from NCBI taxonomy rather
        # than from guessing at the first word of the organism name. If this
        # lookup fails, pick_genus_representatives falls back to that guess.
        grep -oE '"taxId": ?[0-9]+' {output.genomes} | grep -oE '[0-9]+$' | sort -u > {params.taxids}
        if ! datasets summary taxonomy taxon --inputfile {params.taxids} \
                --as-json-lines > {output.taxonomy} 2>> {log}; then
            echo "WARNING: taxonomy lookup failed; genus will be taken from organism names" >> {log}
            : > {output.taxonomy}
        fi
        rm -f {params.taxids}
        """


rule pick_genus_representatives:
    input:
        genomes  = "data/phylogeny/ncbi/genome_summary.jsonl",
        taxonomy = "data/phylogeny/ncbi/taxonomy_summary.jsonl",
    output:
        table = "data/phylogeny/genus_representatives.tsv",
    threads:
        1
    resources:
        mem_mb  = resources["ncbi_datasets"]["mem_mb"],
        runtime = 30,
    log:
        "logs/phylogeny/pick_genus_representatives.log"
    script:
        "../scripts/pick_genus_representatives.py"


rule download_genus_representatives:
    input:
        table = "data/phylogeny/genus_representatives.tsv",
    output:
        genomes = directory("data/phylogeny/reference_genomes"),
        table   = "data/phylogeny/genus_representatives_downloaded.tsv",
    params:
        tmp = "data/phylogeny/ncbi/download",
    threads:
        8
    resources:
        mem_mb  = resources["ncbi_datasets"]["mem_mb"],
        runtime = resources["ncbi_datasets"]["runtime"],
    log:
        "logs/phylogeny/download_genus_representatives.log"
    container:
        "docker://staphb/ncbi-datasets:18.37.0"
    conda:
        "../envs/ncbi_datasets.yml"
    shell:
        """
        rm -rf {params.tmp}
        mkdir -p {params.tmp} {output.genomes}

        # Column 2 of the table is the accession (column 1 the tip label).
        tail -n +2 {input.table} | cut -f2 > {params.tmp}/accessions.txt
        echo "Downloading $(wc -l < {params.tmp}/accessions.txt) genomes" > {log}

        # Dehydrated download + rehydrate is NCBI's recommended route for
        # large sets: the zip only holds the file list, rehydrate then
        # fetches the sequences in parallel and can resume if interrupted.
        datasets download genome accession \
            --inputfile {params.tmp}/accessions.txt \
            --include genome \
            --dehydrated \
            --filename {params.tmp}/representatives.zip 2>> {log}
        unzip -q -o {params.tmp}/representatives.zip -d {params.tmp}
        datasets rehydrate --directory {params.tmp} --max-workers {threads} 2>> {log}

        # Rename <accession>/*.fna to <tip label>.fna -- mashtree uses the file
        # name as the tip label. Genomes NCBI couldn't serve are skipped with
        # a warning instead of failing the whole set.
        head -n 1 {input.table} > {output.table}
        tail -n +2 {input.table} | while IFS=$'\\t' read -r label acc rest; do
            fna=$(find {params.tmp}/ncbi_dataset/data/"$acc" -name '*.fna' 2>/dev/null | head -n 1 || true)
            if [ -n "$fna" ] && [ -s "$fna" ]; then
                mv "$fna" {output.genomes}/"$label".fna
                printf '%s\\t%s\\t%s\\n' "$label" "$acc" "$rest" >> {output.table}
            else
                echo "WARNING: no genome downloaded for $acc ($label) -- left out of the tree" >> {log}
            fi
        done
        echo "$(tail -n +2 {output.table} | wc -l) genomes ready in {output.genomes}" >> {log}
        rm -rf {params.tmp}
        """
