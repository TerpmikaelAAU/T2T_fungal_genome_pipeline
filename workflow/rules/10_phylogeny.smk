# ============================================================================
#  Phylogeny: this run's final assemblies among one genome per genus
# ============================================================================
# The last step of the pipeline. Waits for every sample's two deliverables
# (results/<sample>/<selector>/<sample>_<selector>_final.fasta, one per
# SELECTOR) and for the per-genus NCBI reference set (0_4_ncbi_genus_
# representatives.smk, downloaded in parallel with everything else), then
# builds one alignment-free Mash distance tree over all of them with
# mashtree -- the same approach as the BAGS pipeline's phylogeny.
#
# With one reference per genus this places each sample at genus level
# (which genus/family it sits in); it is not meant to resolve species.
#
# Tip labels: references are <Organism_name>_<accession>, samples are
# <sample>_<selector>. results/phylogeny/tip_labels.tsv maps every tip to
# its source and NCBI lineage, e.g. for colouring the tree in iTOL.

def phylogeny_sample_genomes(wildcards):
    return [f"results/{n}/{sel}/{n}_{sel}_final.fasta"
            for n in SAMPLES for sel in SELECTORS]

rule phylogeny_tree:
    input:
        references = "data/phylogeny/reference_genomes",
        ref_table  = "data/phylogeny/genus_representatives_downloaded.tsv",
        samples    = phylogeny_sample_genomes,
    output:
        tree   = "results/phylogeny/genus_tree.nwk",
        matrix = "results/phylogeny/genus_tree_distances.tsv",
        tips   = "results/phylogeny/tip_labels.tsv",
    params:
        stage = "data/phylogeny/tree_input",
    threads:
        32
    resources:
        mem_mb  = resources["mashtree"]["mem_mb"],
        runtime = resources["mashtree"]["runtime"],
    log:
        "logs/phylogeny/phylogeny_tree.log"
    container:
        "docker://quay.io/biocontainers/mashtree:1.4.6--pl5321h87e0c26_4"
    conda:
        "../envs/mashtree.yml"
    shell:
        """
        rm -rf {params.stage}
        mkdir -p {params.stage}/tmp $(dirname {output.tree})
        cp {input.ref_table} {output.tips}

        # Samples: copied in as <sample>_<selector>.fasta so the tip label
        # drops the "_final" suffix. Empty assemblies (possible on tiny or
        # failed inputs) would crash mash, so they're left out with a warning.
        : > {log}
        for f in {input.samples}; do
            label=$(basename "$f" _final.fasta)
            if [ -s "$f" ]; then
                cp "$f" {params.stage}/"$label".fasta
                # Same 12 columns as ref_table: label, (no accession/genus), organism, ...
                printf '%s\\t-\\t-\\t%s (this study)' "$label" "$label" >> {output.tips}
                printf '\\t-%.0s' 1 2 3 4 5 6 7 8 >> {output.tips}
                echo >> {output.tips}
            else
                echo "WARNING: $f is empty -- left out of the tree" >> {log}
            fi
        done

        genomes=$(ls {input.references}/*.fna {params.stage}/*.fasta 2>/dev/null || true)
        count=$(echo $genomes | wc -w)
        echo "Building tree from $count genomes" >> {log}
        if [ "$count" -lt 3 ]; then
            echo "Not enough genomes to build a tree (found $count, need at least 3)." >> {log}
            echo "(Insufficient_Genomes_$count:0.0);" > {output.tree}
            : > {output.matrix}
        else
            mashtree --numcpus {threads} \
                --tempdir {params.stage}/tmp \
                --outmatrix {output.matrix} \
                $genomes > {output.tree} 2>> {log}
        fi
        rm -rf {params.stage}
        """
