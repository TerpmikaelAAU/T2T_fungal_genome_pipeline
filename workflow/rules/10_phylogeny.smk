# ============================================================================
#  Phylogeny: this run's final assemblies among one genome per genus (UFCG)
# ============================================================================
# A protein phylogeny from UFCG's universal fungal core genes
# (https://ufcg.steineggerlab.com, Kim et al. 2023, NAR) -- the same tool
# used for the paper (scripts/UFCG_Phylogeny/). Distance-based trees (Mash)
# saturate between genera, so they can't resolve a tree spanning all Fungi.
#
#   1. ufcg_profile_reference -- core genes of each NCBI genus representative
#                                (0_4_ncbi_genus_representatives.smk). One job
#                                per genome, no sample dependency: these run
#                                in parallel with the assemblies, and are kept
#                                so they only ever run once.
#   2. ufcg_profile_sample    -- the same for each sample's two final
#                                assemblies (one per SELECTOR), once the
#                                pipeline has picked and polished them.
#   3. ufcg_align             -- align every core gene (MAFFT) and concatenate
#   4. phylogeny_tree         -- FastTree (LG+gamma) on the concatenation
#
# Tip labels: references are <Organism_name>_<accession>, samples are
# <sample>_<selector>. results/phylogeny/tip_labels.tsv maps every tip to
# its NCBI lineage, e.g. for colouring the tree in iTOL.
#
# A genome UFCG can't profile (or an empty assembly) gets an empty .ucg
# placeholder and is left out of the tree with a warning, instead of one bad
# genome out of thousands blocking the whole tree.

# Official UFCG image: unlike the bioconda one it ships UFCG's core gene
# database, which can't be downloaded into a read-only container at runtime.
UFCG_CONTAINER = "docker://endix1029/ufcg:v1.0.6"

def ufcg_reference_profiles(wildcards):
    """One profile per reference genome that actually downloaded -- only
    known once the download checkpoint has run."""
    table = checkpoints.download_genus_representatives.get().output.table
    with open(table) as fh:
        next(fh)  # header
        refs = [line.split("\t", 1)[0] for line in fh if line.strip()]
    return expand("data/phylogeny/ufcg_profiles/references/{ref}.ucg", ref=refs)

def ufcg_sample_profiles(wildcards):
    return [f"data/phylogeny/ufcg_profiles/samples/{n}_{sel}.ucg"
            for n in SAMPLES for sel in SELECTORS]


rule ufcg_profile_reference:
    input:
        "data/phylogeny/reference_genomes/{ref}.fna"
    output:
        "data/phylogeny/ufcg_profiles/references/{ref}.ucg"
    params:
        work = "data/phylogeny/ufcg_work/references/{ref}",
    threads:
        8
    resources:
        mem_mb  = resources["ufcg_profile"]["mem_mb"],
        runtime = resources["ufcg_profile"]["runtime"],
    log:
        "logs/phylogeny/ufcg_profile/{ref}.log"
    container:
        UFCG_CONTAINER
    conda:
        "../envs/ufcg.yml"
    shell:
        """
        rm -rf {params.work}
        mkdir -p {params.work}/out {params.work}/tmp
        # Label = file name minus .fna, so the .ucg is named after it too.
        if ufcg profile -i {input} -o {params.work}/out -w {params.work}/tmp \
                -t {threads} -f --nocolor > {log} 2>&1 \
                && [ -s {params.work}/out/{wildcards.ref}.ucg ]; then
            mv {params.work}/out/{wildcards.ref}.ucg {output}
        else
            echo "WARNING: UFCG profile failed -- {wildcards.ref} left out of the tree" >> {log}
            : > {output}
        fi
        rm -rf {params.work}
        """


rule ufcg_profile_sample:
    input:
        "results/{input}/{selector}/{input}_{selector}_final.fasta"
    output:
        "data/phylogeny/ufcg_profiles/samples/{input}_{selector}.ucg"
    params:
        work = "data/phylogeny/ufcg_work/samples/{input}_{selector}",
    threads:
        8
    resources:
        mem_mb  = resources["ufcg_profile"]["mem_mb"],
        runtime = resources["ufcg_profile"]["runtime"],
    log:
        "logs/phylogeny/ufcg_profile/{input}_{selector}.log"
    container:
        UFCG_CONTAINER
    conda:
        "../envs/ufcg.yml"
    shell:
        """
        rm -rf {params.work}
        mkdir -p {params.work}/in {params.work}/out {params.work}/tmp
        : > {log}
        if [ ! -s {input} ]; then
            echo "WARNING: {input} is empty -- left out of the tree" >> {log}
            : > {output}
        else
            # Copied in as <sample>_<selector>.fasta: that becomes the label.
            cp {input} {params.work}/in/{wildcards.input}_{wildcards.selector}.fasta
            if ufcg profile -i {params.work}/in/{wildcards.input}_{wildcards.selector}.fasta \
                    -o {params.work}/out -w {params.work}/tmp \
                    -t {threads} -f --nocolor >> {log} 2>&1 \
                    && [ -s {params.work}/out/{wildcards.input}_{wildcards.selector}.ucg ]; then
                mv {params.work}/out/{wildcards.input}_{wildcards.selector}.ucg {output}
            else
                echo "WARNING: UFCG profile failed -- left out of the tree" >> {log}
                : > {output}
            fi
        fi
        rm -rf {params.work}
        """


rule ufcg_align:
    input:
        references = ufcg_reference_profiles,
        samples    = ufcg_sample_profiles,
        ref_table  = "data/phylogeny/genus_representatives_downloaded.tsv",
    output:
        alignment = "results/phylogeny/ufcg_concatenated_alignment.fasta",
        tips      = "results/phylogeny/tip_labels.tsv",
    params:
        work = "data/phylogeny/ufcg_work/align",
    threads:
        32
    resources:
        mem_mb  = resources["ufcg_align"]["mem_mb"],
        runtime = resources["ufcg_align"]["runtime"],
    log:
        "logs/phylogeny/ufcg_align.log"
    container:
        UFCG_CONTAINER
    conda:
        "../envs/ufcg.yml"
    shell:
        """
        rm -rf {params.work}
        mkdir -p {params.work}/profiles $(dirname {output.alignment})
        : > {log}

        # Only non-empty profiles go in (empty = placeholder, see header).
        for f in {input.references} {input.samples}; do
            if [ -s "$f" ]; then
                cp "$f" {params.work}/profiles/
            else
                echo "WARNING: $(basename "$f" .ucg) has no UFCG profile -- left out" >> {log}
            fi
        done

        # tip_labels.tsv: the reference table, then one row per sample tip,
        # for exactly the genomes that made it into the alignment.
        ls {params.work}/profiles | sed 's/\\.ucg$//' > {params.work}/labels.txt
        head -n 1 {input.ref_table} > {output.tips}
        awk -F'\t' 'NR == FNR {{ keep[$1]; next }} FNR > 1 && ($1 in keep)' \
            {params.work}/labels.txt {input.ref_table} >> {output.tips}
        for f in {input.samples}; do
            label=$(basename "$f" .ucg)
            if [ -s "$f" ]; then
                # Same 12 columns as ref_table: label, (no accession/genus), organism, ...
                printf '%s\\t-\\t-\\t%s (this study)' "$label" "$label" >> {output.tips}
                printf '\\t-%.0s' 1 2 3 4 5 6 7 8 >> {output.tips}
                echo >> {output.tips}
            fi
        done

        count=$(wc -l < {params.work}/labels.txt)
        echo "Aligning core genes of $count genomes" >> {log}
        if [ "$count" -lt 4 ]; then
            echo "Not enough genomes for a tree (found $count, need at least 4)." >> {log}
            : > {output.alignment}
        else
            # Plain `mafft` (FFT-NS-2) rather than UFCG's default mafft-linsi:
            # L-INS-i is far too slow for thousands of sequences per gene.
            ufcg align -i {params.work}/profiles -o {params.work}/align \
                -a protein -t {threads} -m mafft --nocolor >> {log} 2>&1
            mv {params.work}/align/aligned_concatenated.fasta {output.alignment}
        fi
        # Every profile job is done by now; drop their (empty) parent dirs too.
        rm -rf $(dirname {params.work})
        """


rule phylogeny_tree:
    input:
        alignment = "results/phylogeny/ufcg_concatenated_alignment.fasta",
    output:
        tree = "results/phylogeny/genus_tree.nwk",
    threads:
        32
    resources:
        mem_mb  = resources["fasttree"]["mem_mb"],
        runtime = resources["fasttree"]["runtime"],
    log:
        "logs/phylogeny/phylogeny_tree.log"
    container:
        "docker://quay.io/biocontainers/fasttree:2.2.0--h7b50bb2_1"
    conda:
        "../envs/fasttree.yml"
    shell:
        """
        if [ ! -s {input.alignment} ]; then
            echo "Empty alignment -- see logs/phylogeny/ufcg_align.log" > {log}
            echo "(Insufficient_Genomes:0.0);" > {output.tree}
        else
            # FastTreeMP is the multithreaded build; fall back to FastTree.
            FT=$(command -v FastTreeMP || command -v FastTree)
            OMP_NUM_THREADS={threads} "$FT" -lg -gamma {input.alignment} > {output.tree} 2> {log}
        fi
        """
