# ============================================================================
#  dorado correct, split into blocks and across CPU/GPU
# ============================================================================
# Replaces the old single monolithic `dorado correct` job, which needed ~200 GB
# on the one GPU node and could not be resumed if it died.
#
#   correct_num_blocks : asks dorado how many blocks the input splits into.
#                        Also BUILDS THE FASTQ INDEX -- this must happen once,
#                        before any block job starts, or the block jobs race
#                        each other trying to create it.
#   correct_overlap    : all-vs-all overlaps -> PAF. CPU only, big memory.
#                        Runs on the fat zen3x/zen5x nodes, in parallel.
#   correct_inference  : reads the PAF back and corrects. GPU (or CPU).
#   correct_merge      : concatenates the per-block FASTA.

checkpoint correct_num_blocks:
    input:
        a = get_correct_input
    output:
        a = "data/dorado/{input}/num_blocks.txt"
    threads:
        16
    resources:
        mem_mb = resources["correct_blocks"]["mem_mb"],
        runtime = resources["correct_blocks"]["runtime"],
    log:
        "logs/correct_num_blocks/{input}.log"
    shell:
        """
        "{config[dorado]}" correct {input.a} \
            --index-size {config[dorado_correct][index_size]} \
            --compute-num-blocks > {output.a} 2> {log}
        """


rule correct_overlap:
    input:
        a  = get_correct_input,
        nb = "data/dorado/{input}/num_blocks.txt",
    output:
        a = temp("data/dorado/{input}/block_{block}.paf")
    threads:
        64
    resources:
        mem_mb = scaled_mem(2.5, 64000),
        runtime = scaled_time(0.08, 720),
    log:
        "logs/correct_overlap/{input}_block{block}.log"
    shell:
        """
        "{config[dorado]}" correct {input.a} \
            --index-size {config[dorado_correct][index_size]} \
            --run-block-id {wildcards.block} \
            --to-paf \
            --device cpu \
            --threads {threads} > {output.a} 2> {log}
        """


rule correct_inference:
    input:
        a   = get_correct_input,
        paf = "data/dorado/{input}/block_{block}.paf",
    output:
        a = temp("data/dorado/{input}/block_{block}.fasta")
    threads:
        16
    resources:
        mem_mb = scaled_mem(0.5, 32000, GPU_MEM_CAP),
        runtime = scaled_time(0.05, 480),
        **({"gres": INFER_GRES} if INFER_GRES else {}),
    log:
        "logs/correct_inference/{input}_block{block}.log"
    shell:
        """
        "{config[dorado]}" correct {input.a} \
            --from-paf {input.paf} \
            --device {INFER_DEVICE} > {output.a} 2> {log}
        """


def correct_blocks(wildcards):
    """Expand to one FASTA per block, once the checkpoint has told us how many."""
    nb_file = checkpoints.correct_num_blocks.get(**wildcards).output.a
    with open(nb_file) as fh:
        n = int(fh.read().strip())
    if n > MAX_BLOCKS:
        raise WorkflowError(
            f"dorado reported {n} blocks for '{wildcards.input}', above "
            f"max_blocks={MAX_BLOCKS}. Increase dorado_correct.index_size "
            f"(fewer, larger blocks) or raise max_blocks if this is expected."
        )
    return expand("data/dorado/{sample}/block_{block}.fasta",
                  sample=wildcards.input, block=range(n))


rule correct_merge:
    input:
        correct_blocks
    output:
        a = temp("data/dorado/{input}.fasta")
    threads:
        2
    resources:
        mem_mb = resources["correct_merge"]["mem_mb"],
        runtime = resources["correct_merge"]["runtime"],
    shell:
        """
        cat {input} > {output.a}
        """