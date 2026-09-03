import os
from snakemake.utils import min_version

min_version("9.0")

configfile: "config/config.yaml"

# ============================================================================
#  Sample handling / entry points
# ============================================================================
SAMPLES = config["samples"]
VALID_TYPES = {"pod5", "bam", "fastq"}

# Fail loudly and early on a malformed config rather than deep in the DAG.
for _name, _s in SAMPLES.items():
    if "type" not in _s or "path" not in _s:
        raise WorkflowError(f"Sample '{_name}': both 'type' and 'path' are required.")
    if _s["type"] not in VALID_TYPES:
        raise WorkflowError(
            f"Sample '{_name}': type '{_s['type']}' is not one of {sorted(VALID_TYPES)}."
        )
    if not os.path.isabs(_s["path"]):
        raise WorkflowError(f"Sample '{_name}': 'path' must be absolute.")

def stype(name):
    return SAMPLES[name]["type"]

def has_bam(name):
    """Polishing needs a BAM: available for pod5 (we make one) and bam entries."""
    return stype(name) in ("pod5", "bam")

def subsampled(name):
    return bool(SAMPLES[name].get("subsample"))

def busco_lineage(name):
    return SAMPLES[name].get("busco_lineage", "fungi_odb10")

def organelle_db(name):
    return SAMPLES[name].get("organelle", "fungus_mt")

def wants_organelle(name):
    return str(organelle_db(name)).lower() != "none"

# --- input resolvers used by the rules -------------------------------------
def get_bam(wildcards):
    """BAM for a sample: basecalled by us, or supplied by the user."""
    s = SAMPLES[wildcards.input]
    return (f"data/dorado_basecall/{wildcards.input}.bam"
            if s["type"] == "pod5" else s["path"])

def get_raw_fastq(wildcards):
    """FASTQ entering porechop: converted from BAM, or supplied by the user."""
    s = SAMPLES[wildcards.input]
    return (s["path"] if s["type"] == "fastq"
            else f"data/samtools/Fastq/{wildcards.input}.fastq")

def get_correct_input(wildcards):
    """Reads fed to dorado correct: coverage-capped only if the user asked."""
    n = wildcards.input
    return (f"data/rasusa/Coverage/{n}.fastq" if subsampled(n)
            else f"data/chopper/L10kbQ10/{n}.fastq")


# ============================================================================
#  Adaptive resources
# ============================================================================
# Two problems these solve:
#   1. Fixed budgets tuned for 40 Mb fungal genomes are wrong for a 240 Mb
#      oomycete with 100+ GB of reads.
#   2. A job killed by the OOM reaper should come back with MORE memory,
#      not the same amount again.
#
# `input.size_mb` is the on-disk size of the rule's inputs, so memory tracks
# the actual data. `attempt` starts at 1 and increments on each retry, so a
# killed job doubles its request. Set `retries:` in the profile to enable.
#
# CAVEAT: for gzipped inputs, size_mb is the COMPRESSED size. The factors
# below already account for roughly 3x expansion on .gz reads.

CPU_MEM_CAP = 1_500_000   # zen5x has 2.3 TB; stay well under
GPU_MEM_CAP = 240_000     # bio-node10 has 256 TB total
MAX_RUNTIME = 10080       # 7 days, BioCloud's ceiling

def scaled_mem(factor, floor_mb, cap_mb=CPU_MEM_CAP):
    """mem_mb = max(floor, input_size * factor), doubling on each retry."""
    def _f(wildcards, input, attempt):
        try:
            base = max(floor_mb, int(input.size_mb * factor))
        except (AttributeError, TypeError):
            base = floor_mb
        return int(min(cap_mb, base * (2 ** (attempt - 1))))
    return _f

def scaled_time(factor, floor_min, cap_min=MAX_RUNTIME):
    """runtime in minutes, growing 1.5x on each retry (walltime kills too)."""
    def _f(wildcards, input, attempt):
        try:
            base = max(floor_min, int(input.size_mb * factor))
        except (AttributeError, TypeError):
            base = floor_min
        return int(min(cap_min, base * (1.5 ** (attempt - 1))))
    return _f

# ============================================================================
#  Resources (fixed floors; see scaled_mem/scaled_time above)
#  runtime is an INTEGER NUMBER OF MINUTES (Snakemake >=8).
#  Partitions are NOT set: BioCloud assigns them from the mem-per-CPU ratio.
# ============================================================================
resources = {
    "hifiasm":          {"mem_mb": 30000,  "runtime": 240},
    "flye":             {"mem_mb": 30000,  "runtime": 240},
    "busco":            {"mem_mb": 10000,  "runtime": 480},
    "rasusa":           {"mem_mb": 5000,   "runtime": 40},
    "chopper":          {"mem_mb": 5000,   "runtime": 40},
    "fga":              {"mem_mb": 5000,   "runtime": 60},
    "porechop_api":     {"mem_mb": 100000, "runtime": 600},
    "dorado_basecall":  {"mem_mb": 100000, "runtime": 10080},
    "seqkit":           {"mem_mb": 10000,  "runtime": 60},
    "dorado_align":     {"mem_mb": 75000,  "runtime": 720},
    "dorado_polish":    {"mem_mb": 100000, "runtime": 600},
    # dorado correct, split into its three stages
    "correct_blocks":   {"mem_mb": 32000,  "runtime": 120},
    "correct_overlap":  {"mem_mb": 250000, "runtime": 720},
    "correct_infer":    {"mem_mb": 60000,  "runtime": 480},
    "correct_merge":    {"mem_mb": 5000,   "runtime": 60},
}

# GPU request for BioCloud's single A10 node (bio-node10).
GPU_GRES = "gpu:a10:1"

# Inference may be pinned to CPU via config, in which case it must NOT ask
# for the GPU -- otherwise it queues behind GPU work for no reason.
INFER_DEVICE = config["dorado_correct"].get("inference_device", "cuda:all")
INFER_GRES = "" if INFER_DEVICE == "cpu" else GPU_GRES

# Safety valve: a typo in index_size (e.g. "4" instead of "4G") could make
# dorado report a huge block count and flood the scheduler. Snakemake never
# submits more than `jobs:` at once, but this fails fast and loudly instead.
MAX_BLOCKS = config["dorado_correct"].get("max_blocks", 200)

include: "workflow/rules/0_1_basecall.smk"
include: "workflow/rules/0_2_bam_to_fastq.smk"
include: "workflow/rules/1_porechop_api.smk"
include: "workflow/rules/2_chopper_dorado_correction.smk"
include: "workflow/rules/2_chopper_Flye_mitochondria.smk"
include: "workflow/rules/2_chopper_ultralong.smk"
include: "workflow/rules/2_flye.smk"
include: "workflow/rules/2_rasusa.smk"
include: "workflow/rules/3_dorado_correct.smk"
include: "workflow/rules/3_seqtk_fasta_to_fastq.smk"
include: "workflow/rules/4_seqkit_10kb_50kb.smk"
include: "workflow/rules/5_hifiasm.smk"
include: "workflow/rules/6_0_fga_to_fa.smk"
include: "workflow/rules/6_00_lowest_contig_count.smk"
include: "workflow/rules/6_1_dorado_align.smk"
include: "workflow/rules/6_3_dorado_polish.smk"
include: "workflow/rules/7_BUSCO.smk"
include: "workflow/rules/7_getorganelle_database.smk"
include: "workflow/rules/7_getorganelle.smk"
include: "workflow/rules/9_stats.smk"

# ============================================================================
#  Targets -- built per sample according to its entry point
# ============================================================================
def final_targets():
    t = []
    for n in SAMPLES:
        # Final assembly: polished where a BAM exists, otherwise the best raw one
        if has_bam(n):
            t.append(f"data/dorado_polish/{n}.fasta")
        else:
            t.append(f"data/contig/{n}_lowest_contig_file.fa")
        t.append(f"data/busco/{n}/BUSCO")
        t.append(f"results/{n}/summary.txt")
        t.append(f"results/{n}/{n}_final.fasta")
        if wants_organelle(n):
            t.append(f"data/getorganelle/{n}/Mitochondria")
    return t

rule all:
    input:
        final_targets()