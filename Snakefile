import os
from snakemake.utils import min_version

min_version("9.0")

configfile: "config/config.yaml"

# ============================================================================
#  T2T fungal (and, as of v2, non-fungal) genome assembly pipeline
# ============================================================================
# Structure of this file, top to bottom:
#   1. Sample handling       -- validates config['samples'], entry-point
#                                helpers (pod5/bam/fastq), per-sample option
#                                accessors (busco_lineage, organelle, filter
#                                grid override, ...)
#   2. Input resolvers        -- functions that pick the right upstream file
#                                for a rule given a sample's entry point and
#                                config (get_raw_fastq, get_assembly_input,
#                                hifiasm_inputs, ...)
#   3. Dorado binary          -- resolves a fetched-or-user-supplied dorado
#   4. Adaptive resources     -- scaled_mem/scaled_time size requests off
#                                actual input size and escalate on retry
#   5. Rule includes          -- one file per pipeline stage, under
#                                workflow/rules/, numbered in DAG order
#   6. Targets (`rule all`)   -- assembled per sample according to its entry
#                                point and options (see final_targets())
#
# See config/config.yaml for what's configurable and README.md for how to
# run this.

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
        # A relative path (e.g. "example_data/...") resolves against this
        # repo's root, so the bundled example data works after any clone.
        _s["path"] = os.path.join(workflow.basedir, _s["path"])
    if _s.get("organelle") and not _s.get("organelle_type"):
        raise WorkflowError(
            f"Sample '{_name}': organelle: true requires 'organelle_type' "
            f"(e.g. fungus_mt) -- see config/config.yaml."
        )

# minq/minlen identify one cell of the read-filtering grid (see config.yaml
# `filter:`); block identifies one dorado-correct chunk. Constraining them to
# digits keeps "{input}_q10_l5000" parsing as input="{input}", minq="10",
# minlen="5000" instead of Snakemake trying every split. `selector` picks
# which grid-winner strategy a downstream result came from (see SELECTORS
# below and 6_01_highest_busco.smk).
wildcard_constraints:
    minq     = r"\d+",
    minlen   = r"\d+",
    block    = r"\d+",
    selector = r"lowest_contig|highest_busco",

# Two independent ways to pick a winner across a sample's (min_q, min_len)
# grid: fewest contigs (contig_count, 6_00) or highest BUSCO completeness
# (highest_busco, 6_01). Both go through the same rest of the pipeline
# (polish if possible, BUSCO on the final assembly, stats, results/) --
# see final_targets() below -- so every sample ends up with two parallel
# results/{sample}/<selector>/ deliverables to compare.
SELECTORS = ["lowest_contig", "highest_busco"]

def stype(name):
    return SAMPLES[name]["type"]

def has_bam(name):
    """Polishing needs a BAM: available for pod5 (we make one) and bam entries."""
    return stype(name) in ("pod5", "bam")

def subsampled(name):
    return bool(SAMPLES[name].get("subsample"))

def busco_lineage(name):
    return SAMPLES[name].get("busco_lineage", "fungi_odb12")

def wants_organelle(name):
    return bool(SAMPLES[name].get("organelle", False))

def organelle_type(name):
    """GetOrganelle's -F organelle-type flag; required when organelle: true."""
    return SAMPLES[name]["organelle_type"]

def wants_dorado_correct(name):
    """Per-sample override of the global dorado_correct.enabled default."""
    return bool(SAMPLES[name].get("dorado_correct", config["dorado_correct"]["enabled"]))

def wants_ultralong(name):
    """Per-sample override of the global ultralong.enabled default."""
    return bool(SAMPLES[name].get("ultralong", config["ultralong"]["enabled"]))

def filter_grid(name):
    """This sample's (min_q, min_len) grid: a per-sample override if given,
    else the global config['filter'] default.

    config['filter'] is shared by every sample unless overridden here, so
    widening the grid for one sample used to silently change what
    contig_count/assembly_stats required from every other sample too --
    including combinations never actually run for them.

    When dorado_correct is on for this sample, dorado correct itself always
    runs once on a single fixed (min_q, min_len) cutoff -- see
    config.yaml `dorado_correct.min_q`/`min_len` -- independent of this
    grid. AFTER correction, min_len is swept again (min_q isn't: correction
    discards real quality, so there's nothing left to filter on); min_q here
    is just dorado_correct's fixed value, kept as a folder-naming tag."""
    grid = SAMPLES[name].get("filter", config["filter"])
    if wants_dorado_correct(name):
        return {"min_q": [config["dorado_correct"]["min_q"]], "min_len": grid["min_len"]}
    return grid

# --- input resolvers used by the rules -------------------------------------
def get_bam(wildcards):
    """BAM for a sample: basecalled by us, or supplied by the user."""
    s = SAMPLES[wildcards.input]
    return (f"data/dorado_basecall/{wildcards.input}.bam"
            if s["type"] == "pod5" else s["path"])

def is_gz_fastq(name):
    s = SAMPLES[name]
    return s["type"] == "fastq" and s["path"].endswith(".gz")

def get_raw_fastq(wildcards):
    """Plain-text FASTQ entering porechop/chopper: converted from BAM, decompressed
    once if the user's fastq.gz would otherwise be re-zcat'd by every grid point, or
    the raw path as-is when it's already plain text."""
    n = wildcards.input
    s = SAMPLES[n]
    if is_gz_fastq(n):
        return f"data/decompressed/{n}.fastq"
    return s["path"] if s["type"] == "fastq" else f"data/samtools/Fastq/{n}.fastq"

def trim_adapters(name):
    """Adapter trimming is opt-out per sample (default: on)."""
    return SAMPLES[name].get("porechop", True)

def get_trimmed_fastq(wildcards):
    """Porechop output, or the raw reads when trimming is disabled."""
    n = wildcards.input
    return (f"data/porechopped/{n}.fastq" if trim_adapters(n)
            else get_raw_fastq(wildcards))

def get_correct_input(wildcards):
    """Reads for one (sample, min_q, min_len) grid cell, coverage-capped
    only if the user asked. Used directly as hifiasm's input when
    dorado_correct is off for this sample."""
    n, q, l = wildcards.input, wildcards.minq, wildcards.minlen
    return (f"data/rasusa/Coverage/{n}_q{q}_l{l}.fastq" if subsampled(n)
            else f"data/chopper/{n}_q{q}_l{l}.fastq")

def get_dorado_correct_input(wildcards):
    """Reads fed to dorado correct: this sample's full read set, filtered
    ONCE at dorado_correct's own fixed (min_q, min_len) cutoff -- see
    config.yaml -- independent of the assembly grid. Coverage-capped only
    if the user asked."""
    n = wildcards.input
    dc = config["dorado_correct"]
    q, l = dc["min_q"], dc["min_len"]
    return (f"data/rasusa/Coverage/{n}_q{q}_l{l}.fastq" if subsampled(n)
            else f"data/chopper/{n}_q{q}_l{l}.fastq")

def get_assembly_input(wildcards):
    """FASTQ fed to hifiasm for one grid cell: hifiasm --ont does its own
    ONT-specific correction, and on at least one large low-N50 dataset
    dorado correct discarded the large majority of reads and made the
    assembly worse (see config.yaml `dorado_correct`), so that's off by
    default.

    When dorado_correct IS on for this sample, hifiasm instead gets the
    single dorado-corrected read set, filtered again by min_len only (see
    3_2_length_filter_corrected.smk) -- min_q plays no further part once
    correction has discarded the reads' real quality."""
    n = wildcards.input
    if wants_dorado_correct(n):
        return f"data/dorado_filtered/{n}_l{wildcards.minlen}.fastq"
    return get_correct_input(wildcards)

def get_ultralong_input(wildcards):
    """hifiasm's --ul input: a fixed cutoff independent of the filter grid,
    reusing the same chopper rule. Only built when ultralong is enabled for
    this sample."""
    u = config["ultralong"]
    return f"data/chopper/{wildcards.input}_q{u['min_q']}_l{u['min_len']}.fastq"

def hifiasm_inputs(wildcards):
    d = {"a": get_assembly_input(wildcards)}
    if wants_ultralong(wildcards.input):
        d["b"] = get_ultralong_input(wildcards)
    return d



# ============================================================================
#  Dorado binary
# ============================================================================
# config may give either:
#   dorado: {version: "2.1.2"}          -> fetched by rule get_dorado
#   dorado: {path: /abs/path/dorado}    -> use an existing install, no download
# A bare string (the old format) is still accepted and treated as a path.
_dorado_cfg = config["dorado"]
if isinstance(_dorado_cfg, str):
    DORADO_VERSION = None
    DORADO_BIN = _dorado_cfg
else:
    DORADO_VERSION = _dorado_cfg.get("version", "2.1.2")
    DORADO_BIN = (_dorado_cfg.get("path")
                  or f"resources/dorado-{DORADO_VERSION}-linux-x64/bin/dorado")

def dorado_bin(wildcards):
    """Declared as a rule input so dorado jobs wait for the download."""
    return DORADO_BIN

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
    "hifiasm":          {"mem_mb": 30000,  "runtime": 440},
    "flye":             {"mem_mb": 30000,  "runtime": 440},
    "busco":            {"mem_mb": 14000,  "runtime": 880},  # was 10000; peaked at 9950 (99.5%) on a large test sample
    "rasusa":           {"mem_mb": 5000,   "runtime": 400},
    "chopper":          {"mem_mb": 5000,   "runtime": 400},
    "fga":              {"mem_mb": 5000,   "runtime": 600},
    "porechop_api":     {"mem_mb": 100000, "runtime": 6000},
    "dorado_basecall":  {"mem_mb": 100000, "runtime": 100800},
    "seqkit":           {"mem_mb": 10000,  "runtime": 6000},
    "dorado_align":     {"mem_mb": 75000,  "runtime": 72000},
    "dorado_polish":    {"mem_mb": 100000, "runtime": 60000},
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

include: "workflow/rules/0_0_get_dorado.smk"
include: "workflow/rules/0_1_basecall.smk"
include: "workflow/rules/0_2_bam_to_fastq.smk"
include: "workflow/rules/0_3_decompress.smk"
include: "workflow/rules/1_porechop_api.smk"
include: "workflow/rules/2_chopper.smk"
include: "workflow/rules/2_chopper_Flye_mitochondria.smk"
include: "workflow/rules/2_flye.smk"
include: "workflow/rules/2_rasusa.smk"
include: "workflow/rules/3_dorado_correct.smk"
include: "workflow/rules/3_seqtk_fasta_to_fastq.smk"
include: "workflow/rules/3_2_length_filter_corrected.smk"
include: "workflow/rules/5_hifiasm.smk"
include: "workflow/rules/6_0_fga_to_fa.smk"
include: "workflow/rules/6_00_lowest_contig_count.smk"
include: "workflow/rules/6_01_highest_busco.smk"
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
        for sel in SELECTORS:
            # Final assembly: polished where a BAM exists, otherwise the best raw one
            if has_bam(n):
                t.append(f"data/dorado_polish/{n}_{sel}.fasta")
            else:
                t.append(f"data/contig/{n}_{sel}_file.fa")
            t.append(f"data/busco/{n}_{sel}/BUSCO")
            t.append(f"results/{n}/{sel}/summary.txt")
            t.append(f"results/{n}/{sel}/{n}_{sel}_final.fasta")
        if wants_organelle(n):
            t.append(f"data/getorganelle/{n}/Mitochondria")
    return t

rule all:
    input:
        final_targets()