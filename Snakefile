import os
from snakemake.utils import min_version
 
min_version("7.18.2")
 
configfile: "config/config.yaml"
 
# Resources per tool.
# runtime is an INTEGER NUMBER OF MINUTES (Snakemake >=8 requirement).
# gres is passed straight through to `sbatch --gres=...` by the SLURM plugin.
# Partitions are NOT set: BioCloud assigns them automatically.
resources = {
    "hifiasm":         {"mem_mb": 30000,  "runtime": 240},
    "flye":            {"mem_mb": 30000,  "runtime": 240},
    "busco":           {"mem_mb": 10000,  "runtime": 480},
    "rasusa":          {"mem_mb": 5000,   "runtime": 40},
    "chopper":         {"mem_mb": 5000,   "runtime": 40},
    "fga":             {"mem_mb": 5000,   "runtime": 60},
    "porechop_api":    {"mem_mb": 100000, "runtime": 600},
    "dorado":          {"mem_mb": 200000, "runtime": 2880},
    "dorado_basecall": {"mem_mb": 100000, "runtime": 10080},
    "seqkit":          {"mem_mb": 10000,  "runtime": 60},
    "dorado_align":    {"mem_mb": 75000,  "runtime": 720},
    "dorado_polish":   {"mem_mb": 100000, "runtime": 600},
}
 
# GPU request string for BioCloud's single A10 node (bio-node10).
GPU_GRES = "gpu:a10:1"
 
include: "workflow/rules/0_1_basecall.smk"
include: "workflow/rules/0_2_bam_to_fastq.smk"
include: "workflow/rules/1_porechop_api.smk"
include: "workflow/rules/2_chopper_dorado_correction.smk"
include: "workflow/rules/2_chopper_Flye_mitochondria.smk"
include: "workflow/rules/2_chopper_ultralong.smk"
include: "workflow/rules/2_flye.smk"
include: "workflow/rules/2_rasusa.smk"
include: "workflow/rules/3_dorado_correction.smk"
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

#include: "workflow/rules/stats_seqkit.smk"

rule all:
    input:
        expand("data/dorado_basecall/{input}.bam", input=config["input"]),
        expand("data/samtools/Fastq/{input}.fastq", input=config["input"]),
        expand("0.0.1"),
        expand("data/getorganelle/{input}/Mitochondria", input=config["input"]),
        expand("data/porechopped/{input}.fastq",  input=config["input"]),
        expand("data/chopper/L10kbQ10/{input}.fastq",  input=config["input"]),
        expand("data/chopper/Flye/{input}.fastq",  input=config["input"]),
        expand("data/chopper/ultralong/{input}.fastq", input=config["input"]),
        expand("data/flye/{input}/flye_assembly/assembly.fasta",  input=config["input"]),
        expand("data/rasusa/Coverage/{input}.fastq",  input=config["input"]),
        expand("data/dorado/{input}.fasta",  input=config["input"]),
        expand("data/seqtk/fasta_to_fastq/{input}.fastq",  input=config["input"]),
        expand("data/seqkit/dorado/{input}_{length}.fastq", input=config["input"], length=config["length"]),
        expand("data/hifiasm/{input}_{length}/prefix.p_ctg.gfa",  input=config["input"], length=config["length"]),
        expand("data/hifiasm/{input}_{length}/{input}_{length}.fa", input=config["input"], length=config["length"]),
        expand("data/contig/{input}_lowest_contig_file.fa",  input=config["input"]),
        expand("data/dorado_align/{input}.bam", input=config["input"]),
        expand("data/dorado_polish/{input}.fasta", input=config["input"]),
        expand("data/busco/{input}/BUSCO", input=config["input"]),
        


