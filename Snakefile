import os
from snakemake.utils import min_version

min_version("7.18.2")

configfile: "config/config.yaml"

# Define resources
resources = {
    "hifiasm": {"mem_mb": 30000, "time": "04:00:00"},
    "flye": {"mem_mb": 30000, "time": "04:00:00"},
    "busco": {"mem_mb": 10000, "time": "08:00:00"},
    "rasusa": {"mem_mb": 5000, "time": "00:40:00"},
    "chopper": {"mem_mb": 5000, "time": "00:40:00"},
    "fga": {"mem_mb": 5000, "time": "01:00:00"},
    "porechop_api": {"mem_mb": 100000, "time": "10:00:00"},
    "dorado": {"mem_mb": 200000 , "time": "2-00:00:00"},
    "dorado_basecall": {"mem_mb": 100000,"time": "7-00:00:00"},
    "seqkit": {"mem_mb": 10000, "time": "01:00:00"},
    "rasusa": {"mem_mb": 5000, "time": "00:40:00"},
    "dorado_align": {"mem_mb": 75000, "time": "12:00:00"},
}

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
        


