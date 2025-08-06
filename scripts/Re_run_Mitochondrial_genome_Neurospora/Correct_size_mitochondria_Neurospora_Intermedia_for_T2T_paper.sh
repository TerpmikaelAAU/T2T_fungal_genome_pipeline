#!/usr/bin/bash -l
#SBATCH --job-name=Extract_data
#SBATCH --output=extractdata_log/Extract_data%j.out
#SBATCH --error=extractdata_log/Extract_data%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
# Version 16-04-2025 - Mikael Terp
# Exit on first error and if any variables are unset
set -eu

## README ##
#

max_threads="$(nproc)"

mkdir -p test_assembly

INPUT_FASTQ="Path/to/input/.fastq"
output_dir="test_assembly"

conda activate chopper
chopper -q 25 -l 10000 --threads $(nproc) < $INPUT_FASTQ > ${output_dir}/test.fastq
conda deactivate
conda activate flye
flye -t ${max_threads} --nano-hq ${output_dir}/test.fastq -o "$output_dir" --read-error 0.03
conda deactivate



conda activate getorganelle
Input_graph="Input/assembly_graph.gfa"

get_organelle_from_assembly.py -F fungus_mt --config-dir get/organelle/path/0.0.1 -g $Input_graph -o Output/folder -t $max_threads --overwrite
conda deactivate