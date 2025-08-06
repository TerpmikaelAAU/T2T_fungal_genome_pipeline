#!/usr/bin/bash -l
#SBATCH --job-name=NCBi_datasets
#SBATCH --output=UFCG_joblog/NCBi_datasets%j.out
#SBATCH --error=UFCG_joblog/NCBi_datasets%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=shared
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mter@bio.aau.dk
#Version 06-12-2024 - Mikael Terp
# Exit on first error and if any variables are unset
set -eu
max_threads="$(nproc)"


# Define the input and output directories
#input_dir="Genomes/Neurospora/ncbi_dataset/data"
#output_dir="Genomes/Neurospora/ncbi_dataset/data/All"
#input_dir="NCBI_and_own_genomes/Genomes/Aspergillus/ncbi_dataset/data"
#output_dir="NCBI_and_own_genomes/Genomes/Aspergillus/ncbi_dataset/data/All"
input_dir="NCBI_and_own_genomes/Genomes/Trichoderma/ncbi_dataset/data"
output_dir="NCBI_and_own_genomes/Genomes/Trichoderma/ncbi_dataset/data/All"


# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Find and copy all .fna files from directories starting with GCA_ or GCF_ in the input directory
for dir in "$input_dir"/GCA_* "$input_dir"/GCF_*; do
  if [ -d "$dir" ]; then
    cp "$dir"/*.fna "$output_dir"/
  fi
done