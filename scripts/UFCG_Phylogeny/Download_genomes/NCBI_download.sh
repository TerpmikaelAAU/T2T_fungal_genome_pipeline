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

#Refrence Genus, used in directory creation and in NCBI download
genus=("Neurospora" "Aspergillus" "Trichoderma")
conda activate ncbi_datasets

#make directory  
create_directories() {
    for g in "${genus[@]}"; do
        mkdir -p NCBI_and_own_genomes/Genomes/${g}
    done
}

#download reference dataset, based on Genus
download() {
    for g in "${genus[@]}"; do
        datasets download genome taxon "${g}" --reference --include genome
        unzip ncbi_dataset.zip -d  NCBI_and_own_genomes/Genomes/${g}
        rm -r ncbi_dataset.zip
    done
}
# Main execution logic
main() {
    create_directories
    download
}

main
