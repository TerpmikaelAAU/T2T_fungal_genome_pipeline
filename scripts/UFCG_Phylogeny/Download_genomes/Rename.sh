#!/usr/bin/bash -l
#SBATCH --job-name=NCBi_datasets
#SBATCH --output=UFCG/NCBi_datasets_rename%j.out
#SBATCH --error=UFCG/NCBi_datasets_rename%j.err
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

genus=("Neurospora" "Aspergillus" "Trichoderma")

# Rename fasta files based on organism name and move to new directory
rename_and_move_fasta() {
    for g in "${genus[@]}"; do
        metadata_file="Genomes/${g}/ncbi_dataset/data/assembly_data_report.jsonl"
        if [ -f "${metadata_file}" ]; then
            mkdir -p Genomes/${g}/organism_name
            while IFS= read -r line; do
                accession=$(echo "$line" | jq -r '.accession')
                organism=$(echo "$line" | jq -r '.organism.organismName' | tr ' ' '_')
                fasta_file="Genomes/${g}/ncbi_dataset/data/${accession}/*.fna"
                if ls ${fasta_file} 1> /dev/null 2>&1; then
                    new_name="Genomes/${g}/organism_name/${organism}.fna"
                    mv ${fasta_file} "${new_name}"
                fi
            done < "${metadata_file}"
        fi
    done
}

# Main execution logic
main() {
    rename_and_move_fasta
}

main