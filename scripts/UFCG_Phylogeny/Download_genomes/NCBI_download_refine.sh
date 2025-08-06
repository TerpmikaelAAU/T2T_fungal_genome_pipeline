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

#genus=("Neurospora" "Aspergillus" "Trichoderma")
#genus=("Neurospora")
#Refrence Genus, used in directory creation and in NCBI download
#genus=("Aspergillus oryzae" "Aspergillus flavus")
genus=("Trichoderma asperelloides" "Trichoderma asperellum" "Trichoderma yunnanense")
Final_folder_Name="Trichoderma"
conda activate ncbi_datasets

##make directory  
#create_directories() {
#    for g in "${genus[@]}"; do
#        dir_name=$(echo "$g" | tr ' ' '_')  # Replace spaces with underscores
#        mkdir -p NCBI_and_own_genomes/Genomes/Refine/${dir_name}
#    done
#}#

#download reference dataset, based on Genus
#download() {
#    for g in "${genus[@]}"; do
#        dir_name=$(echo "$g" | tr ' ' '_')
#        datasets download genome taxon "${g}" --include genome 
#        unzip ncbi_dataset.zip -d  NCBI_and_own_genomes/Genomes/Refine/${dir_name}
#        rm -r ncbi_dataset.zip
#    done
#}

#--assembly-level "chromosome" --assembly-level "complete"



# Find and copy all .fna files from directories starting with GCA_ or GCF_ in the input directory
#copy() {
#   for g in "${genus[@]}"; do
#       dir_name=$(echo "$g" | tr ' ' '_')
#       input_dir=NCBI_and_own_genomes/Genomes/Refine/${dir_name}/ncbi_dataset/data
#       output_dir=NCBI_and_own_genomes/Genomes/Refine/${dir_name}/ncbi_dataset/data/All
#       mkdir -p NCBI_and_own_genomes/Genomes/Refine/${dir_name}/ncbi_dataset/data/All
#       for dir in "$input_dir"/GCA_* "$input_dir"/GCF_*; do
#            if [ -d "$dir" ]; then
#            cp "$dir"/*.fna "$output_dir"/
#            fi
#        done
#    done
#}

move_fasta() {
   for g in "${genus[@]}"; do
       dir_name=$(echo "$g" | tr ' ' '_')
       input_dir=NCBI_and_own_genomes/Genomes/Refine/${dir_name}/ncbi_dataset/data/All
       mkdir -p NCBI_and_own_genomes/Genomes/Refine/${Final_folder_Name}
       output_dir=NCBI_and_own_genomes/Genomes/Refine/${Final_folder_Name}
       cp "$input_dir"/* "$output_dir"/
   done
}



 #Main execution logic
main() {
    #create_directories
    #download
    #copy
    move_fasta
}

main
