#!/usr/bin/bash -l
#SBATCH --job-name=Telomere
#SBATCH --output=oblog/Telomere/Telomere%j.out
#SBATCH --error=joblog/Telomere/Telomere%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=8G
#SBATCH --partition=shared
#SBATCH --time=00:05:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mail@asdf.com
# Version 21-10-2024 - Mikael Terp
# Exit on first error and if any variables are unset
set -eu

BASE_DIR="Input_a_name_for_a_base_directory"
INPUT_FASTA="Input_fasta_file"
Name="Input_a_name"
telomeric_repeat_string="Input_a_telomeric_repeat_strin_fx:TTAGGG"

conda activate tidk

mkdir -p ${BASE_DIR}/$Name
tidk search --string $telomeric_repeat_string --output ${Name} --dir ${BASE_DIR}/$Name $INPUT_FASTA
tidk plot --tsv ${BASE_DIR}/$Name/${Name}_telomeric_repeat_windows.tsv -o ${BASE_DIR}/${Name}/${Name}_Telomere_plot
conda deactivate




