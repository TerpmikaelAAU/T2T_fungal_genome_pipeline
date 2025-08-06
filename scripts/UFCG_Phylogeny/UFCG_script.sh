#!/usr/bin/bash -l
#SBATCH --job-name=UFCG
#SBATCH --output=UFCG_joblog/UFCG/UFCG%j.out
#SBATCH --error=UFCG_joblog/UFCG/UFCG%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=30G
#SBATCH --partition=shared
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#Version 31-03-2025 - Mikael Terp
# Exit on first error and if any variables are unset
set -eu
max_threads="$(nproc)"

## README ##
#


#Run UFCG code #The input is genomes from NCBI and inhouse genomes

base_dir="UFCG"
Genus=("Genus_name")
export JAVA_HOME=/home/bio.aau.dk/zl01hh/Java/jdk-21.0.2
export PATH=$JAVA_HOME/bin:$PATH

#Do this for refine.
for genus in "${Genus[@]}"; do
    input_genomes="NCBI_and_own_genomes/Genomes/Refine/Trichoderma"
    output_dir="${base_dir}/${genus}"
    mkdir -p "$output_dir"
    conda activate UFCG
    ufcg profile -i "$input_genomes" -o "${output_dir}/Profile" -t ${max_threads} 
    ufcg tree -i "${output_dir}/Profile" -o "${output_dir}/Tree" -t ${max_threads} 
    #ufcg align -i "${output_dir}/Profile" -o "${output_dir}/Align" -t ${max_threads}
    conda deactivate
done


