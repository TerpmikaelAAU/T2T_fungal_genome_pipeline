#!/usr/bin/bash -l
#SBATCH --job-name=snake_assembly_sub
#SBATCH --output=snake_log/job_%j_%x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=shared
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=email@addrese.com

# Exit on first error and if any variables are unset
set -eu
# Activate conda environment with only snakemake
conda activate snakemake_run
###First time running the workflow###
snakemake --profile profile -s Snakefile --conda-create-envs-only
##After you run the above line to create all the envs, comment the above line, uncomment the line below #Main workflow. and change time=7-00:00:00, cpus-per-task=1, --mem=1G

# Start workflow using resources defined in the profile. Snakemake itself 
# requires nothing, 1CPU + 1G mem is enough (not for making conda envs)


# Main workflow
#snakemake --profile profile -s Snakefile 
#If the workflow crashes unlock with "--unlock", add to the line ^. Run it, then remove the unlock flag and rerun. 




