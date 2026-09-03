#!/usr/bin/bash -l
#SBATCH --job-name=snake_assembly_sub
#SBATCH --output=snake_log/job_%j_%x.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_90
#SBATCH --mail-user=email@addrese.com

# Exit on first error, unset variable, or failed pipe
set -euo pipefail

mkdir -p snake_log

conda activate snakemake_run_v9

# Snakemake itself needs almost nothing: it submits each rule as its own
# SLURM job via the executor plugin. 1 CPU + 1G is plenty.
#
# NOTE: --nodes/--ntasks are gone; they are meaningless for the launcher job.
# NOTE: environments are now created with --conda-create-envs-only below.

### FIRST RUN ONLY: build all conda envs, then comment this out ###
snakemake --profile profile --conda-create-envs-only

### MAIN WORKFLOW: uncomment after the envs exist ###
#snakemake --profile profile

# If the workflow crashed and the directory is locked, run once with --unlock,
# then remove the flag and rerun:
#snakemake --profile profile --unlock
