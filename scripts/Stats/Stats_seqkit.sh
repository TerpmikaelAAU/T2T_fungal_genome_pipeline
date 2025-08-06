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


conda activate seqkit
seqkit stats path/to/fasta -Ta > output.csv

