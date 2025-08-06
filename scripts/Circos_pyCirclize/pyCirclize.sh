#!/usr/bin/bash -l
#SBATCH --job-name=pyCirclize
#SBATCH --output=_joblog/pyCirclize/pyCirclize%j.out
#SBATCH --error=_joblog/pyCirclize/pyCirclize%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=5G
#SBATCH --partition=shared
#SBATCH --time=00:10:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
# Version 24-03-2025 - Mikael Terp
# Exit on first error and if any variables are unset
set -eu
max_threads="$(nproc)"
## Readme
# pyCirclize (v1.9.0) + mummer
conda activate pycirclize
python3 mummer.py
conda deactivate