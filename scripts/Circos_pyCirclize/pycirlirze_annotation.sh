#!/usr/bin/bash -l
#SBATCH --job-name=pyCirclize_anno
#SBATCH --output=_joblog/pyCirclize/pyCirclize_anno%j.out
#SBATCH --error=_joblog/pyCirclize/pyCirclize_anno%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
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
python3 Circos_GBK_plot.py
#python3 Circos_GBK_plot_Test.py
conda deactivate