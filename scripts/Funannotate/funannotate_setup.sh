#!/usr/bin/bash -l
#SBATCH --job-name=funannotate
#SBATCH --output=joblog/funannotate/funannotate_setup%j.out
#SBATCH --error=joblog/funannotate/funannotate_setup%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=15
#SBATCH --mem=22G
#SBATCH --partition=shared
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
# Version 31-03-2025 - Mikael Terp
# Exit on first error and if any variables are unset
#set -eu
max_threads="$(nproc)"

## README ##
#


#start up conda ENV
conda activate funannotate

#check that all modules are installed
#funannotate check --show-versions

#download/setup databases to a writable/readable location
#funannotate setup -d $HOME/funannotate_db

#set ENV variable for $FUNANNOTATE_DB
echo "export FUNANNOTATE_DB=$HOME/funannotate_db" > /conda/installation/path/envs/funannotate/etc/conda/activate.d/funannotate.sh
echo "unset FUNANNOTATE_DB" > /conda/installation/path/envs/funannotate/etc/conda/deactivate.d/funannotate.sh

#run tests -- requires internet connection to download data
#funannotate test -t all --cpus $max_threads