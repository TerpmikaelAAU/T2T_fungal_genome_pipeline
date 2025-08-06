05/08/2025 - Mikael Terp
## Scripts
The scripts in the scripts folder, are formatted as bash code to be submitted to a slurm HPC, so to run them you have to change the different lines in the start of the script.
Here the most important to change is probally: #SBATCH --partition=shared
I.e. change it to a partition paticular to your HPC
## General configuration
Make conda environments with the specific package (and version, some of them might work with the newest versions)
Change the input files -> then you can hopefully run them
 
