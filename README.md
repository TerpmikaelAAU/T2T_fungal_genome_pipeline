## Article
Paper: "Gapless near Telomer-to-Telomer Assembly of Neurospora intermedia, Aspergillus oryzae, and Trichoderma asperellum from Nanopore Simplex Reads"
https://doi.org/10.3390/jof11100701

## Raw data requirements
Nanopore fungal DNA sequenced with kit 10.4 chemistry -> Demultiplexed pod5 data
At least 100 x coverage Q10< (Less can work as show in the paper, but sequence to a high depth just to be sure)
## Instructions
# 0.1
Download this repository to you PC by clicking on "<> Code" and then "Download zip". Unzip the file and upload the folder to the server.
Then we need to make a snakemake envrionment. 
Either open the repository in Vscode or cd Repository_path. For these examples we always want to be inside the folder we are working in.
# 0.2
Change partition in profile/config.yaml this is the default ressources, then change it in all rules, and then in slurm slurm_submit.sbatch. gpu is for gpu partition, shared is for jobs that run at 1-2 gb/thread, and general is jobs that run >5 gb/per thread. 
# 0.3
Create a snakemake environment with the file `Snakemake_env.yml`
`conda env create -n "Environment name" -f "Path to yml file"`.
Run the below code, in the terminal
`conda env create -n snakemake_run -f Snakemake_env.yml`.
After this finishes you should have a snakemake env called "snakemake_run". 
# 0.4
Now we need to download dorado which we use for basecalling, read error correction and polishing.
Download dorado from https://github.com/nanoporetech/dorado , the file called "dorado-0.9.6-linux-x64", upload to server into a folder for example called dorado. 
cd into the folder, i.e. cd "path to dorado folder"
Unzip the dorado file:
`tar -xvzf <file_name>.tar.gz`

After this finishes we should be able to configure our run. 

## First time running the workflow
# 1 - General configuration
The `config.yaml`, is where we add our data, and can do some slight configuration of the run. (There are examples in the file)
In the `config.yaml` in the folder `config/`.
add the absolute path to unzipped dorado folder from before. 
add the absolute path to the pod5_pass folder. 
add the the names of your samples, only what they are called in the pod_pass folder. So if you used Barcode 10 and 11, they should be named accordingly. 
# 2 - Create all envrionments
Snakemake can run with just 1G mem and 1 cpu, but not when we make envrionments. So do the following: 
`sbatch slurm_submit.sbatch`
This will create all the environments, but will not run the workflow. 
# 3 - Running the workflow 
Change the Mem, Cpu, Time, and comment/uncomment as described in the `slurm_submit.sbatch` file. Save the file. 
And now we should just be able to run the pipeline with:
`sbatch slurm_submit.sbatch`
easiest way to see if it is running is `squeue --me`
# 4 - Output
snake_log - snakemake log file.
logs - logs specific for the particular step in the pipline.
data - Here is where all the data will be proccesed and the final assembly will be. 
The most important folders in data.
dorado_polish - The final polished assembly.
busco - fungi_odb10 results of the final assembly.
getorganelle - sometimes the mitochondrial genome is lost in this pipline so in getorganelle if there is a fungus_mt.complete.graph this might be a good mitochondrial genome. But should be checked.
# Notes
All rules and environments are under .workflow/rules. 
If you do not have Raw sequencing data, and only basecalled data in BAM format you can still run the workflow you just have to # the basecalling rule in the snakefile, and then change your data into the correct filename.
If you only have FASTQ data, you also have to skip dorado polishing (Rule 6_1,6_2,6_3), i.e # those rules and change the input for rule 7. 

In the /scripts folder, the code to make the figures presented in the paper is available. They 'relatively barebones, but with some manual work they should be able to run for you. 



