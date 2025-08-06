#!/usr/bin/bash -l
#SBATCH --job-name=funannotate
#SBATCH --output=Joblogs/funannotate/funannotate%j.out
#SBATCH --error=Joblogs/funannotate/funannotate%j.err
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


##
#fasta="funannotate_results/barcode10.fasta"
#fungi="aspergillus oryzae"
#outdir="funannotate_results/aspergillus_oryzae_funannotate"
#Name="aspergillus_oryzae"
##

####
fasta="Input fasta"
#Augustus name
fungi="input fungi name fx. fusarium graminearum"
outdir="output/directory"
Name="same as input name just with a underscore fx. fusarium_graminearum"
#####


### Remeber to have a look at the GBK file before submission, so it actually says the correct organism!!

mkdir -p $outdir
mkdir -p ${outdir}/predict
mkdir -p ${outdir}/annotate
mkdir -p ${outdir}/antismash
mkdir -p ${outdir}/interproscan
#
conda activate funannotate
export FUNANNOTATE_DB=path/funannotate_db
funannotate clean -i $fasta -o ${outdir}/${Name}_clean.fasta
funannotate sort -i ${outdir}/${Name}_clean.fasta -o ${outdir}/${Name}_sorted.fasta -b chr --minlen 1
funannotate mask -i ${outdir}/${Name}_sorted.fasta -o ${outdir}/${Name}_masked.fasta --cpus $max_threads 
conda deactivate

conda activate funannotate
export FUNANNOTATE_DB=path/funannotate_db
funannotate predict -i ${outdir}/${Name}_masked.fasta -o ${outdir}/predict -s "$fungi" --busco_seed_species $Name --cpus $max_threads 
conda deactivate
conda activate antismash
antismash -t fungi -c $max_threads --output-dir ${outdir}/antismash --output-basename $Name --asf --cc-mibig --tigrfam --cb-general ${outdir}/predict/predict_results/${Name}.gbk
conda deactivate
conda activate interproscan
export Java path
/home/bio.aau.dk/zl01hh/my_interproscan/interproscan-5.73-104.0/interproscan.sh -f XML -i ${outdir}/predict/predict_results/${Name}.proteins.fa -cpu $max_threads --outfile ${outdir}/interproscan/${Name}.XML
conda deactivate
conda activate funannotate
export FUNANNOTATE_DB=path/funannotate_db
funannotate annotate -i ${outdir}/predict -o ${outdir}/annotate -s "$fungi" --antismash ${outdir}/antismash/${Name}.gbk --iprscan ${outdir}/interproscan/${Name}.XML --fasta ${outdir}/${Name}_masked.fasta
conda deactivate


