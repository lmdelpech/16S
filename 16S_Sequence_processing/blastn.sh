#!/bin/bash -l

##############################
#   Taxonomic clustering     #
##############################

#SBATCH --job-name=blastn

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-02:00:00
#SBATCH --partition=normal
#SBATCH --mem=32GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=blastn.%J.out
#SBATCH --error=blastn.%J.err

conda activate bio

blastn -num_threads 16 -task blastn -query ~/amplicon_analysis/Swarm/amplicons_seeds_uc.fasta -outfmt "7 qseqid sseqid pident length mismatch gapopen qtart qsen sstart send evalue bitscore qcovs qlen gaps ssciname scomname" -out blastn.csv -db ~/amplicon_analysis/16S_ribosomal_RNA

exit 0
