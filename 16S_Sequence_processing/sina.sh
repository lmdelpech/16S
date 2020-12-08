#!/bin/bash -l

##############################
#   Taxonomic clustering     #
##############################

#SBATCH --job-name=sina_clustering

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=0-10:00:00
#SBATCH --partition=normal
#SBATCH --mem=32GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=sina.%J.out
#SBATCH --error=sina.%J.err

conda activate bio

SINA_DB="SILVA_138_SSURef_NR99.arb"

for jobname in $(cat ~/amplicon_analysis/Swarm/jobnames)
do
sina -i ~/amplicon_analysis/Swarm/prefix.${jobname} --intype fasta -o sina_out.${jobname} --outtype fasta --search --meta-fmt csv --overhang remove --insertion forbid --filter none --fs-kmer-no-fast --fs-kmer-len 10 --fs-req 2 --fs-req-full 1 --fs-min 40 --fs-max 40 --fs-weight 1 --fs-full-len 1400 --fs-msc 0.7 --match-score 1 --mismatch-score -1 --pen-gap 5 --pen-gapext 2 --search-cover query --search-iupac optimistic --search-min-sim 0.9 --turn all --lca-quorum 0.7 --search-db ${SINA_DB} --db ${SINA_DB} --lca-fields tax_slv
done

exit 0
