#!/bin/bash -l

##############################
#     Study dereplication    #
##############################

#SBATCH --job-name=study_dereplication

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-01:00:00
#SBATCH --partition=normal
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=study_dereplication.%J.out
#SBATCH --error=study_dereplication.%J.err


export LC_ALL=C
cat ~/amplicon_analysis/Dereplicated/*_dereplicated.fasta | \
awk 'BEGIN {RS = ">" ; FS = "[_\n]"}
     {if (NR != 1) {abundances[$1] += $2 ; sequences[$1] = $3}}
     END {for (amplicon in sequences) {
         print ">" amplicon "_" abundances[amplicon] "_" sequences[amplicon]}}' | \
sort --temporary-directory=$(pwd) -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > all_samples.fasta

exit 0

