#!/bin/bash -l

##############################
#         Reformat           #
##############################

#SBATCH --job-name=reformat

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-00:30:00
#SBATCH --partition=normal
#SBATCH --mem=2GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=reformat.%J.out
#SBATCH --error=reformat.%J.err


#bbmap only needs a little amount of memory even with a lot of reads

conda activate bio

for sample in $(cat ~/amplicon_analysis/samples)
do
reformat.sh in=~/amplicon_analysis/Merged/${sample}.assembled.fastq  out=~/amplicon_analysis/Swarm/${sample}_good.fasta fastawrap=1000
done

exit 0
