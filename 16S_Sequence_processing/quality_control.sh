#!/bin/bash -l

##############################
#      Quality control       #
##############################

#SBATCH --job-name=quality_control

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-01:00:00
#SBATCH --partition=normal
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=fastqc.%J.out
#SBATCH --error=fastqc.%J.err

conda activate bio

for sample in $(cat ~/amplicon_analysis/samples)
do
fastqc -o ~/amplicon_analysis ~/amplicon_analysis/Merged/${sample}.assembled.fastq
done

exit 0


