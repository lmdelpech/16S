#!/bin/bash -l

##############################
#          Clipping          #
##############################

#SBATCH --job-name=primer_clipping

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-01:00:00
#SBATCH --partition=normal
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=cutadapt.%J.out
#SBATCH --error=cutadapt.%J.err

conda activate bio


for sample in $(cat ~/amplicon_analysis/samples)
do
    cutadapt -a ^GTGCCAGCMGCCGCGGTAA...ATTAGAWACCCBDGTAGTCC \
    -A ^GGACTACHVGGGTWTCTAAT...TTACCGCGGCKGCTGGCAC \
    -m 215 -M 285 --discard-untrimmed \
    -o ${sample}_R1_clipped.fq -p ${sample}_R2_clipped.fq \
    ~/amplicon_analysis/Original/${sample}_R1.fq ~/amplicon_analysis/Original/${sample}_R2.fq \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done

exit 0
