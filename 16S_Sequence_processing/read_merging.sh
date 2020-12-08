#!/bin/bash -l

##############################
#          Merging           #
##############################

#SBATCH --job-name=read_merging

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-02:00:00
#SBATCH --partition=normal
#SBATCH --mem=8GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=pear.%J.out
#SBATCH --error=pear.%J.err

conda activate bio

#-j: threads
#-v: minimum overlap size
#-n: min insert length
#-m: max insert length
#-o: output just needs basename, -the rest is added by pear
#-f: forward input reads, from the trimming step, paired reads with sufficient quality to be merged

for sample in $(cat ~/amplicon_analysis/samples)
do
pear -j 4 -v 10 -n 200 -m 500 -f ~/amplicon_analysis/Trimmed/${sample}_R1_paired.fq -r ~/amplicon_analysis/Trimmed/${sample}_R2_paired.fq -o ./${sample} > ./${sample}_merged.log
done 

exit 0
