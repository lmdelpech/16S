#!/bin/bash -l

##############################
#          Trimming          #
##############################

#SBATCH --job-name=quality_trimming

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0-10:00:00
#SBATCH --partition=normal
#SBATCH --mem=8GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err

conda activate bio

#Job
#Paired-End Mode
#java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog <logFile>] >] [-basein <inputBase> | <input 1> <input 2>] [-baseout <outputBase> | <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1> ...


for sample in $(cat ~/amplicon_analysis/samples)
do
trimmomatic PE -threads 4 -trimlog ./${sample}_trim.log \
~/amplicon_analysis/Clipped/${sample}_R1_clipped.fq ~/amplicon_analysis/Clipped/${sample}_R2_clipped.fq \
./${sample}_R1_paired.fq ./${sample}_R1_unpaired.fq ./${sample}_R2_paired.fq ./${sample}_R2_unpaired.fq \
SLIDINGWINDOW:4:15 MINLEN:100
done

exit 0
