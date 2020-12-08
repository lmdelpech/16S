#!/bin/bash -l

##############################
#   Build contingency table  #
##############################

#SBATCH --job-name=contingency_table

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-04:00:00
#SBATCH --partition=normal
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=python.%J.out
#SBATCH --error=python.%J.err

#script python2 only compatible with python2, python3 is default in Stallo
#module load Python/2.7.15-intel-2018b

python ~/amplicon_analysis/amplicon_contingency_table_python3_revised.py ~/amplicon_analysis/Dereplicated/*_dereplicated.fasta > amplicon_table.csv


exit 0
