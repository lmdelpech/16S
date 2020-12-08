#!/bin/bash -l

##############################
#         Dereplicate        #
##############################

#SBATCH --job-name=dereplicate_reads

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-08:00:00
#SBATCH --partition=normal
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=dereplicate.%J.out
#SBATCH --error=dereplicate.%J.err

conda activate bio

for sample in $(cat ~/amplicon_analysis/samples)
do
grep -v "^>" ~/amplicon_analysis/Swarm/${sample}_good.fasta | \
grep -v [^ACGTacgt] | sort -d | uniq -c | \
while read abundance sequence ; do
    hash=$(printf "${sequence}" | sha1sum)
    hash=${hash:0:40}
    printf ">%s_%d_%s\n" "${hash}" "${abundance}" "${sequence}"
done | sort -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > ${sample}_dereplicated.fasta
done

exit 0

