#!/bin/bash -l

##############################
#          Swarming          #
##############################

#SBATCH --job-name=swarming

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0-01:00:00
#SBATCH --partition=normal
#SBATCH --mem=16GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=swarm.%J.out
#SBATCH --error=swarm.%J.err

conda activate bio

swarm -b 3 -d 1 -f -t 4 -l swarm.log -o amplicons.swarms -s amplicons_stats.txt -w amplicons_seeds.fasta ~/amplicon_analysis/all_samples.fasta



# Defining where the files are

STATS="amplicons_stats.txt"
SWARMS="amplicons.swarms"
AMPLICON_TABLE="amplicon_table.csv"
OTU_TABLE="OTU_contingency_table.csv"

# Header
echo -e "OTU\t$(head -n 1 "${AMPLICON_TABLE}")" > "${OTU_TABLE}"

# Compute "per sample abundance" for each OTU
awk -v SWARM="${SWARMS}" -v TABLE="${AMPLICON_TABLE}"  'BEGIN {FS = " "
            while ((getline < SWARM) > 0) {
                swarms[$1] = $0
            }
            FS = "\t"
            while ((getline < TABLE) > 0) {
                table[$1] = $0
            }
           }

     {# Parse the stat file (OTUs sorted by decreasing abundance)
      seed = $3 "_" $4
      n = split(swarms[seed], OTU, "[ _]")
      for (i = 1; i < n; i = i + 2) {
          s = split(table[OTU[i]], abundances, "\t")
          for (j = 1; j < s; j++) {
              samples[j] += abundances[j+1]
          }
      }
      printf "%s\t%s", NR, $3
      for (j = 1; j < s; j++) {
          printf "\t%s", samples[j]
      }
     printf "\n"
     delete samples
     }' "${STATS}" >> "${OTU_TABLE}"


exit 0

