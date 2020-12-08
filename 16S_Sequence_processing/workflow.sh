# Analysis of 16S amplicons
# Author: Lisa-Marie Delpech

# This script is adapted from the 16S workflow authored by Christiane Hassenrück
# Original script: please visit https://github.com/chassenr/NGS/blob/master/AMPLICON/illumina16Sworkflow_master.txt
# Script adapted to run on the Stallo HPC, using Slurm queuing system
# Input sequences do NOT contain barcode and adapter sequences anymore, but still contain primer sequences


#######################################################
#                    General steps                    #
#######################################################

# step 0: shorten file names, rename to running base name
# step 1: primer clipping with Cutadapt
# step 2: trimming with Trimmomatic
# step 3: read merging with Pear
# step 4: quality control with FastQC
# step 5: swarm OTU clustering
# step 6: taxonomic classification
# step 7: further analysis with R

#######################################################
#                  Analysis pipeline                  #
#######################################################


### Step 0: shorten file names, rename to running base name

## Create a home working directory

mkdir amplicon_analysis # this will be our main directory
cd amplicon_analysis

## Create a file containing all base names

ls *_R1.fq | cut -f1-3 -d "_" > samples #depending on the structure of the name of the fq files
# this file allows to make for loops for all samples, using command for sample in $(cat ~/amplicon_analysis/samples); do

## Moving Original files to separate directory

mkdir Original
mv *.fastq ~/amplicon_analysis/Original

# If necessary rename the original files



### STEP 1: primer clipping

# cutadapt v2.8

# /!\ this step was done simultaneously with demultiplexing using cutadapt, running it again would discard all the reads
# Just rename the files with _clipped so that the rest of the pipeline still works

# for each sample the following commands will remove the primer sequence when found (or discard the read)
# it will search both R1 and R2 for both the forward and the reverse primer sequence

cd ~/amplicon_analysis
sbatch ./Jobscripts/primer_clipping.sh

# Change the forward and reverse primer sequences of course...
# It is good to have a primers.fasta file to check the primers
# Check the names of the input and output files
# Input files should be ${sample}_R1.fq or ${sample}_R1.fastq and ${sample}_R2.fq
# Output files should be ${sample}_R1_clipped.fq and ${sample}_R2_clipped.fq

# We can look through the output of the cutadapt stats file we made (“cutadapt_primer_trimming_stats.txt”)
# Here is a line to look at what fractions of reads were retained in each sample (column 2)
# and at what fraction of bps were retained (column 3)

paste ~/amplicon_analysis/samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")")

## Ordering the new files in the folder Clipped

mkdir Clipped
mv *_clipped.fq ~/amplicon_analysis/Clipped
cd Clipped
mkdir Clipped_logs
mv ~/amplicon_analysis/cutadapt.* ~/amplicon_analysis/Clipped/Clipped_logs
mv ~/amplicon_analysis/cutadapt_primer_trimming_stats.txt ~/amplicon_analysis/Clipped/



### STEP 2: quality trimming

# trimmomatic v0.39

# For long inserts (450ish and onwards) recommended after merging, but this is not our case here (292 bp)
# We need all the bases we can get, additionally two identical bases with low quality which are merged will generally have a higher quality score
# Shorter inserts have more overlap and can afford some loss
# For the standard bacterial illumina insert we can do it before merging

# ptrim/paired = identical sequence headers with R1 and R2 strim/unpaired = single output, complementary reads removed
# SLIDINGWINDOW:4:10 is the absolute minimum! 4:15 recommended
# argument order matters! MINLEN should come after trimming.
# add a phred score threshold?

# For each sample, trim the clipped reads with a sliding window of 4 and a quality threshold of 15
# Discard reads less than 100 bps

## Running trimmomatic on the cluster

cd ~/amplicon_analysis
sbatch ./Jobscripts/quality_trimming.sh   # Here remember to change the jobscript according to base name of the samples

## Cleaning up directories

mkdir Trimmed
mv *paired.fq ~/amplicon_analysis/Trimmed
cd Trimmed
mkdir Trimmed_logs
mv ~/amplicon_analysis/*.log ~/amplicon_analysis/Trimmed/Trimmed_logs
mv ~/amplicon_analysis/trimmomatic.* ~/amplicon_analysis/Trimmed/Trimmed_logs

cd ~/amplicon_analysis


### STEP 3: read merging with pear

# pear  v0.9.6

# This will merge reads with a minimum overlap of 10 (-v), has to be changed depending on the amplicons
# The minimum length of the merged reads is 350 bp (-n), has to be changed depending on the amplicons
# For short insert sizes it might be recommended to set a maximum length for the merged reads (-m).
# freakishly long reads generally indicate an error...

# j: threads
# v: overlap
# n: min insert length
# m: max insert length (expect ...)
# o: output just needs basename, other stuff is added by pear
# no trimming (q) enabled as trimmomatic did the work here

## Running pear on the cluster

cd ~/amplicon_analysis
sbatch ./Jobscripts/read_merging.sh

# Cleaning up directories

mkdir Merged
mv *.fastq ~/amplicon_analysis/Merged
cd Merged
mkdir Merged_logs
mv ~/amplicon_analysis/*.log ~/amplicon_analysis/Merged/Merged_logs
mv ~/amplicon_analysis/pear.* ~/amplicon_analysis/Merged/Merged_logs
cd ~/amplicon_analysis



### STEP 4: quality control with FastQC

# FastQC v0.11.9

# Input files of fastQC are the .assembled files from the merging step since those are the only ones with enough quality to be merged
# NB. The output files are not needed for further steps, but only to check quality.

## Running FastQC on the cluster

cd ~/amplicon_analysis
sbatch ./Jobscripts/quality_control.sh

# Output files are named according to pear renaming from the previous steps
# They have .html and .zip extensions
# Unzip step will be necessary for further steps

## Cleaning up directories

mkdir FastQC
mv *.assembled_* ~/amplicon_analysis/FastQC
cd FastQC
mkdir FastQC_logs
mv ~/amplicon_analysis/fastqc.* ./FastQC_logs

## Unzip the output files

cd FastQC
for sample in $(cat ~/amplicon_analysis/samples)
do
unzip ${sample}.assembled_fastqc.zip
done

# Or unzip *.zip

## Copy files on current computer to have graphic overview

#On the terminal of the current computer
cd ~/Desktop
mkdir FastQC
scp lmdelpech@stallo.uit.no:~/amplicon_analysis/FastQC/* ~/Desktop/FastQC

# In each folder associated with each sample, we can find a fastqc_data.txt file

# Some of the quality scoring done by FASTQC results in inappropriate WARN/FAIL statuses:
# summary stats other than the mean are set to 0, probably due to a lack of sequence representation
# thus the following code will parse through the fastqc output re-assign PASS/WARN/FAIL statuses
# based on per base sequence quality...

# Edit: FastQC was originally conceived for genomic data, so some stats can have inappropriate statuses, when still good for amplicon data


################################################### This is no longer needed with this FastQC version #############################################

cd ~/amplicon_analysis

for i in $(cat ~/amplicon_analysis/samples)
do
# pull out the section of the fastqc output concerning the per-base quality
awk '/^>>Per base sequence quality/,/^>>END_MODULE/' ~/amplicon_analysis/FastQC/$i.assembled_fastqc/fastqc_data.txt | grep "^[0-9]" > ./FastQC/$i.assembled_fastqc/fastqc_PBSQ.txt


# check if the number of bins with a median quality value (column 3, $3 in awk) not equal to zero is equal
# to the number of bins with a median quality of at least 15 and a 10th percentile (column 6, $6 in awk) of at least quality 5.
# if that's true, check if the quality has a median qual of 25 and 10th percentile of at least 10 - assign PASS
# else WARN
# else FAIL
# Based on these results, modify the per base sequence quality (row 2) of the summary.txt file and capture
# the results in a new summary file: summary1.txt.

 if [ $(awk '$3!=0' ./FastQC/$i".assembled_fastqc"/fastqc_PBSQ.txt | wc -l) -eq $(awk '$3>=15 && $6>=5' ./FastQC/$i".assembled_fastqc"/fastqc_PBSQ.txt | wc -l) ];
  then
    if [ $(awk '$3!=0' ./FastQC/$i".assembled_fastqc"/fastqc_PBSQ.txt | wc -l) -eq $(awk '$3>=25 && $6>=10' ./FastQC/$i".assembled_fastqc"/fastqc_PBSQ.txt | wc -l) ];
    then
      awk '{if(NR==2)$1="PASS\t"; print $0}' ./FastQC/$i".assembled_fastqc"/summary.txt > ./FastQC/$i".assembled_fastqc"/summary1.txt
    else
      awk '{if(NR==2)$1="WARN\t"; print $0}' ./FastQC/$i".assembled_fastqc"/summary.txt > ./FastQC/$i".assembled_fastqc"/summary1.txt
    fi
  else
    awk '{if(NR==2)$1="FAIL\t"; print $0}' ./FastQC/$i".assembled_fastqc"/summary.txt > ./FastQC/$i".assembled_fastqc"/summary1.txt
  fi

done

################################################### This is no longer needed with this FastQC version #############################################

## Output some diagnostic files

# Combine flags of 'Per base sequence quality' module for all files
grep "Per base sequence quality" ~/amplicon_analysis/FastQC/*.assembled_fastqc/summary.txt > ~/amplicon_analysis/FastQC/QC_summary.txt

# Range of read lengths
grep "Sequence length" ~/amplicon_analysis/FastQC/*.assembled_fastqc/fastqc_data.txt > ~/amplicon_analysis/FastQC/QC_read_length.txt

# Combine flags of 'Sequence Length Distribution' module for all files including most abundant read lengths
for i in $(cat ~/amplicon_analysis/samples)
do
  awk '/^>>Sequence Length Distribution/,/^>>END_MODULE/' ~/amplicon_analysis/FastQC/$i.assembled_fastqc/fastqc_data.txt |\
  sed -e '1,2d' -e '$d' > ~/amplicon_analysis/FastQC/$i.assembled_fastqc/fastqc_SLD.txt

  sort -t $'\t' -k2nr ~/amplicon_analysis/FastQC/$i.assembled_fastqc/fastqc_SLD.txt |\
  head -1 |\
  paste <(grep "Sequence Length Distribution" ~/amplicon_analysis/FastQC/$i.assembled_fastqc/summary.txt) -
done > ~/amplicon_analysis/FastQC/QC_read_distribution.txt

## Count sequences

# Only counting forward read as representative for PE
# Syntax depending on the sequence identifier @MISEQ, @...
# grep -c '^@M05585' ./Original/*_R1.fq > nSeqs_all.txt
grep -c '^@M05585' ./Clipped/*_R1_clipped.fq >> nSeqs_all.txt
grep -c '^@M05585' ./Trimmed/*_R1_paired.fq >> nSeqs_all.txt
grep -c '^@M05585' ./Merged/*.assembled.fastq >> nSeqs_all.txt


## Cleaning up directories

mv nSeqs_all.txt FastQC


### STEP 5: Swarm OTU clustering

# Swarm v3.0.0

# Extract fasta file from fastq and move to new directory
# Set fastawrap to 1000 to prevent line breaks within sequence

# Reformat is a script included in bbmap program, here the script submitted calls this reformat.sh script
# With input fastq files ".assembled.fastq" from Pear ouput
# Output fasta files, ready for Swarm clustering "_good.fasta"

sbatch ./Jobscripts/reformat.sh

# No need for cleaning up directories since it is already in the script
# It takes files from the Merged folder in input
# And outputs the files in fasta format in the Swarm folder

# Now to dereplicate and rename individual reads to save compute and mental anguish downstream...
# The dereplication code is courtesy of the Swarm developers and can be found here:
# https://github.com/torognes/swarm/wiki/Working-with-several-samples

## Running dereplication script on the cluster

sbatch ./Jobscripts/dereplicate.sh

# input files: Swarm/*_good.fasta
# output files: *_dereplicated.fasta

## Cleaning up directories

mkdir Dereplicated
mv *_dereplicated.fasta ~/amplicon_analysis/Dereplicated
cd Dereplicated
mkdir Dereplicated_logs
mv ~/amplicon_analysis/dereplicate.* Dereplicated_logs
cd ~/amplicon_analysis

## Study level of dereplication

export LC_ALL=C
cat ~/amplicon_analysis/Dereplicated/*_dereplicated.fasta | \
awk 'BEGIN {RS = ">" ; FS = "[_\n]"}
     {if (NR != 1) {abundances[$1] += $2 ; sequences[$1] = $3}}
     END {for (amplicon in sequences) {
         print ">" amplicon "_" abundances[amplicon] "_" sequences[amplicon]}}' | \
sort --temporary-directory=$(pwd) -t "_" -k2,2nr -k1.2,1d | \
sed -e 's/\_/\n/2' > all_samples.fasta

# output file: all_samples.fasta

## Build contingency table with python

# The script to build the contingency table is located in the main directory (~/amplicon_analysis/amplicon_contingency_table.py)
# The python script is written for an older version than python3, therefore we have to load an older version on the cluster (in the script)
# possibility to change this, the problem with the python3 script is basically that the separator is wrong

## Running the script on the cluster

sbatch ./Jobscripts/build_amplicon_contingency_table.sh

# output file: amplicon_table.csv

## Cleaning up directories

cd Swarm
mkdir Swarm_logs
mv ~/amplicon_analysis/python.* Swarm_logs
mv ~/amplicon_analysis/reformat.* Swarm_logs
cd ..


## Swarming

# -b light swarms have less than 3 reads associated with them
# -d 1: local edit distance threshold is 1
# fastidious algorithm (-f): light swarms (amplicon abundance less than 3) will be grafted to heavy swarms
# -t set threads to 4
# -l output a log file
# -o the swarm file itself
# -s output a stats file (needed downstream)
# -w output fasta file with seed sequences

## Running swarm clustering on the cluster

sbatch ./Jobscripts/swarming.sh

# input files: all_samples.fasta, created with the study of dereplication
# output files: output of swarm is amplicons.swarms, but we also ask for amplicons_seeds.fasta (-w) and amplicons_stats.txt (-s)
# and for a log file (-l)

## Cleaning up directories

cd Swarm
mv ~/amplicon_analysis/swarm.* Swarm_logs
cd ~/amplicon_analysis


## Building OTU contingency table for multiple samples
# https://github.com/torognes/swarm/wiki/Working-with-several-samples

## Defining where the files are
cd ~/amplicon_analysis
STATS="amplicons_stats.txt"
SWARMS="amplicons.swarms"
AMPLICON_TABLE="amplicon_table.csv"
OTU_TABLE="OTU_contingency_table.csv"

# NB. Add this to the swarming jobscript

## Header
echo -e "OTU\t$(head -n 1 "${AMPLICON_TABLE}")" > "${OTU_TABLE}"

## Compute "per sample abundance" for each OTU
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

# We may want to check if large swarms are taxonomically consistent
# by classifying more than their seed sequences

## Cleaning up directories

mv amplicons_seeds.fasta amplicons_stats.txt amplicons.swarms amplicon_table.csv OTU_contingency_table.csv ~/amplicon_analysis/Swarm

### STEP 6: taxonomic classification

# At this stage, we could consider removing very rare swarms (less than one or two reads per swarm)
# As a large chunk of the swarms are rare (and will probably be removed from analysis later), we can save compute time here
# As always, whether this is advisable or not depends on our question

# Convert lowercase sequences to uppercase sequences in amplicons_seeds.fasta
cd Swarm
awk '{print /^>/ ? $0 : toupper($0)}' amplicons_seeds.fasta > amplicons_seeds_uc.fasta
cd ~/amplicon_analysis

#################################################### Option with BLASTN ##################################################

# Blastn has multicore support, and also supports multi sequences files
# The amplicons_seeds_uc.fasta thus doesn't need to be split

## Running blastn on the cluster

#num_threads: cpus per task (has multithreads support)
# -query: the input file
# -task: which blast to perform
# -outfmt: output format, 6 for tabular, 7 for tabular with annotations, precise the columns
# -out: name of the output file
# -db: compulsory, path and root name of the database
#
# -outfmt options:
# 	 qseqid means Query Seq-id
#    	       qgi means Query GI
#    	      qacc means Query accesion
#    	   qaccver means Query accesion.version
#    	      qlen means Query sequence length
#    	    sseqid means Subject Seq-id
#    	 sallseqid means All subject Seq-id(s), separated by a ';'
#    	       sgi means Subject GI
#    	    sallgi means All subject GIs
#    	      sacc means Subject accession
#    	   saccver means Subject accession.version
#    	   sallacc means All subject accessions
#    	      slen means Subject sequence length
#    	    qstart means Start of alignment in query
#    	      qend means End of alignment in query
#    	    sstart means Start of alignment in subject
#    	      send means End of alignment in subject
#    	      qseq means Aligned part of query sequence
#    	      sseq means Aligned part of subject sequence
#    	    evalue means Expect value
#    	  bitscore means Bit score
#    	     score means Raw score
#    	    length means Alignment length
#    	    pident means Percentage of identical matches
#    	    nident means Number of identical matches
#    	  mismatch means Number of mismatches
#    	  positive means Number of positive-scoring matches
#    	   gapopen means Number of gap openings
#    	      gaps means Total number of gaps
#    	      ppos means Percentage of positive-scoring matches
#    	    frames means Query and subject frames separated by a '/'
#    	    qframe means Query frame
#    	    sframe means Subject frame
#    	      btop means Blast traceback operations (BTOP)
#    	    staxid means Subject Taxonomy ID
#    	  ssciname means Subject Scientific Name
#    	  scomname means Subject Common Name
#    	sblastname means Subject Blast Name
#    	 sskingdom means Subject Super Kingdom
#    	   staxids means unique Subject Taxonomy ID(s), separated by a ';'
#    			 (in numerical order)
#    	 sscinames means unique Subject Scientific Name(s), separated by a ';'
#    	 scomnames means unique Subject Common Name(s), separated by a ';'
#    	sblastnames means unique Subject Blast Name(s), separated by a ';'
#    			 (in alphabetical order)
#    	sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
#    			 (in alphabetical order)
#    	    stitle means Subject Title
#    	salltitles means All Subject Title(s), separated by a '<>'
#    	   sstrand means Subject Strand
#    	     qcovs means Query Coverage Per Subject
#    	   qcovhsp means Query Coverage Per HSP
#    	    qcovus means Query Coverage Per Unique Subject (blastn only)
#    When not provided, the default value is:
#    'qaccver saccver pident length mismatch gapopen qstart qend sstart send
#    evalue bitscore', which is equivalent to the keyword 'std'
#    The supported format specifier for option 17 is:
#    	        SQ means Include Sequence Data
#    	        SR means Subject as Reference Seq


cd ~/amplicon_analysis
sbatch ~/amplicon_analysis/Jobscripts/blastn.sh

## Cleaning up directories

mkdir Blastn
mv blastn.* Blastn

#################################################### Option with BLASTN ##################################################

## Splitting seed sequence fasta file for parallel processing

# This step might not be necessary, Sina doesn't have multi core support, but we can try running it on the entire file
# the following command splits fasta file in fasta files with 500 sequences each

#asplit '^>' 500 < ~/amplicon_analysis/Swarm/amplicons_seeds_uc.fasta
awk -v size=500 -v pre=prefix -v pad=5 '
/^>/ { n++; if (n % size == 1) { close (fname); fname = sprintf("%s.%0" pad "d", pre, n) } } { print >> fname }
' amplicons_seeds_uc.fasta

## Determine how many chunks there are
JOBCOUNT=$(ls -1 out* | wc -l) #specify the number of files in array job

## Create a file with the names of the chunks
cd Swarm
ls prefix.* | -f2 -d "." > jobnames

## Assign Silva database to a variable
SINA_DB="~/amplicon_analysis/SILVA_138_SSURef_NR99.arb"
# Already in the script

## Running sina on the cluster
sbatch ./Jobscripts/sina_2.sh


## Cleaning up directories

mkdir Sina
mv sina_out.* ~/amplicon_analysis/Sina
cd Sina
mkdir Sina_logs
mv ~/amplicon_analysis/sina.* Sina_logs

# Rename sina_out.csv to avoid confusion in the next steps
mv sina_out.csv sina_table.csv
cd ..


########### This (below) is not necessary anymore with sina_out as a csv :) ############


# Time to gather up the useful info from the split output
# In grep, -h suppresses printing of filenames for results

# Get all the swarm seed hashes (sort of like accessions)

#grep -h '^' $(ls -1v ~/amplicon_analysis/Sina/sina_out.*) | sed 's/^sequence_identifier: //' > amplicons_seeds.accnos
grep -h '>' $(ls -1v ~/amplicon_analysis/Sina/sina_out.fasta) | sed 's/>//' > amplicons_seeds.accnos
# This command gets all the sequence identifiers from the sina_out. alignment file (not from the csv table)
# To get the sequence identifiers from the sina_table.csv:
cut -f1 -d "," ~/amplicon_analysis/Sina/sina_table.csv > amplicons_seeds_table.accnos

# Check if the order is the same as in amplicons_seeds_uc.fasta
#grep '^>' ~/amplicon_analysis/Swarm/amplicons_seeds_uc.fasta | sed 's/^>//' | diff - amplicons_seeds.accnos
grep '>' ~/amplicon_analysis/Swarm/amplicons_seeds_uc.fasta | sed 's/>//' > amplicons_seeds_swarm.accnos
diff amplicons_seeds_swarm.accnos amplicons_seeds_table.accnos > sequence_identifiers_diff.txt

# Get all corresponding taxonomic paths (note the same order as the accnos)
# grep -h '^lca_tax_slv' $(ls -1v ~/amplicon_analysis/Sina/sina_out.*) | sed 's/^lca_tax_slv: //' > amplicons_seeds.tax_slv
# This doesn't work since we're searching for a column in a table

cut -f7 -d "," ~/amplicon_analysis/Sina/sina_table.csv > amplicons_seeds.tax_slv
# to get the 7th column

# Get all alignment qualities (for filtering later)
#grep -h '^align_quality_slv' $(ls -1v ~/amplicon_analysis/Sina/sina_out.*) | sed 's/^align_quality_slv: //' > amplicons_seeds.align_quality_slv

cut -f5 -d "," ~/amplicon_analysis/Sina/sina_table.csv > amplicons_seeds.align_quality_slv

# Merge these output files...
paste amplicons_seeds_table.accnos amplicons_seeds.align_quality_slv amplicons_seeds.tax_slv > amplicons_seeds_taxonomy.txt

# Cleaning up directories
mv amplicons_seeds* sequence_identifiers* ~/amplicon_analysis/Sina

############# This (above) is not necessary anymore with sina_out as a csv :) ############


# Copy final output files to working directory
cp ~/amplicon_analysis/Swarm/OTU_contingency_table.csv ~/amplicon_analysis/Sina/amplicons_seeds_taxonomy.txt ~/amplicon_analysis


### STEP 7: further analysis with R

# Run ReadAmplicon.R
# Not necessarily on the cluster
