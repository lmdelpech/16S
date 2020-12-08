# R and bash scripts for 16S amplicon data

This repository contains R and bash scripts for processing and analysis of 16S amplicon sequencing data.               

## 16S Sequence processing

This folder contains bash scripts adapted for the SLURM queuing system, from the 16S workflow written by Christiane Hassenrück (please visit: <https://github.com/chassenr/NGS>)
The script **workflow.sh** calls the other .sh scripts
The output of this pipeline is an OTU contingency table which contains the OTU identifiers as rows and the samples as columns; and a table that contains the taxonomic paths associated to the OTU identifiers. 
These files should be then processed in R using the **ReadAmplicon.R** script (Sources)

## 16S R Scripts

This folder contains R markdown script files, which are customized for the output of the **ReadAmplicon.R** outputs (Sources).
Each of them provides a short documentation on the required input files and the output.

* **16S_Data_formatting**: script that should be run first, provides important processing for OTU tables and other useful tables. Tables generated in this script are used downstream in all R scripts.  
* **16S_Alpha_diversity**: script to compute diversity indices, and plot according to groups of samples.
* **16S_Beta_diversity**: compute multivariate statistic and visualize ordinations with ggplot2.
* **16S_Taxonomic_composition**: investigate taxonomic composition and provides script to visualize data, as well as output of abundance tables.
* **16S_Indicator_species**: statistical method to identify indicator species, outputs abundance tables of the indicators per group of samples.
* **16S_Environmental_drivers**: compute constrained multivariate analyses using environmental data, visualize ordinations with ggplot2, correlation of indicators with environmental variables.   
* **16S_Functional_analysis**: based on the R package Tax4Fun to infer functions from taxonomy, visualize data with ordinations and heatmap.  

  
## Sources

R scripts called by the R scripts. Most of them were authored by Christiane Hassenrück.  

## Other

**silva138_tax_ssu_curated.tsv** is the file used for parsing taxonomy in the ReadAmplicon.R script. It was generated based on the  tax_slv_ssu.txt file from the SILVA archive by Christiane Hassenrück (<https://gitlab.leibniz-zmt.de/chh/bioinf/blob/master/silva138_tax_ssu_curated.tsv>)

:exclamation: These files are optimized for the studied environmental system and dependant on metadata of the system. 

:dna: The original sequencing data have been archived at the European Nucleotide Archive. I can provide a sample dataset to run the scripts.

:question: Please let me know if any problem or question.

