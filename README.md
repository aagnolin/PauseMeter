# Introduction
Suite of scripts to perform qualitative and quantitative analyses on ribosome profiling data.
NOTE: these scripts are meant to be used in conjunction with alt_predict (https://github.com/BiosystemsDataAnalysis/PausePredictionTools).

# PauseMeter
## Description
This script performs multiple actions on the alt_predict_v2 output files (normalized with Normalize_alt_predict.R) supplied in the input folder. It first creates a table containing only the highest peak(s) per gene and after applying the user-defined count threshold (or not if not specified) it outputs an additional table indicating the pause codon of all peaks present in the input file.
The script then uses the output file containing the extracted pause codons and an additional excel file (supplied by the user) containing the codon usage values of each codon of the model organism and generates plots comparing codon usage vs codon occupancy for each file and saves them as SVG images. The data used for producing the plots is also exported as a .csv file.
Finally, a density plot is generated to show which is the fraction of peaks that has been used in the analysis based on the filter_threshold.
## Calculation of pause scores
Codon pause score is calculated by first dividing the count of each peak in the A-site by the ribosomal gene coverage of the gene on which the peak was mapped, obtaining the pause score for each peak. Secondly, the pause scores of peaks mapped on the same codon (intended as trinucleotide sequence) are summed to obtain the codon pause score. Finally, the codon pause score is divided by the total number of RPF reads in ORFs to obtain the normalised codon pause score. 
## Usage
Required libraries:
- dplyr
- magrittr
- stringr
- ggplot2
- readxl
- ggpubr
- tools
- tidyr
- ggbreak

Usage: Rscript PauseMeter.R <input_folder> <output_folder> <codon_table> [filter_threshold]

Arguments:

<input_folder>     : Path to the input folder containing CSV files

<output_folder>    : Path to the output folder where files will be saved

<codon_table>      : Path to the codon usage table Excel file

[filter_threshold] : Optional filter threshold for data (only include data with count >= filter_threshold) [DEFAULT = 1]

# RiboScout
## Description
RiboScout is an alt_predict-based R script that subsets ribosome peaks mapped within defined coordinates relative to the start or end of genes to perform further calculations. 
This script can calculate initiation and ORF density ratios, log2-asymmetry score, and plot metagene profiles.
The coordinates can be defined by the user when selecting the "ratios" method or are fixed when selecting the "asymmetry" or "metagene" methods.
## Method "ratios"
If method "ratios" is entered, the user can choose two values (value_1 and value_2) representing the start and the end of the first range. The values will be summed (or subtracted if the values are < 0) to the start coordinate of each gene (i.e. position 0). The user will choose another two values (value_3 and value_4) representing the start and the end of the second range, with value_3 being added or subtracted to the start position of the gene, and value_4 being added or subtracted to the end position of the gene. Importantly, value_2 must be >= value_1, value_3 must be >= value_2 and > 0, and value_4 must be > value_3.
The script will then perform a two-sided search of ribosomes mapped within the first and then the second range and will combine the two data frames of subsets of peaks and include the calculations of Translation_initiatio_ratio and the ORF_translation_ratio, i.e. the average ribosome density in each range. These two metrics are calculated as follows:
- Translation_initiation_ratio: Sum of counts of all peaks in the first coordinate range / (absolute value of value_1 + absolute value of value_2)
- ORF_translation_ratio: Sum of counts of all peaks in the second coordinate range / (gene_length - value_3 + value_4)

Suggested values are:
- value_1 = -15
- value_2 = 5
- value_3 = 5
- value_4 = 0

With this selection, the first range will include the Ribosome Binding Site (RBS) - for most bacteria - and the start and second codon, while the second range will start from the end of the second codon and stop at the end of each gene.
## Method "asymmetry"
Asymmetry score has been used as an important metric for ribosome profiling methodology testing (Gerashchenko and Gladyshev, 2014; Mohammad et al., 2019).
If method "asymmetry" is entered, the script will use the information on genes coordinates and gene lengths stored in the supplied gene_info_df .csv file to locate the first and second half of each gene, and it will sum the counts of ribosomal peaks on each gene for both halves. Finally, the script will calculate the log2-asymmetry score for each gene as follows:
- log2_asymmetry_score =  log2(Sum_Norm_count_second_half / Sum_Norm_count_first_half)

Additionally, the density_per_nt is calculated as (Sum_Norm_count_first_half + Sum_Norm_count_second_half)/gene_length.

A box plot combining the log2-asymmetry scores of all genes is then generated.
By default, all genes are used for plotting, but the user can generate box plots applying filtering of genes based on their density per nucleotide which is included in the asymmetry score output .csv file. 
## Method "metagene"
Finally, if method "metagene" is entered, the script will subset all ribosome peaks in the range of coordinates between -20 and 200 (with 0 being the start of the gene) and will assign a position to each ribosomal peak relative to the start coordinate of the gene on which it was mapped (relative_position).
The resulting data is used to generate a metagene profile plot.
## Usage
Required libraries:
- dplyr
- readr
- ggplot2

Usage: Rscript RiboScout.R Normalized_alt_predict_file.csv gene_info_df.csv <method> <subtract/add_value_1> <subtract/add_value_2> <subtract/add_value_3> <subtract/add_value_4>

Arguments:

  Normalized_alt_predict_file.csv: output file from alt_predict_v2.py normalized with Normalize_alt_predict.R
  
  gene_info_df.csv: file containing gene information generated with CreateGeneInfo.R
  
  <method>: Analysis method ('ranges', 'asymmetry', or 'metagene')
  
  <subtract/add_value_1>: Value to subtract or add relative to start position of each gene (first range) [Only for 'ranges' method]
  
  <subtract/add_value_2>: Value to subtract or add relative to start position of each gene (first range) [Only for 'ranges' method]
  
  <subtract/add_value_3>: Value to subtract or add relative to start position of each gene (second range) [Only for 'ranges' method]
  
  <subtract/add_value_4>: Value to subtract or add relative to end position of each gene (second range) [Only for 'ranges' method]

# Normalize_alt_predict
## Description
This script normalizes ribosome profiling counts in alt_predict output files by multiplying the counts of each mapped ribosome profiling peak by the ratio 
between the sum of counts of the file with the highest number of reads and the sum of counts of the one that is going to be normalized. In addition, it filters out non coding RNAs that might still be present after mapping. This allows the different samples' output files and plots obtained with the other scripts of this repository to be comparable. 
## Usage
Required libraries:
- dplyr
- stringr

Usage: Rscript normalize_riboseq.R <input_directory> <output_directory>

Arguments:

  <input_directory>   The path to the directory containing input CSV files.
  
  <output_directory>  The path to the directory where normalized CSV files will be saved.

# CreateGeneInfo
## Description
This script combines the sequence information of a reference genome .fasta file and the gene annotations, including gene and locus name, strand direction, etc. into a single data frame that can be used with other scripts of this repository to generate pause score calculations, metagene profiles and more.
The script has been optimized for the reference .fasta and .gff files of B. subtilis NC_000964.3 (files test_annotation.gff and test_genome.fa are included in the PauseMeter repository for testing).
Given that especially .gff files can vary in their format and descriptions depending on the model organism, this script is meant as help to obtain a final dataset that has a format which is compatible with the other scripts of the PauseMeter repository. It is therefore possible that substantial modifications to this script may be required in case a different model organism is used.
## Usage
Required libraries:
- tidyverse
- Biostrings

This script can to be run once in a R environment to create a final gene_info_df .csv file.

# Additional information
A .zip folder containing an Excel file named Test_codon_usage_table.xlsx is provided together with the PauseMeter.R script for testing purposes. Users can create a file with the same layout containing the codon usage values of their model organism.
Additionally, a .csv file named Test_gene_info_df to test RiboScout.R, and the Test_annotation.gff and Test_genome.fa to help with CreateGeneInfo.R are provided.
