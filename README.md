# Introduction
This script performs multiple actions on the alt_predict_v2 output files supplied in the input folder. It first creates a table containing only the highest peak(s) per gene and after applying the user-defined count threshold (or not if not specified) it outputs an additional table indicating the pause codon of all peaks present in the input file.
The script then uses the output file containing the extracted pause codons and an additional excel file (supplied by the user) containing the codon usage values of each codon of the model organism and generates plots comparing codon usage vs codon occupancy for each file and saves them as SVG images. The data used for producing the plots is also exported as a .csv file.
Finally, a density plot is generated to show which is the fraction of peaks that was used in the analysis based on the filter_threshold.
# Calculation of pause scores
Codon pause score is calculated by first dividing the count of each peak in the A-site by the ribosomal gene coverage of the gene on which the peak was mapped, obtaining the pause score for each peak. Secondly, the pause scores of peaks mapped on the same codon (intended as trinucleotide sequence) are summed to obtain the codon pause score. Finally, the codon pause score is divided by the total number of RPF reads in ORFs to obtain the normalised codon pause score. 
# PauseMeter
Usage: Rscript PauseMeter.R <input_folder> <output_folder> <codon_table> [filter_threshold]

Arguments:

<input_folder>     : Path to the input folder containing CSV files

<output_folder>    : Path to the output folder where files will be saved

<codon_table>      : Path to the codon usage table Excel file

[filter_threshold] : Optional filter threshold for data (only include data with count >= filter_threshold) [DEFAULT = 1]
# Additional information
An Excel file named Test_codon_usage_table.xlsx is provided together with the script for testing purposes. Users can create a file with the same layout containing the codon usage values of their model organism.
