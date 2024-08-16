##

# GitHub: aagnolin (https://github.com/aagnolin/PauseMeter)

# Author:  Alberto Agnolin (alberto.agnolin.1@gmail.com)
# Name of script: CreateGeneInfo.R
# Summary: 

# NOTE: This script combines the sequence information of a reference genome .fasta file and the gene annotations, 
# including gene and locus name, strand direction, etc. into a single data frame that can be used to generate
# pause score calculations, metagene profiles and more using scripts PauseMeter.R and RiboScout.R.
# The script has been optimized for the reference .fasta and .gff files of B. subtilis NC_000964.3 
# (files test_annotation.gff and test_genome.fa are included in the PauseMeter repository for testing).
# Given that especially .gff files can vary in their format and descriptions depending on the model organism, 
# this script is meant as a help to obtain a final dataset that has a format which is compatible with the other
# scripts of the PauseMeter repository. It is therefore possible that substantial modifications to this script 
# may be required in case a different model organism is used.

# Read the README.md file for more information

##


library(tidyverse)
library(Biostrings)

# Read the reference genome in FASTA format
genome <- readDNAStringSet("PATH/genome.fasta")

# Read the GFF file
gff_lines <- readLines("PATH/annotation.gff")

# Extract information for genes 
locus_tags <- list()
start_positions <- list()
end_positions <- list()
sequences <- DNAStringSet()
strands <- character()
gene_names <- list()  # Initialize list for gene names

# Process each line of the GFF file
for (line in gff_lines) {
  # Check if the line contains gene information
  if (grepl("gene\\s", line)) {
    # Split the line by tabs
    line_parts <- unlist(strsplit(line, "\t"))
    
    # Extract the relevant information
    locus_tag <- gsub(".*locus_tag=(\\S+).*", "\\1", line_parts[9])
    
    # Remove semicolon and text after it
    locus_tag <- gsub(";.*", "", locus_tag)
    
    start_pos <- as.numeric(line_parts[4])
    end_pos <- as.numeric(line_parts[5])
    
    # Extract the strand information from the 7th column
    strand <- line_parts[7]
    
    # Extract gene name and remove unwanted string after it
    gene_info <- gsub(".*gene=(\\S+);.*", "\\1", line_parts[9])
    gene_name <- gsub(";.*", "", gene_info)
    
    # Check if the LocusTag starts with "B" (specific for B. subtilis annotations)
    if (substring(locus_tag, 1, 1) == "B") {
      # Append the values to the lists
      locus_tags <- append(locus_tags, locus_tag)
      start_positions <- append(start_positions, start_pos)
      end_positions <- append(end_positions, end_pos)
      strands <- append(strands, strand)
      gene_names <- append(gene_names, gene_name)  # Append gene name
      
      # Extract the DNA sequence for the current locus tag and adjust for the strand direction
      if (strand == "+") {
        seq <- subseq(genome, start_pos, end_pos)
      } else if (strand == "-") {
        seq <- reverseComplement(subseq(genome, start_pos, end_pos))
      } else {
        # Handle unrecognized strand information
        seq <- DNAString("")
      }
      
      sequences <- c(sequences, seq)
    }
  }
}

# Create a data frame with all combined information for genes
gene_info_df <- data.frame(
  locus_tag = unlist(locus_tags),
  StartPosition = unlist(start_positions),
  EndPosition = unlist(end_positions),
  Sequence = as.character(sequences),
  Strand = unlist(strands),
  gene = unlist(gene_names) 
) %>% mutate(gene_length = EndPosition - StartPosition)

# Write gene_info_df
write_csv(gene_info_df, "PATH/gene_info_df.csv")