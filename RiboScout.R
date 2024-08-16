##

# GitHub: aagnolin (https://github.com/aagnolin/PauseMeter)

# Author:  Alberto Agnolin (alberto.agnolin.1@gmail.com)
# Name of script: RiboScout
# Summary: alt_predict-based tool that subsets ribosome pausing peaks based on user-defined coordinates to calculate initiation and ORF translation/pausing ratios.
# In addition, the script can calculate log2 asymmetry score and output metagene profile plots

# NOTE: this script is meant to be used in conjunction with alt_predict (https://github.com/BiosystemsDataAnalysis/PausePredictionTools) and Normalize_alt_predict.R.
# Read the README.md file for more information

##

# Define usage function
usage <- function() {
  cat("Usage: Rscript PeakFinder_two_sided.R Normalized_alt_predict_file.csv gene_info_df.csv <method> <subtract/add_value_1> <subtract/add_value_2> <subtract/add_value_3> <subtract/add_value_4>\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  Normalized_alt_predict_file.csv: output file from alt_predict_v2.py normalized with Normalize_alt_predict.R\n")
  cat("  gene_info_df.csv: file containing gene information generated with CreateGeneInfo.R\n")
  cat("  <method>: Analysis method ('ranges' or 'asymmetry')\n")
  cat("  <subtract/add_value_1>: Value to subtract or add relative to start position of each gene (first range) [Only for 'ranges' method]\n")
  cat("  <subtract/add_value_2>: Value to subtract or add relative to start position of each gene (first range) [Only for 'ranges' method]\n")
  cat("  <subtract/add_value_3>: Value to subtract or add relative to start position of each gene (second range) [Only for 'ranges' method]\n")
  cat("  <subtract/add_value_4>: Value to subtract or add relative to end position of each gene (second range) [Only for 'ranges' method]\n")
}

# Extract command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments is provided
if (length(args) < 3 || length(args) > 7) {
  cat("Incorrect number of arguments. Check usage below.\n\n")
  usage()
  quit(status = 1, save = "no")
}

# Assign arguments to variables
alt_predict_file <- args[1]
gene_info_df_file <- args[2]
method <- args[3]  # Last argument is method

# Initialize value arguments
value_1 <- value_2 <- value_3 <- value_4 <- NA

if (method == "ranges") {
  if (length(args) != 7) {
    cat("Incorrect number of arguments. Check usage below.\n\n")
    usage()
    quit(status = 1, save = "no")
  }
  value_1 <- as.numeric(args[4])
  value_2 <- as.numeric(args[5])
  value_3 <- as.numeric(args[6])
  value_4 <- as.numeric(args[7])
} else if (method == "asymmetry") {
  if (length(args) != 3) {
    cat("Incorrect number of arguments. Check usage below.\n\n")
    usage()
    quit(status = 1, save = "no")
  }
} else if (method == "metagene") {
  if (length(args) != 3) {
    cat("Incorrect number of arguments. Check usage below.\n\n")
    usage()
    quit(status = 1, save = "no")
  }
} else {
  usage()
  quit(status = 1, save = "no")
}

# Load necessary libraries
suppressPackageStartupMessages(library(dplyr))
library(readr)

# Load input files
input_data <- read_csv(alt_predict_file, show_col_types = F)

# Load gene info
gene_info_df <- read_csv(gene_info_df_file, show_col_types = F)

# Assign names to output files
input_file_name <- tools::file_path_sans_ext(basename(alt_predict_file))
output_left_side <- paste0(input_file_name, "_output_left_side.csv")
output_metagene <- paste0(input_file_name, "_metagene.csv")
output_right_side <- paste0(input_file_name, "_output_right_side.csv")
merged_ratios <- paste0(input_file_name, "_merged_ratios.csv")
merged_halves <- paste0(input_file_name, "_merged_asymmetry.csv")
plot_name <- paste0(input_file_name, "_metagene_plot.pdf")


# First range
# Choose target positions based on the method (assumes that the provided gene_info_df takes strand direction of genes into consideration)
if (method == "ranges") {
  target_sequences <- gene_info_df %>%
    mutate(target_start = StartPosition + value_1,
           target_end = StartPosition + value_2) %>%
    select(locus_tag, target_start, target_end, StartPosition)
} else if (method == "asymmetry") {
  target_sequences <- gene_info_df %>%
    mutate(target_start = StartPosition,
           target_end = floor(StartPosition + (gene_length / 2))) %>%
    select(locus_tag, target_start, target_end, gene_length, StartPosition)
} else if (method == "metagene") {
  target_sequences <- gene_info_df %>%
    mutate(target_start = StartPosition - 20,
           target_end = StartPosition + 200) %>%
    select(locus_tag, target_start, target_end, StartPosition)
}  
# Merge input data and target sequences
input_data_merge <- merge(input_data, target_sequences, by = "locus_tag", all = TRUE)

# Sequences outside the ORFs in alt_predict output files do not have a locus tag assigned, but this is necessary if we want to analyze peaks in the upstream region of the ORFs
## Sort target_sequences by target_start for binary search
target_sequences_sorted <- target_sequences[order(target_sequences$target_start), ]

if (method == "ranges") {
  cat("Finding peaks in first position range...\n")
} else if (method == "asymmetry") {
  cat("Finding peaks in first half of genes...\n")
} else if (method == "metagene") {
  cat("Finding peaks in metagene profile range...")
}

## Binary search to assign locus_tag based on position
assign_locus_tag <- function(position) {
  left <- 1
  right <- nrow(target_sequences_sorted)
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (position >= target_sequences_sorted[mid, "target_start"] && position <= target_sequences_sorted[mid, "target_end"]) {
      return(target_sequences_sorted[mid, c("locus_tag", "target_start", "target_end", "StartPosition")])
    } else if (position < target_sequences_sorted[mid, "target_start"]) {
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }
  
  return(NA)
}

## Assign locus_tag, target_start, and target_end for rows with NA locus_tag that are in the target position range
for (i in 1:nrow(input_data_merge)) {
  if (is.na(input_data_merge[i, "locus_tag"])) {
    result <- assign_locus_tag(input_data_merge[i, "position"])
    input_data_merge[i, c("locus_tag", "target_start", "target_end", "StartPosition")] <- result
  }
}

# Filter data that is between the provided positions and calculate position of peak relative to start nucleotide 
input_data_target <- group_by(input_data_merge, locus_tag) %>% filter(position >= target_start & position <= target_end)
input_data_target <- mutate(input_data_target, relative_position = position - StartPosition) %>%  select(-sequence)

# Add gene length to filtered data
input_data_target <- merge(input_data_target, gene_info_df, by = "locus_tag", all = TRUE) %>% select(-"gene_length.x")

# Filter out ncRNA that may have been found outside their loci as they were removed in the normalized alt_predict files
input_data_target <- input_data_target %>% filter(!grepl("^BSU_", locus_tag))

# Write output file of first range for methods "ranges" and "asymmetry"
if (method == "ranges" | method == "asymmetry") {
write_csv(input_data_target, output_left_side)
}

# Calculate initiation ratio if the method is "ranges", 
# only sum Norm_count if the method is "asymmetry", 
# or sum Norm_count and include relative position if the method is "metagene"
if (method == "ranges") {
  # Sum Norm_count in target positions for each gene, then divide by the target sequence length
  df_Ribo_reads_target_1 <- input_data_target %>% 
    group_by(locus_tag) %>% 
    summarize(Sum_Norm_count_initiation = sum(Norm_count)) %>% 
    mutate(Translation_initiation_Ratio = Sum_Norm_count_initiation / ((abs(value_1) + abs(value_2))))
  
} else if (method == "asymmetry") {
  # Sum Norm_count in the first half of the gene if the method is "asymmetry"
  df_Ribo_reads_target_1 <- input_data_target %>% 
    group_by(locus_tag) %>% 
    summarize(Sum_Norm_count_first_half = sum(Norm_count))
} else if (method == "metagene") {
  # Sum Norm_count for metagene target region and include relative position
  df_Ribo_reads_target_1 <- input_data_target %>% 
    group_by(locus_tag) %>% 
    summarize(Sum_Norm_count = sum(Norm_count), relative_position = relative_position)
  write_csv(df_Ribo_reads_target_1, output_metagene)
  # Create metagene profile plot
  cat("generating metagene profile plot...\n")
  library(ggplot2)
  p <- ggplot(data = df_Ribo_reads_target_1,
              mapping = aes(x = relative_position)) +
    geom_freqpoly(bins = 220, linewidth = 0.8) +
    geom_bar(alpha = 0.5, width = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "black", alpha = 0.6, linewidth = 0.7) +
    scale_x_continuous(limits = c(-20, 200)) +
    scale_y_continuous(expand = c(0,0)) +
    theme_bw() +
    labs(x = "Position") +
    theme(plot.subtitle = element_text(face = "bold"),
          plot.caption = element_text(face = "bold"),
          axis.title = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          legend.title = element_text(face = "bold"),
          legend.key = element_rect(linetype = "solid")) +
    theme(axis.ticks = element_line(linewidth = 0.5)) + 
    theme(panel.background = element_rect(fill = NA)) + 
    theme(axis.line = element_line(linetype = "solid")) + 
    theme(axis.line = element_line(linetype = "blank"), 
          panel.grid.major = element_line(colour = "gray89",
                                          linetype = "blank"),
          panel.grid.minor = element_line(linetype = "blank"),
          panel.background = element_rect(linetype = "solid"))
  
  # Generate plot
  ggsave(plot = p, filename = plot_name, device = "pdf")
  
  cat("complete")
  # Stop the script execution if the "metagene" method is used
  quit(status = 0, save = "no")
}

# Second side (only for method "ranges" and "asymmetry")
if (method == "ranges") {
  # Second range
  target_sequences <- gene_info_df %>%
    mutate(target_start = StartPosition + value_3,
           target_end = EndPosition + value_4) %>%
    select(locus_tag, target_start, target_end, StartPosition)
} else if (method == "asymmetry") {
  # Second half
  target_sequences <- gene_info_df %>%
    mutate(target_start = ceiling(StartPosition + (gene_length / 2)),
           target_end = EndPosition) %>%
    select(locus_tag, target_start, target_end, gene_length, StartPosition)
}

# Merge input data and target sequences
input_data_merge <- merge(input_data, target_sequences, by = "locus_tag", all = TRUE)

# Sequences outside the ORFs in alt_predict output files do not have a locus tag assigned, but this is necessary if we want to analyze peaks in the upstream region of the ORFs
## Sort target_sequences by target_start for binary search
target_sequences_sorted <- target_sequences[order(target_sequences$target_start), ]

if (method == "ranges") {
  cat("Finding peaks in second position range...\n")
} else if (method == "asymmetry") {
  cat("Finding peaks in second half of genes...\n")
}

## Binary search to assign locus_tag based on position
assign_locus_tag <- function(position) {
  left <- 1
  right <- nrow(target_sequences_sorted)
  
  while (left <= right) {
    mid <- floor((left + right) / 2)
    
    if (position >= target_sequences_sorted[mid, "target_start"] && position <= target_sequences_sorted[mid, "target_end"]) {
      return(target_sequences_sorted[mid, c("locus_tag", "target_start", "target_end", "StartPosition")])
    } else if (position < target_sequences_sorted[mid, "target_start"]) {
      right <- mid - 1
    } else {
      left <- mid + 1
    }
  }
  
  return(NA)
}

## Assign locus_tag, target_start, and target_end for rows with NA locus_tag that are in the target position range
for (i in 1:nrow(input_data_merge)) {
  if (is.na(input_data_merge[i, "locus_tag"])) {
    result <- assign_locus_tag(input_data_merge[i, "position"])
    input_data_merge[i, c("locus_tag", "target_start", "target_end", "StartPosition")] <- result
  }
}

# Filter data that is between the provided positions and calculate position of peak relative to start nucleotide
input_data_target <- group_by(input_data_merge, locus_tag) %>% filter(position >= target_start & position <= target_end)
input_data_target <- mutate(input_data_target, relative_position = position - StartPosition) %>%  select(-sequence)

# Add gene length to filtered data
input_data_target <- merge(input_data_target, gene_info_df, by = "locus_tag", all = TRUE) %>% select(-"gene_length.x")

# Filter out ncRNA that may have been found outside their loci as they were removed in the normalized alt_predict files
input_data_target <- input_data_target %>% filter(!grepl("^BSU_", locus_tag))

# Write output file 2
write_csv(input_data_target, output_right_side)

# Calculate ORF translation ratio and merge the two sides if the method is "ranges", 
# only sum Norm_count of second half and merge the two halves if the method is "asymmetry"
if (method == "ranges") {
  # Sum Norm_count in target positions for each gene, then divide by the target sequence length
  df_Ribo_reads_target_2 <- input_data_target %>% 
    group_by(locus_tag) %>% 
    summarize(Sum_Norm_count_ORF = sum(Norm_count)) %>% 
    merge(gene_info_df, by = "locus_tag") %>% 
    mutate(ORF_translation_Ratio = Sum_Norm_count_ORF / (gene_length - value_2))
  # Merge the output data frames of the two ranges
  Merged_ratios <- merge(df_Ribo_reads_target_1, df_Ribo_reads_target_2, by = "locus_tag", all = TRUE) %>% 
    select(c(-"StartPosition", -"EndPosition", -"Sequence", -"Strand"))
  # Write final output
  write_csv(Merged_ratios, merged_ratios)
} else if (method == "asymmetry") {
  # Sum Norm_count in the second half of the gene if the method is "asymmetry"
  df_Ribo_reads_target_2 <- input_data_target %>% 
    group_by(locus_tag) %>% 
    summarize(Sum_Norm_count_second_half = sum(Norm_count))
  # Merge the output data frames of the two halves
  merged_halves_df <- merge(df_Ribo_reads_target_1, df_Ribo_reads_target_2, by = "locus_tag", all = TRUE)
  # Calculate log2 asymmetry score
  merged_halves_df <- merged_halves_df %>% 
    mutate(log2_asymmetry_score = log2(Sum_Norm_count_second_half / Sum_Norm_count_first_half))

  # Write final output
  write_csv(merged_halves_df, merged_halves)
}
cat("complete\n")