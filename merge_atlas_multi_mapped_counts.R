#!/usr/bin/env Rscript
# Clean raw (e.g., EXO) sonde data
# Copyright Jackson M. Tsuji (Neufeld Lab PhD student), 2018
# Created Mar. 7, 2018
# Description: Adds featureCounts information for additional read mapped samples onto the ATLAS annotations table. Built as part of the custom ATLAS coassembly module in atlas-extensions.
              # Command line compatible

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("/Users/JTsuji/Research_General/Bioinformatics/02_git/") # your working directory where files are stored
  ATLAS_TABLE <- "L227-2013-6m_annotations.txt"
  FEATURECOUNTS_TABLE <- "L227-2013-6m_counts.txt"
  OUTPUT_FILE <- "L227-2013-6m_annotations_merged.txt"
}
#####################################################

#####################################################
## Load required packages: ##########################
library(getopt)
library(plyr)
suppressMessages(library(dplyr))
#####################################################

SCRIPT_VERSION <- "1.0.22-coassembly-r3"

parse_command_line_input <- function() {
  ### Grab arguments
  # Arguments required:
  # -a input ATLAS annotations TSV table
  # -f input featureCounts TSV table
  # -o output TSV filename
  params <- matrix(c('atlas_table', 'a', 1, "character",
                     'featureCounts_table', 'f', 1, "character",
                     'output_file', 'o', 1, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("merge_atlas_multi_mapped_counts.R: merges additional featureCounts tables to the ATLAS annotations table.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n")
    cat(paste("Version: ", SCRIPT_VERSION, "\n\n", sep = ""))
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    cat("Details:\n", "-atlas_table\t\t\tFilepath for TSV-format ATLAS annotations table. [Required]\n",
        "-featureCounts_table\t\tFilepath for TSV-format featureCounts count file. [Required]\n",
        "-output_file\t\t\tFilepath for output TSV-format annotations file. [Required]\n\n")
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$atlas_table) ) {
    stop("Input ATLAS annotations table filepath required. Try -h for help message.")
  }
  if ( is.null(opt$featureCounts_table) ) {
    stop("Input featureCounts count table filepath required. Try -h for help message.")
  }
  if ( is.null(opt$output_file) ) {
    stop("Output filepath required. Try -h for help message.")
  }
  
  # Make variables from provided input and save as global variables (<<-)
  ATLAS_TABLE <<- opt$atlas_table
  FEATURECOUNTS_TABLE <<- opt$featureCounts_table
  OUTPUT_FILE <<- opt$output_file
  
}

merge_tables <- function(atlas_table_filepath, featureCounts_table_filepath) {
  # Quick description: loads, modifies, and merges tables
  # Input the filepaths for the two tables (supplied by user)
  
  # Read tables
  atlas_tbl <- read.table(file = atlas_table_filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "", quote = "")
  counts_tbl <- read.table(file = featureCounts_table_filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  # Change counts column names
  colnames(counts_tbl)[1:2] <- c("locus_tag", "contig_id")
  
  # Simplify sample names
  num_mapped_samples <- length(colnames(counts_tbl)) - 6
  for (i in 1:num_mapped_samples) {
    sample_col <- i + 6
    sample_name_old <- colnames(counts_tbl)[sample_col]
    
    sample_name_new <- strsplit(sample_name_old, split = "sequence_alignment")[[1]][1]
    sample_name_new <- substr(sample_name_new, start = 1, stop = nchar(sample_name_new) - 1)
    sample_name_new <- paste(sample_name_new, "_count", sep = "")
    
    colnames(counts_tbl)[sample_col] <- sample_name_new
  }
  
  # Remove unneeded columns
  counts_tbl <- counts_tbl[,-c(3:6)]
  
  # Merge tables
  atlas_tbl_merged <- dplyr::left_join(atlas_tbl, counts_tbl, by = c("contig_id", "locus_tag"))
  
  return(atlas_tbl_merged)

}

main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  cat(paste("Running merge_atlas_multi_mapped_counts.R, version: ", SCRIPT_VERSION, "\n\n", sep = ""))
  cat(paste("ATLAS table filepath:", ATLAS_TABLE, "\n"))
  cat(paste("featureCounts table filepath:", FEATURECOUNTS_TABLE, "\n"))
  cat(paste("Output filepath:", OUTPUT_FILE, "\n"))
  cat("\n")
  
  # Import data
  cat("Loading and merging data...\n")
  merged_annotations_table <- merge_tables(ATLAS_TABLE, FEATURECOUNTS_TABLE)
  
  # Write output
  cat("Writing output...\n")
  write.table(merged_annotations_table, file = OUTPUT_FILE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  cat("\nmerge_atlas_multi_mapped_counts.R: finished.\n")
}

main()
