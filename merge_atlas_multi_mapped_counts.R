#!/usr/bin/env Rscript
# Clean raw (e.g., EXO) sonde data
# Copyright Jackson M. Tsuji (Neufeld Lab PhD student), 2018
# Created Jan. 17, 2018
# Description: Converts sonde data of various formats (EXO raw, Jared Wolfe, ELA RBR) to a common format for downstream analysis
              # Command line compatible

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("/Users/JTsuji/Research_General/Bioinformatics/02_git/ELA-Sep2017-sonde-data/") # your working directory where files are stored
  SONDE_FILENAME <- "raw-UW/170909_L239.tsv"
  ENCODING <- "UTF-8" # Generally correct
  DATA_TYPE <- "EXO" # EXO, ELA-RBR, or Jared-Wolfe
  LAKE_ID <- "L239" # or state NA to not add this; or state 'auto' to detect from SITE parameter in EXO file
  OUTPUT_FILENAME <- "cleaned/170909_L239_cleaned.tsv"
}
#####################################################

#####################################################
## Load required packages: ##########################
library(getopt)
#####################################################

SCRIPT_VERSION <- "1.0.0"

parse_command_line_input <- function() {
  ### Grab arguments
  # Arguments required:
  # -i input TSV sonde data filename
  # -o output filename
  # -t sonde data type
  # -l lake ID
  # -e file encoding
  params <- matrix(c('input', 'i', 1, "character",
                     'output', 'o', 1, "character",
                     'data_type', 't', 2, "character",
                     'lake_ID', 'l', 2, "character",
                     'file_encoding', 'e', 2, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("clean-sonde-data.R: script for basic parsing of raw sonde data.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n")
    cat(paste("Version: ", SCRIPT_VERSION, "\n\n", sep = ""))
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    cat("Details:\n", "-input\t\t\tFilename for TSV-format raw sonde data. [Required]\n",
        "-output\t\tFilename for output TSV-format file. [Required]\n",
        "-data_type\t\tFormat of input sonde data. Must be 'EXO', 'ELA-RBR', or 'Jared-Wolfe'. [Default: EXO]\n",
        "-lake_ID\t\tOptional: lake ID to add as first column to file. **Set 'auto' to find lake ID based on 'Site' parameter for EXO files only. [Default: ID not added]\n",
        "-file_encoding\t\tText encoding; usually you don't need to worry about this. [Default: UTF-8]\n\n")
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$input) ) {
    stop("Input sonde data filename required. Try -h for help message.")
  }
  if ( is.null(opt$output) ) {
    stop("Output data filename required. Try -h for help message.")
  }
  
  # Add default setting if not provided for other inputs
  if ( is.null(opt$data_type) ) {
    opt$data_type <- "EXO"
  }
  if ( is.null(opt$lake_ID) ) {
    opt$lake_ID <- NA
  }
  if ( is.null(opt$file_encoding) ) {
    opt$file_encoding <- "UTF-8"
  }
  
  # Make variables from provided input and save as global variables (<<-)
  SONDE_FILENAME <<- opt$input
  ENCODING <<- opt$file_encoding
  DATA_TYPE <<- opt$data_type
  LAKE_ID <<- opt$lake_ID
  OUTPUT_FILENAME <<- opt$output
  
}



# Define file locations
atlas_path="L227-2013-6m_annotations.txt"
counts_path="L227-2013-6m_counts.txt"
output_path="L227-2013-6m_annotations_merged.txt"

# Read tables
atlas_tbl <- read.table(file = atlas_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "", quote = "")
counts_tbl <- read.table(file = counts_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

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
library(plyr)
library(dplyr)

atlas_tbl_merged <- dplyr::left_join(atlas_tbl, counts_tbl, by = c("contig_id", "locus_tag"))


# Write output
write.table(atlas_tbl_merged, file = output_path, sep = "\t", col.names = TRUE, row.names = FALSE)




main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  cat(paste("Running clean-sonde-data.R, version: ", SCRIPT_VERSION, "\n\n", sep = ""))
  cat(paste("Input filename:", SONDE_FILENAME, "\n"))
  cat(paste("Output filename:", OUTPUT_FILENAME, "\n"))
  cat(paste("Sonde data format:", DATA_TYPE, "\n"))
  cat(paste("Lake ID to append:", LAKE_ID, "\n"))
  cat(paste("Input file encoding:", ENCODING, "\n"))
  cat("\n")
  
  # Import data
  cat("Loading data...\n")
  raw_data <- load_raw_file(SONDE_FILENAME, ENCODING)
  
  cat(paste("Parsing data as ", DATA_TYPE, "...\n", sep = ""))
  parsed_data <- parse_sonde_data(raw_data, DATA_TYPE)
  
  # Assign lake ID if a value was given
  if ( is.na(LAKE_ID) == FALSE ) {
    # Grab lake ID from site parameter if "auto" is selected for EXO data
    if (LAKE_ID == "auto") {
      if (DATA_TYPE == "EXO") {
        LAKE_ID <- get_exo_site_ID(raw_data)
        cat(paste("Found lake ID '", LAKE_ID, "' in EXO file parameters.\n", sep = ""))
      } else {
        stop("Cannot use 'auto' function for LAKE_ID except for EXO data type.")
      }
    }
    
    cat(paste("Assigning lake ID '", LAKE_ID, "'...\n", sep = ""))
    parsed_data <- assign_lake_ID(parsed_data, LAKE_ID)
  }
  
  # Write output file
  cat("Writing output...\n")
  write.table(parsed_data, file = OUTPUT_FILENAME, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  cat("\nclean-sonde-data.R: finished.\n")
}

main()
