#!/usr/bin/env Rscript
# Copyright Jackson M. Tsuji (Neufeld Lab PhD student), 2018
# Created Mar. 7, 2018
# Description: Adds contig_length information onto the ATLAS annotations table. Built as part of the custom ATLAS coassembly module in atlas-extensions.
              # Command line compatible

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- FALSE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("C://Users/Jackson Tsuji/Downloads/") # your working directory where files are stored
  ATLAS_TABLE <- "vs2_p2/02_amoA/CA-NE8-2010_annotations_amoA.tsv"
  CONTIG_LENGTH_TABLE <- "contigs/postfilter_coverage_stats_NE8.txt"
  OUTPUT_FILE <- "vs2_p2/02_amoA/CA-NE8-2010_annotations_amoA_with_contig_lengths.tsv"
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
  # -f input contig_length TSV table
  # -o output TSV filename
  params <- matrix(c('atlas_table', 'a', 1, "character",
                     'contig_length_table', 'l', 1, "character",
                     'output_file', 'o', 1, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("add_contig_lengths_to_annotations_table.R: merges additional contig_length tables to the ATLAS annotations table.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n")
    cat(paste("Version: ", SCRIPT_VERSION, "\n\n", sep = ""))
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    cat("Details:\n", "-atlas_table\t\t\tFilepath for TSV-format ATLAS annotations table. [Required]\n",
        "-contig_length_table\t\tFilepath for TSV-format postfilter_coverage_stats.txt file. [Required]\n",
        "-output_file\t\t\tFilepath for output TSV-format annotations file. [Required]\n\n")
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$atlas_table) ) {
    stop("Input ATLAS annotations table filepath required. Try -h for help message.")
  }
  if ( is.null(opt$contig_length_table) ) {
    stop("Input postfilter_coverage_stats.txt contig length table filepath required. Try -h for help message.")
  }
  if ( is.null(opt$output_file) ) {
    stop("Output filepath required. Try -h for help message.")
  }
  
  # Make variables from provided input and save as global variables (<<-)
  ATLAS_TABLE <<- opt$atlas_table
  CONTIG_LENGTH_TABLE <<- opt$contig_length_table
  OUTPUT_FILE <<- opt$output_file
  
}

merge_tables <- function(atlas_table_filepath, contig_length_table_filepath) {
  # Quick description: loads, modifies, and merges tables
  # Input the filepaths for the two tables (supplied by user)
  
  # Read tables
  atlas_tbl <- read.table(file = atlas_table_filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
  contig_length_tbl <- read.table(file = contig_length_table_filepath, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  
  # Reduce contig length table to needed fields, and rename columns
  # TODO: add sanity check that the provided table is indeed postfilter_coverage_stats.txt
  contig_length_tbl <- contig_length_tbl[,c(1,3)]
  colnames(contig_length_tbl) <- c("contig_id", "contig_length")
  
  # Merge tables
  atlas_tbl_merged <- dplyr::left_join(atlas_tbl, contig_length_tbl, by = "contig_id")
  
  return(atlas_tbl_merged)

}

main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  cat(paste("Running add_contig_lengths_to_annotations_table.R, version: ", SCRIPT_VERSION, "\n\n", sep = ""))
  cat(paste("ATLAS table filepath:", ATLAS_TABLE, "\n"))
  cat(paste("Contig length table filepath:", CONTIG_LENGTH_TABLE, "\n"))
  cat(paste("Output filepath:", OUTPUT_FILE, "\n"))
  cat("\n")
  
  # Import data
  cat("Loading and merging data...\n")
  merged_annotations_table <- merge_tables(ATLAS_TABLE, CONTIG_LENGTH_TABLE)
  
  # Write output
  cat("Writing output...\n")
  write.table(merged_annotations_table, file = OUTPUT_FILE, sep = "\t", col.names = TRUE, row.names = FALSE)
  
  cat("\nadd_contig_lengths_to_annotations_table.R: finished.\n")
}

main()
