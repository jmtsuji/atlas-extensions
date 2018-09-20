#!/usr/bin/env Rscript
# calculate_coverage_stats.R
# Copyright Jackson M. Tsuji, 2018
# Neufeld Research Group
# See required R packages below.

### Script logic
# Input: read_length and samtools_coverage_stats.tsv (run with samtools depth -aa), and base_output_name
# Mask bases at the start and end of contigs up to length read_length, because mapping stats will be distorted (low) here:
# Summarise average/sd depth per UNmasked contig, plus occurences of 0's per contig
# Summarise average/sd depth per masked contig, plus occurences of 0's per contig
# Summarise overall average/sd coverage for both masked and unmasked, plus % of contigs with 0's above some threshold
# Report per-contig stats as a supplementary table file (written to file - OUTPUT)
# Report overall stats summary as a table column as STDOUT.

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("/home/jmtsuji/Downloads") # your working directory where files are stored
  
  # Input file (should be the output of 'samtools depth -aa' when run on a single BAM file)
  input_coverage_table_filename <- c("example_samtools_coverage_depth.tsv")
  
  # Output file for supplementary per-contig stats. Optional. Set to NULL for the script to determine an output name based on the input table name.
  output_contig_stats_table_filename <- NULL
  
  # Expected read length for the metagenomic dataset
  read_length <- 250
  
  # ID of input bin and metagenome
  bin_ID <- "bin1"
  metagenome_ID <= "metagenome1"
  
  # OR advanced feature: set the above two to 'NULL' and use the info from the input filename to tease out the metagenome_ID and bin_ID, separated by a delimiter defined here.
  name_delimiter <- "_to_"
  
  # Advanced feature: how many 0's across the positions of a contig does it take for the script to consider that contig to have a zero-coverage event?
  zero_coverage_threshold <- 1
  
}
#####################################################

#####################################################
## Load required packages: ##########################
library(getopt)
library(tools)
library(plyr)
suppressMessages(library(dplyr))
#####################################################

# Function: assign command line input to variables, or throw help message
parse_command_line_input <- function() {
  
  # Define flags
  params <- matrix(c('samtools_coverage_table', 'i', 1, "character",
                     'output_contig_stats_table', 'o', 2, "character",
                     'read_length', 'l', 1, "numeric",
                     'bin_ID', 'b', 2, "character",
                     'metagenome_ID', 'm', 2, "character",
                     'ID_delimiter', 'd', 2, "logical",
                     'zero_coverage_threshold', 'z', 2, "numeric",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("calculate_coverage_stats.R: Summarizes 'samtools depth' coverage information forunassembled reads mapped to genome bins.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n\n")
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    cat("Details:\n", "--samtools_coverage_table\tFilepath for TSV-format coverage table output by 'samtools depth -aa' for a single BAM file. [Required]\n",
        "--output_contig_stats_table\tFilepath for the TSV-format per-contig coverage stats table that this script will write [Optional - can predict based on input filename]\n",
        "--read_length\t\t\tExpected length of individual reads for the unassembled reads (e.g., 200). [Required]\n",
        "--bin_ID\t\t\tID of the input genome bin (e.g., 'Bin001') [Optional]\n",
        "--metagenome_ID\t\tID of the metagenome used for mapping reads (e.g., 'Metagenome001'). [Optional]\n",
        "--ID_delimiter\t\t\tAdvanced feature: if you want to parse out the bin_ID and metagenome_ID from the filename of the coverage_table, then input the string delimiter between the bin_ID and metagenome_ID in the filename (e.g., '_to_' for 'Metagenome001_to_Bin001.tsv'). Filename MUST contain 'metagenome_ID[delimiter]bin_ID.extension' and NOTHING ELSE for this to work. [Optional]\n",
        "--zero_coverage_threshold\tAdvanced feature: how many 0's across the positions of a contig does it take for the script to consider that contig to have a zero-coverage event? [Optional; default 1]\n\n")
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$samtools_coverage_table) ) {
    stop("Coverage table from 'samtools depth -aa' required. Try -h for help message.")
  }
  if ( is.null(opt$read_length) ) {
    stop("Expected read length required. Try -h for help message.")
  }
  
  # Determine the ID input mode or throw an error
  too_much_input <- !(is.null(opt$bin_ID)) && !(is.null(opt$metagenome_ID)) && !(is.null(opt$ID_delimiter))
  if ( too_much_input == TRUE ) {
    stop("Too many input options provided for sample ID. Try -h for help message.")
  } else if ( !(is.null(opt$bin_ID)) && !(is.null(opt$metagenome_ID)) ) {
    parse_from_input_filename <<- FALSE
  } else if ( !(is.null(opt$ID_delimiter)) ) {
    parse_from_input_filename <<- TRUE
  } else {
    stop("Not enough input given to determine the bin_ID and metagenome_ID. Try -h for help message.")
  }
  
  # Set defaults
  if ( is.null(opt$zero_coverage_threshold) == TRUE ) {
    opt$zero_coverage_threshold <- 1
  }
  if ( is.null(opt$output_contig_stats_table) == TRUE ) {
    opt$output_contig_stats_table <- paste(file_path_sans_ext(opt$samtools_coverage_table), 
                                           "_per_contig_stats.tsv", sep = "")
  }
  
  # Make variables from provided input and save as global variables (<<-)
  input_coverage_table_filename <<- opt$samtools_coverage_table
  output_contig_stats_table_filename <<- opt$output_contig_stats_table
  read_length <<- opt$read_length
  zero_coverage_threshold <<- opt$zero_coverage_threshold
  if ( parse_from_input_filename == FALSE ) {
    bin_ID <<- opt$bin_ID
    metagenome_ID <<- opt$metagenome_ID
  } else {
    name_delimiter <<- opt$ID_delimiter
  }
  
}


# Function: helper function for timestamps. Returns a length 1 character vector.
ts <- function() {
  
  datetime_utc <- format(as.POSIXlt(Sys.time(), tz = "UTC"), "%a %d %b %Y %H:%M:%S %Z")
  
  date_message <- paste("[ ", datetime_utc, " ]: ", sep = "")
  
  return(date_message)
}


# Function: read samtools coverage table. Returns a data frame (with headers)
read_coverage_table <- function(table_filename) {
  coverage_table <- read.table(file = table_filename, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  
  # Assign header names (HARD-CODED)
  colnames(coverage_table) <- c("contig_ID", "position", "coverage")
  
  return(coverage_table)
}


# Function: calculate the length of each contig. Returns a data frame.
calculate_contig_lengths <- function(coverage_table) {
  
  contig_lengths <- dplyr::group_by(coverage_table, contig_ID)
  contig_lengths <- dplyr::summarise(contig_lengths, contig_length = max(position))
  
  return(contig_lengths)
  
}


# Function: remove the ends of contigs up to the read_length, because coverage stats can be skewed in these regions. Returns a data frame of same format as input.
filter_contig_ends <- function(coverage_table, contig_lengths, read_length) {
  # Map the length onto the coverage table
  coverage_table <- dplyr::left_join(coverage_table, contig_lengths, by = c("contig_ID"))
  
  # Calculate the reverse positioning numbers
  coverage_table$reverse_position <- coverage_table$contig_length - coverage_table$position + 1
  # This starts at 1 for the back nucleotide in the contig and counts up from there
  
  # Filter out contigs shorter than 2*read_length
  coverage_table_filtered <- dplyr::filter(coverage_table, contig_length > (2 * read_length))
  
  # Filter out positions at the flanking ends of the contigs
  coverage_table_filtered <- dplyr::filter(coverage_table_filtered, position > read_length & reverse_position > read_length)
  
  # Remove unnecessary columns
  coverage_table_filtered$contig_length <- NULL
  coverage_table_filtered$reverse_position <- NULL
  
  return(coverage_table_filtered)
  
}


# Function: calculates mean and sd of coverage on a per-contig basis, and adds contig length info. Returns a data frame.
calculate_depth_per_contig <- function(coverage_table, contig_lengths) {
  
  # Calculate depth
  depth_per_contig <- dplyr::group_by(coverage_table, contig_ID)
  depth_per_contig <- dplyr::summarise(depth_per_contig, coverage_mean = mean(coverage), coverage_sd = sd(coverage))
  
  # Join on contig length information
  depth_per_contig <- dplyr::left_join(depth_per_contig, contig_lengths, by = "contig_ID")
  
  # Re-order columns
  depth_per_contig <- dplyr::select(depth_per_contig, contig_ID, contig_length, everything())
  
  return(depth_per_contig)
  
}


# Funtion: determine the number of 0-coverage sites across each contig. Returns a data frame.
count_zero_coverage_events <- function(coverage_table) {
  
  # Get all zero coverage events
  zero_coverage_events <- dplyr::filter(coverage_table, coverage == 0)
  
  # Count the totals per contig
  zero_coverage_count <- dplyr::group_by(zero_coverage_events, contig_ID)
  zero_coverage_count <- dplyr::summarise(zero_coverage_count, number_of_zero_coverage_sites = n())
  
  return(zero_coverage_count)
  # Will return an empty tibble if there are no zero-coverage events
  
}


# Function: wraps the above two functions to produce a single output data frame. Returns a data frame.
calculate_per_contig_stats <- function(coverage_table, contig_lengths) {
  
  depth_per_contig <- calculate_depth_per_contig(coverage_table, contig_lengths)
  zero_coverage_count <- count_zero_coverage_events(coverage_table)
  
  per_contig_stats <- dplyr::left_join(depth_per_contig, zero_coverage_count, by = "contig_ID")
  
  return(per_contig_stats)
  
}


# Function: determines overall stats for the input sample for coverage and 0-coverage events. Returns a list of three length 1 numeric vectors.
calculate_overall_stats <- function(coverage_table, per_contig_stats, zero_coverage_threshold) {
  
  # Calculate mean and std. dev. of coverage across all entries
  coverage_mean <- mean(coverage_table$coverage)
  coverage_sd <- sd(coverage_table$coverage)
  
  ### Calculate percent of contigs with a zero-coverage event
  # First, filiter out contigs with a number of zero-coverage sites below some arbitrary threshold (e.g., if 1, means that contigs with even 1 zero-coverage site will be counted)
  per_contig_stats <- dplyr::filter(per_contig_stats, number_of_zero_coverage_sites >= zero_coverage_threshold)
  # Then calculate how many contigs have at least 'threshold' zero-coverage sites
  number_of_zero_coverage_contigs <- nrow(per_contig_stats)
  total_contigs <- length(unique(coverage_table$contig_ID))
  percent_contigs_with_zero_coverage_event <- number_of_zero_coverage_contigs / total_contigs * 100
  
  coverage_list <- list(coverage_mean, coverage_sd, percent_contigs_with_zero_coverage_event)
  names(coverage_list) <- c("coverage_mean", "coverage_sd", "percent_contigs_with_zero_coverage_event")
  
  return(coverage_list)
  
}


# Function: binds unmasked and masked per-contig stats tables into a single table. Returns data frame.
bind_supplemental_table <- function(bin_ID, metagenome_ID, per_contig_stats, per_contig_stats_filtered) {
  
  # Join table info together
  colnames(per_contig_stats_filtered) <- c("contig_ID", "contig_length", "coverage_mean_filtered", "coverage_mean_sd", 
                                           "number_of_zero_coverage_sites_filtered")
  contig_stats_joined <- dplyr::left_join(per_contig_stats, per_contig_stats_filtered, by = c("contig_ID", "contig_length"))
  
  # Add bin and metagenome ID and move to front
  contig_stats_joined$bin_ID <- bin_ID
  contig_stats_joined$metagenome_ID <- metagenome_ID
  contig_stats_joined <- dplyr::select(contig_stats_joined, bin_ID, metagenome_ID, everything())
  
  return(contig_stats_joined)
  
}


# Function: prints tab-separated table row with whole-sample coverage stats (unmasked and masked) to STDOUT. Two lines: header, content.
summarize_stats_for_printing <- function(bin_ID, metagenome_ID, overall_stats, overall_stats_filtered) {
  
  headers <- (paste("bin_ID", "metagenome_ID", 
                    "coverage_mean", "coverage_sd", 
                    "percent_contigs_with_zero_coverage_event",
                    "coverage_mean_filtered", "coverage_sd_filtered", 
                    "percent_contigs_with_zero_coverage_event_filtered", sep = "\t"))
  
  content <- paste(bin_ID, metagenome_ID, 
                         overall_stats$coverage_mean, overall_stats$coverage_sd, 
                         overall_stats$percent_contigs_with_zero_coverage_event,
                         overall_stats_filtered$coverage_mean, overall_stats_filtered$coverage_sd, 
                         overall_stats_filtered$percent_contigs_with_zero_coverage_event, sep = "\t")
  
  write(headers, file = "")
  write(content, file = "")
  
}


main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  message(ts(), "Running calculate_coverage_stats.R")
  message(ts(), "Input samtools depth coverage table: ", input_coverage_table_filename)
  message(ts(), "Output supplementary table (main results to STDOUT): ", output_contig_stats_table_filename)
  message(ts(), "Expected read length: ", read_length)
  message(ts(), "Threshold for searching for zero-coverage sites on a contig: ", zero_coverage_threshold)
  message(ts(), "bin_ID and metagenome_ID explicitly provided: ", !(parse_from_input_filename))
  if (parse_from_input_filename == FALSE) {
    
    message(ts(), "bin_ID: ", bin_ID)
    message(ts(), "metagenome_ID: ", metagenome_ID)
    
  } else {
    
    # Parse out sample names
    message(ts(), "Using name_delimiter to extract bin_ID and metagenome_ID from input filename: ", name_delimiter)
    message(ts(), "NOTE: Assuming filename structure is 'metagenome_ID[delimiter]bin_ID.extension'")
    
    # Remove path and extension
    filename_base <- file_path_sans_ext(basename(input_coverage_table_filename))
    
    # Split by the user-provided delimiter
    filename_split <- strsplit(filename_base, split = name_delimiter)[[1]]
    
    # Quick check to make sure the delimiter worked
    if (length(filename_split) != 2) {
      stop(paste(ts(), "After parsing filename based on delimiter '", name_delimiter, "', there were ", 
                 length(filename_split), " entries (expecting 2). Exiting...", sep = ""))
    }
    
    metagenome_ID <- filename_split[1]
    bin_ID <- filename_split[2]
    
    message(ts(), "bin_ID: ", bin_ID)
    message(ts(), "metagenome_ID: ", metagenome_ID)
  }
  
  # Read table and get lengths
  message(ts(), "Reading input table")
  coverage_table <- read_coverage_table(input_coverage_table_filename)
  message(ts(), "Calculating contig lengths")
  contig_lengths <- calculate_contig_lengths(coverage_table)
  
  # Calculate stats on unmasked data (keeping contig ends)
  message(ts(), "Calculating stats (non-masked)")
  per_contig_stats <- calculate_per_contig_stats(coverage_table, contig_lengths)
  overall_stats <- calculate_overall_stats(coverage_table, per_contig_stats, 
                                           zero_coverage_threshold = zero_coverage_threshold)
  
  # Calculate stats on masked data (removing contig ends)
  message(ts(), "Calculating stats (masked based on read_length = ", read_length, ")")
  coverage_table_filtered <- filter_contig_ends(coverage_table, contig_lengths, read_length)
  per_contig_stats_filtered <- calculate_per_contig_stats(coverage_table_filtered, contig_lengths)
  overall_stats_filtered <- calculate_overall_stats(coverage_table_filtered, per_contig_stats_filtered, 
                                                    zero_coverage_threshold = zero_coverage_threshold)
  
  ### Write supplemental table with per-contig stats
  message(ts(), "Binding supplementary table of per-contig coverage statsand writing to file")
  supplemental_table <- bind_supplemental_table(bin_ID, metagenome_ID, per_contig_stats, per_contig_stats_filtered)
  write.table(x = supplemental_table, file = output_contig_stats_table_filename, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # Print whole-sample stats to STDOUT
  message(ts(), "Summarizing overall stats (writing to STDOUT)")
  summarize_stats_for_printing(bin_ID, metagenome_ID, overall_stats, overall_stats_filtered)
  
  message(ts(), "Done.")
  
}

main()

