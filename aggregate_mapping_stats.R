#!/usr/bin/env Rscript
# aggregate_mapping_stats.R
# Copyright Jackson M. Tsuji, 2018
# Neufeld Research Group
# Aggregates the stats tables created by 'calculate_bin_abundance_in_metagenome.sh' into a final summary file
# See required R packages below.

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("/home/jmtsuji/Downloads") # your working directory where files are stored
  genome_stats_filename <- 
  read_count_summary_filename <- ""
  bin_mapping_summary_filename <- ""
  coverage_summary_filename <- ""
  output_final_stats_summary_filename <- ""
  
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
  params <- matrix(c('genome_stats_table', 'g', 1, "character",
                     'read_count_summary_table', 'r', 1, "character",
                     'bin_mapping_summary_table', 'm', 1, "character",
                     'coverage_summary_table', 'c', 1, "character",
                     'output_summary_table', 'o', 1, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("aggregate_mapping_stats.R: Internal script in 'calculate_bin_abundance_in_metagenome.sh' procedure for aggregating stats tables.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n\n")
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$genome_stats_table) ) {
    stop("Genome bin stats table required (-g). Try -h for help message.")
  }
  if ( is.null(opt$read_count_summary_table) ) {
    stop("Metagenome read count summary table required (-r). Try -h for help message.")
  }
  if ( is.null(opt$bin_mapping_summary_table) ) {
    stop("Bin mapping summary table required (-m). Try -h for help message.")
  }
  if ( is.null(opt$coverage_summary_table) ) {
    stop("Coverage summary table required (-c). Try -h for help message.")
  }
  if ( is.null(opt$output_summary_table) ) {
    stop("Output filename required (-o). Try -h for help message.")
  }
  
  # Make variables from provided input and save as global variables (<<-)
  genome_stats_filename <<- opt$genome_stats_table
  read_count_summary_filename <<- opt$read_count_summary_table
  bin_mapping_summary_filename <<- opt$bin_mapping_summary_table
  coverage_summary_filename <<- opt$coverage_summary_table
  output_final_stats_summary_filename <<- opt$output_summary_table
  
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


# Function: reads TSV table with desired input options, given the filename
read_tsv_table <- function(table_filename) {
  
  tsv_table <- read.table(table_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  return(tsv_table)
  
}


main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  message(ts(), "Running aggregate_mapping_stats.R")
  message(ts(), "Genome stats table: ", genome_stats_filename)
  message(ts(), "Metagenome read count table): ", read_count_summary_filename)
  message(ts(), "Bin mapping stats table: ", bin_mapping_summary_filename)
  message(ts(), "Coverage summary table: ", coverage_summary_filename)
  message(ts(), "Output summary stats table: ", output_final_stats_summary_filename)
  
  # Read each table and truncate to the desired columns
  message(ts(), "Reading and reformatting input tables")
  
  
  ##### 1. Genome stats
  genome_stats <- read_tsv_table(genome_stats_filename)
  ### Example
  # n_contigs	contig_bp	gap_pct	ctg_N50	ctg_L50	ctg_N90	ctg_L90	ctg_max	gc_avg	gc_std	filename
  # 208	614011	0.000	41	3843	155	1340	15350	0.55948	0.05012	/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/L227-2013-6m/genomic_bins/L227-2013-6m.001.fasta
  # 1056	4068125	0.000	208	5328	746	1715	66968	0.44385	0.03202	/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/L227-2013-6m/genomic_bins/L227-2013-6m.002.fasta
  # 1026	2235255	0.000	281	2367	820	1180	21504	0.67816	0.02542	/Analysis/jmtsuji/Hippodrome_legacy/180123_ELA111314_atlas_r1.0.22_plus_full/L227-2013-6m/genomic_bins/L227-2013-6m.003.fasta
  
  genome_stats <- dplyr::select(genome_stats, filename, contig_bp, n_contigs)
  genome_stats <- dplyr::rename(genome_stats, genome = filename, genome_num_contigs = n_contigs, genome_length_bp = contig_bp)
  genome_stats$genome <- tools::file_path_sans_ext(basename(genome_stats$genome))
  
  
  ##### 2. Metagenome stats
  metagenome_stats <- read_tsv_table(read_count_summary_filename)
  ### Example
  # metagenome	total_reads	R1_reads	R2_reads	se_reads
  # L227-2013-6m  50392710  23400501  23400501  3591708
  
  # Check that read stats are okay
  metagenome_stats$total_reads_check <- metagenome_stats$R1_reads + metagenome_stats$R2_reads + metagenome_stats$se_reads
  if ( identical(metagenome_stats$total_reads, metagenome_stats$total_reads_check) == FALSE ) {
    stop("Problem with read counts for the metagenome stats file provided. R1 + R2 + se does not always equal the total! Exiting...")
  }
  
  metagenome_stats <- dplyr::select(metagenome_stats, metagenome, total_reads)
  metagenome_stats <- dplyr::rename(metagenome_stats, metagenome_total_reads = total_reads)
  
  
  ##### 3. Bin mapping stats
  bin_mapping_stats <- read_tsv_table(bin_mapping_summary_filename)
  ### Example
  # genome	metagenome	mapped_reads
  # Bin001  L227-2013-6m  100572
  
  # No need to modify.
  
  ##### 4. Coverage stats
  coverage_stats <- read_tsv_table(coverage_summary_filename)
  ### Example
  # genome  metagenome_   coverage_mean   coverage_sd     percent_contigs_with_zero_coverage_event        coverage_mean_filtered  coverage_sd_filtered    percent_contigs_with_ze$
  # Bin001        L227-2013-6m    7.46117798098696        9.32951522614479        89.3238434163701        8.46412909077082        9.73329436732834        75.6227758007117
  # Bin001        L227-2013-6m    3.57531510053049        4.71256012846482        94.8473282442748        3.81342539634244        4.82750775745679        81.1068702290076
  # Bin001        L227-2013-6m    1.68890355501351        3.2299455150493 98.6149584487535        1.81322588275404        3.32302874589562        90.3047091412742
  
  coverage_stats <- dplyr::select(coverage_stats, genome, metagenome, coverage_mean, coverage_sd)
  
  
  ##### Combine tables
  message(ts(), "Combining tables")
  output_table <- dplyr::left_join(bin_mapping_stats, coverage_stats, by = c("genome", "metagenome"))
  output_table <- dplyr::left_join(output_table, genome_stats, by = "genome")
  output_table <- dplyr::left_join(output_table, metagenome_stats, by = "metagenome")
  
  
  # Write output table
  message(ts(), "Writing output table")
  write.table(x = output_table, file = output_final_stats_summary_filename, sep = "\t",
              row.names = FALSE, col.names = TRUE)
  
  message(ts(), "Done.")
  
}

main()

