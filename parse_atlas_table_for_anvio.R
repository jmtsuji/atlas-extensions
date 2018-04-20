#!/usr/bin/env Rscript
# Extract information from the ATLAS table for use in Anvi'o
# Copyright Neufeld lab, 2018
# Description: Generates three summary files from the ATLAS table needed for import into Anvi'o.
# See required R packages below.

#####################################################
## User variables: #################################
RUN_COMMAND_LINE <- TRUE # If selected, all user input here is ignored, and terminal-based input is expected instead.

# Set other user variables here
if (RUN_COMMAND_LINE == FALSE) {
  setwd("/Users/JTsuji/Downloads") # your working directory where files are stored
  atlas_table_filename <- "CA-L227-2014_annotations_test_50.txt"
  output_tax_filename <- "CA-L227-2014_taxonomy.tsv"
  output_contig_bin_filename <- "CA-L227-2014_binning_results.tsv"
  output_bin_info_filename <- "CA-L227-2014_bins_info.tsv"
}
#####################################################

#####################################################
## Load required packages: ##########################
library(getopt)
library(plyr)
suppressMessages(library(dplyr))
#####################################################

SCRIPT_VERSION <- "1.0.22-coassembly-r4-dev"

parse_command_line_input <- function() {
  ### Grab arguments
  # Arguments required:
  # -a input ATLAS annotations TSV table
  # -f input featureCounts TSV table
  # -o output TSV filename
  params <- matrix(c('atlas_table', 'a', 1, "character",
                     'output_tax_file', 't', 1, "character",
                     'output_contig_bin_file', 'c', 1, "character",
                     'output_bin_info_file', 'b', 1, "character",
                     'help', 'h', 2, "character"), byrow=TRUE, ncol=4)
  
  opt <- getopt(params)
  
  # If help was called, print help message and exit
  if ( !is.null(opt$help) ) {
    
    cat("parse_atlas_table_for_anvio.R: Creates taxonomy and bin info for Anvi'o from the ATLAS annotations table.\n")
    cat("Copyright Neufeld Lab, 2018\n")
    cat("Contact Jackson M. Tsuji (jackson.@uwaterloo.ca) for bug reports or feature requests.\n")
    cat(paste("Version: ", SCRIPT_VERSION, "\n\n", sep = ""))
    
    cat(getopt(params, usage = TRUE))
    
    cat("\n")
    
    cat("Details:\n", "--atlas_table\t\t\tFilepath for TSV-format ATLAS annotations table. [Required]\n",
        "--output_tax_file\t\tFilepath for output TSV-format gene-specific taxonomy file. [Required]\n",
        "--output_contig_bin_file\t\tFilepath for output TSV-format summary file of contig assignmen to bins. [Required]\n",
        "--output_bin_info_file\t\t\tFilepath for output TSV-format summary file of bin IDs with CheckM taxonomy - colours can be modified by user. [Required]\n\n")
    
    quit(status = 1)
  }
  
  # Exit if required inputs are not provided
  if ( is.null(opt$atlas_table) ) {
    stop("Input ATLAS annotations table filepath required. Try -h for help message.")
  }
  if ( is.null(opt$output_tax_file) ) {
    stop("Output taxonomy filepath required. Try -h for help message.")
  }
  if ( is.null(opt$output_contig_bin_file) ) {
    stop("Output contig-bin filepath required. Try -h for help message.")
  }
  if ( is.null(opt$output_bin_info_file) ) {
    stop("Output bin info filepath required. Try -h for help message.")
  }
  
  # Make variables from provided input and save as global variables (<<-)
  atlas_table_filename <<- opt$atlas_table
  output_tax_filename <<- opt$output_tax_file
  output_contig_bin_filename <<- opt$output_contig_bin_file
  output_bin_info_filename <<- opt$output_bin_info_file
  
}

parse_taxonomy_preface <- function(taxonomy_rank_entry) {
  # e.g., taxonomy_rank_entry <- "k__Bacteria"
  # Want to convert to "Bacteria"
  split <- strsplit(taxonomy_rank_entry, split = "__")[[1]][2]
  
  return(split)
}

parse_greengenes_taxonomy <- function(greengenes_entry_vector, contig_id, locus_tag) { 
  # E.g., greengenes_entry_vector <- "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Methylococcales;f__Methylococcaceae;g__Methylobacter;s__?"
  # locus_tag <- "CA-L227-2014_00003"
  # Output is that single line formatted as a data frame matching Anvi'o's specifications.
  
  # Stop early if empty
  if (greengenes_entry_vector == "") {
    return(NULL)
  } else {
    
    # Parse by semicolon
    separated <- strsplit(greengenes_entry_vector, ";")[[1]]
    
    # Get rid of double underscore
    refined <- unlist(lapply(separated, parse_taxonomy_preface))
    
    df <- matrix(refined, ncol = 7)
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    colnames(df) <- c("t_kingdom", "t_phylum", "t_class", "t_order", "t_family", "t_genus", "t_species")
    
    # Add locus_tag and contig ID for santity checking (not yet added!)
    df$locus_tag <- locus_tag
    df$contig_id <- contig_id
    
    # Renumber numerically for gene_callers_id
    # Note: NOT exactly the same as the locus_tag ID!
    # TODO - add a check that the ATLAS table is sorted ahead of time.
    gene_callers_id <- seq(1:nrow(df))
    
    # # Parse locus_tag to get gene_callers_id
    # gene_callers_id <- strsplit(locus_tag, split = "_")[[1]]
    # gene_callers_id <- as.numeric(gene_callers_id[length(gene_callers_id)])
    
    # Add to table and re-arrange for clarity
    df$gene_callers_id <- gene_callers_id
    df <- df[,c(ncol(df), 1:(ncol(df) - 1))]
    
    # Remove "?" and replace with blank
    for (i in 1:ncol(df)) {
      df[1,i] <- gsub(pattern = "?", replacement = "", fixed = TRUE, x = df[1,i])
    }
    
    # Bring column number down to Anvi'o requirements
    df_trunc <- df
    df_trunc$t_kingdom <- NULL
    df_trunc$locus_tag <- NULL
    df_trunc$contig_id <- NULL
    
    # Return the full table for testing and the truncated table for anvio
    return(list(df_trunc, df))
  }
  
}

remove_special_characters <- function(input_vector) {
  removed <- gsub(pattern = "[[:punct:]]", replacement = "_", x = input_vector, fixed = FALSE)
  
  return(removed)
}

summarize_contig_bin_mapping <- function(atlas_table) {
  # Get required columns
  cols_to_grab <- c("contig_id", "bin_id")
  mapping <- atlas_table[,cols_to_grab]
  
  # Collapse to unique entries
  mapping <- mapping[!duplicated(mapping),]
  
  # # Remove special characters
  # mapping_cleaned <- lapply(1:ncol(mapping), function(x) { remove_special_characters(mapping[,x]) })
  # names(mapping_cleaned) <- cols_to_grab
  # mapping <- as.data.frame(dplyr::bind_cols(mapping_cleaned))
  
  mapping$bin_id <- remove_special_characters(mapping$bin_id)
  
  return(mapping)
  
  # Note: remember to export with no column names as required by anvi'o.
}

summarize_bin_info <- function(atlas_table) {
  # Get required data
  bin_info <- unique(atlas_table$bin_id)
  
  # Sort by name
  bin_info <- bin_info[order(bin_info)]
  
  # Remove special characters
  bin_info <- remove_special_characters(bin_info)
  
  # Make data frame template
  bin_df <- data.frame("bin_id" = bin_info,
                       "source" = "metabat2",
                       "colour" = "#000000")
  
  return(bin_df)
  
  # Note: remember to export with no column names as required by anvi'o.
}

main <- function() {
  # Run command line version if requested
  if (RUN_COMMAND_LINE == TRUE) {
    parse_command_line_input()
  }
  
  # Startup messages
  cat(paste("Running parse_atlas_table_for_anvio.R, version: ", SCRIPT_VERSION, "\n\n", sep = ""))
  cat(paste("ATLAS table filepath:", atlas_table_filename, "\n"))
  cat(paste("Output taxonomy filepath:", output_tax_filename, "\n"))
  cat(paste("Output contig/bin summary filepath:", output_contig_bin_filename, "\n"))
  cat(paste("Output bin info template filepath:", output_bin_info_filename, "\n"))
  cat("\n")
  
  # Read the table
  cat("Reading ATLAS table\n")
  atlas_table <- read.table(atlas_table_filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  
  # Make the taxonomy table
  cat("Parsing taxonomy\n")
  parsed_tax <- lapply(1:nrow(atlas_table), function(x) { parse_greengenes_taxonomy(greengenes_entry_vector = atlas_table[x,9], 
                                                                                    contig_id = atlas_table[x,1], locus_tag = atlas_table[x,2]) })
  parsed_tax <- dplyr::bind_rows(parsed_tax)
  
  # Make the contig-bin mapping file
  cat("Summarizing contig-to-bin mapping\n")
  contig_bin_mapping <- summarize_contig_bin_mapping(atlas_table)
  
  # Make the bin info file template
  cat("Generating bin info template file\n")
  bin_info <- summarize_bin_info(atlas_table)
  
  # Export
  cat("Writing output\n")
  write.table(parsed_tax, file = output_tax_filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write.table(contig_bin_mapping, file = output_contig_bin_filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(bin_info, file = output_bin_info_filename, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  cat("\nparse_atlas_table_for_anvio.R: finished.\n")
}

main()
