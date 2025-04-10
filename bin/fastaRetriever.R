#!/usr/bin/env Rscript

# Sequence retrieval script for in silico vaccine design pipeline
# This script retrieves protein sequences from NCBI using rentrez

# Check and install missing packages
packages_to_install <- c("rentrez", "Biostrings")
new_packages <- packages_to_install[!(packages_to_install %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("Installing missing packages:", paste(new_packages, collapse=", "), "\n")
  
  # Try CRAN first
  if("rentrez" %in% new_packages) {
    install.packages("rentrez", repos = "http://cran.us.r-project.org")
  }
  
  # Use BiocManager for Bioconductor packages
  if("Biostrings" %in% new_packages) {
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings")
  }
}

# Suppress startup messages
suppressPackageStartupMessages({
  library(rentrez)
  library(Biostrings)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: fastaRetriever.R <accession> [outfile]\n")
  cat("  <accession>: NCBI protein accession number\n")
  cat("  [outfile]: Optional output file name (default: accession_sequence.fasta)\n")
  quit(status = 1)
}

accession <- args[1]
outfile <- ifelse(length(args) >= 2, args[2], paste0(accession, "_sequence.fasta"))

# Set up error handling
options(warn = 2)

# Improved logging function
log_message <- function(msg, is_error = FALSE) {
  if (is_error) {
    cat(paste("ERROR:", msg, "\n"), file = stderr())
  } else {
    cat(paste(msg, "\n"))
  }
}

# Function to retrieve sequence by accession number
retrieve_specific_sequence <- function(accession_number, output_file = NULL) {
  log_message(paste("Retrieving sequence for accession number:", accession_number))
  
  if (is.null(output_file)) {
    output_file <- paste0(accession_number, "_sequence.fasta")
  }
  
  # Handle special case for testing
  if (accession_number == "test") {
    test_seq <- "MVSLVKSDQIGTSTLNQRMEKIVLLLAAKCQTPMGAICPYLGSPSFNPNQKIITI"
    writeLines(c(">Test_sequence", test_seq), output_file)
    log_message(paste("Test sequence created and saved to", output_file))
    return(output_file)
  }
  
  # Fetch the sequence using rentrez
  tryCatch({
    # First, check if the accession exists
    search_results <- entrez_search(db = "protein", term = accession_number)
    
    if (search_results$count == 0) {
      log_message(paste("No sequence found for accession:", accession_number), is_error = TRUE)
      quit(status = 1)
    }
    
    sequence_data <- entrez_fetch(
      db = "protein", 
      id = accession_number, 
      rettype = "fasta", 
      retmode = "text"
    )
    
    if (nchar(sequence_data) < 10) {
      log_message("Retrieved sequence is too short or empty. Check the accession number.", is_error = TRUE)
      quit(status = 1)
    }
    
    # Save the sequence to a FASTA file
    writeLines(sequence_data, output_file)
    log_message(paste("Sequence downloaded and saved to", output_file))
    return(output_file)
  }, error = function(e) {
    log_message(paste("Failed to retrieve sequence:", e$message), is_error = TRUE)
    quit(status = 1)
  })
}

# Execute the retrieval
retrieve_specific_sequence(accession, outfile)