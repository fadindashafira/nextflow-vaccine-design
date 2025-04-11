#!/usr/bin/env Rscript

# Install packages if not available
if (!requireNamespace("rentrez", quietly = TRUE)) {
    install.packages("rentrez", repos = "http://cran.us.r-project.org")
}

library(rentrez)

# Set up error handling
options(warn = 2)

# Retrieve sequence for given accession number
accession <- "BAL61222.1"

cat("Retrieving sequence for accession number:", accession, "\n")

# Fetch the sequence using rentrez
tryCatch({
    # First, check if the accession exists
    search_results <- entrez_search(db = "protein", term = accession)

    if (search_results$count == 0) {
        stop("No sequence found for accession: ", accession)
    }

    sequence_data <- entrez_fetch(
        db = "protein", 
        id = accession, 
        rettype = "fasta", 
        retmode = "text"
    )

    if (nchar(sequence_data) < 10) {
        stop("Retrieved sequence is too short or empty")
    }

    # Save the sequence to a FASTA file
    fasta_file <- paste0(accession, "_sequence.fasta")
    writeLines(sequence_data, fasta_file)

    cat("Sequence downloaded and saved to", fasta_file, "\n")
}, error = function(e) {
    cat("ERROR: Failed to retrieve sequence. ", e$message, "\n", file = stderr())
    quit(status = 1)
})
