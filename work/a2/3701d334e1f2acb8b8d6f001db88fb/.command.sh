#!/usr/bin/env Rscript

library(dplyr)

read_and_standardize_csv <- function(file_path) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    standard_columns <- c("sequence", "start", "end", "score", "type", "method", "source", "hla", "ic50")
    result_df <- data.frame(
        sequence = df$sequence %||% NA_character_,
        start = df$start %||% NA_integer_,
        end = df$end %||% NA_integer_,
        score = df$score %||% NA_real_,
        type = df$type %||% NA_character_,
        method = df$method %||% NA_character_,
        source = df$source %||% NA_character_,
        hla = df$hla %||% NA_character_,
        ic50 = df$ic50 %||% NA_real_,
        stringsAsFactors = FALSE
    )
    return(result_df)
}

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || any(is.na(x))) y else x
}

bcell_df <- read_and_standardize_csv("BAL61222.1_sequence_bcell_epitopes.csv")
tcell_i_df <- read_and_standardize_csv("BAL61222.1_sequence_tcell_i_epitopes.csv")
tcell_ii_df <- read_and_standardize_csv("BAL61222.1_sequence_tcell_ii_epitopes.csv")

combined_df <- rbind(bcell_df, tcell_i_df, tcell_ii_df)

combined_df <- combined_df %>%
    mutate(
        consensus_score = case_when(
            type == "B-cell" ~ score,
            !is.na(ic50) ~ 1 - (pmin(ic50, 5000) / 5000),
            TRUE ~ score
        )
    ) %>%
    arrange(desc(consensus_score))

write.csv(combined_df, "combined_epitopes.csv", row.names = FALSE)

cat("Epitopes combined from multiple prediction methods.\n")
cat("Total epitopes:", nrow(combined_df), "\n")
cat("Epitope types:\n")
print(table(combined_df$type))
