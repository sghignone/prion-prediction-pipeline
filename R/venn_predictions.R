#!/usr/bin/env Rscript

################################################################################
# Venn Diagram: PrionScan vs PAPA Predictions
################################################################################
# Compares prediction outputs from PrionScan and PAPA to visualize
# shared and unique protein predictions.
#
# Usage:
#   Rscript R/venn_predictions.R --prionscan <file> --papa <file>
#   Rscript R/venn_predictions.R --db sprot  # Auto-detect latest files
#
# Output:
#   - Venn diagram PNG in reports/
#   - Summary statistics to console
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggVennDiagram)
  library(optparse)
})

################################################################################
# CONFIGURATION
################################################################################

DATA_PROCESSED_DIR <- "data/processed"
REPORTS_DIR <- "reports"
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

################################################################################
# COMMAND-LINE ARGUMENTS
################################################################################

option_list <- list(
  make_option(c("--prionscan", "-p"), type = "character", default = NULL,
              help = "Path to PrionScan predictions TSV"),
  make_option(c("--papa", "-a"), type = "character", default = NULL,
              help = "Path to PAPA predictions TSV"),
  make_option(c("--db", "-d"), type = "character", default = NULL,
              help = "Database name (sprot/trembl) to auto-detect latest files"),
  make_option(c("--output", "-o"), type = "character", default = NULL,
              help = "Output PNG filename (default: auto-generated)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

################################################################################
# HELPER FUNCTIONS
################################################################################

#' Find most recent file matching pattern
find_latest_file <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)

  # Sort by modification time, return most recent
  file_info <- file.info(files)
  files[order(file_info$mtime, decreasing = TRUE)][1]
}

################################################################################
# MAIN
################################################################################

main <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("VENN DIAGRAM: PrionScan vs PAPA Predictions\n")
  cat("================================================================================\n")

  # Resolve input files
  prionscan_file <- opt$prionscan
  papa_file <- opt$papa

  # Auto-detect if --db provided

if (!is.null(opt$db)) {
    db <- tolower(opt$db)
    cat(sprintf("Auto-detecting latest files for database: %s\n", db))

    prionscan_file <- find_latest_file(
      DATA_PROCESSED_DIR,
      paste0("^", db, "_prion_predictions_.*\\.tsv$")
    )
    papa_file <- find_latest_file(
      DATA_PROCESSED_DIR,
      paste0("^", db, "_papa_predictions_.*\\.tsv$")
    )
  }

  # Validate inputs
  if (is.null(prionscan_file) || !file.exists(prionscan_file)) {
    stop("PrionScan file not found. Use --prionscan <file> or --db <sprot|trembl>")
  }
  if (is.null(papa_file) || !file.exists(papa_file)) {
    stop("PAPA file not found. Use --papa <file> or --db <sprot|trembl>")
  }

  cat(sprintf("PrionScan file: %s\n", basename(prionscan_file)))
  cat(sprintf("PAPA file:      %s\n", basename(papa_file)))

  # Read prediction files
  cat("\nReading prediction files...\n")

  prionscan_data <- read_tsv(prionscan_file, show_col_types = FALSE)
  papa_data <- read_tsv(papa_file, show_col_types = FALSE)

  # Extract positive predictions
  # PrionScan: all rows in output are positive (Score >= 50)
  prionscan_ids <- unique(prionscan_data$Protein_ID)

  # PAPA: filter by Above_Threshold == TRUE
  papa_positive <- papa_data %>% filter(Above_Threshold == TRUE)
  papa_ids <- unique(papa_positive$Protein_ID)

  cat(sprintf(" ↳ PrionScan positive predictions: %d proteins\n", length(prionscan_ids)))
  cat(sprintf(" ↳ PAPA positive predictions:      %d proteins\n", length(papa_ids)))

  # Compute set operations
  both <- intersect(prionscan_ids, papa_ids)
  prionscan_only <- setdiff(prionscan_ids, papa_ids)
  papa_only <- setdiff(papa_ids, prionscan_ids)
  all_unique <- union(prionscan_ids, papa_ids)

  cat("\n--- SET ANALYSIS ---\n")
  cat(sprintf("PrionScan only:  %d proteins\n", length(prionscan_only)))
  cat(sprintf("PAPA only:       %d proteins\n", length(papa_only)))
  cat(sprintf("Both methods:    %d proteins\n", length(both)))
  cat(sprintf("Total unique:    %d proteins\n", length(all_unique)))

  # Calculate percentages
  if (length(all_unique) > 0) {
    cat(sprintf("\nOverlap: %.1f%% of all predictions found by both methods\n",
                100 * length(both) / length(all_unique)))
  }

  # Create Venn diagram
  cat("\nGenerating Venn diagram...\n")

  venn_list <- list(
    PrionScan = prionscan_ids,
    PAPA = papa_ids
  )

  # Determine output filename
  output_file <- opt$output
  if (is.null(output_file)) {
    db_prefix <- if (!is.null(opt$db)) paste0(opt$db, "_") else ""
    output_file <- file.path(REPORTS_DIR,
                             paste0(db_prefix, "venn_predictions_", timestamp, ".png"))
  }

  # Ensure reports directory exists
  if (!dir.exists(REPORTS_DIR)) {
    dir.create(REPORTS_DIR, recursive = TRUE)
  }

  # Generate plot
  p <- ggVennDiagram(venn_list,
                     label_alpha = 0,
                     category.names = c("PrionScan", "PAPA")) +
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    scale_color_manual(values = c("PrionScan" = "#2C3E50", "PAPA" = "#E74C3C")) +
    labs(
      title = "Prion-Like Domain Predictions: PrionScan vs PAPA",
      subtitle = sprintf("Total: %d unique proteins | Overlap: %d (%.1f%%)",
                         length(all_unique), length(both),
                         ifelse(length(all_unique) > 0,
                                100 * length(both) / length(all_unique), 0)),
      caption = sprintf("Generated: %s", Sys.time())
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      plot.caption = element_text(size = 8, color = "gray50"),
      legend.position = "none"
    )

  ggsave(output_file, p, width = 8, height = 6, dpi = 300)

  cat(sprintf(" ↳ Saved to: %s\n", output_file))

  # Final summary
  cat("\n================================================================================\n")
  cat("SUMMARY\n")
  cat("================================================================================\n")
  cat(sprintf("PrionScan only:  %5d proteins (%.1f%%)\n",
              length(prionscan_only),
              ifelse(length(all_unique) > 0, 100 * length(prionscan_only) / length(all_unique), 0)))
  cat(sprintf("PAPA only:       %5d proteins (%.1f%%)\n",
              length(papa_only),
              ifelse(length(all_unique) > 0, 100 * length(papa_only) / length(all_unique), 0)))
  cat(sprintf("Both methods:    %5d proteins (%.1f%%)\n",
              length(both),
              ifelse(length(all_unique) > 0, 100 * length(both) / length(all_unique), 0)))
  cat(sprintf("─────────────────────────────────\n"))
  cat(sprintf("Total unique:    %5d proteins\n", length(all_unique)))
  cat("\n")
  cat(sprintf("Output: %s\n", output_file))
  cat("================================================================================\n")
}

################################################################################
# ENTRY POINT
################################################################################

if (!interactive()) {
  main()
}
