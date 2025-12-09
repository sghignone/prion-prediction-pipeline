#!/usr/bin/env Rscript

################################################################################
# Prion Domain Prediction Script (v1.0.1)
################################################################################
# Predicts prion domains in protein sequences using sliding window analysis
#
# See docs/changelog_prion_parser.md for version history.
#
# This script:
# 1. Reads tabular cache from data_preanalysis.R (or raw .dat files)
# 2. Applies sliding window prion domain scoring algorithm
# 3. Outputs predictions in original and tabular formats
#
# Original algorithm: Vladimir Espinosa Angarica (2012)
# R translation and optimization: 2024-2025
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
  library(stringr)
  library(readr)
})

################################################################################
# CONFIGURATION
################################################################################

# Directory structure (matches data_preanalysis.R)
DATA_RAW_DIR <- "data/raw"
DATA_CACHE_DIR <- "data/cache"
DATA_PROCESSED_DIR <- "data/processed"
REPORTS_DIR <- "reports"

# Analysis parameters
WINDOW <- 60      # Size of sliding window
CUTOFF <- 50      # Minimum score threshold for predictions
CORES <- detectCores() - 1
if (CORES < 1) CORES <- 1

################################################################################
# CONSTANTS AND DATA STRUCTURES
################################################################################

# Amino acid scores for prion domain prediction
aa_scores <- list(
  PrD = c(
    A = -0.568, C = -3.807, D = -1.507, E = -2.766, F = -0.478,
    G = 0.040, H = -0.131, I = -1.515, K = -1.883, L = -1.556,
    M = 0.170, N = 2.511, P = 0.227, Q = 2.044, R = -1.196,
    S = 0.733, T = -0.268, V = -1.716, W = -3.459, Y = 0.786,
    O = 0, U = 0, B = 0, Z = 0, X = 0
  ),
  Core = c(
    A = -1.053, C = -3.807, D = -2.054, E = -3.766, F = -0.893,
    G = -0.219, H = -1.064, I = -1.737, K = -2.561, L = -2.352,
    M = 0.115, N = 2.959, P = -0.467, Q = 2.385, R = -1.322,
    S = 0.532, T = -0.921, V = -2.301, W = -3.459, Y = 1.025,
    O = 0, U = 0, B = 0, Z = 0, X = 0
  )
)

# Log-odds scores for proline pair distances (0-60 amino acids)
proline_pair_logodd <- c(
  -5.009, -5.002, -5.173, -5.366, -5.566, -5.735, -5.880, -6.015,
  -6.144, -6.206, -6.313, -6.419, -6.503, -6.576, -6.663, -6.757,
  -6.833, -6.898, -6.976, -7.066, -7.136, -7.242, -7.300, -7.379,
  -7.485, -7.570, -7.643, -7.676, -7.788, -7.862, -7.933, -8.019,
  -8.090, -8.137, -8.226, -8.330, -8.403, -8.482, -8.568, -8.637,
  -8.701, -8.778, -8.864, -8.932, -8.999, -9.086, -9.180, -9.238,
  -9.307, -9.350, -9.434, -9.500, -9.575, -9.637, -9.733, -9.784,
  -9.848, -9.919, -10.003, -10.062, -10.132
)
names(proline_pair_logodd) <- 0:60

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Clean organism name from SwissProt OS line format
os_name_cleaning <- function(osline) {
  # Remove square brackets around genus names
  cleaned <- gsub("\\[([^\\]]+)\\]", "\\1", osline)
  # Remove single quotes
  cleaned <- gsub("'", "", cleaned)
  # Remove trailing punctuation
  cleaned <- gsub("(,|, and|\\.)$", "", cleaned)

  if (grepl("[^\\(\\)]+\\(.+\\)[^\\(\\)]+$", cleaned)) {
    sci_name <- gsub("\\.$", "", cleaned)
  } else {
    match <- regmatches(cleaned, regexpr("\\S[^\\(]+", cleaned))
    sci_name <- ifelse(length(match) > 0, match[1], cleaned)
  }

  gsub("\\s+$", "", sci_name)
}

################################################################################
# CACHE DISCOVERY
################################################################################

#' Find most recent tabular cache file
#' @param db_type Database type: "sprot" or "trembl"
#' @return Path to cache file or NULL
find_tabular_cache <- function(db_type = NULL) {
  pattern <- if (!is.null(db_type)) {
    paste0("^", db_type, "_fungi_tabular_\\d{8}_\\d{6}\\.tsv$")
  } else {
    "^(sprot|trembl)_fungi_tabular_\\d{8}_\\d{6}\\.tsv$"
  }

  cache_files <- list.files(DATA_CACHE_DIR, pattern = pattern, full.names = TRUE)

  if (length(cache_files) == 0) {
    return(NULL)
  }

  # Get file info with modification times
  file_info <- file.info(cache_files)
  file_info$path <- cache_files

  # Return most recent
  most_recent <- file_info[which.max(file_info$mtime), "path"]
  return(most_recent)
}

################################################################################
# PRION DOMAIN PREDICTION ALGORITHM
################################################################################

#' Calculate prion domain score for a sequence window
calculate_window_score <- function(domain, aa_scores_prd, proline_logodd) {
  domain_chars <- strsplit(domain, "")[[1]]
  domain_score <- sum(sapply(domain_chars, function(aa) {
    ifelse(aa %in% names(aa_scores_prd), aa_scores_prd[aa], 0)
  }))

  proline_positions <- gregexpr("P", domain)[[1]]

  correction_factor <- 0
  if (length(proline_positions) >= 2 && proline_positions[1] != -1) {
    for (i in 1:(length(proline_positions) - 1)) {
      dist <- proline_positions[i + 1] - proline_positions[i]
      if (dist <= 60) {
        correction_factor <- correction_factor +
          proline_logodd[as.character(dist)]
      }
    }
  }

  round(domain_score + correction_factor, 3)
}

#' Analyze a single protein for prion domains
#' @param protein_id Protein identifier
#' @param organism Organism name
#' @param taxonomy Full taxonomy string
#' @param taxid NCBI Taxonomy ID
#' @param sequence Amino acid sequence
#' @param window Sliding window size
#' @param cutoff Score cutoff
#' @return Prediction result or NULL
analyze_protein <- function(protein_id, organism, taxonomy, taxid, sequence,
                            window, cutoff, aa_scores_prd, proline_logodd) {
  # Skip sequences shorter than window
  if (nchar(sequence) < window) return(NULL)

  # Clean organism name
  org_clean <- os_name_cleaning(organism)

  best_result <- list(Score = -Inf, Seq = "", Window = 0)

  # Slide window across sequence
  for (i in 0:(nchar(sequence) - window)) {
    domain <- substr(sequence, i + 1, i + window)
    score <- calculate_window_score(domain, aa_scores_prd, proline_logodd)

    if (score > best_result$Score) {
      best_result <- list(Score = score, Seq = domain, Window = i)
    }
  }

  # Only return results that meet cutoff
  if (best_result$Score >= cutoff) {
    return(list(
      seq_id = protein_id,
      organism = org_clean,
      taxonomy = taxonomy,
      taxid = taxid,
      score = best_result$Score,
      window = best_result$Window,
      domain = best_result$Seq
    ))
  }

  NULL
}

################################################################################
# PARALLEL PROCESSING FROM TABULAR CACHE
################################################################################

#' Process tabular cache with parallel workers
#' @param cache_file Path to tabular TSV file
#' @param num_cores Number of parallel workers
#' @return List of predictions and organism counts
process_tabular_cache <- function(cache_file, num_cores) {
  cat("\n=== Loading Tabular Cache ===\n")
  cat(sprintf("File: %s\n", basename(cache_file)))

  # Read cache file
  start_time <- Sys.time()
  cache_data <- read_tsv(
    cache_file,
    col_types = cols(
      Protein_ID = col_character(),
      Organism = col_character(),
      Full_Taxonomy = col_character(),
      NCBI_TaxID = col_character(),
      Sequence = col_character()
    ),
    show_col_types = FALSE,
    progress = TRUE
  )

  load_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf(" ↳ Loaded %s proteins in %.1f sec\n",
              format(nrow(cache_data), big.mark = ","),
              as.numeric(load_time)))

  # Filter out entries without sequences
  cache_data <- cache_data[!is.na(cache_data$Sequence) & nchar(cache_data$Sequence) >= WINDOW, ]
  cat(sprintf(" ↳ %s proteins with sequences >= %d aa\n",
              format(nrow(cache_data), big.mark = ","), WINDOW))

  cat("\n=== Processing in Parallel ===\n")
  cat(sprintf("Using %d CPU cores\n", num_cores))
  cat(sprintf("Window size: %d amino acids\n", WINDOW))
  cat(sprintf("Score cutoff: %d\n\n", CUTOFF))

  # Split data for parallel processing
  n_rows <- nrow(cache_data)
  chunk_size <- ceiling(n_rows / num_cores)

  process_start <- Sys.time()

  # Parallel processing
  results_list <- mclapply(1:num_cores, function(worker_id) {
    start_idx <- (worker_id - 1) * chunk_size + 1
    end_idx <- min(worker_id * chunk_size, n_rows)

    if (start_idx > n_rows) return(list(predictions = list(), counts = list()))

    local_predictions <- list()
    local_counts <- list()

    for (i in start_idx:end_idx) {
      row <- cache_data[i, ]

      # Count organism
      org_name <- os_name_cleaning(row$Organism)
      if (is.null(local_counts[[org_name]])) {
        local_counts[[org_name]] <- 0
      }
      local_counts[[org_name]] <- local_counts[[org_name]] + 1

      # Analyze for prion domain
      pred <- analyze_protein(
        protein_id = row$Protein_ID,
        organism = row$Organism,
        taxonomy = row$Full_Taxonomy,
        taxid = row$NCBI_TaxID,
        sequence = row$Sequence,
        window = WINDOW,
        cutoff = CUTOFF,
        aa_scores_prd = aa_scores$PrD,
        proline_logodd = proline_pair_logodd
      )

      if (!is.null(pred)) {
        if (is.null(local_predictions[[pred$organism]])) {
          local_predictions[[pred$organism]] <- list()
        }
        local_predictions[[pred$organism]][[length(local_predictions[[pred$organism]]) + 1]] <- pred
      }
    }

    list(predictions = local_predictions, counts = local_counts)
  }, mc.cores = num_cores)

  process_time <- difftime(Sys.time(), process_start, units = "secs")

  # Merge results from all workers
  cat("Merging results from workers...\n")
  all_predictions <- list()
  all_counts <- list()

  for (result in results_list) {
    # Merge predictions
    for (organism in names(result$predictions)) {
      if (is.null(all_predictions[[organism]])) {
        all_predictions[[organism]] <- list()
      }
      all_predictions[[organism]] <- c(all_predictions[[organism]],
                                       result$predictions[[organism]])
    }

    # Merge counts
    for (organism in names(result$counts)) {
      if (is.null(all_counts[[organism]])) {
        all_counts[[organism]] <- 0
      }
      all_counts[[organism]] <- all_counts[[organism]] + result$counts[[organism]]
    }
  }

  cat(sprintf(" ↳ Processing complete in %.1f sec\n", as.numeric(process_time)))
  cat(sprintf(" ↳ Found predictions in %d organisms\n", length(all_predictions)))

  list(predictions = all_predictions, counts = all_counts)
}

################################################################################
# STREAMING FILE PARSER (FALLBACK FOR RAW .dat FILES)
################################################################################

#' Stream SwissProt entries from file (fallback when no cache available)
#' Retained for compatibility with direct .dat file processing
stream_swissprot_file <- function(dat_file, num_cores) {
  cat("\n=== Streaming from Raw .dat File ===\n")
  cat("Note: Consider running data_preanalysis.R first for faster processing\n\n")

  # Validate file exists
  if (!file.exists(dat_file)) {
    stop("Not a valid PATH to the Sequence file!\n")
  }

  # Count entries first
  cat("Counting entries...\n")
  if (grepl("\\.gz$", dat_file, ignore.case = TRUE)) {
    con <- gzfile(dat_file, "rt")
  } else {
    con <- file(dat_file, "rt")
  }

  total_entries <- 0
  while (length(line <- readLines(con, n = 1000, warn = FALSE)) > 0) {
    total_entries <- total_entries + sum(grepl("^//", line))
  }
  close(con)

  cat(sprintf("Found %s entries\n", format(total_entries, big.mark = ",")))

  # Process with parallel workers (simplified - each reads full file)
  entries_per_worker <- ceiling(total_entries / num_cores)

  cat(sprintf("Processing with %d workers (~%s entries each)...\n",
              num_cores, format(entries_per_worker, big.mark = ",")))

  results_list <- mclapply(1:num_cores, function(worker_id) {
    start_entry <- (worker_id - 1) * entries_per_worker + 1
    end_entry <- min(worker_id * entries_per_worker, total_entries)

    local_predictions <- list()
    local_counts <- list()

    if (grepl("\\.gz$", dat_file, ignore.case = TRUE)) {
      con <- gzfile(dat_file, "rt")
    } else {
      con <- file(dat_file, "rt")
    }

    entry_num <- 0
    current_id <- ""
    current_os <- ""
    current_oc <- ""
    current_ox <- ""
    current_seq <- ""
    in_sequence <- FALSE

    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      if (startsWith(line, "ID ")) {
        current_id <- sub("^ID\\s+(\\w+)\\s+.*", "\\1", line)
      } else if (startsWith(line, "OS ")) {
        os_part <- sub("^OS\\s+", "", line)
        current_os <- if (current_os == "") os_part else paste(current_os, os_part)
      } else if (startsWith(line, "OC ")) {
        oc_part <- sub("^OC\\s+", "", line)
        oc_part <- sub(";$", "", oc_part)
        current_oc <- if (current_oc == "") oc_part else paste(current_oc, oc_part, sep = "; ")
      } else if (startsWith(line, "OX ")) {
        ox_match <- regmatches(line, regexpr("NCBI_TaxID=(\\d+)", line))
        if (length(ox_match) > 0) {
          current_ox <- sub("NCBI_TaxID=", "", ox_match)
        }
      } else if (startsWith(line, "SQ ")) {
        in_sequence <- TRUE
      } else if (in_sequence && !startsWith(line, "//") && grepl("[A-Z]", line)) {
        current_seq <- paste0(current_seq, gsub("[^A-Z]", "", line))
      } else if (startsWith(line, "//")) {
        entry_num <- entry_num + 1

        if (entry_num >= start_entry && entry_num <= end_entry) {
          org_name <- os_name_cleaning(current_os)

          if (is.null(local_counts[[org_name]])) {
            local_counts[[org_name]] <- 0
          }
          local_counts[[org_name]] <- local_counts[[org_name]] + 1

          if (nchar(current_seq) >= WINDOW) {
            pred <- analyze_protein(
              protein_id = current_id,
              organism = current_os,
              taxonomy = current_oc,
              taxid = current_ox,
              sequence = current_seq,
              window = WINDOW,
              cutoff = CUTOFF,
              aa_scores_prd = aa_scores$PrD,
              proline_logodd = proline_pair_logodd
            )

            if (!is.null(pred)) {
              if (is.null(local_predictions[[pred$organism]])) {
                local_predictions[[pred$organism]] <- list()
              }
              local_predictions[[pred$organism]][[length(local_predictions[[pred$organism]]) + 1]] <- pred
            }
          }
        }

        if (entry_num > end_entry) break

        current_id <- ""
        current_os <- ""
        current_oc <- ""
        current_ox <- ""
        current_seq <- ""
        in_sequence <- FALSE
      }
    }

    close(con)
    list(predictions = local_predictions, counts = local_counts)
  }, mc.cores = num_cores)

  # Merge results
  all_predictions <- list()
  all_counts <- list()

  for (result in results_list) {
    for (organism in names(result$predictions)) {
      if (is.null(all_predictions[[organism]])) {
        all_predictions[[organism]] <- list()
      }
      all_predictions[[organism]] <- c(all_predictions[[organism]],
                                       result$predictions[[organism]])
    }
    for (organism in names(result$counts)) {
      if (is.null(all_counts[[organism]])) {
        all_counts[[organism]] <- 0
      }
      all_counts[[organism]] <- all_counts[[organism]] + result$counts[[organism]]
    }
  }

  list(predictions = all_predictions, counts = all_counts)
}

################################################################################
# OUTPUT WRITING
################################################################################

#' Write prediction outputs
write_outputs <- function(results, output_file, tabular_file, db_name) {
  all_predictions <- results$predictions

  cat("\n=== Writing Output Files ===\n")

  # OUTPUT 1: Original format
  cat(sprintf("Writing: %s\n", output_file))
  sink(output_file)

  for (organism in sort(names(all_predictions))) {
    predictions <- all_predictions[[organism]]
    cat(sprintf(">%s: Total=%d\n", organism, length(predictions)))

    for (pred in predictions) {
      cat(sprintf("\t%s\tWindow Position=%d; Score=%.3f | Prion Domain: %s\n",
                  pred$seq_id, pred$window, pred$score, pred$domain))
    }
    cat("\n")
  }

  sink()

  # OUTPUT 2: Tabular format
  cat(sprintf("Writing: %s\n", tabular_file))

  tabular_rows <- list()
  idx <- 1

  for (organism in names(all_predictions)) {
    for (pred in all_predictions[[organism]]) {
      tabular_rows[[idx]] <- data.frame(
        Protein_ID = pred$seq_id,
        Organism = pred$organism,
        Taxonomy = pred$taxonomy,
        NCBI_TaxID = ifelse(is.null(pred$taxid) || is.na(pred$taxid), "", pred$taxid),
        Window_Position = pred$window,
        Score = pred$score,
        Prion_Domain = pred$domain,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }

  if (length(tabular_rows) > 0) {
    tabular_data <- do.call(rbind, tabular_rows)
    tabular_data <- tabular_data[order(tabular_data$Organism, -tabular_data$Score), ]
    write_tsv(tabular_data, tabular_file)
  } else {
    # Empty file with headers
    write_tsv(
      data.frame(
        Protein_ID = character(),
        Organism = character(),
        Taxonomy = character(),
        NCBI_TaxID = character(),
        Window_Position = integer(),
        Score = numeric(),
        Prion_Domain = character()
      ),
      tabular_file
    )
  }

  cat(" ↳ Done\n")
}

################################################################################
# MAIN EXECUTION
################################################################################

main <- function() {
  option_list <- list(
    make_option(
      c("-s", "--seq-file"),
      type = "character",
      default = NULL,
      help = "PATH to SwissProt format file (.dat or .dat.gz). Optional if cache exists.",
      metavar = "FILE"
    ),
    make_option(
      c("-d", "--db"),
      type = "character",
      default = NULL,
      help = "Database type: 'sprot' or 'trembl'. Auto-detects from cache if not specified.",
      metavar = "TYPE"
    ),
    make_option(
      c("-o", "--output-file"),
      type = "character",
      default = NULL,
      help = "Output file for original format predictions (default: auto-generated)",
      metavar = "FILE"
    ),
    make_option(
      c("-t", "--tabular-output"),
      type = "character",
      default = NULL,
      help = "Output file for tab-separated format (default: auto-generated)",
      metavar = "FILE"
    ),
    make_option(
      c("-c", "--cores"),
      type = "integer",
      default = CORES,
      help = paste0("Number of CPU cores to use (default: ", CORES, ")"),
      metavar = "NUMBER"
    )
  )

  opt_parser <- OptionParser(
    option_list = option_list,
    description = paste("\n=== Prion Domain Prediction v1.0 ===\n",
                        "Predicts prion domains using sliding window analysis.\n",
                        "Uses cached tabular files from data_preanalysis.R for speed."),
    epilogue = paste("\nReferences:",
                     "  Angarica et al. (2013) BMC Genomics",
                     "  Alberti et al. (2009) Cell 137(1):146-158",
                     "\nExamples:",
                     "  Rscript R/prion_parser.R --db sprot",
                     "  Rscript R/prion_parser.R -s data/raw/uniprot_trembl_fungi.dat.gz",
                     "  Rscript R/prion_parser.R --db sprot -o results.dat -t results.tsv",
                     sep = "\n")
  )

  opt <- parse_args(opt_parser)
  num_cores <- opt$cores

  cat("\n")
  cat("================================================================================\n")
  cat("PRION DOMAIN PREDICTION PIPELINE (v1.0)\n")
  cat("================================================================================\n")
  cat("See docs/changelog_prion_parser.md for version history.\n")
  cat("================================================================================\n")

  # Determine input source
  cache_file <- NULL
  dat_file <- NULL
  db_name <- opt$db

  # Try to find cache first
  if (!is.null(db_name)) {
    cache_file <- find_tabular_cache(db_name)
  } else {
    cache_file <- find_tabular_cache(NULL)
  }

  if (!is.null(cache_file)) {
    cat(sprintf("\nFound tabular cache: %s\n", basename(cache_file)))
    # Extract db_name from cache filename
    if (is.null(db_name)) {
      if (grepl("sprot", basename(cache_file))) db_name <- "sprot"
      else if (grepl("trembl", basename(cache_file))) db_name <- "trembl"
      else db_name <- "uniprot"
    }
  } else if (!is.null(opt$`seq-file`)) {
    dat_file <- opt$`seq-file`
    cat(sprintf("\nNo cache found, using raw file: %s\n", basename(dat_file)))
    if (is.null(db_name)) {
      if (grepl("sprot", basename(dat_file), ignore.case = TRUE)) db_name <- "sprot"
      else if (grepl("trembl", basename(dat_file), ignore.case = TRUE)) db_name <- "trembl"
      else db_name <- "uniprot"
    }
  } else {
    print_help(opt_parser)
    stop("\nError: No cache found and no --seq-file specified.\n",
         "Run data_preanalysis.R first or provide a .dat file.\n",
         call. = FALSE)
  }

  # Generate output filenames
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  if (is.null(opt$`output-file`)) {
    opt$`output-file` <- file.path(DATA_PROCESSED_DIR,
                                   paste0(db_name, "_prion_predictions_", timestamp, ".dat"))
  }

  if (is.null(opt$`tabular-output`)) {
    opt$`tabular-output` <- file.path(DATA_PROCESSED_DIR,
                                      paste0(db_name, "_prion_predictions_", timestamp, ".tsv"))
  }

  cat(sprintf("Output: %s\n", opt$`output-file`))
  cat(sprintf("Tabular: %s\n", opt$`tabular-output`))

  # Process
  if (!is.null(cache_file)) {
    results <- process_tabular_cache(cache_file, num_cores)
  } else {
    results <- stream_swissprot_file(dat_file, num_cores)
  }

  # Write outputs
  write_outputs(results, opt$`output-file`, opt$`tabular-output`, db_name)

  # Summary
  total_predictions <- sum(sapply(results$predictions, length))

  cat("\n")
  cat("================================================================================\n")
  cat("ANALYSIS COMPLETE\n")
  cat("================================================================================\n")
  cat(sprintf("Total organisms with predictions: %d\n", length(results$predictions)))
  cat(sprintf("Total predictions: %s\n", format(total_predictions, big.mark = ",")))
  cat("\nOutput files:\n")
  cat(sprintf("  1. %s (original format)\n", opt$`output-file`))
  cat(sprintf("  2. %s (tabular format)\n", opt$`tabular-output`))
  cat("================================================================================\n")
}

################################################################################
# SCRIPT ENTRY POINT
################################################################################

if (!interactive()) {
  main()
}
