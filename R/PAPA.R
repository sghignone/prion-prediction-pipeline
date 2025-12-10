#!/usr/bin/env Rscript

################################################################################
# PAPA - Prion Aggregation Prediction Algorithm (v1.0)
################################################################################
# Predicts prion-forming propensity using sliding windows and FoldIndex filtering
#
# See docs/changelog_PAPA.md for version history.
#
# This script:
# 1. Reads tabular cache from data_preanalysis.R (or raw .dat files)
# 2. Applies PAPA sliding window algorithm with FoldIndex filtering
# 3. Outputs predictions in tabular format
#
# Original Python implementation: http://www.prion.bcm.tmc.edu/
# R translation and pipeline integration: 2025
#
# Citation:
#   Toombs JA, McCarty BR, Ross ED. (2010)
#   "Compositional determinants of prion formation in yeast."
#   Mol Cell Biol 30(1):319-332
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
  library(readr)
  library(stringr)
})

################################################################################
# CONFIGURATION
################################################################################

# Directory structure (matches data_preanalysis.R and PrionScan.R)
DATA_RAW_DIR <- "data/raw"
DATA_CACHE_DIR <- "data/cache"
DATA_PROCESSED_DIR <- "data/processed"
REPORTS_DIR <- "reports"

# PAPA algorithm parameters
WINDOW_SIZE <- 41       # Default super-window size
HALF_WINDOW <- 20       # Half-window for local averaging
SCORE_THRESHOLD <- 0.05 # Default threshold for positive prediction

# Parallel processing
CORES <- detectCores() - 1
if (CORES < 1) CORES <- 1

################################################################################
# CONSTANTS AND DATA STRUCTURES
################################################################################

# Valid amino acids
AMINO_ACIDS <- c('A','C','D','E','F','G','H','I','K','L',
                 'M','N','P','Q','R','S','T','V','W','Y')

# Amino acid prion propensities (from training data)
PROPENSITIES <- c(
  A = -0.396490246, C =  0.415164505, D = -1.276997939, E = -0.605023827,
  F =  0.838732498, G = -0.039220713, H = -0.278573356, I =  0.813697862,
  K = -1.576748587, L = -0.040005335, M =  0.673729095, N =  0.080295334,
  P = -1.197447496, Q =  0.069168387, R = -0.405858577, S =  0.133912418,
  T = -0.11457038,  V =  0.813697862, W =  0.666735081, Y =  0.77865336
)

# Amino acid hydrophobicity values (for FoldIndex calculation)
HYDROPHOBICITY <- c(
  A = 0.7,         C = 0.777777778, D = 0.111111111, E = 0.111111111,
  F = 0.811111111, G = 0.455555556, H = 0.144444444, I = 1,
  K = 0.066666667, L = 0.922222222, M = 0.711111111, N = 0.111111111,
  P = 0.322222222, Q = 0.111111111, R = 0,           S = 0.411111111,
  T = 0.422222222, V = 0.966666667, W = 0.4,         Y = 0.355555556
)

# Amino acid charge values (for FoldIndex calculation)
CHARGE <- c(
  A = 0, C = 0, D = -1, E = -1, F = 0,
  G = 0, H = 0, I = 0,  K = 1,  L = 0,
  M = 0, N = 0, P = 0,  Q = 0,  R = 1,
  S = 0, T = 0, V = 0,  W = 0,  Y = 0
)

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Find most recent tabular cache file
#' @param db_type Database type: "sprot" or "trembl" (NULL for any)
#' @return Path to most recent cache file, or NULL if none found
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

#' Clean organism name from SwissProt OS line format
os_name_cleaning <- function(osline) {
  cleaned <- gsub("\\[([^\\]]+)\\]", "\\1", osline)
  cleaned <- gsub("'", "", cleaned)
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
# PAPA ALGORITHM FUNCTIONS
################################################################################

#' Get window boundaries for a given position
get_window <- function(sequence, position, half_window_size) {
  seq_len <- length(sequence)
  start <- max(position - half_window_size + 1, 1)
  end <- min(seq_len, position + half_window_size)
  return(list(start = start, end = end))
}

#' Calculate score for a single window position
window_score <- function(sequence, position, aa_dict, half_window_size,
                        ignore_consecutive_prolines = FALSE) {
  window <- get_window(sequence, position, half_window_size)
  start <- window$start
  end <- window$end

  score <- 0.0
  count <- 0

  for (i in start:end) {
    aa <- sequence[i]
    if (!(aa %in% AMINO_ACIDS)) next

    if (!ignore_consecutive_prolines) {
      score <- score + aa_dict[aa]
      count <- count + 1
    } else {
      if (aa != 'P') {
        score <- score + aa_dict[aa]
        count <- count + 1
      } else if ((i > 1 && sequence[i-1] == 'P') ||
                 (i > 2 && sequence[i-2] == 'P')) {
        # Skip consecutive prolines
      } else {
        score <- score + aa_dict[aa]
        count <- count + 1
      }
    }
  }

  if (count == 0) return(0)
  return(score / count)
}

#' Calculate window scores for all positions in sequence
window_scores <- function(sequence, aa_dict, half_window_size,
                         ignore_consecutive_prolines = FALSE) {
  seq_len <- length(sequence)
  scores <- numeric(seq_len)

  for (i in 1:seq_len) {
    scores[i] <- window_score(sequence, i, aa_dict, half_window_size,
                              ignore_consecutive_prolines)
  }

  return(scores)
}

#' Calculate FoldIndex for sequence
#' FoldIndex formula: 2.785(H) - |R| - 1.151
fold_index <- function(sequence, half_window_size, window_size) {
  charges <- window_scores(sequence, CHARGE, half_window_size, FALSE)
  hydrophobicities <- window_scores(sequence, HYDROPHOBICITY, half_window_size, FALSE)

  fold_index_list <- 2.785 * hydrophobicities - abs(charges) - 1.151

  return(super_window_scores(sequence, fold_index_list, half_window_size,
                            window_size, NULL))
}

#' Calculate super-window averaged scores
super_window_scores <- function(sequence, window_scores, half_window_size,
                               window_size, fold_index_scores = NULL) {
  seq_len <- length(sequence)
  max_pos <- seq_len - window_size + 1

  if (max_pos <= 0) return(numeric(0))

  scores <- numeric(max_pos)

  for (i in 1:max_pos) {
    if (!is.null(fold_index_scores) &&
        i <= length(fold_index_scores) &&
        fold_index_scores[i] > 0) {
      scores[i] <- NA
      next
    }

    score <- 0.0
    weights <- 0.0

    for (j in i:(i + window_size - 1)) {
      window <- get_window(sequence, j, half_window_size)
      window_len <- window$end - window$start + 1
      score <- score + window_len * window_scores[j]
      weights <- weights + window_len
    }

    scores[i] <- score / weights
  }

  return(scores)
}

#' Main PAPA classification function
classify_papa <- function(sequence, window_size, half_window_size,
                         ignore_fold_index = FALSE) {

  fold_index_list <- fold_index(sequence, half_window_size, window_size)

  window_propensities <- window_scores(sequence, PROPENSITIES,
                                       half_window_size, TRUE)

  if (ignore_fold_index) {
    scores <- super_window_scores(sequence, window_propensities,
                                  half_window_size, window_size, NULL)
  } else {
    scores <- super_window_scores(sequence, window_propensities,
                                  half_window_size, window_size,
                                  fold_index_list)
  }

  valid_scores <- scores[!is.na(scores)]

  if (length(valid_scores) == 0) {
    max_score <- -1.0
    max_position <- -1
  } else {
    max_score <- max(valid_scores)
    max_position <- which.max(scores) - 1  # 0-indexed
  }

  return(list(
    max_score = max_score,
    max_position = max_position,
    scores = scores,
    fold_index_list = fold_index_list
  ))
}

################################################################################
# PROTEIN ANALYSIS WRAPPER
################################################################################

#' Analyze single protein for PAPA prediction
#' @return List with prediction results or NULL if no significant hit
analyze_protein_papa <- function(protein_id, organism, taxonomy, taxid, sequence,
                                 window_size, half_window_size, threshold,
                                 ignore_fold_index) {

  sequence_chars <- strsplit(toupper(sequence), "")[[1]]

  if (length(sequence_chars) <= window_size) {
    return(NULL)  # Sequence too short
  }

  result <- classify_papa(sequence_chars, window_size, half_window_size,
                         ignore_fold_index)

  org_clean <- os_name_cleaning(organism)

  # Return result (even below threshold, for completeness)
  list(
    protein_id = protein_id,
    organism = org_clean,
    taxonomy = taxonomy,
    taxid = taxid,
    score = result$max_score,
    position = result$max_position,
    above_threshold = result$max_score >= threshold
  )
}

################################################################################
# PARALLEL PROCESSING FROM TABULAR CACHE
################################################################################

#' Process tabular cache with parallel workers
process_tabular_cache <- function(cache_file, num_cores, window_size,
                                  half_window_size, threshold, ignore_fold_index) {
  cat("\n=== Loading Tabular Cache ===\n")
  cat(sprintf("File: %s\n", basename(cache_file)))

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
  cat(sprintf(" -> Loaded %s proteins in %.1f sec\n",
              format(nrow(cache_data), big.mark = ","),
              as.numeric(load_time)))

  # Filter sequences that are too short
  cache_data <- cache_data[!is.na(cache_data$Sequence) &
                           nchar(cache_data$Sequence) > window_size, ]
  cat(sprintf(" -> %s proteins with sequences > %d aa\n",
              format(nrow(cache_data), big.mark = ","), window_size))

  cat("\n=== Processing with PAPA Algorithm ===\n")
  cat(sprintf("Using %d CPU cores\n", num_cores))
  cat(sprintf("Window size: %d\n", window_size))
  cat(sprintf("Score threshold: %.2f\n", threshold))
  cat(sprintf("Ignore FoldIndex: %s\n\n", ignore_fold_index))

  n_rows <- nrow(cache_data)
  chunk_size <- ceiling(n_rows / num_cores)

  process_start <- Sys.time()

  results_list <- mclapply(1:num_cores, function(worker_id) {
    start_idx <- (worker_id - 1) * chunk_size + 1
    end_idx <- min(worker_id * chunk_size, n_rows)

    if (start_idx > n_rows) return(list())

    local_results <- list()

    for (i in start_idx:end_idx) {
      row <- cache_data[i, ]

      pred <- analyze_protein_papa(
        protein_id = row$Protein_ID,
        organism = row$Organism,
        taxonomy = row$Full_Taxonomy,
        taxid = row$NCBI_TaxID,
        sequence = row$Sequence,
        window_size = window_size,
        half_window_size = half_window_size,
        threshold = threshold,
        ignore_fold_index = ignore_fold_index
      )

      if (!is.null(pred)) {
        local_results[[length(local_results) + 1]] <- pred
      }
    }

    local_results
  }, mc.cores = num_cores)

  process_time <- difftime(Sys.time(), process_start, units = "secs")

  # Flatten results
  all_results <- do.call(c, results_list)

  cat(sprintf(" -> Processing complete in %.1f sec\n", as.numeric(process_time)))
  cat(sprintf(" -> Analyzed %s proteins\n", format(length(all_results), big.mark = ",")))

  # Count above threshold
  above_threshold <- sum(sapply(all_results, function(x) x$above_threshold))
  cat(sprintf(" -> %s proteins with PAPA score >= %.2f\n",
              format(above_threshold, big.mark = ","), threshold))

  all_results
}

################################################################################
# OUTPUT FUNCTIONS
################################################################################

#' Write results to TSV file
write_papa_results <- function(results, output_file, threshold) {
  cat(sprintf("\nWriting results to: %s\n", output_file))

  # Create data frame from results
  df <- data.frame(
    Protein_ID = sapply(results, function(x) x$protein_id),
    Organism = sapply(results, function(x) x$organism),
    Taxonomy = sapply(results, function(x) x$taxonomy),
    NCBI_TaxID = sapply(results, function(x) x$taxid),
    PAPA_Score = sapply(results, function(x) round(x$score, 6)),
    PAPA_Position = sapply(results, function(x) x$position),
    Above_Threshold = sapply(results, function(x) x$above_threshold),
    stringsAsFactors = FALSE
  )

  # Sort by score descending
  df <- df[order(df$PAPA_Score, decreasing = TRUE), ]

  write_tsv(df, output_file)

  cat(sprintf(" -> Wrote %s entries\n", format(nrow(df), big.mark = ",")))

  # Summary stats
  cat("\n=== Summary Statistics ===\n")
  cat(sprintf("Total proteins analyzed: %s\n", format(nrow(df), big.mark = ",")))
  cat(sprintf("Proteins above threshold (%.2f): %s\n",
              threshold, format(sum(df$Above_Threshold), big.mark = ",")))
  cat(sprintf("Max PAPA score: %.4f\n", max(df$PAPA_Score)))
  cat(sprintf("Mean PAPA score: %.4f\n", mean(df$PAPA_Score)))
  cat(sprintf("Median PAPA score: %.4f\n", median(df$PAPA_Score)))
}

################################################################################
# STREAMING FALLBACK (FOR RAW .dat FILES)
################################################################################

#' Parse SwissProt .dat file (streaming fallback)
parse_swissprot_streaming <- function(file_path) {
  if (grepl("\\.gz$", file_path)) {
    con <- gzfile(file_path, "r")
  } else {
    con <- file(file_path, "r")
  }
  on.exit(close(con))

  records <- list()
  current_id <- NULL
  current_org <- ""
  current_oc <- ""
  current_ox <- ""
  current_seq <- c()

  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    if (startsWith(line, "//")) {
      if (!is.null(current_id) && length(current_seq) > 0) {
        records[[length(records) + 1]] <- list(
          id = current_id,
          organism = current_org,
          taxonomy = current_oc,
          taxid = current_ox,
          sequence = paste(current_seq, collapse = "")
        )
      }
      current_id <- NULL
      current_org <- ""
      current_oc <- ""
      current_ox <- ""
      current_seq <- c()
      next
    }

    if (startsWith(line, "ID ")) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      current_id <- parts[2]
    } else if (startsWith(line, "OS ")) {
      current_org <- paste0(current_org, sub("^OS\\s+", "", line))
    } else if (startsWith(line, "OC ")) {
      current_oc <- paste0(current_oc, sub("^OC\\s+", "", line))
    } else if (startsWith(line, "OX ")) {
      ox_match <- regmatches(line, regexpr("NCBI_TaxID=(\\d+)", line))
      if (length(ox_match) > 0) {
        current_ox <- sub("NCBI_TaxID=", "", ox_match)
      }
    } else if (grepl("^\\s+[A-Z]", line)) {
      seq_line <- gsub("[^A-Z]", "", line)
      if (nchar(seq_line) > 0) {
        current_seq <- c(current_seq, seq_line)
      }
    }
  }

  if (!is.null(current_id) && length(current_seq) > 0) {
    records[[length(records) + 1]] <- list(
      id = current_id,
      organism = current_org,
      taxonomy = current_oc,
      taxid = current_ox,
      sequence = paste(current_seq, collapse = "")
    )
  }

  records
}

#' Process raw .dat file (fallback when no cache available)
process_raw_file <- function(input_file, num_cores, window_size,
                            half_window_size, threshold, ignore_fold_index) {
  cat("\n=== Streaming from Raw .dat File ===\n")
  cat("Note: Consider running data_preanalysis.R first for faster processing\n\n")

  cat(sprintf("Parsing: %s\n", basename(input_file)))
  records <- parse_swissprot_streaming(input_file)
  cat(sprintf(" -> Loaded %s proteins\n", format(length(records), big.mark = ",")))

  cat("\n=== Processing with PAPA Algorithm ===\n")
  cat(sprintf("Using %d CPU cores\n", num_cores))

  results_list <- mclapply(records, function(rec) {
    analyze_protein_papa(
      protein_id = rec$id,
      organism = rec$organism,
      taxonomy = rec$taxonomy,
      taxid = rec$taxid,
      sequence = rec$sequence,
      window_size = window_size,
      half_window_size = half_window_size,
      threshold = threshold,
      ignore_fold_index = ignore_fold_index
    )
  }, mc.cores = num_cores)

  # Remove NULLs
  results_list <- results_list[!sapply(results_list, is.null)]

  cat(sprintf(" -> Analyzed %s proteins\n", format(length(results_list), big.mark = ",")))

  results_list
}

################################################################################
# MAIN FUNCTION
################################################################################

main <- function() {
  option_list <- list(
    make_option(
      c("-d", "--db"),
      type = "character",
      default = NULL,
      help = "Database type: 'sprot' or 'trembl'. Auto-detects from cache if not specified.",
      metavar = "TYPE"
    ),
    make_option(
      c("-i", "--input"),
      type = "character",
      default = NULL,
      help = "Direct input file (.dat or .dat.gz). Use this for fallback when no cache available.",
      metavar = "FILE"
    ),
    make_option(
      c("-o", "--output"),
      type = "character",
      default = NULL,
      help = "Output file path (default: auto-generated in data/processed/)",
      metavar = "FILE"
    ),
    make_option(
      c("-c", "--cores"),
      type = "integer",
      default = CORES,
      help = sprintf("Number of CPU cores to use [default: %d]", CORES),
      metavar = "NUMBER"
    ),
    make_option(
      c("-w", "--window_size"),
      type = "integer",
      default = WINDOW_SIZE,
      help = sprintf("Window size for PAPA algorithm [default: %d]", WINDOW_SIZE),
      metavar = "INT"
    ),
    make_option(
      c("-t", "--threshold"),
      type = "double",
      default = SCORE_THRESHOLD,
      help = sprintf("Score threshold for positive prediction [default: %.2f]", SCORE_THRESHOLD),
      metavar = "NUM"
    ),
    make_option(
      c("--ignore_fold_index"),
      action = "store_true",
      default = FALSE,
      help = "Ignore FoldIndex filtering (predict on all regions)"
    )
  )

  opt_parser <- OptionParser(
    option_list = option_list,
    description = paste("\n=== PAPA v1.0 ===\n",
                        "Prion Aggregation Prediction Algorithm.\n",
                        "Uses cached tabular files from data_preanalysis.R for speed."),
    epilogue = paste("\nReferences:",
                     "  Toombs et al. (2010) Mol Cell Biol 30(1):319-332",
                     "\nExamples:",
                     "  Rscript R/PAPA.R --db sprot",
                     "  Rscript R/PAPA.R --db trembl --cores 8",
                     "  Rscript R/PAPA.R --input data/raw/uniprot_sprot_fungi.dat.gz",
                     sep = "\n")
  )

  opt <- parse_args(opt_parser)
  num_cores <- opt$cores
  window_size <- opt$window_size
  half_window_size <- as.integer(window_size / 2)
  threshold <- opt$threshold
  ignore_fold_index <- opt$ignore_fold_index

  cat("\n")
  cat("================================================================================\n")
  cat("PAPA - Prion Aggregation Prediction Algorithm (v1.0)\n")
  cat("================================================================================\n")
  cat("See docs/changelog_PAPA.md for version history.\n")
  cat("================================================================================\n")

  # Determine input source
  cache_file <- NULL
  input_file <- NULL
  db_name <- opt$db

  # Try to find cache first
  if (!is.null(db_name)) {
    cache_file <- find_tabular_cache(db_name)
  } else if (is.null(opt$input)) {
    cache_file <- find_tabular_cache(NULL)
  }

  # Fallback to direct input
  if (is.null(cache_file) && !is.null(opt$input)) {
    input_file <- opt$input
    if (!file.exists(input_file)) {
      stop("Input file not found: ", input_file)
    }
  }

  # Validate we have some input

  if (is.null(cache_file) && is.null(input_file)) {
    cat("\nERROR: No input available.\n")
    cat("Either:\n")
    cat("  1. Run data_preanalysis.R first to create tabular cache, OR\n")
    cat("  2. Use --input to specify a .dat file directly\n\n")
    stop("No valid input source found", call. = FALSE)
  }

  # Detect db_name from cache file if not specified
  if (is.null(db_name) && !is.null(cache_file)) {
    if (grepl("sprot", cache_file)) {
      db_name <- "sprot"
    } else if (grepl("trembl", cache_file)) {
      db_name <- "trembl"
    } else {
      db_name <- "unknown"
    }
  } else if (is.null(db_name)) {
    db_name <- "custom"
  }

  # Generate output filename
  if (is.null(opt$output)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file <- file.path(DATA_PROCESSED_DIR,
                             sprintf("%s_papa_predictions_%s.tsv", db_name, timestamp))
  } else {
    output_file <- opt$output
  }

  # Ensure output directory exists
  if (!dir.exists(dirname(output_file))) {
    dir.create(dirname(output_file), recursive = TRUE)
  }

  cat(sprintf("\nInput source: %s\n",
              ifelse(!is.null(cache_file), basename(cache_file), basename(input_file))))
  cat(sprintf("Output file: %s\n", output_file))
  cat(sprintf("Database: %s\n", db_name))
  cat(sprintf("Cores: %d\n", num_cores))

  # Process
  if (!is.null(cache_file)) {
    results <- process_tabular_cache(cache_file, num_cores, window_size,
                                     half_window_size, threshold, ignore_fold_index)
  } else {
    results <- process_raw_file(input_file, num_cores, window_size,
                               half_window_size, threshold, ignore_fold_index)
  }

  # Write output
  write_papa_results(results, output_file, threshold)

  cat("\n================================================================================\n")
  cat("PAPA analysis complete!\n")
  cat(sprintf("Results saved to: %s\n", output_file))
  cat("================================================================================\n")
}

# Execute main if script is run directly
if (!interactive()) {
  main()
}
