#!/usr/bin/env Rscript

################################################################################
# Prion Domain Prediction Script (R version v13 - Streaming Architecture)
################################################################################
# An R script to read proteomes in SwissProt format and predict prion domains
# Optimized for very large files (100GB+) using streaming processing
# Translated from Perl to R - 2024
# Original Copyright (C) Vladimir Espinosa Angarica 2012
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program comes with absolutely NO WARRANTY
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
  library(stringr)
})

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

# Analysis parameters
WINDOW <- 60      # Size of sliding window
CUTOFF <- 50      # Minimum score threshold for predictions
CORES <- detectCores() - 1  # Number of CPU cores to use
if (CORES < 1) CORES <- 1

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Clean organism name from SwissProt OS line format
os_name_cleaning <- function(osline) {
  osline <- gsub("(,|, and|\\.)$", "", osline)
  
  if (grepl("[^\\(\\)]+\\(.+\\)[^\\(\\)]+$", osline)) {
    sci_name <- gsub("\\.$", "", osline)
  } else {
    match <- regmatches(osline, regexpr("\\S[^\\(]+", osline))
    sci_name <- ifelse(length(match) > 0, match[1], osline)
  }
  
  sci_name <- gsub("\\s+$", "", sci_name)
  return(sci_name)
}

#' Extract organism name from full sequence identifier
get_organism <- function(line) {
  patterns <- c(
    "^>(\\w+);\\sRecName: Full=(.+);SEQUENCE.*\\[(.*)\\]$",
    "^>(\\w+);\\sRecName: Full=(.+);AltName.*\\[(.*)\\]$",
    "^>(\\w+);\\sSubName: Full=(.+);SEQUENCE.*\\[(.*)\\]$",
    "^>(\\w+);\\sSubName: Full=(.+);AltName.*\\[(.*)\\]$"
  )
  
  for (pattern in patterns) {
    if (grepl(pattern, line)) {
      matches <- str_match(line, pattern)
      return(matches[4])
    }
  }
  return(NA)
}

################################################################################
# STREAMING FILE PARSER
################################################################################

#' Stream SwissProt entries from file
#'
#' Reads SwissProt file entry-by-entry without loading entire file into memory.
#' Yields one entry at a time as a list containing sequence data.
#'
#' @param dat_file Path to SwissProt format file (.dat or .dat.gz)
#' @param callback Function to call for each parsed entry
#' @return Total number of entries processed
stream_swissprot_file <- function(dat_file, callback) {
  # Validate file exists
  if (!file.exists(dat_file)) {
    stop("Not a valid PATH to the Sequence file!!!!!\n")
  }
  
  # Open file connection (handles .gz automatically)
  if (grepl("\\.gz$", dat_file, ignore.case = TRUE)) {
    cat("Opening gzipped file stream...\n")
    con <- gzfile(dat_file, "rt")
  } else {
    con <- file(dat_file, "rt")
  }
  
  # Initialize entry variables
  seq_id <- NULL
  description <- ""
  organism <- ""
  sci_name <- ""
  taxonomy_lineage <- ""
  ncbi_taxid <- NA
  sequence <- ""
  entry_count <- 0
  
  # Stream through file line by line
  while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
    
    if (grepl("^ID\\s+(\\w+)\\s+", line)) {
      seq_id <- str_match(line, "^ID\\s+(\\w+)\\s+")[2]
      
    } else if (grepl("^DE\\s+(.*)", line) || grepl("^SQ\\s+(.*)", line)) {
      match <- str_match(line, "^(?:DE|SQ)\\s+(.*)")
      description <- paste0(description, match[2])
      
    } else if (grepl("^OS\\s+(.*)", line)) {
      match <- str_match(line, "^OS\\s+(.*)")
      organism <- ifelse(organism == "", match[2], paste(organism, match[2]))
      sci_name <- os_name_cleaning(organism)
      
    } else if (grepl("^OC\\s+(.*)", line)) {
      match <- str_match(line, "^OC\\s+(.*)")
      taxonomy_part <- gsub(";$", "", match[2])
      taxonomy_lineage <- ifelse(taxonomy_lineage == "", 
                                 taxonomy_part, 
                                 paste(taxonomy_lineage, taxonomy_part, sep = "; "))
      
    } else if (grepl("^OX\\s+NCBI_TaxID=(\\d+)", line)) {
      match <- str_match(line, "^OX\\s+NCBI_TaxID=(\\d+)")
      ncbi_taxid <- match[2]
      
    } else if (grepl("^\\s+(.*)", line)) {
      match <- str_match(line, "^\\s+(.*)")
      seq_content <- gsub("\\s", "", match[2])
      sequence <- paste0(sequence, seq_content)
      
    } else if (grepl("^//", line)) {
      # End of entry - process it
      if (!is.null(seq_id) && nchar(sequence) > 0) {
        entry <- list(
          id = paste0(">", seq_id, "; ", description, "  [", sci_name, "]"),
          seq_id = seq_id,
          organism = sci_name,
          taxonomy = taxonomy_lineage,
          taxid = ncbi_taxid,
          sequence = sequence
        )
        
        # Call the callback function with this entry
        callback(entry)
        entry_count <- entry_count + 1
        
        # Progress indicator every 1000 entries
        if (entry_count %% 1000 == 0) {
          cat(sprintf("\rProcessed %d entries...", entry_count), file = stderr())
        }
      }
      
      # Reset for next entry
      seq_id <- NULL
      description <- ""
      organism <- ""
      sci_name <- ""
      taxonomy_lineage <- ""
      ncbi_taxid <- NA
      sequence <- ""
    }
  }
  
  close(con)
  cat(sprintf("\rProcessed %d entries total.\n", entry_count), file = stderr())
  return(entry_count)
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
  
  corrected_score <- domain_score + correction_factor
  return(round(corrected_score, 3))
}

#' Analyze a single entry for prion domains
#'
#' @param entry List containing sequence information
#' @param window Integer window size
#' @param cutoff Numeric score cutoff
#' @param aa_scores_prd Named vector of amino acid scores
#' @param proline_logodd Named vector of proline pair log-odds
#' @return Prediction result or NULL if no prediction meets cutoff
analyze_entry <- function(entry, window, cutoff, aa_scores_prd, proline_logodd) {
  seq <- entry$sequence
  
  # Skip sequences shorter than window
  if (nchar(seq) < window) return(NULL)
  
  # Extract organism
  organism <- ifelse(is.na(entry$organism) || entry$organism == "", 
                     "Unknown", 
                     entry$organism)
  
  best_result <- list(Score = -Inf, Seq = "", Window = 0)
  
  # Slide window across sequence
  for (i in 0:(nchar(seq) - window)) {
    domain <- substr(seq, i + 1, i + window)
    corrected_score <- calculate_window_score(domain, aa_scores_prd, 
                                              proline_logodd)
    
    if (corrected_score > best_result$Score) {
      best_result <- list(
        Score = corrected_score,
        Seq = domain,
        Window = i
      )
    }
  }
  
  # Only return results that meet cutoff
  if (best_result$Score >= cutoff) {
    return(list(
      seq_id = entry$seq_id,
      full_id = entry$id,
      organism = organism,
      taxonomy = entry$taxonomy,
      taxid = entry$taxid,
      score = best_result$Score,
      window = best_result$Window,
      domain = best_result$Seq
    ))
  }
  
  return(NULL)
}

################################################################################
# PARALLEL STREAMING PROCESSOR
################################################################################

#' Process file with parallel workers using streaming
#'
#' Divides file into chunks by entry count and processes in parallel
#'
#' @param dat_file Path to input file
#' @param num_cores Number of parallel workers
#' @return List of all predictions
parallel_stream_process <- function(dat_file, num_cores) {
  
  # PASS 1: Quick scan to count entries (fast - only reads ID lines)
  cat("\n=== Counting entries for workload distribution ===\n")
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
  
  cat(sprintf("Found %d entries in file\n", total_entries))
  
  # Calculate entries per worker
  entries_per_worker <- ceiling(total_entries / num_cores)
  cat(sprintf("Each of %d workers will process ~%d entries\n", 
              num_cores, entries_per_worker))
  
  # PASS 2: Process with parallel workers
  cat("\n=== Processing entries in parallel ===\n")
  
  # Shared data structures (using environment for pass-by-reference)
  results <- new.env()
  results$predictions <- list()
  results$organism_counts <- list()
  results$lock <- 0  # Simple counter for progress
  
  # Worker function - processes a range of entries
  worker_process <- function(worker_id, start_entry, end_entry, dat_file) {
    local_predictions <- list()
    local_counts <- list()
    
    if (grepl("\\.gz$", dat_file, ignore.case = TRUE)) {
      con <- gzfile(dat_file, "rt")
    } else {
      con <- file(dat_file, "rt")
    }
    
    # Stream to starting position
    entry_num <- 0
    seq_id <- NULL
    description <- ""
    organism <- ""
    sci_name <- ""
    taxonomy_lineage <- ""
    ncbi_taxid <- NA
    sequence <- ""
    
    while (length(line <- readLines(con, n = 1, warn = FALSE)) > 0) {
      
      if (grepl("^ID\\s+(\\w+)\\s+", line)) {
        seq_id <- str_match(line, "^ID\\s+(\\w+)\\s+")[2]
        
      } else if (grepl("^DE\\s+(.*)", line) || grepl("^SQ\\s+(.*)", line)) {
        match <- str_match(line, "^(?:DE|SQ)\\s+(.*)")
        description <- paste0(description, match[2])
        
      } else if (grepl("^OS\\s+(.*)", line)) {
        match <- str_match(line, "^OS\\s+(.*)")
        organism <- ifelse(organism == "", match[2], paste(organism, match[2]))
        sci_name <- os_name_cleaning(organism)
        
      } else if (grepl("^OC\\s+(.*)", line)) {
        match <- str_match(line, "^OC\\s+(.*)")
        taxonomy_part <- gsub(";$", "", match[2])
        taxonomy_lineage <- ifelse(taxonomy_lineage == "", 
                                   taxonomy_part, 
                                   paste(taxonomy_lineage, taxonomy_part, sep = "; "))
        
      } else if (grepl("^OX\\s+NCBI_TaxID=(\\d+)", line)) {
        match <- str_match(line, "^OX\\s+NCBI_TaxID=(\\d+)")
        ncbi_taxid <- match[2]
        
      } else if (grepl("^\\s+(.*)", line)) {
        match <- str_match(line, "^\\s+(.*)")
        seq_content <- gsub("\\s", "", match[2])
        sequence <- paste0(sequence, seq_content)
        
      } else if (grepl("^//", line)) {
        entry_num <- entry_num + 1
        
        # Process if within this worker's range
        if (entry_num >= start_entry && entry_num <= end_entry) {
          if (!is.null(seq_id) && nchar(sequence) > 0) {
            entry <- list(
              id = paste0(">", seq_id, "; ", description, "  [", sci_name, "]"),
              seq_id = seq_id,
              organism = sci_name,
              taxonomy = taxonomy_lineage,
              taxid = ncbi_taxid,
              sequence = sequence
            )
            
            # Count organism
            org_name <- ifelse(is.na(sci_name) || sci_name == "", "Unknown", sci_name)
            if (is.null(local_counts[[org_name]])) {
              local_counts[[org_name]] <- 0
            }
            local_counts[[org_name]] <- local_counts[[org_name]] + 1
            
            # Analyze entry
            pred <- analyze_entry(entry, WINDOW, CUTOFF, 
                                  aa_scores$PrD, proline_pair_logodd)
            
            if (!is.null(pred)) {
              if (is.null(local_predictions[[pred$organism]])) {
                local_predictions[[pred$organism]] <- list()
              }
              local_predictions[[pred$organism]][[pred$full_id]] <- pred
            }
          }
        }
        
        # Stop if past this worker's range
        if (entry_num > end_entry) break
        
        # Reset for next entry
        seq_id <- NULL
        description <- ""
        organism <- ""
        sci_name <- ""
        taxonomy_lineage <- ""
        ncbi_taxid <- NA
        sequence <- ""
      }
    }
    
    close(con)
    
    return(list(predictions = local_predictions, counts = local_counts))
  }
  
  # Launch parallel workers
  worker_ranges <- list()
  for (i in 1:num_cores) {
    start <- (i - 1) * entries_per_worker + 1
    end <- min(i * entries_per_worker, total_entries)
    worker_ranges[[i]] <- c(start, end)
  }
  
  cat(sprintf("Launching %d parallel workers...\n", num_cores))
  
  results_list <- mclapply(1:num_cores, function(i) {
    worker_process(i, worker_ranges[[i]][1], worker_ranges[[i]][2], dat_file)
  }, mc.cores = num_cores)
  
  # Merge results from all workers
  cat("\n=== Merging results from workers ===\n")
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
  
  return(list(predictions = all_predictions, counts = all_counts))
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
      help = "PATH to SwissProt format file (MANDATORY). Can be .dat or .dat.gz",
      metavar = "FILE"
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
      help = "Output file for tab-separated format with taxonomy (default: auto-generated)",
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
    description = paste("\n=== Prion Domain Prediction v13 - Streaming Architecture ===\n",
                        "Optimized for very large files (100GB+) using streaming processing"),
    epilogue = paste("\nNew in v13: Streaming architecture for memory-efficient processing",
                     "\nReferences:",
                     "Angarica et al. (2013) BMC Genomics",
                     "Alberti et al. (2009) Cell 137(1):146-158",
                     "\nExample usage:",
                     "  ./prion_parser.R -s uniprot_trembl_fungi.dat.gz",
                     "  ./prion_parser.R -s input.dat -o results.dat -t results.tsv -c 16",
                     sep = "\n")
  )
  
  opt <- parse_args(opt_parser)
  
  if (is.null(opt$`seq-file`)) {
    print_help(opt_parser)
    stop("\nError: --seq-file is required!\n", call. = FALSE)
  }
  
  # Auto-generate output filenames
  if (is.null(opt$`output-file`)) {
    input_basename <- basename(opt$`seq-file`)
    output_name <- gsub("\\.dat(\\.gz)?$", "_predictions.dat", 
                        input_basename, ignore.case = TRUE)
    opt$`output-file` <- output_name
    cat(paste0("Output file: ", output_name, " (auto-generated)\n"))
  }
  
  if (is.null(opt$`tabular-output`)) {
    input_basename <- basename(opt$`seq-file`)
    tabular_name <- gsub("\\.dat(\\.gz)?$", "_predictions_table.tsv", 
                         input_basename, ignore.case = TRUE)
    opt$`tabular-output` <- tabular_name
    cat(paste0("Tabular output: ", tabular_name, " (auto-generated)\n"))
  }
  
  num_cores <- opt$cores
  
  # Process file with streaming
  cat("\n=== Starting Streaming Analysis ===\n")
  cat(paste0("Using ", num_cores, " CPU cores\n"))
  cat(paste0("Window size: ", WINDOW, " amino acids\n"))
  cat(paste0("Score cutoff: ", CUTOFF, "\n\n"))
  
  results <- parallel_stream_process(opt$`seq-file`, num_cores)
  all_predictions <- results$predictions
  all_counts <- results$counts
  
  cat(sprintf("\n✓ Found predictions in %d organisms\n", length(all_predictions)))
  
  # Write outputs
  cat("\n=== Writing output files ===\n")
  
  # OUTPUT 1: Original format
  cat(paste0("Writing original format: ", opt$`output-file`, "\n"))
  sink(opt$`output-file`)
  
  for (organism in names(all_predictions)) {
    predictions <- all_predictions[[organism]]
    total_predictions <- length(predictions)
    
    cat(paste0(">", organism, ": Total=", total_predictions, "\n"))
    
    for (pred in predictions) {
      cat(paste0("\t", pred$seq_id, 
                 "\tWindow Position=", pred$window, 
                 "; Score=", pred$score, 
                 " | Prion Domain: ", pred$domain, "\n"))
    }
    cat("\n")
  }
  
  sink()
  
  # OUTPUT 2: Tabular format
  cat(paste0("Writing tabular format: ", opt$`tabular-output`, "\n"))
  
  tabular_data <- data.frame(
    Protein_ID = character(),
    Organism = character(),
    Taxonomy = character(),
    NCBI_TaxID = character(),
    Window_Position = integer(),
    Score = numeric(),
    Prion_Domain = character(),
    stringsAsFactors = FALSE
  )
  
  for (organism in names(all_predictions)) {
    predictions <- all_predictions[[organism]]
    
    for (pred in predictions) {
      taxid_value <- ifelse(is.na(pred$taxid), "NA", pred$taxid)
      
      tabular_data <- rbind(tabular_data, data.frame(
        Protein_ID = pred$seq_id,
        Organism = pred$organism,
        Taxonomy = pred$taxonomy,
        NCBI_TaxID = taxid_value,
        Window_Position = pred$window,
        Score = pred$score,
        Prion_Domain = pred$domain,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  tabular_data <- tabular_data[order(tabular_data$Organism, -tabular_data$Score), ]
  
  write.table(tabular_data, file = opt$`tabular-output`, 
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Summary
  cat("\n=== Analysis Complete ===\n")
  cat(paste0("Total organisms: ", length(all_predictions), "\n"))
  cat(paste0("Total predictions: ", sum(sapply(all_predictions, length)), "\n"))
  cat("\nOutput files:\n")
  cat(paste0("  1. ", opt$`output-file`, " (original format)\n"))
  cat(paste0("  2. ", opt$`tabular-output`, " (tabular with taxonomy & TaxID)\n"))
  cat("\n")
}

################################################################################
# SCRIPT ENTRY POINT
################################################################################

if (!interactive()) {
  main()
}