#!/usr/bin/env Rscript

################################################################################
# FunGuild Pre-Analysis Pipeline for UniProt Fungi Database (v1.0.1)
################################################################################
# Prepares taxonomic-guild mappings from complete UniProt fungi proteomes
#
# See docs/changelog_data_preanalysis.md for version history.
#
# This script:
# 1. Downloads FunGuild database (timestamped, cached for 14 days)
# 2. Locates UniProt fungi proteome file
# 2.5 Converts .dat to tabular TSV format (single-pass, cached)
# 3. Parses tabular TSV (vectorized, ultra-fast)
# 4. Queries FunGuild for each unique genus
# 5. Generates diagnostic reports, unmatched genera log, and visualizations
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(FUNGuildR)
  library(patchwork)
})

################################################################################
# CONFIGURATION
################################################################################

# Directory structure (3-tier)
# data/raw/       - User-provided input files (.dat, .dat.gz)
# data/cache/     - Intermediate/cached files (converted TSV, FunGuild DB)
# data/processed/ - Final outputs for downstream scripts (taxonomy, guild mapping)
DATA_RAW_DIR <- "data/raw"
DATA_CACHE_DIR <- "data/cache"
DATA_PROCESSED_DIR <- "data/processed"
REPORTS_DIR <- "reports"

# Database URLs
FUNGUILD_DB_URL <- "http://www.stbates.org/funguild_db_2.php"

# Timestamp for this run
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

# Cache validity (days)
CACHE_VALIDITY_DAYS <- 14

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Clean organism name from SwissProt OS line format
#' Vectorized version for use with mutate()
#'
#' Handles special cases:
#' - Square brackets: [Bisifusarium] -> Bisifusarium
#' - Single quotes: 'Aporospora -> Aporospora
#' - Trailing punctuation
#' - Parenthetical information
os_name_cleaning_vec <- function(osline) {
  # Remove square brackets around genus names (e.g., [Bisifusarium] -> Bisifusarium)
  cleaned <- str_replace_all(osline, "\\[([^\\]]+)\\]", "\\1")


  # Remove single quotes from names (e.g., 'Aporospora -> Aporospora)
  cleaned <- str_replace_all(cleaned, "'", "")

  # Remove trailing punctuation
  cleaned <- str_replace_all(cleaned, "(,|, and|\\.)$", "")

  # Handle parenthetical patterns
  # If pattern is: text (stuff) more_text -> keep as is
  # Otherwise extract text before parenthesis
  result <- case_when(
    str_detect(cleaned, "[^\\(\\)]+\\(.+\\)[^\\(\\)]+$") ~ str_replace(cleaned, "\\.$", ""),
    TRUE ~ str_extract(cleaned, "^[^\\(]+")
  )

  # Trim whitespace
  str_trim(result)
}

#' Parse scientific name into genus and species (vectorized)
#' @param sci_name Character vector of scientific names
#' @return List with genus and species character vectors
parse_scientific_name_vec <- function(sci_name) {
  # Remove strain/isolate information in parentheses
  cleaned <- str_replace_all(sci_name, "\\s*\\([^\\)]+\\)", "")
  cleaned <- str_trim(cleaned)

  # Extract genus (first word)
  genus <- str_extract(cleaned, "^\\S+")

  # Extract species (second word, if present)
  # Match: first_word whitespace second_word
  species <- str_match(cleaned, "^\\S+\\s+(\\S+)")[, 2]

  list(genus = genus, species = species)
}

#' Safe FunGuild query with error handling
safe_funguild_query <- function(genus, funguild_db) {
  tryCatch({
    result <- FUNGuildR::funguild_query(genus, field = "taxon", funguild_db)

    if (is.null(result) || nrow(result) == 0) {
      return(NULL)
    }

    result$query_genus <- genus
    result$query_date <- Sys.Date()

    return(result)
  }, error = function(e) {
    return(NULL)
  })
}

################################################################################
# STEP 0: CREATE DIRECTORY STRUCTURE
################################################################################

create_directories <- function() {
  cat("\n=== Step 0: Creating Directory Structure ===\n")

  dirs <- c(DATA_RAW_DIR, DATA_CACHE_DIR, DATA_PROCESSED_DIR, REPORTS_DIR)

  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      cat(sprintf(" ↳ Created directory: %s\n", d))
    } else {
      cat(sprintf(" ↳ Directory exists: %s\n", d))
    }
  }
}

################################################################################
# STEP 1: DOWNLOAD FUNGUILD DATABASE
################################################################################

check_existing_funguild_db <- function() {
  existing_dbs <- list.files(DATA_CACHE_DIR,
                             pattern = "^funguild_db_\\d{8}_\\d{6}\\.rds$",
                             full.names = TRUE)

  if (length(existing_dbs) == 0) {
    return(NULL)
  }

  db_info <- data.frame(
    file = existing_dbs,
    basename = basename(existing_dbs),
    stringsAsFactors = FALSE
  )

  db_info$timestamp_str <- gsub("funguild_db_(\\d{8}_\\d{6})\\.rds", "\\1", db_info$basename)
  db_info$timestamp <- as.POSIXct(db_info$timestamp_str, format = "%Y%m%d_%H%M%S")
  db_info <- db_info[order(db_info$timestamp, decreasing = TRUE), ]

  most_recent <- db_info[1, ]
  age_days <- as.numeric(difftime(Sys.time(), most_recent$timestamp, units = "days"))

  return(list(
    file = most_recent$file,
    age_days = age_days,
    timestamp = most_recent$timestamp
  ))
}

download_funguild_db <- function() {
  cat("\n=== Step 1: FunGuild Database ===\n")

  existing <- check_existing_funguild_db()

  if (!is.null(existing)) {
    cat("Found existing FunGuild database:\n")
    cat(sprintf("  File: %s\n", basename(existing$file)))
    cat(sprintf("  Age: %.1f days\n", existing$age_days))

    if (existing$age_days < CACHE_VALIDITY_DAYS) {
      cat(sprintf(" ↳ Database is recent (< %d days old), using cached version\n",
                  CACHE_VALIDITY_DAYS))
      funguild_db <- readRDS(existing$file)
      cat(sprintf(" ↳ Loaded database with %d records\n", nrow(funguild_db)))
      return(list(db = funguild_db, file = existing$file, is_new = FALSE))
    } else {
      cat(sprintf(" ↳ Database is old (≥ %d days), downloading fresh version\n",
                  CACHE_VALIDITY_DAYS))
    }
  } else {
    cat("No existing FunGuild database found.\n")
  }

  cat(sprintf("Downloading from: %s\n", FUNGUILD_DB_URL))

  tryCatch({
    funguild_db <- get_funguild_db(db = FUNGUILD_DB_URL)

    if (is.null(funguild_db) || nrow(funguild_db) == 0) {
      stop("Failed to download FunGuild database or database is empty")
    }

    cat(sprintf(" ↳ Downloaded database with %d records\n", nrow(funguild_db)))

    new_file <- file.path(DATA_CACHE_DIR, paste0("funguild_db_", timestamp, ".rds"))

    return(list(db = funguild_db, file = new_file, is_new = TRUE))

  }, error = function(e) {
    stop(sprintf("Error downloading FunGuild database: %s", e$message))
  })
}

################################################################################
# STEP 2: LOCATE UNIPROT FUNGI FILE
################################################################################

find_uniprot_databases <- function() {
  cat("\n=== Step 2: Locating UniProt Fungi Databases ===\n")

  existing_files <- list.files(
    DATA_RAW_DIR,
    pattern = "^uniprot_(sprot|trembl)_fungi\\.dat(\\.gz)?$",
    full.names = TRUE,
    ignore.case = FALSE
  )

  if (length(existing_files) > 0) {
    cat(sprintf("Found %d existing database(s):\n", length(existing_files)))
    for (i in seq_along(existing_files)) {
      file_size <- file.size(existing_files[i]) / (1024^3)
      cat(sprintf("  [%d] %s (%.2f GB)\n", i, basename(existing_files[i]), file_size))
    }
    return(existing_files)
  }

  return(character(0))
}

################################################################################
# STEP 2.5: CONVERT UNIPROT .dat TO TABULAR FORMAT
################################################################################

#' Check for existing tabular conversion
#' @param dat_file Path to the .dat file
#' @return Path to existing tabular file or NULL
check_existing_tabular <- function(dat_file) {
  db_name <- if (grepl("sprot", basename(dat_file), ignore.case = TRUE)) {
    "sprot"
  } else if (grepl("trembl", basename(dat_file), ignore.case = TRUE)) {
    "trembl"
  } else {
    "uniprot"
  }

  # Look for existing tabular files (same db type)
  existing_tabular <- list.files(
    DATA_CACHE_DIR,
    pattern = paste0("^", db_name, "_fungi_tabular_\\d{8}_\\d{6}\\.tsv$"),
    full.names = TRUE
  )

  if (length(existing_tabular) == 0) {
    return(NULL)
  }

  # Extract timestamps and find most recent
  tab_info <- data.frame(
    file = existing_tabular,
    basename = basename(existing_tabular),
    stringsAsFactors = FALSE
  )

  tab_info$timestamp_str <- gsub(
    paste0(db_name, "_fungi_tabular_(\\d{8}_\\d{6})\\.tsv"),
    "\\1",
    tab_info$basename
  )
  tab_info$timestamp <- as.POSIXct(tab_info$timestamp_str, format = "%Y%m%d_%H%M%S")
  tab_info <- tab_info[order(tab_info$timestamp, decreasing = TRUE), ]

  most_recent <- tab_info[1, ]
  age_days <- as.numeric(difftime(Sys.time(), most_recent$timestamp, units = "days"))

  cat(sprintf(" ↳ Found existing tabular conversion:\n"))
  cat(sprintf("     %s (%.1f days old)\n", most_recent$basename, age_days))

  if (age_days < CACHE_VALIDITY_DAYS) {
    return(most_recent$file)
  } else {
    cat(sprintf(" ↳ Tabular file is old (≥ %d days), will re-convert\n",
                CACHE_VALIDITY_DAYS))
    return(NULL)
  }
}

#' Single-pass streaming converter: .dat -> TSV with one row per entry
#' Includes protein sequences for downstream  analysis
#' @param dat_file Path to UniProt .dat or .dat.gz file
#' @return Path to created TSV file
convert_uniprot_to_tabular <- function(dat_file) {
  cat("\n=== Step 2.5: Converting UniProt .dat to Tabular Format ===\n")
  cat(sprintf("Input: %s\n", basename(dat_file)))

  db_name <- if (grepl("sprot", basename(dat_file), ignore.case = TRUE)) {
    "sprot"
  } else if (grepl("trembl", basename(dat_file), ignore.case = TRUE)) {
    "trembl"
  } else {
    "uniprot"
  }

  tabular_file <- file.path(DATA_CACHE_DIR, paste0(db_name, "_fungi_tabular_", timestamp, ".tsv"))

  cat(sprintf(" ↳ Creating: %s\n", basename(tabular_file)))
  cat(" ↳ Single-pass streaming conversion (this may take a while for large files)...\n")
  cat(" ↳ Includes protein sequences for downstream analyses\n")

  start_time <- Sys.time()

  # Open input connection
  con <- if (grepl("\\.gz$", dat_file, ignore.case = TRUE)) {
    gzfile(dat_file, "rt")
  } else {
    file(dat_file, "rt")
  }
  on.exit(close(con), add = TRUE)

  # Pre-allocate vectors (estimate based on file size)
  file_size_mb <- file.size(dat_file) / (1024^2)
  # Rough estimate: ~5KB per entry compressed, ~50KB uncompressed
  if (grepl("\\.gz$", dat_file, ignore.case = TRUE)) {
    estimated_entries <- min(ceiling(file_size_mb * 1000 / 5), 5e7)
  } else {
    estimated_entries <- min(ceiling(file_size_mb * 1000 / 50), 5e7)
  }

  cat(sprintf(" ↳ Estimated entries: ~%s\n", format(estimated_entries, big.mark = ",")))

  ID_vec <- character(estimated_entries)
  OS_vec <- character(estimated_entries)
  OC_vec <- character(estimated_entries)
  SEQ_vec <- character(estimated_entries)
  entry_idx <- 1L

  # Current entry accumulators
  current_id <- ""
  current_os <- ""
  current_oc <- ""
  current_seq <- ""
  in_sequence <- FALSE

  # Progress tracking
  lines_read <- 0L
  last_progress <- 0L

  while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
    lines_read <- lines_read + 1L

    # Progress every 1M lines
    if (lines_read - last_progress >= 1000000L) {
      cat(sprintf("\r ↳ Read %s lines, %s entries...",
                  format(lines_read, big.mark = ","),
                  format(entry_idx - 1L, big.mark = ",")))
      last_progress <- lines_read
    }

    if (startsWith(line, "ID ")) {
      # Extract ID (format: "ID   Q8XXX   ...")
      current_id <- sub("^ID\\s+(\\w+)\\s+.*", "\\1", line)

    } else if (startsWith(line, "OS ")) {
      # Accumulate OS lines (can span multiple lines)
      os_part <- sub("^OS\\s+", "", line)
      if (current_os == "") {
        current_os <- os_part
      } else {
        current_os <- paste0(current_os, " ", os_part)
      }

    } else if (startsWith(line, "OC ")) {
      # Accumulate OC lines (can span multiple lines)
      oc_part <- sub("^OC\\s+", "", line)
      # Remove trailing semicolon if present
      oc_part <- sub(";$", "", oc_part)
      if (current_oc == "") {
        current_oc <- oc_part
      } else {
        current_oc <- paste0(current_oc, "; ", oc_part)
      }

    } else if (startsWith(line, "SQ ")) {
      # Start of sequence section
      in_sequence <- TRUE

    } else if (in_sequence && !startsWith(line, "//") && !startsWith(line, "ID ") && 
               grepl("[A-Z]", line)) {
      # Sequence line: any line with uppercase letters (while in sequence section)
      # Stops when we hit "//" or start of next entry
      seq_part <- gsub("[^A-Z]", "", line)
      if (seq_part != "") {
        current_seq <- paste0(current_seq, seq_part)
      }

    } else if (startsWith(line, "//")) {
      # End of entry - save if we have data
      if (current_id != "") {
        # Expand vectors if needed
        if (entry_idx > length(ID_vec)) {
          new_size <- length(ID_vec) * 2L
          ID_vec <- c(ID_vec, character(new_size - length(ID_vec)))
          OS_vec <- c(OS_vec, character(new_size - length(OS_vec)))
          OC_vec <- c(OC_vec, character(new_size - length(OC_vec)))
          SEQ_vec <- c(SEQ_vec, character(new_size - length(SEQ_vec)))
        }

        ID_vec[entry_idx] <- current_id
        OS_vec[entry_idx] <- current_os
        OC_vec[entry_idx] <- current_oc
        SEQ_vec[entry_idx] <- current_seq
        entry_idx <- entry_idx + 1L
      }

      # Reset for next entry
      current_id <- ""
      current_os <- ""
      current_oc <- ""
      current_seq <- ""
      in_sequence <- FALSE
    }
  }

  cat("\n")

  # Trim vectors to actual length
  actual_entries <- entry_idx - 1L
  ID_vec <- ID_vec[1:actual_entries]
  OS_vec <- OS_vec[1:actual_entries]
  OC_vec <- OC_vec[1:actual_entries]
  SEQ_vec <- SEQ_vec[1:actual_entries]

  # Create data frame
  tabular_data <- data.frame(
    Protein_ID = ID_vec,
    Organism = OS_vec,
    Full_Taxonomy = OC_vec,
    Sequence = SEQ_vec,
    stringsAsFactors = FALSE
  )

  # Write TSV
  write_tsv(tabular_data, tabular_file, na = "")

  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  file_size_out <- file.size(tabular_file) / (1024^2)

  cat(sprintf(" ↳ Conversion complete in %.1f minutes\n", as.numeric(elapsed)))
  cat(sprintf(" ↳ Output: %s entries → %s (%.1f MB)\n",
              format(nrow(tabular_data), big.mark = ","),
              basename(tabular_file),
              file_size_out))

  return(tabular_file)
}

#' Select UniProt database and convert to tabular format
#' @param db_choice Optional command-line database choice
#' @return List with dat_file and tabular_file paths
select_or_download_uniprot <- function(db_choice = NULL) {
  existing_files <- find_uniprot_databases()

  if (length(existing_files) == 0) {
    stop(paste(
      "\nError: No UniProt fungi databases found in data/raw/ directory.",
      "\nPlease download one of the following and place in data/raw/:",
      "\n  - uniprot_sprot_fungi.dat or .dat.gz (Swiss-Prot, ~0.06 GB)",
      "\n  - uniprot_trembl_fungi.dat or .dat.gz (TrEMBL, ~72 GB)",
      "\nDownload from:",
      "\n  https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/",
      sep = ""
    ))
  }

  # Select the database file
  dat_file <- NULL

  if (length(existing_files) == 1) {
    cat(sprintf("\n ↳ Found one database: %s\n", basename(existing_files[1])))
    cat(" ↳ Using automatically\n")
    dat_file <- existing_files[1]

  } else {
    # Multiple databases - check if choice was provided
    if (!is.null(db_choice)) {
      if (tolower(db_choice) == "sprot") {
        match_idx <- grep("sprot", basename(existing_files), ignore.case = TRUE)
        if (length(match_idx) > 0) {
          dat_file <- existing_files[match_idx[1]]
          cat(sprintf("\n ↳ Using database: %s (matched 'sprot')\n", basename(dat_file)))
        }
      } else if (tolower(db_choice) == "trembl") {
        match_idx <- grep("trembl", basename(existing_files), ignore.case = TRUE)
        if (length(match_idx) > 0) {
          dat_file <- existing_files[match_idx[1]]
          cat(sprintf("\n ↳ Using database: %s (matched 'trembl')\n", basename(dat_file)))
        }
      }

      if (is.null(dat_file)) {
        choice_num <- suppressWarnings(as.integer(db_choice))
        if (!is.na(choice_num) && choice_num >= 1 && choice_num <= length(existing_files)) {
          dat_file <- existing_files[choice_num]
          cat(sprintf("\n ↳ Using database: %s (choice %d)\n", basename(dat_file), choice_num))
        }
      }

      if (is.null(dat_file)) {
        cat(sprintf("\nError: Invalid database choice '%s'\n", db_choice))
        cat("Available databases:\n")
        for (i in seq_along(existing_files)) {
          cat(sprintf("  %d: %s\n", i, basename(existing_files[i])))
        }
        stop("Please use --db with 'sprot', 'trembl', or a valid number (1, 2, etc.)")
      }

    } else {
      # No command-line choice - interactive selection
      cat("\nMultiple databases found. Please select:\n")
      for (i in seq_along(existing_files)) {
        cat(sprintf("  [%d] Use %s\n", i, basename(existing_files[i])))
      }

      if (!interactive()) {
        stop(paste(
          "\nError: Multiple databases found but no --db argument provided.",
          "\nWhen using Rscript with multiple databases, specify which to use:",
          "\n  Rscript R/data_preanalysis.R --db sprot",
          "\n  Rscript R/data_preanalysis.R --db trembl",
          "\n  Rscript R/data_preanalysis.R --db 1",
          "\nOr run interactively in R console/RStudio to select.",
          sep = ""
        ))
      }

      choice_str <- readline(prompt = "\nSelect option: ")
      choice <- suppressWarnings(as.integer(choice_str))

      if (is.na(choice) || choice < 1 || choice > length(existing_files)) {
        stop("Invalid selection - must be a number between 1 and ", length(existing_files))
      }

      dat_file <- existing_files[choice]
      cat(sprintf(" ↳ Using: %s\n", basename(dat_file)))
    }
  }

  # Check for existing tabular conversion
  existing_tabular <- check_existing_tabular(dat_file)

  if (!is.null(existing_tabular)) {
    cat(sprintf(" ↳ Using cached tabular conversion\n"))
    return(list(
      dat_file = dat_file,
      tabular_file = existing_tabular,
      used_cache = TRUE
    ))
  }

  # Convert to tabular format
  tabular_file <- convert_uniprot_to_tabular(dat_file)

  return(list(
    dat_file = dat_file,
    tabular_file = tabular_file,
    used_cache = FALSE
  ))
}

################################################################################
# STEP 3: PARSE TABULAR TSV (ULTRA-FAST VECTORIZED)
################################################################################

#' Extract taxonomy from tabular TSV file (vectorized, no loops)
#' @param tabular_file Path to TSV file from convert_uniprot_to_tabular()
#' @return List with taxonomy data and metadata
extract_taxonomy_from_tabular <- function(tabular_file) {
  cat("\n=== Step 3: Parsing Tabular Data (Ultra-Fast Vectorized) ===\n")
  cat(sprintf("Loading: %s\n", basename(tabular_file)))

  start_time <- Sys.time()

  # Fast TSV read with readr
  taxonomy_raw <- read_tsv(
    tabular_file,
    col_types = cols(
      Protein_ID = col_character(),
      Organism = col_character(),
      Full_Taxonomy = col_character()
    ),
    show_col_types = FALSE,
    progress = TRUE
  )

  load_time <- difftime(Sys.time(), start_time, units = "secs")
  cat(sprintf(" ↳ Loaded %s rows in %.1f sec\n",
              format(nrow(taxonomy_raw), big.mark = ","),
              as.numeric(load_time)))

  # Determine db_name from filename
  db_name <- if (grepl("sprot", basename(tabular_file), ignore.case = TRUE)) {
    "sprot"
  } else if (grepl("trembl", basename(tabular_file), ignore.case = TRUE)) {
    "trembl"
  } else {
    "uniprot"
  }

  taxonomy_file <- file.path(DATA_PROCESSED_DIR, paste0(db_name, "_fungi_taxonomy_", timestamp, ".tsv"))

  cat(" ↳ Parsing organism names (vectorized)...\n")
  parse_start <- Sys.time()

  # Vectorized parsing using dplyr and stringr
  parsed <- taxonomy_raw %>%
    mutate(
      # Clean organism names (vectorized version of os_name_cleaning)
      # This now handles square brackets and single quotes
      Organism_clean = os_name_cleaning_vec(Organism)
    ) %>%
    mutate(
      # Parse genus and species from cleaned organism name
      Genus = str_extract(Organism_clean, "^\\S+"),
      # Extract species: second word after genus
      Species = str_match(Organism_clean, "^\\S+\\s+(\\S+)")[, 2]
    ) %>%
    select(
      Protein_ID,
      Organism = Organism_clean,
      Genus,
      Species,
      Full_Taxonomy
    )

  parse_time <- difftime(Sys.time(), parse_start, units = "secs")
  total_time <- difftime(Sys.time(), start_time, units = "secs")

  cat(sprintf(" ↳ Parsing complete in %.1f sec\n", as.numeric(parse_time)))
  cat(sprintf(" ↳ Total Step 3 time: %.1f sec\n", as.numeric(total_time)))

  # Save to file
  write_tsv(parsed, taxonomy_file)
  cat(sprintf(" ↳ Saved to: %s\n", basename(taxonomy_file)))

  # Summary stats
  n_unique_genera <- n_distinct(parsed$Genus, na.rm = TRUE)
  n_unique_species <- n_distinct(parsed$Species, na.rm = TRUE)
  cat(sprintf(" ↳ Unique genera: %s\n", format(n_unique_genera, big.mark = ",")))
  cat(sprintf(" ↳ Unique species: %s\n", format(n_unique_species, big.mark = ",")))

  return(list(
    data = parsed,
    db_name = db_name,
    taxonomy_file = taxonomy_file
  ))
}

################################################################################
# STEP 4: QUERY FUNGUILD FOR EACH GENUS
################################################################################

query_funguild_for_genera <- function(taxonomy_result, funguild_db) {
  cat("\n=== Step 4: Querying FunGuild for Genera ===\n")

  taxonomy_data <- taxonomy_result$data
  db_name <- taxonomy_result$db_name

  guild_mapping_file <- file.path(DATA_PROCESSED_DIR,
                                  paste0(db_name, "_genus_guild_mapping_", timestamp, ".tsv"))

  # Get unique genera (excluding NA)
  unique_genera <- taxonomy_data %>%
    filter(!is.na(Genus), Genus != "") %>%
    pull(Genus) %>%
    unique() %>%
    sort()

  cat(sprintf("Found %s unique genera to query\n", format(length(unique_genera), big.mark = ",")))

  # Test query first
  cat("\nTesting with 'Amanita'...\n")
  test_result <- safe_funguild_query("Amanita", funguild_db)
  if (!is.null(test_result)) {
    cat(sprintf(" ↳ Test successful! Found %d records for Amanita\n", nrow(test_result)))
  } else {
    warning("Test query failed - proceeding anyway")
  }

  # Query all genera with progress bar
  cat("\nQuerying FunGuild database...\n")
  pb <- txtProgressBar(min = 0, max = length(unique_genera), style = 3)

  results_list <- list()
  success_count <- 0
  unmatched_genera <- character(0)  # Track genera without guild data

  for (i in seq_along(unique_genera)) {
    genus <- unique_genera[i]
    result <- safe_funguild_query(genus, funguild_db)

    if (!is.null(result)) {
      results_list[[length(results_list) + 1]] <- result
      success_count <- success_count + 1
    } else {
      unmatched_genera <- c(unmatched_genera, genus)
    }

    setTxtProgressBar(pb, i)
  }

  close(pb)

  cat(sprintf("\n\n ↳ Successfully queried %d/%d genera (%.1f%%)\n",
              success_count, length(unique_genera),
              100 * success_count / length(unique_genera)))
  cat(sprintf(" ↳ Genera without guild data: %d\n", length(unmatched_genera)))

  # Combine all results
  if (length(results_list) > 0) {
    guild_data <- bind_rows(results_list)
    cat(sprintf(" ↳ Retrieved %d total guild records\n", nrow(guild_data)))

    write_tsv(guild_data, guild_mapping_file)
    cat(sprintf(" ↳ Saved to: %s\n", guild_mapping_file))

    return(list(
      guild_data = guild_data,
      total_genera = length(unique_genera),
      matched_genera = success_count,
      unmatched_genera = unmatched_genera,
      guild_mapping_file = guild_mapping_file
    ))
  } else {
    stop("No successful queries - check database connection")
  }
}

################################################################################
# STEP 5: GENERATE DIAGNOSTIC REPORT AND PLOTS
################################################################################

generate_diagnostics <- function(taxonomy_result, guild_results) {
  cat("\n=== Step 5: Generating Diagnostic Report and Plots ===\n")

  taxonomy_data <- taxonomy_result$data
  db_name <- taxonomy_result$db_name
  guild_data <- guild_results$guild_data
  total_genera <- guild_results$total_genera
  matched_genera <- guild_results$matched_genera
  unmatched_genera <- guild_results$unmatched_genera

  report_file <- file.path(REPORTS_DIR, paste0(db_name, "_funguild_report_", timestamp, ".txt"))
  plot_file <- file.path(REPORTS_DIR, paste0(db_name, "_funguild_plots_", timestamp, ".png"))
  unmatched_log_file <- file.path(REPORTS_DIR, paste0(db_name, "_unmatched_genera_", timestamp, ".tsv"))

  # Genus statistics
  genus_counts <- taxonomy_data %>%
    filter(!is.na(Genus)) %>%
    count(Genus, name = "protein_count") %>%
    arrange(desc(protein_count))

  coverage_pct <- 100 * matched_genera / total_genera

  # Create unmatched genera log with protein counts
  cat(sprintf("Creating unmatched genera log: %s\n", basename(unmatched_log_file)))

  unmatched_with_counts <- genus_counts %>%
    filter(Genus %in% unmatched_genera) %>%
    arrange(desc(protein_count))

  write_tsv(unmatched_with_counts, unmatched_log_file)
  cat(sprintf(" ↳ Logged %d genera without guild data\n", nrow(unmatched_with_counts)))

  # Guild summary
  guild_summary <- guild_data %>%
    group_by(guild) %>%
    summarise(
      genus_count = n_distinct(query_genus),
      record_count = n(),
      avg_confidence = mean(
        case_when(
          confidenceRanking == "Highly Probable" ~ 3,
          confidenceRanking == "Probable" ~ 2,
          confidenceRanking == "Possible" ~ 1,
          TRUE ~ 0
        ),
        na.rm = TRUE
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(genus_count))

  # Trophic mode distribution
  trophic_summary <- guild_data %>%
    filter(!is.na(trophicMode), trophicMode != "") %>%
    mutate(trophicMode = str_trim(trophicMode)) %>%
    count(trophicMode, name = "count") %>%
    arrange(desc(count))

  # Growth form distribution
  growth_summary <- guild_data %>%
    filter(!is.na(growthForm), growthForm != "") %>%
    mutate(growthForm = str_trim(growthForm)) %>%
    count(growthForm, name = "count") %>%
    arrange(desc(count))

  # Write report
  cat(sprintf("Writing report to: %s\n", report_file))

  sink(report_file)

  cat("================================================================================\n")
  cat("FUNGUILD PRE-ANALYSIS DIAGNOSTIC REPORT (v1.0.1)\n")
  cat("================================================================================\n")
  cat(sprintf("Generated: %s\n", Sys.time()))
  cat(sprintf("Database: %s\n", toupper(db_name)))
  cat(sprintf("Timestamp: %s\n\n", timestamp))

  cat("--- INPUT DATA ---\n")
  cat(sprintf("Total proteins: %s\n", format(nrow(taxonomy_data), big.mark = ",")))
  cat(sprintf("Total unique genera: %s\n", format(total_genera, big.mark = ",")))
  cat(sprintf("Total unique species: %s\n",
              format(n_distinct(taxonomy_data$Species, na.rm = TRUE), big.mark = ",")))
  cat("\n")

  cat("--- FUNGUILD COVERAGE ---\n")
  cat(sprintf("Genera with guild data: %d/%d (%.1f%%)\n", matched_genera, total_genera, coverage_pct))
  cat(sprintf("Genera without guild data: %d (%.1f%%)\n",
              length(unmatched_genera), 100 - coverage_pct))
  cat(sprintf("Total guild records: %d\n", nrow(guild_data)))
  cat(sprintf("Average records per matched genus: %.1f\n", nrow(guild_data) / matched_genera))
  cat("\n")

  cat("--- TOP 20 GENERA BY PROTEIN COUNT ---\n")
  top_genera <- head(genus_counts, 20)
  for (i in 1:nrow(top_genera)) {
    cat(sprintf("%2d. %-30s %8s proteins\n", i, top_genera$Genus[i],
                format(top_genera$protein_count[i], big.mark = ",")))
  }
  cat("\n")

  cat("--- GENERA WITHOUT GUILD DATA (TOP 20 BY PROTEIN COUNT) ---\n")
  top_unmatched <- head(unmatched_with_counts, 20)
  if (nrow(top_unmatched) > 0) {
    for (i in 1:nrow(top_unmatched)) {
      cat(sprintf("%2d. %-30s %8s proteins\n", i, top_unmatched$Genus[i],
                  format(top_unmatched$protein_count[i], big.mark = ",")))
    }
  } else {
    cat("  (All genera have guild data)\n")
  }
  cat(sprintf("\nFull list saved to: %s\n", basename(unmatched_log_file)))
  cat("\n")

  cat("--- TOP 15 GUILDS BY GENUS COUNT ---\n")
  top_guilds <- head(guild_summary, 15)
  for (i in 1:nrow(top_guilds)) {
    cat(sprintf("%2d. %-40s %3d genera, %4d records\n",
                i, top_guilds$guild[i], top_guilds$genus_count[i], top_guilds$record_count[i]))
  }
  cat("\n")

  cat("--- TROPHIC MODE DISTRIBUTION ---\n")
  for (i in 1:min(nrow(trophic_summary), 15)) {
    cat(sprintf("%-20s %5d records\n", trophic_summary$trophicMode[i], trophic_summary$count[i]))
  }
  cat("\n")

  cat("--- GROWTH FORM DISTRIBUTION ---\n")
  for (i in 1:min(nrow(growth_summary), 15)) {
    cat(sprintf("%-20s %5d records\n", growth_summary$growthForm[i], growth_summary$count[i]))
  }
  cat("\n")

  cat("--- OUTPUT FILES ---\n")
  cat(sprintf("FunGuild DB: %s\n", file.path(DATA_CACHE_DIR, paste0("funguild_db_", timestamp, ".rds"))))
  cat(sprintf("Taxonomy: %s\n", taxonomy_result$taxonomy_file))
  cat(sprintf("Guild mapping: %s\n", guild_results$guild_mapping_file))
  cat(sprintf("Unmatched genera: %s\n", unmatched_log_file))
  cat(sprintf("Report: %s\n", report_file))
  cat(sprintf("Plots: %s\n", plot_file))
  cat("\n")
  cat("================================================================================\n")
  cat("END OF REPORT\n")
  cat("================================================================================\n")
  sink()
  cat(" ↳ Report written\n")

  # Create visualizations
  cat(sprintf("Creating plots: %s\n", plot_file))

  # Plot 1: Top 20 Genera by Protein Count
  p1 <- genus_counts %>%
    head(20) %>%
    mutate(Genus = fct_reorder(Genus, protein_count)) %>%
    ggplot(aes(x = protein_count, y = Genus)) +
    geom_col(fill = "#2C3E50") +
    geom_text(aes(label = format(protein_count, big.mark = ",")), hjust = -0.1, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2)), labels = scales::comma) +
    labs(
      title = "Top 20 Genera by Protein Count",
      subtitle = sprintf("From %s total proteins in %s",
                        format(nrow(taxonomy_data), big.mark = ","), toupper(db_name)),
      x = "Number of Proteins",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Plot 2: Top Guilds
  p2 <- guild_summary %>%
    head(15) %>%
    mutate(guild = str_wrap(guild, width = 30), guild = fct_reorder(guild, genus_count)) %>%
    ggplot(aes(x = genus_count, y = guild)) +
    geom_col(fill = "#3498db") +
    geom_text(aes(label = genus_count), hjust = -0.2, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Top 15 Guilds by Genus Count", x = "Number of Genera", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Plot 3: FunGuild Coverage
  coverage_data <- data.frame(
    Category = c("Matched", "Unmatched"),
    Count = c(matched_genera, total_genera - matched_genera),
    Percentage = c(coverage_pct, 100 - coverage_pct)
  )
  p3 <- ggplot(coverage_data, aes(x = "", y = Count, fill = Category)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    geom_text(
      aes(label = sprintf("%d\n(%.1f%%)", Count, Percentage)),
      position = position_stack(vjust = 0.5),
      size = 4, fontface = "bold"
    ) +
    scale_fill_manual(values = c("Matched" = "#27ae60", "Unmatched" = "#95a5a6")) +
    labs(
      title = "FunGuild Coverage",
      subtitle = sprintf("%d of %d genera", matched_genera, total_genera)
    ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )

  # Plot 4: Trophic Modes
  p4 <- trophic_summary %>%
    head(10) %>%
    mutate(trophicMode = str_wrap(trophicMode, width = 25),
           trophicMode = fct_reorder(trophicMode, count)) %>%
    ggplot(aes(x = count, y = trophicMode)) +
    geom_col(fill = "#e74c3c") +
    geom_text(aes(label = count), hjust = -0.2, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Trophic Mode Distribution", subtitle = "Top 10 modes",
         x = "Number of Records", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Plot 5: Growth Forms
  p5 <- growth_summary %>%
    head(10) %>%
    mutate(growthForm = str_wrap(growthForm, width = 25),
           growthForm = fct_reorder(growthForm, count)) %>%
    ggplot(aes(x = count, y = growthForm)) +
    geom_col(fill = "#9b59b6") +
    geom_text(aes(label = count), hjust = -0.2, size = 3) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "Growth Form Distribution", subtitle = "Top 10 forms",
         x = "Number of Records", y = NULL) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  # Combine plots
  combined_plot <- (p1 | p2) / (p3 | p4 | p5) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = sprintf("FunGuild Pre-Analysis: %s Database (v1.0.1)", toupper(db_name)),
      subtitle = sprintf("Generated: %s", Sys.time()),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
    )

  ggsave(plot_file, combined_plot, width = 16, height = 10, dpi = 300)
  cat(" ↳ Plots saved\n")

  return(list(
    total_proteins = nrow(taxonomy_data),
    total_genera = total_genera,
    matched_genera = matched_genera,
    unmatched_genera_count = length(unmatched_genera),
    coverage_pct = coverage_pct,
    total_guild_records = nrow(guild_data),
    report_file = report_file,
    plot_file = plot_file,
    unmatched_log_file = unmatched_log_file
  ))
}

################################################################################
# MAIN EXECUTION
################################################################################

main <- function() {
  cat("\n")
  cat("================================================================================\n")
  cat("FUNGUILD PRE-ANALYSIS PIPELINE (v1.0.1)\n")
  cat("================================================================================\n")
  cat("Prepares genus-guild mapping from complete UniProt fungi database\n")
  cat("See docs/changelog_data_preanalysis.md for version history.\n")
  cat("================================================================================\n")

  # Parse command-line arguments
  args <- commandArgs(trailingOnly = TRUE)
  db_choice <- NULL

  if (length(args) > 0) {
    if (any(args %in% c("--db", "-d"))) {
      db_idx <- which(args %in% c("--db", "-d"))
      if (length(args) > db_idx) {
        db_choice <- args[db_idx + 1]
        cat(sprintf("Command-line database selection: %s\n", db_choice))
      }
    } else if (args[1] %in% c("--help", "-h")) {
      cat("\nUsage:\n")
      cat("  Rscript R/data_preanalysis.R [options]\n\n")
      cat("Options:\n")
      cat("  --db, -d <choice>    Select database: 'sprot', 'trembl', or number\n")
      cat("  --help, -h           Show this help message\n\n")
      cat("Examples:\n")
      cat("  Rscript R/data_preanalysis.R --db sprot\n")
      cat("  Rscript R/data_preanalysis.R --db trembl\n")
      cat("  Rscript R/data_preanalysis.R --db 1\n\n")
      cat("See docs/changelog_data_preanalysis.md for version history.\n\n")
      return(invisible())
    }
  }

  start_time <- Sys.time()

  # Step 0: Create directories
  create_directories()

  # Step 1: Download FunGuild database
  funguild_result <- download_funguild_db()
  funguild_db <- funguild_result$db
  funguild_db_file <- funguild_result$file

  if (funguild_result$is_new) {
    saveRDS(funguild_db, funguild_db_file)
    cat(sprintf(" ↳ Saved FunGuild DB to: %s\n", funguild_db_file))
  }

  # Step 2 + 2.5: Locate UniProt fungi + convert to tabular
  uniprot_result <- select_or_download_uniprot(db_choice)

  # Step 3: Extract taxonomy from tabular format (ultra-fast)
  taxonomy_result <- extract_taxonomy_from_tabular(uniprot_result$tabular_file)

  # Step 4: Query FunGuild
  guild_results <- query_funguild_for_genera(taxonomy_result, funguild_db)

  # Step 5: Generate diagnostics
  summary_stats <- generate_diagnostics(taxonomy_result, guild_results)

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  # Final summary
  cat("\n")
  cat("================================================================================\n")
  cat("PIPELINE COMPLETE\n")
  cat("================================================================================\n")
  cat(sprintf("Total execution time: %.1f minutes\n\n", as.numeric(elapsed)))

  cat("SUMMARY:\n")
  cat(sprintf("  Total proteins analyzed:     %s\n",
              format(summary_stats$total_proteins, big.mark = ",")))
  cat(sprintf("  Total unique genera:         %s\n",
              format(summary_stats$total_genera, big.mark = ",")))
  cat(sprintf("  Genera with guild data:      %d (%.1f%%)\n",
              summary_stats$matched_genera, summary_stats$coverage_pct))
  cat(sprintf("  Genera without guild data:   %d (%.1f%%)\n",
              summary_stats$unmatched_genera_count, 100 - summary_stats$coverage_pct))
  cat(sprintf("  Total guild records:         %d\n", summary_stats$total_guild_records))
  cat("\n")

  if (uniprot_result$used_cache) {
    cat("CACHING: Used cached tabular file (fast re-run)\n\n")
  } else {
    cat("CACHING: Created new tabular file (will be reused on next run)\n\n")
  }

  cat("OUTPUT FILES:\n")
  cat(sprintf("  FunGuild database:   %s\n", funguild_db_file))
  cat(sprintf("  Tabular conversion:  %s\n", uniprot_result$tabular_file))
  cat(sprintf("  Taxonomy table:      %s\n", taxonomy_result$taxonomy_file))
  cat(sprintf("  Guild mapping:       %s\n", guild_results$guild_mapping_file))
  cat(sprintf("  Unmatched genera:    %s\n", summary_stats$unmatched_log_file))
  cat(sprintf("  Diagnostic report:   %s\n", summary_stats$report_file))
  cat(sprintf("  Diagnostic plots:    %s\n", summary_stats$plot_file))
  cat("\n")

  cat("OUTPUT FILE USAGE:\n")
  cat("  - Taxonomy table: Contains Protein_ID, Organism, Genus, Species, Full_Taxonomy\n")
  cat("  - Guild mapping: Contains ecological guild data for each genus\n")
  cat("  - Unmatched genera: List of genera without FunGuild data (sorted by protein count)\n")
  cat("  - Review diagnostic report for coverage analysis\n")
  cat("  - Check plots for taxonomic and ecological distributions\n")
  cat("\n")
  cat("================================================================================\n")
}

################################################################################
# SCRIPT ENTRY POINT
################################################################################

if (!interactive()) {
  main()
}
