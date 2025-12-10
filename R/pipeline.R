#!/usr/bin/env Rscript
# =============================================================================
# Pipeline Orchestrator for Prion Prediction Pipeline
# Version: 1.1
# =============================================================================
#
# Usage:
#   Rscript R/pipeline.R --db sprot              # Run full pipeline
#   Rscript R/pipeline.R --db sprot --step prion_parser  # Run single step
#   Rscript R/pipeline.R --list                  # Show available steps
#   Rscript R/pipeline.R --check --db sprot      # Validate without running
#   Rscript R/pipeline.R --diagram               # Export diagram to docs/
#   Rscript R/pipeline.R --help                  # Show help
#
# =============================================================================
# WORKFLOW DIAGRAM (Mermaid)
# =============================================================================
#
# ```mermaid
# flowchart TD
#     subgraph inputs["Raw Inputs"]
#         DAT[("UniProt .dat files<br/>data/raw/")]
#     end
#
#     subgraph step1["Step 1: Data Pre-Analysis"]
#         S1[data_preanalysis.R]
#     end
#
#     subgraph cache["Cache Layer"]
#         TSV[("Tabular Cache<br/>data/cache/")]
#         FGDB[("FunGuild DB<br/>data/cache/")]
#     end
#
#     subgraph step2["Step 2: Prion Prediction"]
#         S2[prion_parser.R]
#     end
#
#     subgraph outputs["Processed Outputs"]
#         TAX[("Taxonomy TSV<br/>data/processed/")]
#         GUILD[("Guild Mapping<br/>data/processed/")]
#         PRION[("Prion Predictions<br/>data/processed/")]
#     end
#
#     DAT --> S1
#     S1 --> TSV
#     S1 --> FGDB
#     S1 --> TAX
#     S1 --> GUILD
#     TSV --> S2
#     S2 --> PRION
# ```
#
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

SCRIPT_VERSION <- "1.1"

# Directory structure
DATA_RAW_DIR      <- "data/raw"
DATA_CACHE_DIR    <- "data/cache"
DATA_PROCESSED_DIR <- "data/processed"
REPORTS_DIR       <- "reports"
DOCS_DIR          <- "docs"

# Database name mapping
DB_NAMES <- list(
  sprot  = "uniprot_sprot_fungi",
  trembl = "uniprot_trembl_fungi"
)

# =============================================================================
# PIPELINE REGISTRY
# =============================================================================
# Each step defines:
#   - script: Path to R script
#   - description: Human-readable description
#   - inputs: Named list of input file patterns (with {db} placeholder)
#   - outputs: Named list of output file patterns (with {db} placeholder)
#   - depends_on: Vector of step names this step depends on (NULL if none)
#   - args: Function to generate CLI arguments for the script

PIPELINE_STEPS <- list(

  data_preanalysis = list(
    script = "R/data_preanalysis.R",
    description = "Convert .dat to TSV, extract taxonomy, map FunGuild",
    inputs = list(
      raw = "{DATA_RAW_DIR}/{db_name}.dat(.gz)?"
    ),
    outputs = list(
      cache    = "{DATA_CACHE_DIR}/{db_prefix}_fungi_tabular_*.tsv",
      taxonomy = "{DATA_PROCESSED_DIR}/{db_prefix}_fungi_taxonomy_*.tsv",
      guild    = "{DATA_PROCESSED_DIR}/{db_prefix}_genus_guild_mapping_*.tsv"
    ),
    depends_on = NULL,
    args = function(db) c("--db", db)
  ),

  prion_parser = list(
    script = "R/prion_parser.R",
    description = "Predict prion domains from cached sequences",
    inputs = list(
      cache = "{DATA_CACHE_DIR}/{db_prefix}_fungi_tabular_*.tsv"
    ),
    outputs = list(
      predictions = "{DATA_PROCESSED_DIR}/{db_prefix}_prion_predictions_*.tsv"
    ),
    depends_on = c("data_preanalysis"),
    args = function(db) c("--db", db)
  )

)

# Step execution order
STEP_ORDER <- c("data_preanalysis", "prion_parser")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Get database prefix from db type
get_db_prefix <- function(db) {

  switch(db,
    sprot  = "sprot",
    trembl = "trembl",
    stop("Unknown database: ", db)
  )
}

#' Get full database name from db type
get_db_name <- function(db) {
  DB_NAMES[[db]]
}

#' Expand pattern with database-specific values
expand_pattern <- function(pattern, db) {
  db_prefix <- get_db_prefix(db)
  db_name <- get_db_name(db)

  pattern <- gsub("\\{DATA_RAW_DIR\\}", DATA_RAW_DIR, pattern)
  pattern <- gsub("\\{DATA_CACHE_DIR\\}", DATA_CACHE_DIR, pattern)
  pattern <- gsub("\\{DATA_PROCESSED_DIR\\}", DATA_PROCESSED_DIR, pattern)
  pattern <- gsub("\\{db_prefix\\}", db_prefix, pattern)
  pattern <- gsub("\\{db_name\\}", db_name, pattern)

  pattern
}

#' Find files matching a pattern (supports regex and glob-like patterns)
find_files <- function(pattern, dir = ".") {
  # Convert pattern to regex if it contains glob wildcards
  if (grepl("\\*", pattern)) {
    # Extract directory and file pattern
    if (grepl("/", pattern)) {
      dir_part <- dirname(pattern)
      file_part <- basename(pattern)
    } else {
      dir_part <- dir
      file_part <- pattern
    }

    # Convert glob to regex
    regex_pattern <- gsub("\\.", "\\\\.", file_part)
    regex_pattern <- gsub("\\*", ".*", regex_pattern)
    regex_pattern <- paste0("^", regex_pattern, "$")

    if (dir.exists(dir_part)) {
      files <- list.files(dir_part, pattern = regex_pattern, full.names = TRUE)
      return(files)
    } else {
      return(character(0))
    }
  } else {
    # Direct file path - check if exists (with optional .gz)
    if (file.exists(pattern)) {
      return(pattern)
    } else if (file.exists(paste0(pattern, ".gz"))) {
      return(paste0(pattern, ".gz"))
    } else if (grepl("\\.gz\\)\\?$", pattern)) {
      # Pattern ends with (.gz)? - try both
      base_pattern <- sub("\\(\\.gz\\)\\?$", "", pattern)
      if (file.exists(base_pattern)) return(base_pattern)
      if (file.exists(paste0(base_pattern, ".gz"))) return(paste0(base_pattern, ".gz"))
    }
    return(character(0))
  }
}

#' Find the most recent file matching a pattern
find_latest_file <- function(pattern, db) {
  expanded <- expand_pattern(pattern, db)
  files <- find_files(expanded)

  if (length(files) == 0) return(NULL)

  # Sort by modification time (most recent first)
  mtimes <- file.mtime(files)
  files[order(mtimes, decreasing = TRUE)][1]
}

#' Check if inputs exist for a step
check_inputs <- function(step_name, db) {
  step <- PIPELINE_STEPS[[step_name]]
  results <- list()

  for (input_name in names(step$inputs)) {
    pattern <- step$inputs[[input_name]]
    file <- find_latest_file(pattern, db)
    results[[input_name]] <- list(
      pattern = expand_pattern(pattern, db),
      found = !is.null(file),
      file = file
    )
  }

  results
}

#' Check if outputs exist for a step
check_outputs <- function(step_name, db) {
  step <- PIPELINE_STEPS[[step_name]]
  results <- list()

  for (output_name in names(step$outputs)) {
    pattern <- step$outputs[[output_name]]
    file <- find_latest_file(pattern, db)
    results[[output_name]] <- list(
      pattern = expand_pattern(pattern, db),
      found = !is.null(file),
      file = file
    )
  }

  results
}

#' Check if dependencies are satisfied
check_dependencies <- function(step_name, db) {
  step <- PIPELINE_STEPS[[step_name]]

  if (is.null(step$depends_on)) {
    return(list(satisfied = TRUE, missing = character(0)))
  }

  missing <- character(0)
  for (dep in step$depends_on) {
    outputs <- check_outputs(dep, db)
    all_found <- all(sapply(outputs, function(x) x$found))
    if (!all_found) {
      missing <- c(missing, dep)
    }
  }

  list(satisfied = length(missing) == 0, missing = missing)
}

#' Log message with timestamp
log_msg <- function(..., log_file = NULL) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  if (!is.null(log_file)) {
    cat(msg, "\n", file = log_file, append = TRUE)
  }
}

#' Run a single pipeline step
run_step <- function(step_name, db, log_file = NULL) {
  step <- PIPELINE_STEPS[[step_name]]

  log_msg("========================================", log_file = log_file)
  log_msg("Running step: ", step_name, log_file = log_file)
  log_msg("Script: ", step$script, log_file = log_file)
  log_msg("========================================", log_file = log_file)

  # Check dependencies
  dep_check <- check_dependencies(step_name, db)
  if (!dep_check$satisfied) {
    log_msg("ERROR: Missing dependencies: ", paste(dep_check$missing, collapse = ", "), log_file = log_file)
    return(FALSE)
  }

  # Check inputs
  input_check <- check_inputs(step_name, db)
  all_inputs_found <- all(sapply(input_check, function(x) x$found))

  if (!all_inputs_found) {
    missing_inputs <- names(input_check)[!sapply(input_check, function(x) x$found)]
    log_msg("ERROR: Missing inputs: ", paste(missing_inputs, collapse = ", "), log_file = log_file)
    for (name in missing_inputs) {
      log_msg("  - ", name, ": ", input_check[[name]]$pattern, log_file = log_file)
    }
    return(FALSE)
  }

  # Build command
  script_path <- step$script
  args <- step$args(db)

  log_msg("Arguments: ", paste(args, collapse = " "), log_file = log_file)
  log_msg("", log_file = log_file)

  # Execute
  start_time <- Sys.time()
  result <- system2("Rscript", c(script_path, args), stdout = "", stderr = "")
  end_time <- Sys.time()

  elapsed <- difftime(end_time, start_time, units = "mins")

  if (result == 0) {
    log_msg("Step completed successfully in ", round(elapsed, 2), " minutes", log_file = log_file)
    return(TRUE)
  } else {
    log_msg("ERROR: Step failed with exit code ", result, log_file = log_file)
    return(FALSE)
  }
}

#' Print pipeline status
print_status <- function(db = NULL) {
  cat("\n")
  cat("=== Pipeline Status ===\n")
  cat("\n")

  for (step_name in STEP_ORDER) {
    step <- PIPELINE_STEPS[[step_name]]
    cat(sprintf("%-20s %s\n", step_name, step$description))
    cat(sprintf("%-20s Script: %s\n", "", step$script))

    if (!is.null(step$depends_on)) {
      cat(sprintf("%-20s Depends on: %s\n", "", paste(step$depends_on, collapse = ", ")))
    }

    if (!is.null(db)) {
      # Show input/output status
      inputs <- check_inputs(step_name, db)
      outputs <- check_outputs(step_name, db)

      cat(sprintf("%-20s Inputs:\n", ""))
      for (name in names(inputs)) {
        status <- if (inputs[[name]]$found) "[OK]" else "[MISSING]"
        cat(sprintf("%-20s   %s %s\n", "", status, name))
      }

      cat(sprintf("%-20s Outputs:\n", ""))
      for (name in names(outputs)) {
        status <- if (outputs[[name]]$found) "[OK]" else "[--]"
        cat(sprintf("%-20s   %s %s\n", "", status, name))
      }
    }

    cat("\n")
  }
}

#' Get the Mermaid diagram code (without markdown fences)
get_diagram_code <- function() {
'flowchart TD
    subgraph inputs["Raw Inputs"]
        DAT[("UniProt .dat files<br/>data/raw/")]
    end

    subgraph step1["Step 1: Data Pre-Analysis"]
        S1[data_preanalysis.R]
    end

    subgraph cache["Cache Layer"]
        TSV[("Tabular Cache<br/>data/cache/")]
        FGDB[("FunGuild DB<br/>data/cache/")]
    end

    subgraph step2["Step 2: Prion Prediction"]
        S2[prion_parser.R]
    end

    subgraph outputs["Processed Outputs"]
        TAX[("Taxonomy TSV<br/>data/processed/")]
        GUILD[("Guild Mapping<br/>data/processed/")]
        PRION[("Prion Predictions<br/>data/processed/")]
    end

    DAT --> S1
    S1 --> TSV
    S1 --> FGDB
    S1 --> TAX
    S1 --> GUILD
    TSV --> S2
    S2 --> PRION'
}

#' Check if internet is available by testing mermaid.ink
check_internet <- function() {
  tryCatch({
    con <- url("https://mermaid.ink", "r", timeout = 5)
    close(con)
    TRUE
  }, error = function(e) FALSE)
}

#' Generate PNG using mermaid.ink API
generate_png_api <- function(diagram_code, output_file) {
  # Base64 encode the diagram
  encoded <- base64enc::base64encode(charToRaw(diagram_code))

  # URL-safe base64 (replace + with -, / with _)
  encoded <- gsub("\\+", "-", encoded)
  encoded <- gsub("/", "_", encoded)

  # Fetch PNG from mermaid.ink

  url <- paste0("https://mermaid.ink/img/", encoded, "?type=png")

  tryCatch({
    download.file(url, output_file, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) {
    message("  API error: ", e$message)
    FALSE
  })
}

#' Generate PNG using mmdc CLI
generate_png_cli <- function(input_file, output_file) {
  mmdc_path <- Sys.which("mmdc")

  if (mmdc_path == "") {
    return(FALSE)
  }

  result <- system2(mmdc_path, c("-i", input_file, "-o", output_file),
                    stdout = FALSE, stderr = FALSE)
  result == 0
}

#' Export diagram to docs/ (Mermaid code + PNG)
print_diagram <- function() {
  # Ensure docs directory exists
  if (!dir.exists(DOCS_DIR)) {
    dir.create(DOCS_DIR, recursive = TRUE)
  }

  diagram_code <- get_diagram_code()
  md_file <- file.path(DOCS_DIR, "pipeline_diagram.md")
  png_file <- file.path(DOCS_DIR, "pipeline_diagram.png")

  # Write Mermaid markdown file
  md_content <- paste0(
    "# Pipeline Workflow Diagram\n\n",
    "```mermaid\n",
    diagram_code,
    "\n```\n"
  )
  writeLines(md_content, md_file)
  cat("Mermaid code saved to: ", md_file, "\n", sep = "")

  # Also print to console
  cat("\n", md_content, "\n", sep = "")

  # Generate PNG
  cat("Generating PNG...\n")

  # Check for base64enc package (needed for API method)
  has_base64enc <- requireNamespace("base64enc", quietly = TRUE)

  # Try API first if internet available and base64enc exists
  png_generated <- FALSE

  if (has_base64enc && check_internet()) {
    cat("  Using mermaid.ink API...\n")
    png_generated <- generate_png_api(diagram_code, png_file)
    if (png_generated) {
      cat("PNG saved to: ", png_file, "\n", sep = "")
    }
  }

  # Fall back to CLI if API failed or unavailable

if (!png_generated) {
    if (!has_base64enc) {
      cat("  Note: base64enc package not installed (needed for API method)\n")
    }
    cat("  Trying mmdc CLI...\n")

    # mmdc needs a .mmd file, not .md with fences
    mmd_file <- file.path(DOCS_DIR, "pipeline_diagram.mmd")
    writeLines(diagram_code, mmd_file)

    png_generated <- generate_png_cli(mmd_file, png_file)

    if (png_generated) {
      cat("PNG saved to: ", png_file, "\n", sep = "")
      # Clean up temp .mmd file
      unlink(mmd_file)
    } else {
      cat("  mmdc not found. Install with: npm install -g @mermaid-js/mermaid-cli\n")
      cat("  PNG generation skipped.\n")
      # Keep .mmd file for manual conversion
      cat("  Raw diagram saved to: ", mmd_file, " (for manual conversion)\n", sep = "")
    }
  }
}

#' Show help message
show_help <- function() {
  cat("\n")
  cat("Prion Prediction Pipeline Orchestrator v", SCRIPT_VERSION, "\n", sep = "")
  cat("====================================================\n")
  cat("\n")
  cat("Usage:\n")
  cat("  Rscript R/pipeline.R --db <sprot|trembl> [options]\n")
  cat("\n")
  cat("Options:\n")
  cat("  -d, --db <TYPE>      Database to process: sprot or trembl (required for run)\n")
  cat("  -s, --step <NAME>    Run only specified step (default: run all)\n")
  cat("  -f, --from <NAME>    Run from specified step to end\n")
  cat("  -l, --list           List all steps and their status\n")
  cat("  -c, --check          Validate inputs/outputs without running\n")
  cat("  -m, --diagram        Export workflow diagram to docs/ (Mermaid + PNG)\n")
  cat("  -h, --help           Show this help message\n")
  cat("\n")
  cat("Examples:\n")
  cat("  Rscript R/pipeline.R --db sprot              # Run full pipeline\n")
  cat("  Rscript R/pipeline.R --db sprot --step prion_parser  # Single step\n")
  cat("  Rscript R/pipeline.R --list --db sprot       # Show status\n")
  cat("  Rscript R/pipeline.R --check --db sprot      # Dry run\n")
  cat("  Rscript R/pipeline.R --diagram               # Export diagram to docs/\n")
  cat("\n")
}

# =============================================================================
# CLI ARGUMENT PARSING
# =============================================================================

option_list <- list(
  make_option(c("-d", "--db"), type = "character", default = NULL,
              help = "Database type: sprot or trembl"),
  make_option(c("-s", "--step"), type = "character", default = NULL,
              help = "Run only specified step"),
  make_option(c("-f", "--from"), type = "character", default = NULL,
              help = "Run from specified step to end"),
  make_option(c("-l", "--list"), action = "store_true", default = FALSE,
              help = "List all steps"),
  make_option(c("-c", "--check"), action = "store_true", default = FALSE,
              help = "Check inputs/outputs without running"),
  make_option(c("-m", "--diagram"), action = "store_true", default = FALSE,
              help = "Export workflow diagram to docs/ (Mermaid + PNG)"),
  make_option(c("-h", "--help"), action = "store_true", default = FALSE,
              help = "Show help")
)

parser <- OptionParser(option_list = option_list, add_help_option = FALSE)
args <- parse_args(parser)

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Handle --help
if (args$help) {
  show_help()
  quit(status = 0)
}

# Handle --diagram
if (args$diagram) {
  print_diagram()
  quit(status = 0)
}

# Handle --list
if (args$list) {
  print_status(args$db)
  quit(status = 0)
}

# Validate --db for running/checking
if (is.null(args$db) && !args$list && !args$diagram) {
  cat("ERROR: --db is required. Use --db sprot or --db trembl\n")
  cat("Use --help for more information.\n")
  quit(status = 1)
}

if (!is.null(args$db) && !args$db %in% c("sprot", "trembl")) {
  cat("ERROR: --db must be 'sprot' or 'trembl'\n")
  quit(status = 1)
}

# Determine which steps to run
steps_to_run <- STEP_ORDER

if (!is.null(args$step)) {
  if (!args$step %in% names(PIPELINE_STEPS)) {
    cat("ERROR: Unknown step '", args$step, "'. Available steps: ",
        paste(names(PIPELINE_STEPS), collapse = ", "), "\n", sep = "")
    quit(status = 1)
  }
  steps_to_run <- args$step
}

if (!is.null(args$from)) {
  if (!args$from %in% names(PIPELINE_STEPS)) {
    cat("ERROR: Unknown step '", args$from, "'. Available steps: ",
        paste(names(PIPELINE_STEPS), collapse = ", "), "\n", sep = "")
    quit(status = 1)
  }
  from_idx <- which(STEP_ORDER == args$from)
  steps_to_run <- STEP_ORDER[from_idx:length(STEP_ORDER)]
}

# Handle --check (dry run)
if (args$check) {
  cat("\n=== Pipeline Check (", args$db, ") ===\n\n", sep = "")

  all_ok <- TRUE
  for (step_name in steps_to_run) {
    step <- PIPELINE_STEPS[[step_name]]
    cat("Step: ", step_name, "\n", sep = "")

    # Check dependencies
    dep_check <- check_dependencies(step_name, args$db)
    if (!dep_check$satisfied) {
      cat("  [BLOCKED] Missing dependencies: ", paste(dep_check$missing, collapse = ", "), "\n")
      all_ok <- FALSE
      next
    }

    # Check inputs
    inputs <- check_inputs(step_name, args$db)
    inputs_ok <- all(sapply(inputs, function(x) x$found))

    if (inputs_ok) {
      cat("  [READY] All inputs available\n")
    } else {
      cat("  [MISSING] Some inputs missing:\n")
      for (name in names(inputs)) {
        if (!inputs[[name]]$found) {
          cat("    - ", name, ": ", inputs[[name]]$pattern, "\n", sep = "")
        }
      }
      all_ok <- FALSE
    }

    # Check outputs
    outputs <- check_outputs(step_name, args$db)
    outputs_exist <- all(sapply(outputs, function(x) x$found))

    if (outputs_exist) {
      cat("  [EXISTS] Outputs already exist (will be regenerated)\n")
    }

    cat("\n")
  }

  if (all_ok) {
    cat("Pipeline ready to run.\n")
  } else {
    cat("Pipeline has missing inputs or dependencies.\n")
  }

  quit(status = 0)
}

# =============================================================================
# RUN PIPELINE
# =============================================================================

# Ensure reports directory exists
if (!dir.exists(REPORTS_DIR)) {
  dir.create(REPORTS_DIR, recursive = TRUE)
}

# Create log file
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file <- file.path(REPORTS_DIR, paste0("pipeline_run_", args$db, "_", timestamp, ".log"))

log_msg("Pipeline Orchestrator v", SCRIPT_VERSION, log_file = log_file)
log_msg("Database: ", args$db, log_file = log_file)
log_msg("Steps to run: ", paste(steps_to_run, collapse = " -> "), log_file = log_file)
log_msg("Log file: ", log_file, log_file = log_file)
log_msg("", log_file = log_file)

# Run each step
pipeline_start <- Sys.time()
success <- TRUE

for (step_name in steps_to_run) {
  step_success <- run_step(step_name, args$db, log_file)

  if (!step_success) {
    log_msg("Pipeline stopped due to error in step: ", step_name, log_file = log_file)
    success <- FALSE
    break
  }

  log_msg("", log_file = log_file)
}

pipeline_end <- Sys.time()
total_time <- difftime(pipeline_end, pipeline_start, units = "mins")

log_msg("========================================", log_file = log_file)
if (success) {
  log_msg("Pipeline completed successfully!", log_file = log_file)
} else {
  log_msg("Pipeline failed.", log_file = log_file)
}
log_msg("Total time: ", round(total_time, 2), " minutes", log_file = log_file)
log_msg("Log saved to: ", log_file, log_file = log_file)
log_msg("========================================", log_file = log_file)

quit(status = if (success) 0 else 1)
