# Changelog - Data Pre-Analysis

All notable changes to the FunGuild Pre-Analysis Pipeline are documented in this file.

## [v1.0.3] - 2025

### Features
- **NEW**: FASTA sequence export in Step 3
  - Generates `{db}_fungi_sequences_TIMESTAMP.fasta` in `data/processed/`
  - Header format: `>PROTEIN_ID|NCBI_TaxID|Organism name`
  - Uses existing tabular cache (no re-parsing of .dat files)
  - Ready for downstream tools (PLAAC, external predictors)

---

## [v1.0.2] - 2025

### Features
- **NEW**: NCBI_TaxID column added to tabular cache
  - Extracts taxonomy ID from OX lines in UniProt .dat files
  - Available for downstream PrionScan and other analyses

---

## [v1.0.1] - 2025

### Bug Fixes
- **FIXED**: Sequence extraction during .dat to TSV conversion
  - Previous regex `^\\s+[A-Z\\s]+$` was too strict and missed some sequence lines
  - New approach uses robust check for any line containing uppercase letters while in sequence section
  - Extracts only amino acid characters with `gsub("[^A-Z]", "", line)`

---

## [v1.0] - 2025

Complete rewrite with optimized architecture. Reset version numbering for new directory layout.

### Directory Structure
New 3-tier data organization:
- `data/raw/` - User-provided input files (.dat, .dat.gz)
- `data/cache/` - Intermediate/cached files (converted TSV, FunGuild DB)
- `data/processed/` - Final outputs for downstream scripts (taxonomy, guild mapping)

### Tabular Cache
The cached TSV includes all data needed for downstream analysis:
- `Protein_ID` - UniProt accession
- `Organism` - Full organism name from OS line
- `Full_Taxonomy` - Complete OC taxonomy string
- `NCBI_TaxID` - NCBI Taxonomy ID from OX line
- `Sequence` - Full amino acid sequence

### Features
- Single-pass .dat to TSV conversion (eliminates redundant file reads)
- Tabular caching for fast re-runs (14-day validity)
- Ultra-fast vectorized parsing (replaces parallel processing)
- Improved organism name cleaning: handles `[brackets]` and `'quotes`
- Unmatched genera log with protein counts for diagnostic purposes
- Simpler codebase (no parallel cluster management)
- Expected 4-5x faster first run, 15-25x faster re-runs

---

## Historical Versions (Pre-reset)

### v41.3
- **IMPROVED**: Organism name cleaning now handles:
  - Square brackets: `[Bisifusarium]` -> `Bisifusarium`
  - Single quotes: `'Aporospora` -> `Aporospora`
- Better genus/species extraction for edge cases in UniProt naming

### v41.2
- **NEW**: Unmatched genera log - generates TSV of genera without FunGuild data
- **NEW**: Diagnostic report includes top unmatched genera by protein count
- Helps identify coverage gaps for future database improvements

### v41.1
- **STANDALONE**: Removed references to external workflows/pipelines
- Self-contained taxonomic and ecological guild analysis tool

### v41
- **NEW**: Step 2.5 - Single-pass .dat to TSV conversion (eliminates redundant reads)
- **NEW**: Tabular caching - reuses converted TSV files (14-day validity)
- **NEW**: Ultra-fast Step 3 - vectorized parsing of TSV (no parallel overhead)
- **REMOVED**: Parallel file reading (replaced with single-pass approach)
- **SPEED**: Expected 4-5x faster on first run, 15-25x on subsequent runs

### Previous Features (retained through versions)
- Extract full OC taxonomy string for each entry
- FunGuild database caching (14-day validity)
- Windows and Linux compatible
