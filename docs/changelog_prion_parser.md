# Changelog - Prion Domain Prediction

All notable changes to the Prion Domain Prediction script are documented in this file.

## [v1.0.1] - 2025

### Features
- **NEW**: NCBI_TaxID included in prediction output
  - Reads from tabular cache NCBI_TaxID column
  - Included in tabular output format
  - Streaming fallback also extracts OX lines from .dat files

### Bug Fixes
- **FIXED**: Sequence extraction in streaming fallback uses improved regex (matches data_preanalysis.R fix)

---

## [v1.0] - 2025

Complete rewrite leveraging cached tabular files from data_preanalysis.R. Reset version numbering.

### Major Changes
- **Cache-First Architecture**: Reads pre-converted TSV from `data/cache/` instead of streaming .dat files
- **Simplified Parallelism**: Workers process data frame chunks directly (no file seeking)
- **Fallback Mode**: Retains streaming parser for direct .dat file processing when cache unavailable
- **Unified Directory Structure**: Uses same 3-tier layout as data_preanalysis.R

### Performance
- ~10x faster when using cached tabular files
- Parallel processing now operates on in-memory data chunks
- No redundant file reads across workers

### Output
- Predictions saved to `data/processed/` (aligned with pipeline structure)
- Both original (.dat) and tabular (.tsv) formats supported

---

## Historical Versions (Pre-reset)

### v13 (Streaming Architecture)
- Parallel streaming from .dat files
- Each worker reads full file, seeking to assigned entry range
- Memory-efficient for 100GB+ files
- Extracted ID, OS, OC, OX, and sequence per entry

### Original Perl Version
- Original implementation by Vladimir Espinosa Angarica (2012)
- Published in Angarica et al. (2013) BMC Genomics
- Algorithm based on Alberti et al. (2009) Cell 137(1):146-158

## Algorithm

The prion domain prediction uses:
1. **Sliding Window Analysis**: 60 amino acid window scans entire sequence
2. **Amino Acid Scoring**: Each residue has prion-propensity score
3. **Proline Pair Correction**: Log-odds adjustment for proline distances
4. **Cutoff Filtering**: Only predictions scoring >= 50 are reported

### Amino Acid Scores
Residues with highest prion propensity: N (2.511), Q (2.044), Y (0.786), S (0.733)
Residues with lowest prion propensity: C (-3.807), W (-3.459), E (-2.766)
