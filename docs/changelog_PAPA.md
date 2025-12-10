# Changelog - PAPA

All notable changes to PAPA (Prion Aggregation Prediction Algorithm) are documented in this file.

## [v1.1] - 2025

### Output Simplification
- **REMOVED**: Taxonomy and NCBI_TaxID columns from prediction output
  - These fields are redundant with tabular cache
  - Join with `*_fungi_tabular_*.tsv` by Protein_ID when metadata needed
  - Reduces output file size by ~30%
- **Output columns**: Protein_ID, Organism, PAPA_Score, PAPA_Position, Above_Threshold

---

## [v1.0] - 2025

Initial pipeline-integrated version. Major refactor from standalone script.

### Features

- **Cache-First Architecture**: Reads pre-converted TSV from `data/cache/` instead of streaming .dat files
- **Parallel Processing**: Uses `mclapply` with configurable `--cores` option
- **Unified CLI**: Matches PrionScan interface (`--db sprot|trembl`)
- **Auto-Discovery**: Automatically finds most recent tabular cache file
- **Fallback Mode**: Retains streaming parser via `--input` option when cache unavailable
- **Metadata Output**: Includes Organism, Taxonomy, NCBI_TaxID in predictions

### CLI Options

- `--db <sprot|trembl>`: Database selection (auto-discovers cache)
- `--input <FILE>`: Direct .dat file input (fallback)
- `--output <FILE>`: Custom output path (default: auto-generated)
- `--cores <N>`: Number of CPU cores
- `--window_size <N>`: PAPA window size (default: 41)
- `--threshold <N>`: Score threshold (default: 0.05)
- `--ignore_fold_index`: Skip FoldIndex filtering

### Output Format

TSV file with columns:
- `Protein_ID` - UniProt accession
- `Organism` - Cleaned organism name
- `Taxonomy` - Full OC taxonomy string
- `NCBI_TaxID` - NCBI Taxonomy ID
- `PAPA_Score` - Maximum PAPA score
- `PAPA_Position` - Position of max score (0-indexed)
- `Above_Threshold` - Boolean: score >= threshold

### Directory Structure

- Output: `data/processed/{db}_papa_predictions_YYYYMMDD_HHMMSS.tsv`
- Input: `data/cache/{db}_fungi_tabular_*.tsv`

---

## Algorithm

PAPA predicts prion-forming propensity using:

1. **Sliding Window Analysis**: Configurable window (default 41 aa)
2. **Prion Propensity Scores**: Amino acid composition-based scoring
3. **FoldIndex Filtering**: Masks folded regions (optional, enabled by default)
4. **Super-Window Averaging**: Smooths scores across positions

### References

- Toombs JA, McCarty BR, Ross ED. (2010) "Compositional determinants of prion formation in yeast." Mol Cell Biol 30(1):319-332
- Original Python implementation: http://www.prion.bcm.tmc.edu/

### Threshold Interpretation

| Score | Interpretation |
|-------|----------------|
| >= 0.05 | Potential prion-forming propensity |
| >= 0.10 | High confidence candidate |
| < 0.05 | Unlikely prion-forming |
