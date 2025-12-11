# Changelog - Venn Predictions

All notable changes to the Venn Predictions comparison tool are documented in this file.

## [v1.0] - 2025

### Initial Release

First release of `venn_predictions.R` - a tool to compare and visualize overlap between PrionScan and PAPA predictions.

### Features
- **Venn diagram generation**: Creates publication-ready 2-set Venn diagrams using `ggVennDiagram`
- **Auto-detection**: Use `--db sprot` or `--db trembl` to automatically find latest prediction files
- **Manual input**: Specify files explicitly with `--prionscan` and `--papa` arguments
- **Set analysis**: Computes and reports:
  - PrionScan-only predictions
  - PAPA-only predictions
  - Predictions found by both methods
  - Total unique predictions
- **Console summary**: Detailed statistics with counts and percentages
- **PNG output**: High-resolution (300 DPI) diagram saved to `reports/`

### Usage

```bash
# Auto-detect latest files for a database
Rscript R/venn_predictions.R --db sprot

# Specify files explicitly
Rscript R/venn_predictions.R \
  --prionscan data/processed/sprot_prion_predictions_*.tsv \
  --papa data/processed/sprot_papa_predictions_*.tsv

# Custom output filename
Rscript R/venn_predictions.R --db sprot --output my_venn.png
```

### Dependencies
- `tidyverse`
- `ggVennDiagram`
- `optparse`

### Output
- PNG file: `reports/{db}_venn_predictions_TIMESTAMP.png`
- Console: Summary statistics

---

## Notes

This script compares protein IDs between:
- **PrionScan**: All proteins in output file (already filtered at Score >= 50)
- **PAPA**: Proteins where `Above_Threshold == TRUE` (Score >= 0.05)

The overlap indicates proteins predicted to have prion-like domains by both algorithms, representing higher-confidence candidates.
