<!-- badges: start -->
[![R](https://github.com/sghignone/prion-prediction-pipeline/actions/workflows/r.yml/badge.svg)](https://github.com/sghignone/prion-prediction-pipeline/actions/workflows/r.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

# Prion Prediction Pipeline

A comprehensive R-based pipeline for identifying prion-like domains in fungal proteomes. This pipeline analyzes UniProt fungi databases to predict proteins with prion-forming propensity using multiple complementary algorithms.

## Overview

The pipeline processes UniProt SwissProt/TrEMBL fungi databases through a series of analysis steps:

1. **Data Pre-Analysis** - Converts UniProt .dat files to tabular format, extracts taxonomy, and maps ecological guilds via FunGuild
2. **PrionScan** - Predicts Q/N-rich prion domains using probabilistic scoring (Espinosa Angarica algorithm)
3. **PAPA** - Prion Aggregation Prediction Algorithm with FoldIndex filtering (Toombs et al.)
4. **PLAAC** - Prion-Like Amino Acid Composition analysis using log-likelihood scoring (Lancaster et al.)

## Requirements

### R Packages

```r
install.packages(c("tidyverse", "FUNGuildR", "patchwork", "optparse", "parallel", "readr", "stringr", "ggVennDiagram"))
```

### Java (for PLAAC)

```bash
# Option 1: System Java (JRE 8+)
java -version  # Verify installation

# Option 2: Conda environment
conda create -n plaac-java openjdk=11 -c conda-forge
# The script can auto-setup this with: Rscript R/plaac.R --setup-java
```

Git is also required for initial PLAAC repository clone.

### Optional (for diagram generation)

```bash
# For PNG export via mermaid.ink API
install.packages("base64enc")

# OR local CLI (no internet required)
npm install -g @mermaid-js/mermaid-cli
```

## Directory Structure

```
prion-prediction-pipeline/
├── R/                          # Analysis scripts
│   ├── pipeline.R              # Orchestrator (runs all steps)
│   ├── data_preanalysis.R      # Step 1: Data conversion & FunGuild
│   ├── PrionScan.R             # Step 2: Prion domain prediction
│   ├── PAPA.R                  # Step 2: Aggregation propensity
│   ├── plaac.R                 # Step 2: PLAAC prion-like composition
│   └── venn_predictions.R      # Compare all three methods
├── data/
│   ├── raw/                    # Input: UniProt .dat files
│   ├── cache/                  # Intermediate: Tabular cache, FunGuild DB
│   └── processed/              # Output: Predictions, taxonomy, guild mapping
├── docs/                       # Documentation & changelogs
│   ├── workflow.md             # Full workflow description
│   └── changelog_*.md          # Version history per script
└── reports/                    # Pipeline run logs
```

## Quick Start

### 1. Download Input Data

Place UniProt fungi proteome files in `data/raw/`:

```bash
# Swiss-Prot (~50MB, curated)
wget -P data/raw/ https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_fungi.dat.gz

# TrEMBL (~73GB, comprehensive)
wget -P data/raw/ https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_trembl_fungi.dat.gz
```

### 2. Run Full Pipeline

```bash
# Process Swiss-Prot database (recommended for testing)
Rscript R/pipeline.R --db sprot

# Process TrEMBL database (large, requires significant time/memory)
Rscript R/pipeline.R --db trembl
```

### 3. Run Individual Steps

```bash
# Step 1: Data pre-analysis only
Rscript R/data_preanalysis.R --db sprot

# Step 2a: PrionScan only (requires Step 1 cache)
Rscript R/PrionScan.R --db sprot

# Step 2b: PAPA only (requires Step 1 cache)
Rscript R/PAPA.R --db sprot

# Step 2c: PLAAC only (requires Step 1 FASTA output)
Rscript R/plaac.R -i data/processed/sprot_fungi_sequences_*.fasta --db sprot

# Compare predictions (auto-detects all available methods)
Rscript R/venn_predictions.R --db sprot
```

## Pipeline Orchestrator

The `pipeline.R` script coordinates all analysis steps:

```bash
# Show available steps and status
Rscript R/pipeline.R --list --db sprot

# Dry-run (validate inputs without executing)
Rscript R/pipeline.R --check --db sprot

# Run specific step only
Rscript R/pipeline.R --db sprot --step PrionScan

# Export workflow diagram
Rscript R/pipeline.R --diagram
```

## Output Files

### From data_preanalysis.R

| File | Description |
|------|-------------|
| `*_fungi_tabular_*.tsv` | Cached protein data (ID, Organism, Taxonomy, TaxID, Sequence) |
| `*_fungi_taxonomy_*.tsv` | Protein taxonomy mapping |
| `*_fungi_sequences_*.fasta` | FASTA sequences (header: `>ID\|TaxID\|Organism`) |
| `*_genus_guild_mapping_*.tsv` | FunGuild ecological guild annotations |

### From PrionScan.R

| File | Description |
|------|-------------|
| `*_prion_predictions_*.tsv` | Prion domain predictions (Score >= 50) |

Columns: `Protein_ID, Organism, Window_Position, Score, Prion_Domain`

### From PAPA.R

| File | Description |
|------|-------------|
| `*_papa_predictions_*.tsv` | PAPA scores for all proteins |

Columns: `Protein_ID, Organism, PAPA_Score, PAPA_Position, Above_Threshold`

### From plaac.R

| File | Description |
|------|-------------|
| `*_plaac_full_*.tsv` | Full 37-column PLAAC output (all scores and metrics) |
| `*_plaac_predictions_*.tsv` | Standardized predictions for pipeline integration |

Full output includes: COREscore, LLR, PAPA, FoldIndex disorder, MW enrichment, sequence positions.

Standardized columns: `Protein_ID, Organism, Prion_Score, LLR_Score, PAPA_Score, Has_Prion_Domain, Domain_Start, Domain_End, Domain_Length, Domain_Sequence, Protein_Length`

### From venn_predictions.R

| File | Description |
|------|-------------|
| `*_venn_predictions_*.png` | Venn diagram comparing all available prediction methods (2-way or 3-way) |

> **Note**: Taxonomy and NCBI_TaxID are stored in the tabular cache. Join by `Protein_ID` when needed.

## Workflow Diagram

```mermaid
flowchart TD
    DAT[("UniProt .dat")] --> S1[data_preanalysis.R]
    S1 --> TSV[("Tabular Cache")]
    S1 --> TAX[("Taxonomy TSV")]
    S1 --> FASTA[("FASTA Sequences")]
    S1 --> GUILD[("Guild Mapping")]
    TSV --> S2[PrionScan.R]
    TSV --> S3[PAPA.R]
    FASTA --> S4[plaac.R]
    S2 --> PRION[("PrionScan Predictions")]
    S3 --> PAPA_OUT[("PAPA Predictions")]
    S4 --> PLAAC_OUT[("PLAAC Predictions")]
    PRION --> VENN[venn_predictions.R]
    PAPA_OUT --> VENN
    PLAAC_OUT --> VENN
    VENN --> VENN_OUT[("Venn Diagram")]
```

## Algorithm Details

### PrionScan

- **Method**: Sliding window (60 aa) with amino acid propensity scoring
- **Threshold**: Score >= 50 for positive prediction
- **Reference**: Espinosa Angarica et al. (2013) BMC Genomics; Alberti et al. (2009) Cell

### PAPA

- **Method**: Prion propensity with FoldIndex disorder filtering
- **Window**: 41 aa (configurable)
- **Threshold**: Score >= 0.05 for potential prion-forming propensity
- **Reference**: Toombs et al. (2010) Mol Cell Biol 30(1):319-332

### PLAAC

- **Method**: Log-likelihood scoring for Q/N-rich sequences vs background proteome
- **Core length**: 60 aa minimum (configurable)
- **Threshold**: COREscore >= 20 for high-confidence hits
- **Reference**: Lancaster et al. (2014) BMC Bioinformatics; https://github.com/whitehead/plaac

## Performance Notes

- **Tabular caching**: First run converts .dat to TSV (~4-5x faster). Subsequent runs use cache (~15-25x faster)
- **Parallel processing**: PrionScan and PAPA use all available cores by default (`--cores N` to customize)
- **Memory**: TrEMBL processing requires ~16GB+ RAM

## Documentation

See `docs/` for detailed documentation:

- `workflow.md` - Complete analysis workflow and tool descriptions
- `changelog_*.md` - Version history for each script

## License

GNU GPL v3+

## Citations

If you use this pipeline, please cite:

- **PrionScan algorithm**: Espinosa Angarica V, et al. (2013) "PrionScan: an online database of predicted prion domains in complete proteomes." BMC Genomics
- **PAPA algorithm**: Toombs JA, et al. (2010) "Compositional determinants of prion formation in yeast." Mol Cell Biol 30(1):319-332
- **PLAAC algorithm**: Lancaster AK, et al. (2014) "PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition." Bioinformatics
- **FunGuild**: Nguyen NH, et al. (2016) "FUNGuild: An open annotation tool for parsing fungal community datasets by ecological guild." Fungal Ecology
