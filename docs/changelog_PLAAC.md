# Changelog - PLAAC

All notable changes to PLAAC (Prion-Like Amino Acid Composition) launcher are documented here.

## [v1.0] - 2025

Initial pipeline-integrated version.

### Features

- **Java Environment Management**: Auto-detection, Conda integration, setup mode
- **Repository Auto-Setup**: Clones PLAAC from GitHub (`R/tools/plaac/`)
- **Parallel Chunking**: Splits large FASTA for parallel processing (10k seqs/chunk)
- **Dual Output**:
  - `*_plaac_full_*.tsv` - Original 37 PLAAC columns
  - `*_plaac_predictions_*.tsv` - Standardized pipeline format
- **Standalone Mode**: Can run on any FASTA file independently
- **FASTA Validation**: Validates input format before processing
- **Error Handling**: Automatic cleanup of temporary chunk directories

### CLI Options

- `-i, --input <FILE>`: Input FASTA file (REQUIRED)
- `-o, --output-dir <DIR>`: Output directory (default: data/processed/)
- `-d, --db <NAME>`: Database prefix for output naming (sprot/trembl/fungi)
- `-c, --core-length <N>`: Core length (default: 60)
- `-a, --alpha <N>`: Background interpolation (default: 1.0 = S.cerevisiae)
- `-j, --java-env <ENV>`: Conda environment with Java
- `-n, --cores <N>`: CPU cores for parallel processing
- `-p, --per-residue`: Generate per-residue scores
- `--setup-java`: Create Conda environment with Java
- `--setup-plaac`: Download PLAAC repository

### Output Formats

**Full Output** (37 columns):
SEQid, MW, MWstart, MWend, MWlen, LLR, LLRstart, LLRend, LLRlen, NLLR, VITmaxrun, COREscore, COREstart, COREend, CORElen, COREaa, PRDscore, PRDstart, PRDend, PRDlen, PRDaa, PROTlen, HMMall, HMMvit, STARTaa, ENDaa, FInumaa, FImeanhydro, FImeancharge, FImeancombo, FImaxrun, PAPAcombo, PAPAprop, PAPAfi, PAPAllr, PAPAllr2, PAPAcen, PAPAaa

**Pipeline Output** (11 columns):
Protein_ID, Organism, Prion_Score, LLR_Score, PAPA_Score, Has_Prion_Domain, Domain_Start, Domain_End, Domain_Length, Domain_Sequence, Protein_Length

### Directory Structure

- Input: Any FASTA file (standalone) or `data/processed/*_sequences.fasta` (pipeline)
- Output: `data/processed/{db}_plaac_*.tsv`
- Tool: `R/tools/plaac/` (auto-cloned from GitHub)

### Threshold Interpretation

| Score | Interpretation |
|-------|----------------|
| COREscore >= 20 | High confidence prion-like domain |
| COREscore >= 10 | Moderate confidence |
| COREscore < 10 | Low confidence |

### References

- Lancaster AK, et al. (2014) "PLAAC: a web and command-line application to identify proteins with prion-like amino acid composition." Bioinformatics
- https://github.com/whitehead/plaac
