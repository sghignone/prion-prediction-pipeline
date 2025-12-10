# Changelog - Pipeline Orchestrator

All notable changes to the Pipeline Orchestrator are documented in this file.

## [v1.1] - 2025

### Enhanced Diagram Export

- **NEW**: `--diagram` now exports files to `docs/` instead of just printing to console
  - Writes `docs/pipeline_diagram.md` (Mermaid code)
  - Generates `docs/pipeline_diagram.png` (rendered image)
- **PNG Generation**: Dual-method approach
  - Primary: mermaid.ink API (if internet available and base64enc package installed)
  - Fallback: mmdc CLI (if installed via `npm install -g @mermaid-js/mermaid-cli`)

---

## [v1.0] - 2025

Initial release of the pipeline orchestrator.

### Features

- **Registry-Based Architecture**: Each pipeline step defined with inputs, outputs, and dependencies
- **Dependency Resolution**: Automatically validates step dependencies before execution
- **Input/Output Tracking**: Validates file existence before running steps

### CLI Options

- `--db <sprot|trembl>`: Required database selection
- `--step <name>`: Run single step
- `--from <name>`: Run from step to end
- `--list`: Show available steps and status
- `--check`: Dry-run validation
- `--diagram`: Export workflow diagram to docs/

### Logging

- Timestamped run logs saved to `reports/pipeline_run_*.log`
- Console output mirrors log file

### Pipeline Steps (v1.0)

1. **data_preanalysis**: Convert .dat to TSV, extract taxonomy, map FunGuild
2. **prion_parser**: Predict prion domains from cached sequences

### Embedded Workflow Diagram

```mermaid
flowchart TD
    DAT[("UniProt .dat")] --> S1[data_preanalysis.R]
    S1 --> TSV[("Tabular Cache")]
    S1 --> TAX[("Taxonomy TSV")]
    S1 --> GUILD[("Guild Mapping")]
    TSV --> S2[prion_parser.R]
    S2 --> PRION[("Prion Predictions")]
```
