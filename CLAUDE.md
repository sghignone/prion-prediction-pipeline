## Workflow
1.	First, think through the problem. Read the codebase and write a plan in task/todo.md.
2.	The plan should be a checklist of todo items.
3.	Check in with me before starting work-I'll verify the plan.
4.	Then, complete the todos one by one, marking them off as you go.
5.	At every step, give me a high-level explanation of what you changed.
6.	Keep every change simple and minimal. Avoid big rewrites.
7.	At the end, add a review section in todo.md summarizing the changes.

## Required R Packages

- tidyverse
- FUNGuildR
- patchwork


## Input Data

Place UniProt fungi proteome files in `data/raw/`:
- `uniprot_sprot_fungi.dat` or `.dat.gz` (Swiss-Prot, ~50MB)
- `uniprot_trembl_fungi.dat` or `.dat.gz` (TrEMBL, ~73GB)

Download from: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/

## Directory Structure

```
data/
├── raw/          # User-provided input files
├── cache/        # Intermediate/cached files
└── processed/    # Final outputs for downstream scripts
```

