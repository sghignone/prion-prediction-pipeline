# Data Directory

This directory uses a 3-tier structure to organize pipeline data:

## Structure

```
data/
├── raw/          # User-provided input files
├── cache/        # Intermediate/cached files (auto-generated)
└── processed/    # Final outputs for downstream scripts
```

## raw/

Place UniProt fungi database files here:
- `uniprot_sprot_fungi.dat.gz` (Swiss-Prot, ~50MB)
- `uniprot_trembl_fungi.dat.gz` (TrEMBL, ~73GB)

Download from: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/

## cache/

Auto-generated intermediate files (14-day validity):
- `funguild_db_*.rds` - Downloaded FunGuild database
- `*_fungi_tabular_*.tsv` - Converted UniProt .dat to TSV format

These can be safely deleted to force re-generation.

## processed/

Final outputs used by downstream scripts:
- `*_fungi_taxonomy_*.tsv` - Protein ID, Organism, Genus, Species, Full Taxonomy
- `*_genus_guild_mapping_*.tsv` - Ecological guild data for each genus

These files are the primary inputs for subsequent analysis scripts.
