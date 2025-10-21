# AI Assistant Instructions for Somatem

## Project Overview
Somatem is a long-read metagenomics pipeline built with Nextflow, focusing on 4th generation sequencing technologies (Oxford Nanopore and PacBio). The pipeline implements best practices for long-read analysis and is designed to be LLM-accessible.

## Key Architecture Components

### Workflow Structure
- Entry point: `main.nf` orchestrates the pipeline through `ORCHESTRATE_SOMATEM` workflow
- Core subworkflows in `/subworkflows/`:
  - Data preprocessing: Quality visualization (NanoPlot), host removal (hostile), sequence filtering (chopper)
  - Taxonomic classification
  - Genome assembly/MAG analysis
  - Pathogen detection (seqscreen integration)

### Critical Directories
- `/modules/`: Core process definitions
- `/subworkflows/`: Modular analysis components
- `/assets/`: Configuration files and sample sheets
- `/conf/`: Environment and resource configurations
- `/assets/databases/`: Reference database location (configurable)

## Development Workflows

### Environment Setup
```bash
micromamba create -n nf_base_env nextflow
micromamba activate nf_base_env
```

### Running Pipeline
1. Configure metadata:
   ```bash
   cp docs/somatem_docs/metadata_template.yaml assets/custom_metadata.yaml
   ```
2. Execute workflow:
   ```bash
   nextflow run . -param-file assets/custom_metadata.yaml
   ```

### Database Configuration
- Default path: `db_base_dir = "/home/dbs"`
- Large databases (bakta, checkm2, singlem) require ~100GB space
- Update path in `nextflow.config`

## Key Integration Points

### Database Dependencies
- Databases are auto-downloaded when needed
- Key databases:
  - checkm2 (~60GB)
  - gtdbtk (~140GB)
  - singlem metapackage
  - hostile reference databases

### Container & Conda Management
- Uses micromamba for environment management
- Containers managed through nextflow configuration
- Critical for reproducibility across environments

## Project-Specific Conventions

### Parameter Management
- Use YAML configuration files for run parameters
- Store in `assets/` directory
- Follow `metadata_template.yaml` structure

### Resource Handling
- Significant compute requirements for assembly
- Example: 2-sample assembly needs ~128GB RAM, 128 CPUs
- 6TB storage recommended for full pipeline

## Common Operations
- Test with example data: `nextflow run subworkflows/local/get_example_data.nf`
- Database setup: Check `nextflow.config` for `db_base_dir` setting
- Monitor large processes: Assembly can take 6+ hours

## Best Practices
- Always activate `nf_base_env` before running pipeline
- Verify database paths before large runs
- Use shared database directory for multi-user setups
- Monitor storage during MAG analysis workflow