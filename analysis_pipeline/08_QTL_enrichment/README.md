
---

# cis-e/sQTL Annotation and Enrichment Analysis

This repository includes two scripts for annotation and enrichment analysis tasks using GREGOR and the Variant Effect Predictor (VEP).

## Dependencies
- **GREGOR**: Ensure GREGOR is installed and configured.
- **VEP**: Verify that VEP is installed with necessary cache and plugins.

## Scripts and Usage

### 1. GREGOR Analysis Script
Analysis using GREGOR.
#### Parameters
```
- `input_directory`: Directory containing initial data files.
- `output_directory`: Target directory for output files.
- `config_base_directory`: Directory with GREGOR configuration files.
- `work_directory`: Location of GREGOR tool and scripts.
```

### 2. VEP Annotation Script
Annotation using VEP.
#### Parameters
```
- `vcf_file`: Path to the VCF file for annotation.
- `output_file`: Path for the annotated VCF.
- `cache_directory`: Directory for VEP cache and plugins.
```

## Important Notes
- Adjust the resource specifications in Slurm commands based on your system’s capabilities.
- Ensure that the correct paths and modules are specified for both tools in their respective scripts.

---

