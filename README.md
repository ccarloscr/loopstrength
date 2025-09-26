# loopstrength

This repository contains an R script designed to perform statistical analysis of chromatin loops derived from Hi-C data across two experimental conditions. The workflow takes as input previously identified chromatin loops for each condition and evaluates their relative enrichment. The goal is to detect loops that are significantly enriched in one condition compared to the other, providing insights into condition-specific 3D genome organization.


## Requirements
- R (≥ 4.0)
- R packages: ggplot2, ggrepel, dplyr
The script will automatically install missing packages.


## Input files
Place your input files in the data/ folder:
- real_loops.tsv: Sample chromatin loops
- random_loops.tsv: Random chromatin loops (null distribution)
Each file must be tab-delimited and contain the following columns (no header required):


## Configuration file
Edit config_loopstrength.txt to define input and output paths:


## Output
- output_loopstrength.txt: Table with loop statistics, log2 fold changes, p-values, and adjusted p-values.
- volcano_plot.pdf: Volcano plot highlighting significant loops.


## Notes
- Loops smaller than 1 Mb are excluded.
- Empirical p-values are computed using the random loop distribution.
- Adjusted p-values use Benjamini-Hochberg correction.
