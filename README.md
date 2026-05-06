# loopstrength

Quantification of chromatin loop strength using empirical null distributions and BH-corrected p-values.

## Overview

This repository contains an R script for differential chromatin loop strength analysis from Hi-C data across two experimental conditions.

For each previously identified loop, the script computes the log2-transformed fold change (log2FC) of contact strength between conditions, using a pseudocount of 1 to handle zero counts.

Statistical significance is assessed using a non-parametric, randomization-based empirical null approach. A set of randomly shifted loops — preserving the original loop size distribution — is provided alongside the real loops. The log2FC distribution of these random loops serves as the empirical null. Two-sided empirical p-values are computed for each real loop by comparing its absolute log2FC against the null distribution. P-values are then corrected for multiple testing using the Benjamini-Hochberg False Discovery Rate (FDR) method.

Results are exported as a tab-delimited file and visualised as a volcano plot, allowing identification of loops significantly enriched in either condition and providing insight into condition-specific 3D genome organisation.



## Requirements
- R (≥ 4.0)
- R packages: ggplot2, ggrepel, dplyr
  
The script will automatically install missing packages.


## Installation
To clone this repository:
```bash
git clone https://github.com/ccarloscr/loopstrength.git
cd loopstrength
```


## Input files
Place your input files in the [data/](./data/) folder:
- [real_loops.tsv](./data/real_loops.tsv): Sample chromatin loops
- [random_loops.tsv](./data/random_loops.tsv): Random chromatin loops (null distribution)
Each file must be tab-delimited and contain the following columns (no header required):


## Configuration file
Edit [config_loopstrength.txt](./config_loopstrength.txt) to define input and output paths:


## Output
- output_loopstrength.txt: Table with loop statistics, log2 fold changes, p-values, and adjusted p-values.
- volcano_plot.pdf: Volcano plot highlighting significant loops.


## Notes
- Loops smaller than 1 Mb are excluded.
- Empirical p-values are computed using the random loop distribution.
- Adjusted p-values use Benjamini-Hochberg correction.


## License
This project is licensed under the MIT License.
