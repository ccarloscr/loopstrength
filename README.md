# loopstrength

This repository contains an R script designed to perform statistical analysis of chromatin loops derived from Hi-C data across two experimental conditions. The workflow takes as input previously identified chromatin loops for each condition and evaluates their relative enrichment. The goal is to detect loops that are significantly enriched in one condition compared to the other, providing insights into condition-specific 3D genome organization.

Key features:
- Input: list of chromatin loops identified from Hi-C data.
- Calculation of relative loop strength (log2-transformed fold change).
- Bootstrap-based statistical testing to assess differential loop enrichment.
- Output: annotated list of loops enriched in one condition.





