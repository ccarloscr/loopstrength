#!/usr/bin/env bash
# =============================================================================
# config_loopstrength.sh
# Central configuration for the loopstrength pipeline.
# Edit these values before running run_pipeline.sh
# =============================================================================

# --- Input files -------------------------------------------------------------
REAL_LOOPS="data/real_loops.bedpe"              # Real Hi-C loops in BEDPE format
CONTROL_MCOOL="data/control.mcool"              # Control condition mcool file
EXPERIMENTAL_MCOOL="data/experimental.mcool"    # Experimental condition mcool file
CHROMSIZES="data/dm6.chrom.sizes.txt"           # Chromosome sizes file

# --- Random loop generation --------------------------------------------------
RANDOM_PER_REAL=500                             # Random controls per real loop
WARN_GAP=10000                                  # Warn if anchor gap < this (bp)

# --- Hi-C counts parameters ------------------------------------------------
RESOLUTION=4000                                 # Bin resolution in bp (must exist in mcool)
BALANCE=False                                   # Whether to use balanced or raw counts (True/False) 

# --- Analysis parameters -----------------------------------------------------
MIN_LOOP_SIZE=0                                 # Filter loops shorter than this (bp); 0 = no filter
SIG_THRESHOLD=0.05                              # Adjusted p-value significance cutoff
CONDITION_A="Control"                           # Label for condition A (control)
CONDITION_B="Experimental"                      # Label for condition B (experimental)

# --- Output ------------------------------------------------------------------
OUTPUT_DIR="results"                            # Directory for all pipeline outputs

# --- Canonical chromosomes (space-separated) ---------------------------------
CANONICAL=("chr2L" "chr2R" "chr3L" "chr3R" "chr4" "chrX")
