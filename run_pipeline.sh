#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh
#
# Loopstrength pipeline — end-to-end execution
#
# Usage:
#   bash run_pipeline.sh [config_file]
#
# If no config file is provided, defaults to config_loopstrength.sh
# All parameters are defined in the config file.
#
# Pipeline steps:
#   1. Generate matched random control loops (bash)
#   2. Fetch Hi-C counts for real loops (Python)
#   3. Fetch Hi-C counts for random loops (Python)
#   4. Compute loop strength statistics + volcano plot (R)
#
# Output structure:
#   results/
#   ├── random_loops.bedpe
#   ├── sample_loops_counts.tsv
#   ├── random_loops_counts.tsv
#   ├── output_loopstrength.tsv
#   └── volcano_plot.pdf
#
# =============================================================================

# Set strict mode
set -euo pipefail

# --- Load config -------------------------------------------------------------

# Default config file (can be overridden by CLI argument)
CONFIG="${1:-config_loopstrength.sh}"

# Validate config file exists
if [[ ! -f "$CONFIG" ]]; then
  echo "Error: config file not found: $CONFIG" >&2
  exit 1
fi

# Source the config file to load parameters
source "$CONFIG"

# --- Validate required config variables -------------------------------------

# List of required config variables (must be defined in the config file)
required_vars=(
  REAL_LOOPS
  CONTROL_MCOOL
  EXPERIMENTAL_MCOOL
  CHROMSIZES
  RANDOM_PER_REAL
  WARN_GAP
  RESOLUTION
  BALANCE
  MIN_LOOP_SIZE
  SIG_THRESHOLD
  CONDITION_A
  CONDITION_B
  OUTPUT_DIR
)

# Check that all required config variables are set and non-empty
for v in "${required_vars[@]}"; do
  [[ -z "${!v:-}" ]] && {
    echo "Error: required config variable missing: $v"
    exit 1
  }
done

# Additional check for CANONICAL array (must have at least one element)
if [[ ${#CANONICAL[@]} -eq 0 ]]; then
  echo "Error: required config variable missing: CANONICAL"
  exit 1
fi

# --- Locate child scripts ----------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"        # Directory containing this script
SCRIPTS="$SCRIPT_DIR/scripts"                                     # Directory containing child scripts (assumed to be in "scripts/" subdir)

RANDOM_SCRIPT="$SCRIPTS/generate_random_loops_matched.sh"         # Script to generate matched random loops
COUNTS_SCRIPT="$SCRIPTS/counts_from_mcool.py"                     # Script to fetch Hi-C counts from mcool files
ANALYSIS_SCRIPT="$SCRIPTS/loopstrength.R"                         # Script to perform statistical analysis and generate volcano plot

# --- Validate child scripts exist ------------------------------------------------
for f in "$RANDOM_SCRIPT" "$COUNTS_SCRIPT" "$ANALYSIS_SCRIPT"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: required script not found: $f" >&2
    exit 1
  fi
done

# --- Validate required inputs ------------------------------------------------
for f in "$REAL_LOOPS" "$CONTROL_MCOOL" "$EXPERIMENTAL_MCOOL" "$CHROMSIZES"; do
  if [[ ! -f "$f" ]]; then
    echo "Error: required input file not found: $f" >&2
    exit 1
  fi
done

# --- Check required tools ----------------------------------------------------
for tool in python3 Rscript; do
  if ! command -v "$tool" &> /dev/null; then
    echo "Error: required tool not found in PATH: $tool" >&2
    exit 1
  fi
done

# --- Set derived paths -------------------------------------------------------
mkdir -p "$OUTPUT_DIR"                                                 # Create output directory if it doesn't exist

RANDOM_BEDPE="$OUTPUT_DIR/random_loops.bedpe"                          # Output path for random loops BEDPE file
SAMPLE_COUNTS="$OUTPUT_DIR/sample_loops_counts.tsv"                    # Output path for real loop counts
RANDOM_COUNTS="$OUTPUT_DIR/random_loops_counts.tsv"                    # Output path for random loop counts

# --- Logging -----------------------------------------------------------------
LOG="$OUTPUT_DIR/pipeline.log"                                         # Log file for pipeline execution    
echo "=============================================" | tee "$LOG"
echo "  Loopstrength pipeline" | tee -a "$LOG"
echo "  Started: $(date)" | tee -a "$LOG"
echo "  Config:  $CONFIG" | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }               # Timestamped logging function
die() { echo "Error: $*" | tee -a "$LOG" >&2; exit 1; }                # Error handling function

# =============================================================================
# STEP 1 — Generate matched random control loops
# =============================================================================
log "Step 1/4: Generating matched random loops..."

bash "$RANDOM_SCRIPT" \
  --real-loops    "$REAL_LOOPS" \
  --chromsizes    "$CHROMSIZES" \
  --output        "$RANDOM_BEDPE" \
  --per-real      "$RANDOM_PER_REAL" \
  --warn-gap      "$WARN_GAP" \
  --canonical     "${CANONICAL[*]}" \
  2>> "$LOG" || die "Step 1 failed. Check $LOG for details."

log "Step 1 complete. Random loops: $RANDOM_BEDPE"

# Validation to avoid proceeding with downstream steps if the random loop generation failed or produced no output
[[ -s "$RANDOM_BEDPE" ]] || die "Step 1 produced an empty output file. Check $LOG for skipped loop warnings."

# =============================================================================
# STEP 2 — Fetch Hi-C counts for real loops
# =============================================================================
log "Step 2/4: Fetching Hi-C counts for real loops..."

python3 "$COUNTS_SCRIPT" \
  "$REAL_LOOPS" \
  "$CONTROL_MCOOL" \
  "$EXPERIMENTAL_MCOOL" \
  "$RESOLUTION" \
  "$BALANCE" \
  > "$SAMPLE_COUNTS" \
  2>> "$LOG" \
|| die "Step 2 failed. Check $LOG for details."

log "Step 2 complete. Real loop counts: $SAMPLE_COUNTS"

# =============================================================================
# STEP 3 — Fetch Hi-C counts for random loops
# =============================================================================
log "Step 3/4: Fetching Hi-C counts for random loops..."

python3 "$COUNTS_SCRIPT" \
  "$RANDOM_BEDPE" \
  "$CONTROL_MCOOL" \
  "$EXPERIMENTAL_MCOOL" \
  "$RESOLUTION" \
  "$BALANCE" \
  > "$RANDOM_COUNTS" \
  2>> "$LOG" \
|| die "Step 3 failed. Check $LOG for details."

log "Step 3 complete. Random loop counts: $RANDOM_COUNTS"

# =============================================================================
# STEP 4 — Statistical analysis and volcano plot
# =============================================================================
log "Step 4/4: Running loop strength analysis in R..."

# Export variables so R can read them via Sys.getenv()
export SAMPLE_LOOPS="$SAMPLE_COUNTS"
export RANDOM_LOOPS="$RANDOM_COUNTS"
export CONFIG_PATH="$(realpath "$CONFIG")"
export OUTPUT_DIR="$OUTPUT_DIR"
export MIN_LOOP_SIZE="$MIN_LOOP_SIZE"
export SIG_THRESHOLD="$SIG_THRESHOLD"
export CONDITION_A="$CONDITION_A"
export CONDITION_B="$CONDITION_B"

Rscript "$ANALYSIS_SCRIPT" 2>> "$LOG" \
|| die "Step 4 failed. Check $LOG for details."

log "Step 4 complete."

# =============================================================================
# Summary
# =============================================================================
echo "" | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
echo "  Pipeline complete: $(date)" | tee -a "$LOG"
echo "  Outputs in: $OUTPUT_DIR/" | tee -a "$LOG"
echo "    random_loops.bedpe" | tee -a "$LOG"
echo "    sample_loops_counts.tsv" | tee -a "$LOG"
echo "    random_loops_counts.tsv" | tee -a "$LOG"
echo "    output_loopstrength.tsv" | tee -a "$LOG"
echo "    volcano_plot.pdf" | tee -a "$LOG"
echo "    pipeline.log" | tee -a "$LOG"
echo "=============================================" | tee -a "$LOG"
