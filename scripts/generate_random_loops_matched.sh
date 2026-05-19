#!/usr/bin/env bash
# =============================================================================
# generate_random_loops_matched.sh
#
# Generates matched random cis loops for each real Hi-C loop.
#
# For every real loop, generates N random loops preserving:
#   - Chromosome
#   - Anchor sizes
#   - Loop distance (start-to-start)
#
# Usage (CLI) with example values:
#   bash generate_random_loops_matched.sh \
#     --real-loops  data/real_loops.bedpe \
#     --chromsizes  data/dm6.chrom.sizes.txt \
#     --output      results/random_loops.bedpe \
#     --per-real    5 \
#     --warn-gap    10000 \
#     --canonical   "chr2L chr2R chr3L chr3R chr4 chrX"
#
#
# Input BEDPE format (tab-separated, at least 6 columns):
#   chrA  startA  endA  chrB  startB  endB  loop_id  extra_fields...
#
# =============================================================================

# --- Strict mode -------------------------------------------------------------
set -euo pipefail

# --- Default configuration (Initialized empty, must be set via CLI) ----------
REAL_LOOPS=""           # Path to input real loops BEDPE file
CHROMSIZES=""           # Path to input chromosome sizes file
OUTPUT=""               # Path to output random loops BEDPE file
RANDOM_PER_REAL=""      # Number of random loops to generate per real loop
WARN_GAP=""             # Warn if gap between anchors is less than this (bp)
CANONICAL=()            # Array of canonical chromosome names (e.g. "chr2L chr2R chr3L chr3R chr4 chrX")


# --- Parse CLI arguments ------------------------------------------------------
while [[ $# -gt 0 ]]; do                                      # Loop over arguments
  case "$1" in                                                # Match known flags and assign values
    --real-loops)  REAL_LOOPS="$2"; shift 2 ;;                # Path to real loops BEDPE file
    --chromsizes)  CHROMSIZES="$2"; shift 2 ;;                # Path to chromosome sizes file
    --output)      OUTPUT="$2"; shift 2 ;;                    # Path to output random loops BEDPE file
    --per-real)    RANDOM_PER_REAL="$2"; shift 2 ;;           # Number of random loops to generate per real loop
    --warn-gap)    WARN_GAP="$2"; shift 2 ;;                  # Warn if gap between anchors is less than this (bp)
    --canonical)   read -r -a CANONICAL <<< "$2"; shift 2 ;;  # Space-separated list of canonical chromosomes
    *) echo "Error: Unknown argument $1" >&2; exit 1 ;;       # Handle unknown arguments
  esac
done

# --- Strict Parameter Validation ---------------------------------------------
if [[ -z "$REAL_LOOPS" || -z "$CHROMSIZES" || -z "$OUTPUT" || -z "$RANDOM_PER_REAL" || -z "$WARN_GAP" || ${#CANONICAL[@]} -eq 0 ]]; then
  echo "Error: Missing required arguments" >&2
  echo "Usage: $0 --real-loops <file> --chromsizes <file> --output <file> --per-real <int> --warn-gap <int> --canonical <string>" >&2
  exit 1
fi

# --- Create output directory -------------------------------------------------
mkdir -p "$(dirname "$OUTPUT")"

# --- Temporary files ---------------------------------------------------------
TMP_WORK=$(mktemp -d)
trap 'rm -rf "$TMP_WORK"' EXIT

CANONICAL_SIZES="$TMP_WORK/canonical.chrom.sizes"

# --- Filter to canonical chromosome sizes ------------------------------------
# grep options: -F=fixed strings, -f=pattern list from file
# printf "%s\n" outputs each array element on a new line
# > (...) converts printed elements into a temporary virtual file
grep -Ff <(printf "%s\n" "${CANONICAL[@]}") "$CHROMSIZES" > "$CANONICAL_SIZES"


# --- Read chromosome sizes into an associative array for quick lookup --------
# Declare an associative array
declare -A CHROM_SIZES

# Read chromosome sizes into the associative array
# The [[ -n "$chr" ]] condition ensures that we process the last line even if it doesn't end with a newline character
# CANONICAL_SIZES is used here because it contains only the canonical chromosomes
while read -r chr size || [[ -n "$chr" ]]; do
  CHROM_SIZES[$chr]=$size
done < "$CANONICAL_SIZES"

# --- Validate all canonical chromosomes are present --------------------------
for chr in "${CANONICAL[@]}"; do
  if [[ -z "${CHROM_SIZES[$chr]:-}" ]]; then        # The :- syntax prevents an unbound variable error if the chromosome is missing
    echo "Error: chromosome $chr not found in $CHROMSIZES" >&2
    exit 1
  fi
done

# --- Seed for reproducibility ------------------------------------------------------
SEED=123

# --- Build canonical sizes string for awk ------------------------------------
# Pass chromosome sizes as a string "chr:size chr:size ..." to awk
CHROM_SIZE_STR=""
for chr in "${CANONICAL[@]}"; do
  CHROM_SIZE_STR+="${chr}:${CHROM_SIZES[$chr]} "
done

# --- Run everything in a single awk pass -------------------------------------
awk -v seed="$SEED" \
    -v per_real="$RANDOM_PER_REAL" \
    -v warn_gap="$WARN_GAP" \
    -v chrom_size_str="$CHROM_SIZE_STR" \
'
BEGIN {
  srand(seed)

  # Parse chrom sizes string into array
  n = split(chrom_size_str, pairs, " ")
  for (i = 1; i <= n; i++) {
    split(pairs[i], kv, ":")
    chrom_sizes[kv[1]] = kv[2] + 0
  }

  real_used     = 0
  real_skipped  = 0
  rand_generated = 0
  short_gap_warn = 0
  line_no       = 0
}

{
  line_no++

  # Delete any carriage return characters if present
  sub(/\r$/, "")

  # Skip empty, comment, and header lines
  if ($0 == "") next
  if (substr($1, 1, 1) == "#") next
  if ($2 == "startA" || $2 == "start1") next

  chr1=$1; start1=$2+0; end1=$3+0
  chr2=$4; start2=$5+0; end2=$6+0
  loop_id = (NF >= 7) ? $7 : ""

  # Validate required fields
  if (NF < 6) {
    print "Warning: malformed line " line_no | "cat >&2"
    real_skipped++; next
  }

  # Validate numeric coordinates (awk already coerces, check they are > 0 sensibly)
  if (start1 == 0 && $2 != "0") { print "Warning: skipping line " line_no " — coordinates not valid integers." | "cat >&2"; real_skipped++; next }

  # Require cis loop
  if (chr1 != chr2) {
    print "Warning: skipping line " line_no " — trans loop (" chr1 " vs " chr2 ")." | "cat >&2"
    real_skipped++; next
  }

  # Require canonical chromosome
  if (!(chr1 in chrom_sizes)) {
    print "Warning: skipping line " line_no " — non-canonical chromosome: " chr1 "." | "cat >&2"
    real_skipped++; next
  }

  # Require positive anchor sizes
  if (end1 <= start1 || end2 <= start2) {
    print "Warning: skipping line " line_no " — non-positive anchor size." | "cat >&2"
    real_skipped++; next
  }

  # Normalize anchor order (upstream first)
  if (start1 <= start2) {
    up_start=start1; up_end=end1; down_start=start2; down_end=end2
  } else {
    up_start=start2; up_end=end2; down_start=start1; down_end=end1
  }

  sizeA     = up_end   - up_start
  sizeB     = down_end - down_start
  loop_dist = down_start - up_start
  gap       = down_start - up_end

  if (gap < 0) {
    print "Warning: line " line_no " has overlapping anchors (gap=" gap " bp)." | "cat >&2"
  } else if (gap < warn_gap) {
    print "Warning: line " line_no " has short anchor gap: " gap " bp < WARN_GAP=" warn_gap " bp." | "cat >&2"
    short_gap_warn++
  }

  chrom_size  = chrom_sizes[chr1]
  maxStartA_A = chrom_size - sizeA
  maxStartA_B = chrom_size - loop_dist - sizeB
  maxStartA   = (maxStartA_A < maxStartA_B) ? maxStartA_A : maxStartA_B

  if (maxStartA < 0) {
    print "Warning: skipping line " line_no " — loop does not fit on " chr1 "." | "cat >&2"
    real_skipped++; next
  }

  source_id = (loop_id == "") ? "real_" line_no : loop_id

  real_used++

  for (i = 1; i <= per_real; i++) {
    rand_startA = int(rand() * (maxStartA + 1))
    rand_endA   = rand_startA + sizeA
    rand_startB = rand_startA + loop_dist
    rand_endB   = rand_startB + sizeB
    random_id   = source_id "_random_" i

    printf "%s\t%d\t%d\t%s\t%d\t%d\t%s\t%s\n", \
      chr1, rand_startA, rand_endA, \
      chr1, rand_startB, rand_endB, \
      random_id, source_id
    rand_generated++
  }
}

END {
  print "Done."                                                          | "cat >&2"
  print "  Real loops used:    " real_used                              | "cat >&2"
  print "  Real loops skipped: " real_skipped                          | "cat >&2"
  print "  Random loops:       " rand_generated " (" per_real " per real loop)" | "cat >&2"
  print "  Short-gap warnings: " short_gap_warn " (threshold: " warn_gap " bp)" | "cat >&2"
}
' "$REAL_LOOPS" > "$OUTPUT"

echo "  Output: $OUTPUT"
