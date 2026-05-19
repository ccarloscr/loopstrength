#!/usr/bin/env python3
"""
counts_from_mcool.py

Fetches raw Hi-C contact counts for a set of loops from two mcool files.

Usage:
    python counts_from_mcool.py <loops.bedpe> <control.mcool> <experimental.mcool> [resolution] [balance]

Arguments:
    loops.bedpe             Loop anchors in BEDPE format (6–8 columns; see below)
    control.mcool           mcool file for condition A (control)
    experimental.mcool      mcool file for condition B (experimental)
    resolution              Bin resolution in bp (default: 4000). Must exist in mcool.
    balance                 Whether to use balanced counts (True/False, default: False)

BEDPE format (tab-separated):
    chrom1  start1  end1  chrom2  start2  end2  loop_id  source_loop_id

Output (TSV, stdout):
    All input columns + control_count + experimental_count

Author: Carlos Camilleri-Robles
License: MIT
"""

import sys
import pandas as pd
import numpy as np
import cooler


# Column names for 8 column BEDPE files
_BEDPE_COLS = {
    8: ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "loop_id", "source_loop_id"],
}


# Function to load BEDPE file
def load_bedpe(path: str) -> pd.DataFrame:
    """Load a BEDPE file, handling 8 columns and optional headers."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype=str)

    # Drop header row if present (first row has non-numeric start1)
    if not df.iloc[0, 1].lstrip("-").isdigit():
        df = df.iloc[1:].reset_index(drop=True)

    # Warn if input has more than 8 columns
    if df.shape[1] > 8:
        print(
            f"Warning: input has {df.shape[1]} columns; trimming to first 8.",
            file=sys.stderr
        )

    # Trim to supported column count
    ncols = min(df.shape[1], 8)
    df = df.iloc[:, :ncols]

    # Validate column count for fewer than 8 columns
    if ncols not in _BEDPE_COLS:
        raise ValueError(f"Unsupported number of columns: {df.shape[1]}. Expected 8.")

    # Assign column names
    df.columns = _BEDPE_COLS[ncols]

    # Cast coordinate columns to int
    for col in ["start1", "end1", "start2", "end2"]:
        df[col] = df[col].astype(int)

    return df


# Function to fetch counts from cooler
def get_counts(clr: cooler.Cooler, df: pd.DataFrame, balance: bool) -> list:
    """Fetch summed raw or balanced contact counts for each loop anchor pair."""
    counts = []
    for _, row in df.iterrows():
        try:
            matrix = clr.matrix(balance=balance).fetch(
                (row.chrom1, row.start1, row.end1),
                (row.chrom2, row.start2, row.end2),
            )
            counts.append(float(matrix.sum()) if balance else int(matrix.sum()))
        except Exception as e:
            print(
                f"Warning: could not fetch counts for "
                f"{row.chrom1}:{row.start1}-{row.end1} x "
                f"{row.chrom2}:{row.start2}-{row.end2} — {e}",
                file=sys.stderr,
            )
            counts.append(np.nan)
    return counts


# Main function
def main():
    if len(sys.argv) < 4:
        print(__doc__, file=sys.stderr)
        sys.exit(1)

    bedpe_path          = sys.argv[1]
    control_mcool       = sys.argv[2]
    experimental_mcool  = sys.argv[3]
    resolution          = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
    balance_arg         = sys.argv[5] if len(sys.argv) > 5 else "False"
    
    if balance_arg not in ["True", "False"]:
        sys.exit(f"Error: balance argument must be 'True' or 'False', got '{balance_arg}'")
    use_balance = balance_arg == "True"

    print(f"Resolution: {resolution} bp", file=sys.stderr)
    print(f"Balance: {use_balance}", file=sys.stderr)

    # Load coolers
    clr_control = cooler.Cooler(f"{control_mcool}::/resolutions/{resolution}")
    clr_experimental   = cooler.Cooler(f"{experimental_mcool}::/resolutions/{resolution}")

    # Load loops
    loops = load_bedpe(bedpe_path)
    print(f"Loaded {len(loops)} loops from {bedpe_path}", file=sys.stderr)

    # Fetch counts
    print("Fetching control counts...", file=sys.stderr)
    loops["control_count"] = get_counts(clr_control, loops, use_balance)

    print("Fetching experimental counts...", file=sys.stderr)
    loops["experimental_count"] = get_counts(clr_experimental, loops, use_balance)

    # Write TSV to stdout
    loops.to_csv(sys.stdout, sep="\t", index=False)
    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
