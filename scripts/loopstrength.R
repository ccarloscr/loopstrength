###############################################################################
# QUANTIFICATION OF CHROMATIN LOOP STRENGTH
# -----------------------------------------------------------------------------
# Description: Quantifies chromatin loop strength by comparing real loops
#              against matched random controls using log2 fold change and
#              empirical p-values. Outputs a volcano plot and results table.
#
# Input files (defined in config_loopstrength.sh):
#   - sample_loops.tsv : Real loop counts (from counts_from_mcool.py)
#   - random_loops.tsv : Matched random loop counts (null distribution)
#
# Output files:
#   - volcano_plot.pdf        : Volcano plot of loop strength changes
#   - output_loopstrength.tsv : Full results table with statistics
#
# Dependencies: ggplot2, ggrepel, dplyr
#
# Author: Carlos Camilleri-Robles
# License: MIT
###############################################################################

# --- Package setup -----------------------------------------------------------

# Check for required packages and install if missing
cran_packages <- c("ggplot2", "ggrepel", "dplyr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

# Load libraries with suppressed startup messages for cleaner output
invisible(lapply(cran_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))


# --- Load config from environment or config file -----------------------------

# Function to get configuration values with the following priority:
# 1. environment variables (set by run_pipeline.sh)
# 2. config list (parsed from config_loopstrength.sh)
# 3. default value (if provided)
get_config <- function(key, config, default = NULL) {
  env_val <- Sys.getenv(key, unset = NA)                      # Check env var; returns NA if not set
  if (!is.na(env_val) && nchar(env_val) > 0) return(env_val)  # Only use env var if it's set and non-empty
  if (key %in% names(config)) return(config[[key]])           # Return config value if key exists in config list
  if (!is.null(default)) return(default)                      # Return default if env var and config key are both missing
  stop(paste("Missing required config key:", key))            # Throw error if key is missing and no default provided
}

# Initialize empty config list to store values from the config file
config <- list()

# Define path to config file: can be overridden by CONFIG_PATH env var, otherwise defaults to "config_loopstrength.sh"
config_path <- Sys.getenv("CONFIG_PATH", unset = "config_loopstrength.sh")

# Read config file if it exists and parse KEY=value pairs into the initialized config list
if (file.exists(config_path)) {
  lines <- readLines(config_path)
  # Parse KEY=value lines, ignoring comments and export/blank lines
  kv_lines <- grep("^[A-Za-z_][A-Za-z0-9_]*=", lines, value = TRUE)
  for (line in kv_lines) {
    parts <- strsplit(line, "=", fixed = TRUE)[[1]]           # Split on first '=' to separate key and value
    key   <- trimws(parts[1])                                 # Extract key and trim whitespace
    value <- trimws(paste(parts[-1], collapse = "="))         # Rejoin value parts in case value contains '=' and trim whitespace
    value <- gsub('^["\']|["\']$', "", value)                 # Remove surrounding quotes if present
    config[[key]] <- value                                    # Store in config list
  }
}


# --- Extract configuration values with validation and defaults ----------------
sample_loops_path  <- get_config("SAMPLE_LOOPS",    config)
random_loops_path  <- get_config("RANDOM_LOOPS",    config)
output_directory   <- get_config("OUTPUT_DIR",      config, default = "results")
min_loop_size      <- as.numeric(get_config("MIN_LOOP_SIZE", config, default = "0"))
sig_threshold      <- as.numeric(get_config("SIG_THRESHOLD", config, default = "0.05"))
condition_a_label  <- get_config("CONDITION_A",     config, default = "Control")
condition_b_label  <- get_config("CONDITION_B",     config, default = "Experimental")


# --- Validate input files and output directory --------------------------------

# Check that input files exist before proceeding to avoid runtime errors later on
for (f in c(sample_loops_path, random_loops_path)) {
  if (!file.exists(f)) stop(paste("Input file not found:", f))
}

# Create output directory if it doesn't exist, with recursive = TRUE to create any necessary parent directories
if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE)


# --- Load data and validate columns ------------------------------------------------

# Function to read loops from a TSV file and validate column names
read_loops <- function(path) {
  df <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  col_names <- colnames(df)
  n_cols <- length(col_names)
  
  # Extract the last two column names and compare to expected names
  last_two_names <- col_names[(n_cols - 1):n_cols]
  expected_names <- c("control_count", "experimental_count")
  
  # Check if they match exactly
  if (!all(last_two_names == expected_names)) {
    stop(sprintf(
      "Unexpected column names in %s: expected last two columns to be '%s' and '%s', but got '%s' and '%s'.",
      path, expected_names[1], expected_names[2], last_two_names[1], last_two_names[2]
    ))
  }
  
  # If validation passes, return the data frame
  df
}

# Use tryCatch to provide informative error messages if file reading fails
sample_loops <- tryCatch(read_loops(sample_loops_path), 
                        error = function(e) stop("Error reading sample loops: ", e$message))
random_loops <- tryCatch(read_loops(random_loops_path),
                         error = function(e) stop("Error reading random loops: ", e$message))

# Show summary of loaded data to confirm successful loading and check for expected columns
cat(sprintf("Loaded %d real loops and %d random loops.\n",
            nrow(sample_loops), nrow(random_loops)))


# --- Compute loop size and filter --------------------------------------------

# Define loop size as the distance between the two loop anchors
add_loop_size <- function(df) {
  df %>% mutate(loop_size = start2 - start1)
}

# Add loop_size column to both real and random loops
sample_loops <- add_loop_size(sample_loops)
random_loops <- add_loop_size(random_loops)

# Filter out loops smaller than the specified minimum loop size, if applicable
if (min_loop_size > 0) {
  n_before <- nrow(sample_loops)
  sample_loops <- sample_loops %>% filter(loop_size >= min_loop_size)
  random_loops <- random_loops %>% filter(loop_size >= min_loop_size)
  cat(sprintf("Filtered loops < %g bp: %d real loops retained (removed %d).\n",
              min_loop_size, nrow(sample_loops), n_before - nrow(sample_loops)))
}


# --- Compute log2 fold change ------------------------------------------------

# Pseudocount of 1 added to avoid log(0); symmetric for both conditions
compute_logfc <- function(df) {
  df %>% mutate(logFC = log2((experimental_count + 1) / (control_count + 1)))
}

# Compute logFC for both real and random loops
sample_loops <- compute_logfc(sample_loops)
random_loops <- compute_logfc(random_loops)


# --- Empirical p-values (two-sided) ------------------------------------------
# For each real loop, p = proportion of null |logFC| >= observed |logFC|
# Floor at 1/N to avoid p=0 -> Inf on -log10 scale
# Where N is the number of matched null loops per real loop

# Use mapply to vectorize over logFC and loop_id columns of sample_loops,
# calculating p-values based on matched random loops
sample_loops$pval <- mapply(
  function(logfc, lid) {
    # Extract the logFC values of the matched random loops for this real loop ID
    matched_null <- random_loops$logFC[random_loops$source_loop_id == lid]

    # If there are fewer than 50 matched null loops, return NA to indicate insufficient data for p-value calculation
    if (length(matched_null) < 50) {
      return(NA_real_)
    }

    # Calculate the empirical p-value as the proportion of matched null logFC values with absolute value
    # greater than or equal to the observed absolute logFC
    raw <- mean(abs(matched_null) >= abs(logfc))
    max(raw, 1 / length(matched_null))
  },
  # Vectorized over logFC and loop_id columns of the sample_loops data frame
  sample_loops$logFC,
  sample_loops$loop_id
)


# --- Multiple testing correction (Benjamini-Hochberg) -----------------------

# Adjust p-values for multiple testing using the BH method to control FDR
sample_loops$padj <- p.adjust(sample_loops$pval, method = "BH")

# Check for all NA padj values, stop if true to avoid misleading results and empty plots
if (all(is.na(sample_loops$padj))) {
  stop("All padj values are NA — no valid counts were fetched. Check pipeline logs.")
}

# Summarize significant loops based on the specified significance threshold, ignoring NA values
n_sig <- sum(sample_loops$padj < sig_threshold, na.rm = TRUE)
cat(sprintf("Significant loops (padj < %.2f): %d / %d\n",
            sig_threshold, n_sig, nrow(sample_loops)))


# --- Volcano plot ------------------------------------------------------------

# Dynamic y-axis: headroom above max finite -log10(padj)
y_max <- max(-log10(sample_loops$padj[is.finite(-log10(sample_loops$padj))]),
             na.rm = TRUE) * 1.15

# Create volcano plot with ggplot2, coloring points by significance and adding reference lines
p <- ggplot(sample_loops, aes(x = logFC, y = -log10(padj))) +
  geom_point(
    aes(color = padj < sig_threshold),
    alpha = 0.7, size = 1.8) +
  scale_color_manual(
    values = c("TRUE" = "#d62728", "FALSE" = "grey50"),
    labels = c("TRUE" = paste0("padj < ", sig_threshold), "FALSE" = "n.s."),
    name = NULL) +
  geom_hline(yintercept = -log10(sig_threshold),
             linetype = "dashed", color = "#d62728", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey70", linewidth = 0.3) +
  theme_classic(base_size = 12) +
  labs(
    title = paste("Loop Strength:", condition_a_label, "vs", condition_b_label),
    x     = expression(log[2] ~ "Fold Change"),
    y     = expression(-log[10] ~ "(adjusted p-value)")) +
  coord_cartesian(ylim = c(0, y_max))

# Add loop_id labels only if the column exists and there are significant loops to label
if ("loop_id" %in% colnames(sample_loops)) {
  p <- p + geom_text_repel(
    data = subset(sample_loops, padj < sig_threshold),
    aes(label = loop_id), size = 2.5, max.overlaps = 20
  )
}

# Save volcano plot to PDF with specified dimensions, suppressing output to console
pdf(file.path(output_directory, "volcano_plot.pdf"), width = 7, height = 5)
print(p)
invisible(dev.off())


# --- Output table ------------------------------------------------------------
write.table(
  sample_loops,
  file      = file.path(output_directory, "output_loopstrength.tsv"),
  quote     = FALSE,
  row.names = FALSE,
  sep       = "\t"
  # dec defaults to "." — do not override unless specifically needed
)

cat("Done. Results saved in:", output_directory, "\n")