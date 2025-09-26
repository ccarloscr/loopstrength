###############################################################################
# QUANTIFICATION OF CHROMATIN LOOP STRENGTH
# -----------------------------------------------------------------------------
# Description: This script quantifies chromatin loop strength by comparing
#              sample and random loops using log2 fold change and empirical
#              p-values. It generates a volcano plot and outputs a results table.
#
# Input files:
#   - real_loops.tsv: Sample chromatin loops
#   - random_loops.tsv: Random chromatin loops (null distribution)
#
# Output files:
#   - output_loopstrength.txt: Table with loop statistics and p-values
#   - Volcano plot
#
# Dependencies:
#   - R packages: ggplot2, ggrepel, dplyr
#
# Author: Carlos Camilleri-Robles
# Date: 26-09-2025
# License: [MIT]
###############################################################################

# Install packages
cran_packages <- c("ggplot2", "ggrepel", "dplyr")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}


# Load all packages
lapply(cran_packages, library, character.only = TRUE)


# Load configuration from external .txt file
config_path <- "config_loopstrength.txt"
config_lines <- readLines(config_path)
config <- setNames(
  sapply(strsplit(config_lines, "="), `[`, 2),
  sapply(strsplit(config_lines, "="), `[`, 1)
)


# Extract paths from config
sample_loops_path <- config["sample_loops_path"]
random_loops_path <- config["random_loops_path"]
output_directory <- config["output_directory"]


# Validate config
required_keys <- c("sample_loops_path", "random_loops_path", "output_directory")
missing_keys <- setdiff(required_keys, names(config))
if (length(missing_keys) > 0) {
  stop(paste("Missing keys in config file:", paste(missing_keys, collapse = ", ")))
}


# Create output directory if it doesn't exist
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}


# Check input files exist
if (!file.exists(sample_loops_path) || !file.exists(random_loops_path)) {
  stop("One or both input files do not exist. Check your config file.")
}



# Set input files
sample_loops <- tryCatch({
  read.delim(sample_loops_path, header = FALSE)
}, error = function(e) {
  stop("Error reading sample_loops file: ", e$message)
})
random_loops <- tryCatch({
  read.delim(random_loops_path, header = FALSE)
}, error = function(e) {
  stop("Error reading random_loops file: ", e$message)
})


# Set headers for input files
headers <- c("chrA", "startA", "endA", "chrB", "startB", "endB", "loop_id", "strength_control", "strength_mutant")
colnames(sample_loops) <- headers
colnames(random_loops) <- headers


# Define list of loops
looplist <- list(sample_loops = sample_loops,
                 random_loops = random_loops)


# Add loop size
looplist <- lapply(looplist, function(x) {
  x %>% mutate(loop_size = startB - startA)
  })


# Remove loops < 1 Mb
looplist <- lapply(looplist, function(x) {
  x %>% filter(loop_size >= 1e6)
  })


# Calculate log2Fc
looplist <- lapply(looplist, function(x) {
  x$logFC <- log2((x$strength_mutant + 1) / (x$strength_control + 1))
  x
})


# Define the logFC null vector
filtered_random_loops <- looplist$random_loops
null_logFC <- filtered_random_loops$logFC


# Empirical p-value (bilateral)
filtered_sample_loops <- looplist$sample_loops
filtered_sample_loops$pval <- sapply(filtered_sample_loops$logFC, function(x) {
  mean(abs(null_logFC) >= abs(x))
})


# Adjusted p value Bonferroni correction
filtered_sample_loops$padj <- p.adjust(filtered_sample_loops$pval, method = "BH")
significance_threshold <- 0.05


# Volcano plot
pdf(file = file.path(output_directory, "volcano_plot.pdf"), width = 7, height = 5)
ggplot(filtered_sample_loops, aes(x = logFC, y = -log10(padj))) +
  geom_point(alpha = 1) +
  geom_hline(yintercept = -log10(significance_threshold), linetype = "dashed", color = "red") +
  geom_text_repel(
    data = subset(filtered_sample_loops, pval < significance_threshold),
    aes(label = loop_id),
    size = 3
  ) +
  theme_minimal() +
  labs(x = "log2 Fold Change",
       y = "-log10(p-value)") +
  ylim(0, 2.5)
dev.off()


# Output table
write.table(filtered_sample_loops,
            file = file.path(output_directory, "output_loopstrength.txt"),
            quote = FALSE,
            row.names = FALSE,
            sep = "\t",
            dec = ",")



cat("Analysis completed successfully. Results saved in:", output_directory, "\n")
