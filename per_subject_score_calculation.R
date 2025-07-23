# File: calculate_subject_total_burden.R
# Description: This script calculates multiple per-subject total genetic burden scores
#              based on different gene selection criteria (Burden p-value, SKAT-O
#              p-value, and effect size).

# --- Load Libraries ---
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org/")
if (!require("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org/")

library(dplyr)
library(data.table)

# --- Configuration ---
# File containing the per-subject, per-gene burden scores
burden_matrix_file <- "gene_burden_scores_skat_per_subject.tsv"

# File containing the comprehensive gene association results
association_results_file <- "gene_association_results_comprehensive.csv"

# Original phenotype file to get case/control status
pheno_file <- "tb_iris_pheno.fid0.tab.txt"

# Output file for the final per-subject scores
output_file <- "subject_total_burden_scores.csv"

# P-value threshold to define "significant" genes
significance_threshold <- 0.05
or_threshold_upper <- 1.2
or_threshold_lower <- 0.8


# --- 1. Load Association Results and Identify Gene Sets ---
cat("Loading association results from:", association_results_file, "\n")
if (!file.exists(association_results_file)) {
  stop("Association results file not found: ", association_results_file)
}
association_results <- fread(association_results_file)

# Check for required columns
required_cols <- c("Gene", "P.value.Burden", "P.value.SKAT.O", "OR.Burden")
if (!all(required_cols %in% names(association_results))) {
    stop("Association results file is missing required columns: ", paste(setdiff(required_cols, names(association_results)), collapse=", "))
}

# Set 1: Genes significant in the Firth's Burden test
genes_burden_sig <- association_results %>%
  filter(P.value.Burden <= significance_threshold) %>%
  pull(Gene)
cat(sprintf("Identified %d genes from Burden test (p <= %.2f).\n", length(genes_burden_sig), significance_threshold))

# Set 2: Genes from Burden test with notable effect size
genes_burden_effect <- association_results %>%
  filter(P.value.Burden <= significance_threshold & (OR.Burden > or_threshold_upper | OR.Burden < or_threshold_lower)) %>%
  pull(Gene)
cat(sprintf("Identified %d genes from Burden test with OR > %.1f or < %.1f.\n", length(genes_burden_effect), or_threshold_upper, or_threshold_lower))

# Set 3: Genes significant in the SKAT-O test
genes_skato_sig <- association_results %>%
  filter(P.value.SKAT.O <= significance_threshold) %>%
  pull(Gene)
cat(sprintf("Identified %d genes from SKAT-O test (p <= %.2f).\n", length(genes_skato_sig), significance_threshold))

# Set 4: Genes from SKAT-O test with notable Burden effect size
genes_skato_effect <- association_results %>%
  filter(P.value.SKAT.O <= significance_threshold & (OR.Burden > or_threshold_upper | OR.Burden < or_threshold_lower)) %>%
  pull(Gene)
cat(sprintf("Identified %d genes from SKAT-O test with Burden OR > %.1f or < %.1f.\n", length(genes_skato_effect), or_threshold_upper, or_threshold_lower))


# --- 2. Load the Per-Subject Burden Matrix ---
cat("Loading per-subject burden score matrix from:", burden_matrix_file, "\n")
if (!file.exists(burden_matrix_file)) {
  stop("Burden matrix file not found: ", burden_matrix_file)
}
burden_matrix <- fread(burden_matrix_file, header = TRUE)
setnames(burden_matrix, "V1", "IID")

# --- 3. Calculate Total Burden Scores ---
cat("Calculating total burden scores per subject...\n")

# Intersect defined gene sets with columns available in the matrix
all_gene_cols <- setdiff(names(burden_matrix), "IID")
genes_burden_sig_cols <- intersect(genes_burden_sig, all_gene_cols)
genes_burden_effect_cols <- intersect(genes_burden_effect, all_gene_cols)
genes_skato_sig_cols <- intersect(genes_skato_sig, all_gene_cols)
genes_skato_effect_cols <- intersect(genes_skato_effect, all_gene_cols)


# Calculate scores
subject_scores <- burden_matrix %>%
  mutate(
    Total_Burden_All_Genes = rowSums(select(., all_of(all_gene_cols)), na.rm = TRUE),
    Total_Burden_Burden_Sig = rowSums(select(., all_of(genes_burden_sig_cols)), na.rm = TRUE),
    Total_Burden_Burden_Effect = rowSums(select(., all_of(genes_burden_effect_cols)), na.rm = TRUE),
    Total_Burden_SKATO_Sig = rowSums(select(., all_of(genes_skato_sig_cols)), na.rm = TRUE),
    Total_Burden_SKATO_Effect = rowSums(select(., all_of(genes_skato_effect_cols)), na.rm = TRUE)
  ) %>%
  select(IID, Total_Burden_All_Genes, Total_Burden_Burden_Sig, Total_Burden_Burden_Effect, Total_Burden_SKATO_Sig, Total_Burden_SKATO_Effect)


# --- 4. Load Phenotype Data and Merge ---
cat("Loading phenotype data and merging with scores...\n")
if (!file.exists(pheno_file)) {
  stop("Phenotype file not found: ", pheno_file)
}
pheno_data <- fread(pheno_file)

pheno_data <- pheno_data %>%
  mutate(Status = ifelse(PHENOTYPE == 2, "TB-IRIS", "non-IRIS")) %>%
  select(IID, Status)

final_output <- left_join(pheno_data, subject_scores, by = "IID")

# --- 5. Save Final Output ---
cat("Saving final per-subject burden scores to:", output_file, "\n")
fwrite(final_output, file = output_file)

cat("Script finished successfully.\n")
print(head(final_output))
