# Gene-Level Burden, SKAT-O, and Firth's Regression Analysis for TB-IRIS
# This script performs comprehensive gene-based association tests. It:
# 1. Runs a Firth's logistic regression burden test to get Odds Ratios, CIs, and SEs.
# 2. Runs the optimal SKAT-O test to get a robust p-value.
# 3. Creates a per-subject, per-gene burden score matrix.

# --- 1. Installation of Required Packages ---
# Uncomment the following lines to install the necessary packages if you haven't already.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("vcfR")
#
# install.packages(c("dplyr", "tidyr", "purrr", "stringr", "SKAT", "logistf"))

# --- 2. Load Libraries ---
library(vcfR)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(SKAT)
library(logistf)

# --- 3. Define File Paths ---
# Please ensure these file paths are correct.
vcf_file <- "vep_annotated_deg_variants.nochr.vcf.gz"
pheno_covar_file <- "tb_iris_pheno.fid0.tab.txt"
deg_list_file <- "All_DEG_gene_symbol.txt" # File with the list of DEGs

# --- 4. Load and Prepare Phenotype and Covariate Data ---
cat("Loading and preparing phenotype data...\n")
pheno_data <- read.table(pheno_covar_file, header = TRUE, stringsAsFactors = FALSE)

# Prepare data for both logistf (factors) and SKAT (numeric)
pheno_data <- pheno_data %>%
  mutate(
    TB_IRIS_status = PHENOTYPE - 1,
    SEX_factor = as.factor(SEX), # For logistf
    SEX_numeric = as.numeric(as.factor(SEX)) - 1 # For SKAT
  ) %>%
  select(IID, TB_IRIS_status, SEX_factor, SEX_numeric, CD4)

# --- 5. Load VCF and Extract Genotype Data ---
cat("Loading VCF file...\n")
vcf <- read.vcfR(vcf_file)

# Extract genotype (GT) matrix
gt_matrix <- extract.gt(vcf, element = "GT")

# Ensure samples in VCF and phenotype data match
common_samples <- intersect(colnames(gt_matrix), pheno_data$IID)
if (length(common_samples) == 0) {
  stop("No common samples found between VCF and phenotype file. Please check sample IDs (IID).")
}
gt_matrix <- gt_matrix[, common_samples]
pheno_data <- pheno_data %>% filter(IID %in% common_samples)
pheno_data <- pheno_data[match(colnames(gt_matrix), pheno_data$IID), ]

# --- 6. Create Numeric Allele Count Matrix ---
cat("Creating numeric allele count matrix from genotypes...\n")
allele_counts <- matrix(0, nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
allele_counts[gt_matrix == "0/1" | gt_matrix == "1/0"] <- 1
allele_counts[gt_matrix == "1/1"] <- 2
allele_counts[is.na(gt_matrix) | gt_matrix == "./."] <- NA # SKAT handles missing genotypes
colnames(allele_counts) <- colnames(gt_matrix)
allele_counts_transposed <- t(allele_counts) # samples x variants for SKAT

# --- 7. Map Variants to Genes using VEP Annotation ---
cat("Mapping variants to genes...\n")
gene_symbol_index <- 4 # Assume gene symbol is the 4th field in CSQ
info_field <- getINFO(vcf)
csq_annotations <- str_extract(info_field, "CSQ=([^;]+)")
csq_annotations <- gsub("CSQ=", "", csq_annotations)

gene_to_variants_map <- list()
for (i in seq_along(csq_annotations)) {
  if (is.na(csq_annotations[i])) next
  annotations <- str_split(csq_annotations[i], ",")[[1]]
  variant_genes <- unique(na.omit(map_chr(annotations, ~ {
    fields <- str_split(.x, "\\|")[[1]]
    if (length(fields) >= gene_symbol_index) return(fields[gene_symbol_index])
    return(NA)
  })))
  for (gene in variant_genes) {
    if (gene != "" && !is.na(gene) && nchar(gene) > 1) {
      gene_to_variants_map[[gene]] <- c(gene_to_variants_map[[gene]], i)
    }
  }
}
gene_to_variants_map <- lapply(gene_to_variants_map, unique)

# --- 8. Filter Genes against DEG List ---
cat("Filtering genes against the provided DEG list...\n")
deg_list <- read.table(deg_list_file, header = FALSE, stringsAsFactors = FALSE)$V1
genes_from_vcf <- names(gene_to_variants_map)
retained_genes <- intersect(genes_from_vcf, deg_list)
gene_to_variants_map <- gene_to_variants_map[retained_genes]
unique_genes <- names(gene_to_variants_map)

cat(sprintf("%d genes were retained for analysis after filtering.\n", length(unique_genes)))

# --- 9. Prepare SKAT Null Model (accounts for covariates) ---
cat("Fitting SKAT Null Model...\n")
X <- model.matrix(TB_IRIS_status ~ SEX_numeric + CD4, data = pheno_data)
obj_null <- SKAT_Null_Model(pheno_data$TB_IRIS_status ~ X, out_type = "D")

# --- 10. Perform Per-Gene Association Tests ---
cat("Performing per-gene tests...\n")
results_list <- list()
burden_score_matrix <- matrix(NA,
                              nrow = length(common_samples),
                              ncol = length(unique_genes),
                              dimnames = list(pheno_data$IID, unique_genes))

pb <- txtProgressBar(min = 0, max = length(unique_genes), style = 3)
for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  variant_indices <- gene_to_variants_map[[gene]]

  # --- Burden Score Calculation and Firth's Test ---
  burden_alleles <- allele_counts[variant_indices, , drop = FALSE]
  burden_score <- colSums(burden_alleles, na.rm = TRUE)
  burden_score_matrix[, i] <- burden_score
  
  firth_fit <- NULL
  if(length(unique(burden_score)) > 1){
      analysis_df <- pheno_data %>% mutate(burden_score = burden_score)
      firth_fit <- try(logistf(TB_IRIS_status ~ burden_score + SEX_factor + CD4, data = analysis_df))
  }
  
  # --- SKAT-O Test ---
  Z <- allele_counts_transposed[, variant_indices, drop = FALSE]
  if(any(is.na(Z))){
      for(col_idx in 1:ncol(Z)){
          missing_rows <- is.na(Z[, col_idx])
          if(any(missing_rows)){
              mean_val <- mean(Z[!missing_rows, col_idx], na.rm=TRUE)
              Z[missing_rows, col_idx] <- mean_val
          }
      }
  }

  skato_fit <- NULL
  if (ncol(Z) >= 2 && sum(colSums(Z)) > 0) {
      skato_fit <- SKAT(Z, obj_null, method = "SKATO")
  }
  
  # --- Combine and Store Results ---
  res_row <- data.frame(Gene = gene, N.Variants = length(variant_indices))
  
  # Extract Firth results
  if(!is.null(firth_fit) && !inherits(firth_fit, "try-error")){
    burden_idx <- which(names(firth_fit$coefficients) == 'burden_score')
    if(length(burden_idx) > 0){
        res_row$P.value.Burden <- firth_fit$prob[burden_idx]
        res_row$OR.Burden <- exp(firth_fit$coefficients[burden_idx])
        res_row$CI.Lower.Burden <- exp(firth_fit$ci.lower[burden_idx])
        res_row$CI.Upper.Burden <- exp(firth_fit$ci.upper[burden_idx])
        res_row$SE.Burden <- sqrt(diag(firth_fit$var))[burden_idx]
    }
  }

  # Extract SKATO results
  if(!is.null(skato_fit)){
      res_row$P.value.SKAT.O <- skato_fit$p.value
  }

  results_list[[gene]] <- res_row
  
  setTxtProgressBar(pb, i)
}
close(pb)

final_results <- bind_rows(results_list)

# --- 11. Adjust P-values and Finalize ---
cat("\nAdjusting p-values using Benjamini-Hochberg (FDR)...\n")
if(nrow(final_results) > 0){
    if("P.value.Burden" %in% names(final_results)){
        final_results$FDR.Burden <- p.adjust(final_results$P.value.Burden, method = "BH")
    }
    if("P.value.SKAT.O" %in% names(final_results)){
        final_results$FDR.SKAT.O <- p.adjust(final_results$P.value.SKAT.O, method = "BH")
    }
    final_results <- final_results %>% arrange(P.value.SKAT.O)
}

# --- 12. Save Outputs ---
output_csv <- "gene_association_results_comprehensive.csv"
write.csv(final_results, output_csv, row.names = FALSE)
cat(sprintf("\nAnalysis complete. Comprehensive results saved to %s\n", output_csv))
print(head(final_results))

output_tsv <- "gene_burden_scores_skat_per_subject.tsv"
burden_df <- as.data.frame(burden_score_matrix)
write.table(burden_df, file = output_tsv, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("Per-gene burden scores saved to %s\n", output_tsv))
