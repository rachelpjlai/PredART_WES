# File: create_final_summary.R 
# Description: Correctly joins PLINK2 association results with VCF data,
#              calculates genotype counts, maps genes, and filters the final table.

# --- Load Libraries ---
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("vcfR", quietly = TRUE)) install.packages("vcfR")

library(dplyr)
library(tidyr)
library(data.table)
library(vcfR)

# --- Configuration ---
plink_results_file <- "plink_results_ADD_test.tsv"
vcf_file <- "vep_annotated_deg_variants.nochr.vcf.gz"
pheno_file <- "rvtests_pheno.txt" 
deg_list_file <- "All_DEG_gene_symbol.txt"
output_significant_file <- "final_significant_variants_in_DEGs_with_counts.tsv"

# --- Load Data ---
cat("Loading data...\n")
pheno_data <- fread(pheno_file, header = TRUE, sep = "\t") %>%
  select(IID, PHENOTYPE = TB_IRIS_STATUS)

deg_list <- fread(deg_list_file, header = FALSE)$V1

plink_results <- fread(plink_results_file, header = TRUE) %>%
  mutate(VariantID = paste(gsub("chr", "", `#CHROM`), POS, REF, ALT, sep = ":"))

vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

# --- Step 1: Prepare a detailed map of all allele-specific variants from the VCF ---
cat("Preparing allele-specific variant map from VCF...\n")
fix_df_vcf_raw <- as_tibble(vcfR::getFIX(vcf, getINFO = TRUE))
fix_df_vcf_raw$Original_VCF_Row_Index <- 1:nrow(fix_df_vcf_raw)

vcf_allele_lookup <- fix_df_vcf_raw %>%
  mutate(CHROM_no_chr = gsub("chr", "", CHROM)) %>% 
  select(Original_VCF_Row_Index, CHROM_no_chr, POS, REF, ALT, INFO) %>%
  separate_rows(ALT, sep = ",") %>%
  group_by(Original_VCF_Row_Index) %>%
  mutate(Allele_Number_In_VCF = row_number()) %>% 
  ungroup() %>%
  mutate(VariantID = paste(CHROM_no_chr, POS, REF, ALT, sep = ":"))

# --- Step 2: Join PLINK results with the VCF variant map ---
cat("Joining PLINK results with VCF variant map...\n")
plink_results_mapped <- plink_results %>%
  inner_join(vcf_allele_lookup, by = "VariantID")

cat("Successfully mapped", nrow(plink_results_mapped), "PLINK results to VCF data.\n")

# --- Step 3: Extract Genotypes and Calculate Counts ---
cat("Extracting genotypes and calculating counts for", nrow(plink_results_mapped), "variants...\n")
gt_matrix_full <- vcfR::extract.gt(vcf, element = "GT")

classify_genotype <- function(gt_string, alt_allele_number) {
  if (is.na(gt_string) || gt_string %in% c("./.", ".|.")) return("Missing")
  alleles <- suppressWarnings(as.integer(unlist(strsplit(gt_string, "[/|]"))))
  if (any(is.na(alleles))) return("Missing")
  
  count_ref <- sum(alleles == 0)
  count_this_alt <- sum(alleles == alt_allele_number)
  
  if (count_ref == 2) return("RefHom_00")
  if (count_this_alt == 2) return("AltHom_11")
  if (count_ref == 1 && count_this_alt == 1) return("Het_01")
  return("Other")
}

genotype_counts_list <- list()
for (i in 1:nrow(plink_results_mapped)) {
  if (i %% 50000 == 0) cat("  Processing genotype counts for variant", i, "of", nrow(plink_results_mapped), "\n")
  
  vcf_row <- plink_results_mapped$Original_VCF_Row_Index[i]
  allele_num <- plink_results_mapped$Allele_Number_In_VCF[i]
  genotypes <- gt_matrix_full[vcf_row, ]
  
  counts_df <- data.frame(IID = names(genotypes), GT = genotypes) %>%
    inner_join(pheno_data, by = "IID") %>%
    rowwise() %>% 
    mutate(GT_Class = classify_genotype(GT, allele_num)) %>%
    ungroup()

  group_counts <- counts_df %>%
    filter(GT_Class %in% c("RefHom_00", "Het_01", "AltHom_11")) %>%
    group_by(PHENOTYPE) %>%
    count(GT_Class) %>%
    ungroup() %>%
    complete(PHENOTYPE, GT_Class, fill = list(n = 0))

  getCount <- function(df, pheno_val, class_val) {
    val <- df$n[df$PHENOTYPE == pheno_val & df$GT_Class == class_val]
    if (length(val) == 0) return(0)
    return(val)
  }
  
  genotype_counts_list[[i]] <- data.frame(
    VariantID = plink_results_mapped$VariantID[i],
    NonIRIS_RefHom_00 = getCount(group_counts, 0, "RefHom_00"),
    NonIRIS_Het_01    = getCount(group_counts, 0, "Het_01"),
    NonIRIS_AltHom_11 = getCount(group_counts, 0, "AltHom_11"),
    TBIRIS_RefHom_00  = getCount(group_counts, 1, "RefHom_00"),
    TBIRIS_Het_Count  = getCount(group_counts, 1, "Het_01"),
    TBIRIS_AltHom_11  = getCount(group_counts, 1, "AltHom_11")
  )
}
genotype_counts_df <- bind_rows(genotype_counts_list)

# --- Step 4: Final Merge, Annotation, and Filtering ---
cat("Finalizing and filtering summary table...\n")

get_gene_from_csq <- function(info_string) {
  if (is.na(info_string) || !grepl("CSQ=", info_string)) return(NA_character_)
  csq_data_all <- sub(".*CSQ=([^;]+).*", "\\1", info_string)
  first_csq_block <- strsplit(csq_data_all, ",")[[1]][1]
  fields <- strsplit(first_csq_block, "\\|")[[1]]
  if (length(fields) >= 4 && fields[4] != "") { return(fields[4]) } else { return(NA_character_) }
}

plink_results_mapped$Gene <- sapply(plink_results_mapped$INFO, get_gene_from_csq)
plink_results_mapped$IsInDEG <- plink_results_mapped$Gene %in% deg_list

final_table <- inner_join(plink_results_mapped, genotype_counts_df, by = "VariantID")

significant_variants_final <- final_table %>%
  filter(IsInDEG == TRUE) %>%
  filter(P <= 0.05) %>%
  arrange(P)

cat("Found", nrow(significant_variants_final), "nominally significant variant annotations in DEGs.\n")

# Select and rename columns, using the ".x" suffix that dplyr adds during joins
# to get the columns from the original PLINK results data frame.
final_output <- significant_variants_final %>%
  select(
    CHROM = `#CHROM`, 
    POS = POS.x, # Use POS.x from the left table
    PLINK_ID = ID,
    REF = REF.x, 
    ALT = ALT.x, 
    A1_Allele = A1,
    Gene, 
    IsInDEG,
    Odds_Ratio_PLINK = OR, 
    P_Value_PLINK = P,
    NonIRIS_RefHom_00, NonIRIS_Het_01, NonIRIS_AltHom_11,
    TBIRIS_RefHom_00, TBIRIS_Het_Count, TBIRIS_AltHom_11
  )


fwrite(final_output, output_significant_file, sep = "\t", row.names = FALSE, quote=FALSE)
cat("Final summary of significant variants in DEGs saved to:", output_significant_file, "\n")
print(head(final_output))

cat("Script finished.\n")