# File: add_snpeff_annotations.R
# Description: Takes a list of significant variants and adds SnpEff functional
#              annotations by querying the main VCF file.

# --- Load Libraries ---
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("data.table", quietly = TRUE)) install.packages("data.table")
if (!require("stringr", quietly = TRUE)) install.packages("stringr")

library(dplyr)
library(data.table)
library(stringr)

# --- Configuration ---
# Input file: Final table of significant variants
significant_variants_file <- "final_significant_variants_in_DEGs_with_counts_OR.tsv"

# The VCF file that contains the SnpEff ANN annotations (with 'chr' prefix)
vcf_file <- "vep_annotated_deg_variants.vcf.gz"

# Final output file
output_file_annotated <- "final_significant_variants_fully_annotated.tsv"


# --- Step 1: Load significant variants and create a targets file for bcftools ---
cat("Loading significant variants from:", significant_variants_file, "\n")
if (!file.exists(significant_variants_file)) stop("Significant variants file not found:", significant_variants_file)

variants_to_annotate <- fread(significant_variants_file, header = TRUE)

if(nrow(variants_to_annotate) == 0) {
  stop("Input file contains no variants to annotate.")
}
cat("Preparing to annotate", nrow(variants_to_annotate), "variants...\n")

# Create a targets file (CHR POS REF ALT) for bcftools
# Ensures the 'chr' prefix is present to match the VCF
variants_to_annotate %>%
  mutate(CHROM_fixed = if_else(startsWith(CHROM, "chr"), CHROM, paste0("chr", CHROM))) %>%
  select(CHROM_fixed, POS, REF, ALT) %>%
  fwrite("temp_targets_for_annotation.tsv", sep = "\t", col.names = FALSE)


# --- Step 2: Use bcftools to extract SnpEff ANN annotations ---
cat("Extracting SnpEff annotations from VCF using bcftools...\n")
bcftools_cmd <- paste(
  "bcftools query",
  "--targets-file temp_targets_for_annotation.tsv",
  # Extract variant ID and the full ANN field from INFO
  "-f '%CHROM:%POS:%REF:%ALT\t%INFO/ANN\n'",
  vcf_file,
  "> temp_snpeff_annotations.tsv"
)
system(bcftools_cmd)

if (!file.exists("temp_snpeff_annotations.tsv") || file.info("temp_snpeff_annotations.tsv")$size == 0) {
  stop("bcftools query failed to extract any annotations.")
}

# --- Step 3: Load and process the annotations ---
cat("Parsing annotations and merging with results...\n")
snpeff_annotations <- fread("temp_snpeff_annotations.tsv", header = FALSE, 
                            col.names = c("VariantID_VCF", "ANN_string"))

# Helper function to get the first, primary SnpEff annotation term
get_first_snpeff_annotation <- function(ann_string) {
  if (is.na(ann_string) || ann_string == "" || ann_string == ".") return(NA_character_)
  first_ann_block <- str_split(ann_string, ",", simplify = FALSE)[[1]][1]
  fields <- str_split(first_ann_block, "\\|", simplify = FALSE)[[1]]
  if (length(fields) >= 2 && fields[2] != "") {
    # Return the main consequence term, e.g., "missense_variant"
    return(str_split(fields[2], "&", simplify = TRUE)[1,1])
  } else {
    return(NA_character_)
  }
}

# Process the annotations
snpeff_annotations <- snpeff_annotations %>%
  mutate(SnpEff_Annotation = sapply(ANN_string, get_first_snpeff_annotation)) %>%
  # Create a matching VariantID (without 'chr') for joining with the results table
  mutate(VariantID = gsub("^chr", "", VariantID_VCF)) %>%
  select(VariantID, SnpEff_Annotation)


# --- Step 4: Join annotations with main results and save ---
# Create a matching VariantID in your main results table
variants_to_annotate <- variants_to_annotate %>%
  mutate(VariantID = paste(gsub("chr", "", CHROM), POS, REF, ALT, sep = ":"))

# Join the data
final_annotated_table <- left_join(variants_to_annotate, snpeff_annotations, by = "VariantID") %>%
  # Reorder columns to place the new annotation near the gene
  select(CHROM, POS, PLINK_ID, REF, ALT, Gene, SnpEff_Annotation, everything(), -VariantID) %>%
  arrange(P_Value_PLINK)

# Save the final table
fwrite(final_annotated_table, output_file_annotated, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Final annotated summary table saved to:", output_file_annotated, "\n")
print(head(final_annotated_table))

# --- Clean up temporary files ---
file.remove("temp_targets_for_annotation.tsv", "temp_snpeff_annotations.tsv")

cat("Script finished.\n")