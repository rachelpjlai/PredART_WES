# R code to prepare data for the forest plot

library(data.table)
library(dplyr)

# Load curated file
curated_variants <- fread("final_significant_variants_fully_annotated.tsv")

# Load the full, CD4-adjusted PLINK results to get the Standard Error
full_plink_results <- fread("tb_iris_association_results_allchr_CD4.PHENOTYPE.glm.logistic.hybrid")

# Prepare a VariantID for joining
curated_variants <- curated_variants %>%
  mutate(VariantID = paste(gsub("chr", "", CHROM), POS, REF, ALT, sep = ":"))

full_plink_results <- full_plink_results %>%
  filter(TEST == "ADD") %>%
  mutate(VariantID = paste(gsub("chr", "", `#CHROM`), POS, REF, ALT, sep = ":")) %>%
  select(VariantID, SE = `LOG(OR)_SE`) # Select Standard Error of the Log(OR)

# Join the SE to the curated list
data_for_plot <- left_join(curated_variants, full_plink_results, by = "VariantID") %>%
  # Calculate the 95% Confidence Interval for the Odds Ratio
  mutate(
    LogOR = log(Odds_Ratio_ALT),
    CI_Lower = exp(LogOR - 1.96 * SE),
    CI_Upper = exp(LogOR + 1.96 * SE)
  )

# Save this final data table for plotting
fwrite(data_for_plot, "data_for_forest_plot.tsv", sep = "\t")