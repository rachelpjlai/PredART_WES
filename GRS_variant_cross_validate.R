# File: run_variant_grs_cv.R
# Description: This script performs both an in-sample and a rigorous 10-fold
#              cross-validation to evaluate the predictive performance (AUC) of a
#              variant-level Genetic Risk Score (GRS). It exports detailed
#              metrics and generates separate, publication-quality ROC plots.
# Version 10: Final legend ordering adjustment.

# --- 1. Load Libraries ---
# if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
# if (!require("data.table", quietly = TRUE)) install.packages("data.table")
# if (!require("pROC", quietly = TRUE)) install.packages("pROC")
# if (!require("caret", quietly = TRUE)) install.packages("caret")
# if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# if (!require("vcfR", quietly = TRUE)) install.packages("vcfR")
# if (!require("knitr", quietly = TRUE)) install.packages("knitr")


library(dplyr)
library(data.table)
library(pROC)
library(caret)
library(ggplot2)
library(vcfR)
library(knitr)

# --- 2. Configuration ---
vcf_file <- "vep_annotated_deg_variants.nochr.vcf.gz"
pheno_file <- "tb_iris_pheno.fid0.tab.txt"
output_prefix <- "variant_grs"
num_top_variants <- 5
num_folds <- 10

# --- 3. Setup ---
cat("--- Initial Setup ---\n")
plink_path <- Sys.which("plink2")
if (plink_path == "") {
  stop("FATAL ERROR: plink2 executable not found. Please install PLINK 2.0.")
}
pheno_data <- fread(pheno_file) %>% select(IID, PHENOTYPE)
sex_info_for_plink <- fread(pheno_file) %>%
  mutate(SEX_PLINK = case_when(SEX == "M" ~ 1, SEX == "F" ~ 2, TRUE ~ 0)) %>%
  select(FID, IID, SEX_PLINK)
temp_sex_file <- "temp_sex_info.txt"
fwrite(sex_info_for_plink, temp_sex_file, sep = "\t", col.names = FALSE)

cat("Loading full VCF file into memory...\n")
vcf_full <- read.vcfR(vcf_file, verbose = FALSE)
gt_full <- extract.gt(vcf_full, element = "GT")
fix_full <- as.data.frame(getFIX(vcf_full), stringsAsFactors = FALSE) %>%
            mutate(VariantID = paste(CHROM, POS, REF, ALT, sep = ":"))

# Helper function to calculate GRS
calculate_grs <- function(variant_set, subject_ids) {
    if (nrow(variant_set) == 0) return(rep(0, length(subject_ids)))
    variant_ids_to_extract <- paste(variant_set$`#CHROM`, variant_set$POS, variant_set$REF, variant_set$ALT, sep = ":")
    vcf_indices <- which(fix_full$VariantID %in% variant_ids_to_extract)
    if (length(vcf_indices) == 0) return(rep(0, length(subject_ids)))
    
    gt_subset <- gt_full[vcf_indices, subject_ids, drop = FALSE]
    count_alt <- function(gt) { sapply(gt, function(x) if (is.na(x) || x == "./.") 0 else sum(as.integer(unlist(strsplit(x, "[/|]"))))) }
    allele_counts <- t(apply(gt_subset, 1, count_alt))
    
    # Ensure order matches for matrix multiplication
    ordered_variants <- variant_set[match(fix_full$VariantID[vcf_indices], variant_ids_to_extract),]
    weights <- log(ordered_variants$OR)
    grs <- as.vector(weights %*% allele_counts)
    return(grs)
}

# --- 4. In-Sample Analysis (for comparison) ---
cat("\n--- Performing In-Sample Analysis (Full Dataset) ---\n")
plink_full_prefix <- paste0(output_prefix, "_full_dataset")
plink_cmd_full <- paste(
    "plink2", "--vcf", vcf_file, "--pheno", pheno_file, "--update-sex", temp_sex_file, 
    "--split-par hg38", "--glm allow-no-covars firth-fallback", "--out", plink_full_prefix
)
system(plink_cmd_full)
full_assoc_results <- fread(paste0(plink_full_prefix, ".PHENOTYPE.glm.logistic.hybrid")) %>% filter(TEST == "ADD")

top_risk_full <- full_assoc_results %>% filter(OR > 1) %>% arrange(P) %>% head(num_top_variants)
top_protective_full <- full_assoc_results %>% filter(OR < 1) %>% arrange(P) %>% head(num_top_variants)

# Export the top variants from the full analysis
top_variants_export <- bind_rows(top_risk_full, top_protective_full) %>% 
  mutate(LOG_OR = log(OR)) %>%
  select(ID, P, OR, LOG_OR)
fwrite(top_variants_export, paste0(output_prefix, "_top_in_sample_variants.csv"))
cat("Exported top in-sample variants to", paste0(output_prefix, "_top_in_sample_variants.csv"), "\n")

in_sample_scores <- pheno_data %>%
    mutate(
        GRS_Risk = calculate_grs(top_risk_full, IID),
        GRS_Protective = calculate_grs(top_protective_full, IID)
    ) %>%
    mutate(GRS_Combined = GRS_Risk + GRS_Protective)

# --- 5. Cross-Validation ---
cat("\n--- Starting 10-Fold Cross-Validation ---\n")
set.seed(123)
cv_folds <- createFolds(factor(pheno_data$PHENOTYPE), k = num_folds, list = TRUE, returnTrain = FALSE)
out_of_sample_predictions <- data.frame()

for (i in 1:num_folds) {
  cat(sprintf("  Processing Fold %d/%d...\n", i, num_folds))
  test_indices <- cv_folds[[i]]
  train_samples <- pheno_data[-test_indices, ]
  fwrite(train_samples[, .(FID = 0, IID)], "temp_train_samples.txt", sep = "\t", col.names = FALSE)
  
  plink_fold_prefix <- paste0(output_prefix, "_fold", i)
  plink_cmd_fold <- paste(
    "plink2", "--vcf", vcf_file, "--keep temp_train_samples.txt", "--pheno", pheno_file, 
    "--update-sex", temp_sex_file, "--split-par hg38", 
    "--glm allow-no-covars firth-fallback", "--out", plink_fold_prefix
  )
  system(plink_cmd_fold)
  
  plink_results_file <- paste0(plink_fold_prefix, ".PHENOTYPE.glm.logistic.hybrid")
  if (!file.exists(plink_results_file)) { next }
  fold_assoc_results <- fread(plink_results_file) %>% filter(TEST == "ADD")
  
  top_risk_fold <- fold_assoc_results %>% filter(OR > 1) %>% arrange(P) %>% head(num_top_variants)
  top_protective_fold <- fold_assoc_results %>% filter(OR < 1) %>% arrange(P) %>% head(num_top_variants)
  
  test_data_scores <- pheno_data[test_indices, ] %>%
    mutate(
      GRS_Risk = calculate_grs(top_risk_fold, IID),
      GRS_Protective = calculate_grs(top_protective_fold, IID)
    ) %>%
    mutate(GRS_Combined = GRS_Risk + GRS_Protective) %>%
    select(IID, PHENOTYPE, GRS_Risk, GRS_Protective, GRS_Combined)
    
  out_of_sample_predictions <- bind_rows(out_of_sample_predictions, test_data_scores)
}

# --- 6. Final ROC Analysis & Metric Export ---
cat("\n--- Performing Final ROC Analysis and Exporting Metrics ---\n")
all_metrics <- data.frame()
add_metrics <- function(roc_obj, model_name, type) {
    coords_best <- coords(roc_obj, "best", ret=c("threshold", "specificity", "sensitivity"), best.method="closest.topleft")
    data.frame(
        Model = model_name,
        Validation_Type = type,
        AUC = as.numeric(auc(roc_obj)),
        CI_Lower = ci.auc(roc_obj)[1],
        CI_Upper = ci.auc(roc_obj)[3],
        Best_Threshold = coords_best$threshold,
        Sensitivity = coords_best$sensitivity,
        Specificity = coords_best$specificity
    )
}

# Process both in-sample and cross-validated results
all_scores <- list(
    "In-Sample" = in_sample_scores,
    "Cross-Validated" = out_of_sample_predictions
)
roc_list <- list()

for (val_type in names(all_scores)) {
    df <- all_scores[[val_type]]
    if(nrow(df) == 0) next
    df$Status <- as.factor(ifelse(df$PHENOTYPE == 2, "TB_IRIS", "non_IRIS"))
    
    roc_risk <- roc(response = df$Status, predictor = df$GRS_Risk, levels = c("non_IRIS", "TB_IRIS"))
    roc_prot <- roc(response = df$Status, predictor = df$GRS_Protective, levels = c("non_IRIS", "TB_IRIS"))
    roc_comb <- roc(response = df$Status, predictor = df$GRS_Combined, levels = c("non_IRIS", "TB_IRIS"))

    all_metrics <- bind_rows(all_metrics, add_metrics(roc_risk, "Risk", val_type))
    all_metrics <- bind_rows(all_metrics, add_metrics(roc_prot, "Protective", val_type))
    all_metrics <- bind_rows(all_metrics, add_metrics(roc_comb, "Combined", val_type))
    
    roc_list[[val_type]] <- list(Risk=roc_risk, Protective=roc_prot, Combined=roc_comb)
}

fwrite(all_metrics, paste0(output_prefix, "_auc_metrics.csv"))
cat("Exported all AUC metrics to", paste0(output_prefix, "_auc_metrics.csv"), "\n")
print(knitr::kable(all_metrics, digits=3))

# --- 7. Plotting ---
cat("\nGenerating separate ROC plots...\n")

# A helper function to create a single ROC plot to avoid code repetition
create_roc_plot <- function(roc_data_list, plot_title, filename){
    
    # Create the dynamic labels for the legend
    risk_label <- sprintf("Risk (AUC=%.2f)", auc(roc_data_list$Risk))
    prot_label <- sprintf("Protective (AUC=%.2f)", auc(roc_data_list$Protective))
    comb_label <- sprintf("Combined (AUC=%.2f)", auc(roc_data_list$Combined))
    
    # Combine the data for plotting
    plot_data <- bind_rows(
        data.frame(ggroc(roc_data_list$Risk)$data, Model = risk_label),
        data.frame(ggroc(roc_data_list$Protective)$data, Model = prot_label),
        data.frame(ggroc(roc_data_list$Combined)$data, Model = comb_label)
    ) %>%
    # **FIXED**: Convert Model to a factor with the desired order for the legend
    mutate(Model = factor(Model, levels = c(risk_label, prot_label, comb_label)))
    
    
    # Create a named vector for the colors correctly.
    color_values <- c("firebrick", "darkblue", "darkgreen")
    names(color_values) <- c(risk_label, prot_label, comb_label)
    
    roc_plot <- ggplot(plot_data, aes(x = specificity, y = sensitivity, color = Model)) +
        geom_line(linewidth = 1.2) +
        geom_abline(intercept = 1, slope = 1, linetype = "dotted", color = "grey50") +
        scale_x_reverse(name = "Specificity") +
        scale_y_continuous(name = "Sensitivity") +
        # Use the correctly constructed named vector
        scale_color_manual(values = color_values, name = "GRS Model") +
        labs(title = plot_title) +
        coord_equal() +
        theme_bw(base_size = 20) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5), 
            legend.position = c(0.95, 0.05),
            legend.justification = c("right", "bottom"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6),
            legend.background = element_blank(),
            legend.key = element_blank()
        ) +
        guides(color = guide_legend(ncol = 1)) # Stack legend items vertically
    
    ggsave(filename, plot = roc_plot, width = 9, height = 9, device = "png", dpi = 300)
    cat("Saved plot to", filename, "\n")
    return(roc_plot)
}

# Create and save the two plots
if("In-Sample" %in% names(roc_list)){
    print(create_roc_plot(
        roc_list[["In-Sample"]], 
        "In-Sample Performance of GRS", 
        paste0(output_prefix, "_in_sample_roc_plot.png")
    ))
}

if("Cross-Validated" %in% names(roc_list)){
    print(create_roc_plot(
        roc_list[["Cross-Validated"]], 
        "Cross-Validated Performance of GRS", 
        paste0(output_prefix, "_cross_validated_roc_plot.png")
    ))
}


# --- Cleanup ---
cat("Cleaning up temporary files...\n")
file.remove(temp_sex_file)
file.remove("temp_train_samples.txt")
file.remove(list.files(pattern = paste0(output_prefix, "_fold*|", output_prefix, "_full_dataset*")))

cat("Script finished.\n")
