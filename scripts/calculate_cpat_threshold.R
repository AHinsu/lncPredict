#!/usr/bin/env Rscript

################################################################################
# calculate_cpat_threshold.R - Calculate optimal CPAT threshold using ROC
# 
# This script calculates the optimal coding probability threshold for CPAT
# models using 10-fold cross-validation and ROC curve analysis. The threshold 
# is determined by finding the point where sensitivity and specificity curves 
# intersect, following the methodology used by CPAT developers in their 
# original publication.
#
# This implementation follows the CPAT original 10-fold cross-validation
# approach where:
# 1. Data is split into 10 folds
# 2. For each fold, a logistic regression model is trained on 90% of data
# 3. Predictions are made on the held-out 10% using predict(..., type="response")
# 4. All predictions are combined for ROC analysis
#
# Requirements:
#   - R with ROCR package installed
#   - Install via conda: conda install -c conda-forge r-base r-rocr
#   - Or install manually: install.packages("ROCR")
#
# Based on CPAT documentation and original CPAT publication methodology:
# https://cpat.readthedocs.io/en/latest/#how-to-choose-cutoff
#
# Usage: Rscript calculate_cpat_threshold.R <feature_file> <output_prefix>
#
# Arguments:
#   feature_file:    Feature file created by make_logitModel.py containing:
#                    ID, mRNA, ORF, Fickett, Hexamer, Label columns
#                    Label: 1 = coding, 0 = noncoding
#   output_prefix:   Prefix for output files (optional, default: 'threshold')
#
# The script generates:
#   - Optimal threshold value via 10-fold cross-validation
#   - ROC curve plot
#   - Sensitivity vs Specificity plot
#   - Area Under Curve (AUC)
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  if (!require("ROCR", quietly = TRUE)) {
    cat("\n")
    cat("ERROR: ROCR package is not installed.\n")
    cat("\n")
    cat("Please install ROCR using conda:\n")
    cat("  conda install -c conda-forge r-rocr\n")
    cat("\n")
    cat("Or install it manually in R:\n")
    cat("  install.packages('ROCR')\n")
    cat("\n")
    quit(status = 1)
  }
  library(ROCR)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  cat("Usage: Rscript calculate_cpat_threshold.R <feature_file> [output_prefix]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  feature_file:    Feature file created by make_logitModel.py\n")
  cat("                   Must contain: ID, mRNA, ORF, Fickett, Hexamer, Label\n")
  cat("                   Label: 1 = coding, 0 = noncoding\n")
  cat("  output_prefix:   Prefix for output files (optional, default: 'threshold')\n")
  cat("\n")
  cat("Outputs:\n")
  cat("  <prefix>_threshold.txt: Optimal threshold value\n")
  cat("  <prefix>_performance_plots.png: Four-panel performance visualization\n")
  cat("  <prefix>_stats.txt:     Detailed statistics\n")
  quit(status = 1)
}

feature_file <- args[1]
output_prefix <- if (length(args) >= 2) args[2] else "threshold"

cat("\n")
cat("================================================================================\n")
cat("CPAT Threshold Calculation using 10-Fold Cross-Validation\n")
cat("================================================================================\n")
cat("\n")

# Check input file
if (!file.exists(feature_file)) {
  stop("ERROR: Feature file not found: ", feature_file)
}

# Read feature data
cat("Reading feature data...\n")
data <- read.table(feature_file, header = TRUE, sep = "\t", row.names = 1)

cat(sprintf("  Total sequences: %d\n", nrow(data)))
cat("\n")

# Check required columns
required_cols <- c("mRNA", "ORF", "Fickett", "Hexamer", "Label")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("ERROR: Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Check if we have data
if (nrow(data) == 0) {
  stop("ERROR: No data found in feature file")
}

# Count coding and non-coding
n_coding <- sum(data$Label == 1)
n_noncoding <- sum(data$Label == 0)
cat(sprintf("  Coding sequences (Label=1): %d\n", n_coding))
cat(sprintf("  Non-coding sequences (Label=0): %d\n", n_noncoding))
cat("\n")

if (n_coding == 0 || n_noncoding == 0) {
  stop("ERROR: Need both coding and non-coding sequences for training")
}

# ============================================================================
# 10-Fold Cross-Validation (following CPAT original methodology)
# ============================================================================

cat("Performing 10-fold cross-validation...\n")
cat("  This matches the CPAT original publication methodology\n")
cat("\n")

# Total number of sequences
total_n <- nrow(data)

# Create 10 folds
n_folds <- 10
fold_size <- ceiling(total_n / n_folds)

# Initialize vectors to store all predictions and labels
all_predictions <- numeric(total_n)
all_labels <- numeric(total_n)
all_ids <- rownames(data)

# Perform 10-fold cross-validation
for (fold in 1:n_folds) {
  # Define test indices for this fold
  test_start <- (fold - 1) * fold_size + 1
  test_end <- min(fold * fold_size, total_n)
  test_indices <- test_start:test_end
  
  # Training indices are all others
  train_indices <- setdiff(1:total_n, test_indices)
  
  cat(sprintf("  Fold %d: Training on %d sequences, testing on %d sequences\n", 
              fold, length(train_indices), length(test_indices)))
  
  # Extract training data (avoid using attach())
  train_data <- data[train_indices, ]
  test_data <- data[test_indices, ]
  
  # Check for NA values in training data
  if (any(is.na(train_data))) {
    cat(sprintf("    WARNING: Fold %d has NA values in training data, removing them\n", fold))
    train_data <- na.omit(train_data)
  }
  
  # Train logistic regression model with error handling
  tryCatch({
    mylogit <- glm(Label ~ mRNA + ORF + Fickett + Hexamer, 
                   data = train_data,
                   family = binomial(link = "logit"))
    
    # Predict probabilities on test data
    # IMPORTANT: type="response" gives probabilities (0-1 range)
    test_predictions <- predict(mylogit, newdata = test_data, type = "response")
    
    # Store predictions and labels
    all_predictions[test_indices] <- test_predictions
    all_labels[test_indices] <- test_data$Label
    
  }, error = function(e) {
    cat(sprintf("    ERROR: Failed to fit model for fold %d: %s\n", fold, e$message))
    cat(sprintf("    Using fallback: mean coding probability from other folds\n"))
    # Use mean of already computed predictions as fallback
    fallback_prob <- if (fold > 1) mean(all_predictions[1:(test_start-1)], na.rm = TRUE) else 0.5
    all_predictions[test_indices] <- fallback_prob
    all_labels[test_indices] <- test_data$Label
  })
}

cat("\n")
cat("Cross-validation complete!\n")
cat(sprintf("  Total predictions: %d\n", length(all_predictions)))
cat("\n")

# Use the cross-validated predictions for ROC analysis
scores <- all_predictions
labels <- all_labels

# Calculate ROC curve using ROCR
cat("Calculating ROC curve using ROCR...\n")
pred_obj <- prediction(scores, labels)

# Calculate AUC
perf_auc <- performance(pred_obj, "auc")
auc_value <- perf_auc@y.values[[1]]
cat(sprintf("  AUC: %.4f\n", auc_value))
cat("\n")

# Validate AUC is reasonable
if (auc_value < 0.5) {
  cat("WARNING: AUC < 0.5 suggests model is performing worse than random!\n")
  cat("         This may indicate a problem with the training data or features.\n")
  cat("\n")
}

# Calculate sensitivity and specificity at different cutoffs
perf_sens <- performance(pred_obj, measure = "sens")
perf_spec <- performance(pred_obj, measure = "spec")

# Extract cutoffs, sensitivity, and specificity values
cutoffs <- perf_sens@x.values[[1]]
sensitivity <- perf_sens@y.values[[1]]
specificity <- perf_spec@y.values[[1]]

# Find the optimal threshold where sensitivity and specificity are closest
# This follows the CPAT developers' methodology
# Find the point where abs(sensitivity - specificity) is minimized
diff_sens_spec <- abs(sensitivity - specificity)
optimal_idx <- which.min(diff_sens_spec)

optimal_threshold <- cutoffs[optimal_idx]
optimal_sensitivity <- sensitivity[optimal_idx]
optimal_specificity <- specificity[optimal_idx]

cat("Optimal threshold calculation:\n")
cat(sprintf("  Method: Intersection of Sensitivity and Specificity curves (CPAT original methodology)\n"))
cat(sprintf("  Optimal threshold: %.4f\n", optimal_threshold))
cat(sprintf("  Sensitivity: %.4f (%.2f%%)\n", optimal_sensitivity, optimal_sensitivity * 100))
cat(sprintf("  Specificity: %.4f (%.2f%%)\n", optimal_specificity, optimal_specificity * 100))
cat(sprintf("  Difference: %.4f\n", abs(optimal_sensitivity - optimal_specificity)))
cat("\n")

# Calculate statistics at optimal threshold
coding_scores <- scores[labels == 1]
noncoding_scores <- scores[labels == 0]

coding_correctly_classified <- sum(coding_scores >= optimal_threshold)
noncoding_correctly_classified <- sum(noncoding_scores < optimal_threshold)
coding_accuracy <- coding_correctly_classified / length(coding_scores)
noncoding_accuracy <- noncoding_correctly_classified / length(noncoding_scores)
overall_accuracy <- (coding_correctly_classified + noncoding_correctly_classified) / length(scores)

cat("Classification performance at optimal threshold:\n")
cat(sprintf("  Coding accuracy: %.2f%% (%d/%d)\n", 
            coding_accuracy * 100, coding_correctly_classified, length(coding_scores)))
cat(sprintf("  Non-coding accuracy: %.2f%% (%d/%d)\n", 
            noncoding_accuracy * 100, noncoding_correctly_classified, length(noncoding_scores)))
cat(sprintf("  Overall accuracy: %.2f%%\n", overall_accuracy * 100))
cat("\n")

# Save plots
cat("Generating plots...\n")
png_file <- paste0(output_prefix, "_performance_plots.png")
png(png_file, width = 1200, height = 800, res = 120)

par(mfrow = c(2, 2), mar = c(5, 4, 3, 2), cex.axis = 1.2, cex.lab = 1.2)

# Plot 1: ROC curve
perf_roc <- performance(pred_obj, "tpr", "fpr")
plot(perf_roc, 
     col = "blue", 
     lwd = 2,
     xlab = "1-Specificity (FPR)", 
     ylab = "Sensitivity (TPR)",
     main = "ROC Curve",
     cex.axis = 1.2,
     cex.lab = 1.2)
abline(a = 0, b = 1, lty = "dashed", col = "gray50")
text(0.6, 0.2, sprintf("AUC = %.4f", auc_value), cex = 1.2)

# Plot 2: Precision-Recall curve
perf_pr <- performance(pred_obj, measure = "prec", x.measure = "rec")
plot(perf_pr,
     col = "blue",
     lwd = 2,
     xlab = "Recall (TPR)",
     ylab = "Precision (PPV)",
     main = "Precision-Recall Curve",
     xlim = c(0, 1),
     ylim = c(0, 1),
     cex.axis = 1.2,
     cex.lab = 1.2)

# Plot 3: Accuracy vs Cutoff
perf_acc <- performance(pred_obj, "acc")
plot(perf_acc,
     col = "blue",
     lwd = 2,
     xlab = "Coding probability cutoff",
     ylab = "Accuracy",
     main = "Accuracy vs Cutoff",
     ylim = c(0.5, 1),
     cex.axis = 1.2,
     cex.lab = 1.2)
abline(h = 1, lty = "dashed", lwd = 0.5)
abline(v = optimal_threshold, lty = "dashed", col = "red", lwd = 1.5)

# Plot 4: Sensitivity vs Specificity (CPAT original methodology)
plot(cutoffs, sensitivity,
     type = "l",
     col = "blue",
     lwd = 2,
     xlab = "Coding Probability Cutoff",
     ylab = "Performance",
     main = "Sensitivity vs Specificity",
     ylim = c(0, 1),
     cex.axis = 1.2,
     cex.lab = 1.2)
lines(cutoffs, specificity, col = "red", lwd = 2)
# Mark the intersection point
abline(v = optimal_threshold, lty = "dashed", lwd = 1.5, col = "darkgreen")
abline(h = optimal_sensitivity, lty = "dotted", lwd = 0.5, col = "gray50")
abline(h = optimal_specificity, lty = "dotted", lwd = 0.5, col = "gray50")
# Add point at intersection
points(optimal_threshold, optimal_sensitivity, pch = 19, col = "darkgreen", cex = 1.5)
legend("bottomleft",
       col = c("blue", "red", "darkgreen"),
       lwd = 2,
       legend = c("Sensitivity", "Specificity", sprintf("Optimal (%.3f)", optimal_threshold)),
       cex = 1.0)

dev.off()

cat(sprintf("  Performance plots saved to: %s\n", png_file))
cat("\n")

# Save threshold to file
threshold_file <- paste0(output_prefix, "_threshold.txt")
cat("Writing threshold file...\n")
write(sprintf("%.4f", optimal_threshold), threshold_file)
cat(sprintf("  Threshold saved to: %s\n", threshold_file))
cat("\n")

# Save detailed statistics
stats_file <- paste0(output_prefix, "_stats.txt")
cat("Writing statistics file...\n")

sink(stats_file)
cat("CPAT Optimal Threshold Calculation\n")
cat("================================================================================\n")
cat("\n")
cat("Methodology: 10-Fold Cross-Validation (CPAT Original)\n")
cat("  - Data split into 10 folds\n")
cat("  - Each fold: train on 90%, predict on 10%\n")
cat("  - Logistic regression: glm(Label ~ mRNA + ORF + Fickett + Hexamer)\n")
cat("  - Predictions: predict(..., type='response') for probabilities (0-1)\n")
cat("  - All fold predictions combined for ROC analysis\n")
cat("\n")
cat("Input Files:\n")
cat(sprintf("  Feature file: %s\n", feature_file))
cat("\n")
cat("Training Data:\n")
cat(sprintf("  Coding sequences: %d\n", n_coding))
cat(sprintf("  Non-coding sequences: %d\n", n_noncoding))
cat(sprintf("  Total sequences: %d\n", total_n))
cat("\n")
cat("ROC Analysis:\n")
cat(sprintf("  AUC: %.4f\n", auc_value))
cat("\n")
cat("Optimal Threshold (Sensitivity-Specificity Intersection):\n")
cat("  Method: Following CPAT original publication methodology\n")
cat("  The threshold is determined by finding where sensitivity and\n")
cat("  specificity curves intersect (or are closest).\n")
cat("\n")
cat(sprintf("  Threshold: %.4f\n", optimal_threshold))
cat(sprintf("  Sensitivity (True Positive Rate): %.4f (%.2f%%)\n", optimal_sensitivity, optimal_sensitivity * 100))
cat(sprintf("  Specificity (True Negative Rate): %.4f (%.2f%%)\n", optimal_specificity, optimal_specificity * 100))
cat(sprintf("  Absolute difference: %.4f\n", abs(optimal_sensitivity - optimal_specificity)))
cat("\n")
cat("Classification Performance:\n")
cat(sprintf("  Coding correctly classified: %d/%d (%.2f%%)\n", 
            coding_correctly_classified, length(coding_scores), coding_accuracy * 100))
cat(sprintf("  Non-coding correctly classified: %d/%d (%.2f%%)\n", 
            noncoding_correctly_classified, length(noncoding_scores), noncoding_accuracy * 100))
cat(sprintf("  Overall accuracy: %.2f%%\n", overall_accuracy * 100))
cat("\n")
cat("Score Distribution (Cross-Validated Predictions):\n")
cat("  Coding sequences:\n")
cat(sprintf("    Min: %.4f\n", min(coding_scores)))
cat(sprintf("    Max: %.4f\n", max(coding_scores)))
cat(sprintf("    Mean: %.4f\n", mean(coding_scores)))
cat(sprintf("    Median: %.4f\n", median(coding_scores)))
cat("\n")
cat("  Non-coding sequences:\n")
cat(sprintf("    Min: %.4f\n", min(noncoding_scores)))
cat(sprintf("    Max: %.4f\n", max(noncoding_scores)))
cat(sprintf("    Mean: %.4f\n", mean(noncoding_scores)))
cat(sprintf("    Median: %.4f\n", median(noncoding_scores)))
cat("\n")
cat("Interpretation:\n")
cat(sprintf("  Use threshold %.4f for classifying transcripts:\n", optimal_threshold))
cat(sprintf("    - Transcripts with score >= %.4f are classified as CODING\n", optimal_threshold))
cat(sprintf("    - Transcripts with score < %.4f are classified as NON-CODING\n", optimal_threshold))
cat("\n")
sink()

cat(sprintf("  Statistics saved to: %s\n", stats_file))
cat("\n")

cat("================================================================================\n")
cat("Threshold calculation completed successfully!\n")
cat("================================================================================\n")
cat("\n")
cat("Summary:\n")
cat(sprintf("  Method: 10-Fold Cross-Validation (CPAT original)\n"))
cat(sprintf("  Optimal threshold: %.4f\n", optimal_threshold))
cat(sprintf("  AUC: %.4f\n", auc_value))
cat(sprintf("  Sensitivity: %.2f%%\n", optimal_sensitivity * 100))
cat(sprintf("  Specificity: %.2f%%\n", optimal_specificity * 100))
cat("\n")
cat("Output files:\n")
cat(sprintf("  %s (threshold value)\n", threshold_file))
cat(sprintf("  %s (performance plots)\n", png_file))
cat(sprintf("  %s (detailed statistics)\n", stats_file))
cat("\n")
