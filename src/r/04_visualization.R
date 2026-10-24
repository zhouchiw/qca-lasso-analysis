### Load packages

library(pROC)
library(dplyr)
library(grid)

### Output directories

dir.create("outputs/r/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/r/tables", recursive = TRUE, showWarnings = FALSE)

### Define predictors for comparison

comparison_predictors <- list(
  "LASSO model" = list(variable = "predicted_prob", positive_if = "higher"),
  "SFR" = list(variable = "sfr", positive_if = "lower"),
  "Minimal lumen diameter" = list(variable = "minimal_lumen_diameter", positive_if = "lower"),
  "DS (VE)" = list(variable = "ds_ve", positive_if = "higher"),
  "Lesion length" = list(variable = "lesion_length", positive_if = "higher")
)

### Calculate classification metrics from predicted class

metrics_from_class <- function(actual, predicted_class) {

  cm <- table(
    Predicted = factor(predicted_class, levels = c(0, 1)),
    Actual = factor(actual, levels = c(0, 1))
  )
  
  TN <- cm["0", "0"]
  FN <- cm["0", "1"]
  FP <- cm["1", "0"]
  TP <- cm["1", "1"]
  
  data.frame(
    accuracy = (TP + TN) / (TP + TN + FP + FN),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    ppv = TP / (TP + FP),
    npv = TN / (TN + FN)
  )
}

### Evaluate one predictor

evaluate_predictor <- function(
    development_data,
    test_data,
    outcome,
    predictor,
    predictor_name,
    positive_if
) {
  
  # pROC direction
  roc_direction <- ifelse(positive_if == "higher", "<", ">")
  
  # ROC curves
  roc_dev <- roc(
    response = development_data[[outcome]],
    predictor = development_data[[predictor]],
    direction = roc_direction,
    percent = TRUE,
    quiet = TRUE
  )
  
  roc_test <- roc(
    response = test_data[[outcome]],
    predictor = test_data[[predictor]],
    direction = roc_direction,
    percent = TRUE,
    quiet = TRUE
  )
  
  # Cutoff
  cutoff <- as.numeric(
    coords(
      roc_dev,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    )[1]
  )
  
  # Convert predictor into binary predictions
  if (positive_if == "higher") {
    dev_class <- ifelse(development_data[[predictor]] >= cutoff, 1, 0)
    test_class <- ifelse(test_data[[predictor]] >= cutoff, 1, 0)
  } else {
    dev_class <- ifelse(development_data[[predictor]] <= cutoff, 1, 0)
    test_class <- ifelse(test_data[[predictor]] <= cutoff, 1, 0)
  }
  
  # Metrics
  dev_metrics <- metrics_from_class(
    actual = development_data[[outcome]],
    predicted_class = dev_class
  )
  
  test_metrics <- metrics_from_class(
    actual = test_data[[outcome]],
    predicted_class = test_class
  )
  
  # AUC
  dev_auc <- as.numeric(auc(roc_dev)) / 100
  test_auc <- as.numeric(auc(roc_test)) / 100
  dev_ci <- as.numeric(ci.auc(roc_dev)) / 100
  test_ci <- as.numeric(ci.auc(roc_test)) / 100
  
  development_result <- data.frame(
    predictor = predictor_name,
    cutoff = cutoff,
    accuracy = dev_metrics$accuracy,
    sensitivity = dev_metrics$sensitivity,
    specificity = dev_metrics$specificity,
    ppv = dev_metrics$ppv,
    npv = dev_metrics$npv,
    auc = dev_auc,
    auc_ci_lower = dev_ci[1],
    auc_ci_upper = dev_ci[3]
  )
  
  test_result <- data.frame(
    predictor = predictor_name,
    cutoff = cutoff,
    accuracy = test_metrics$accuracy,
    sensitivity = test_metrics$sensitivity,
    specificity = test_metrics$specificity,
    ppv = test_metrics$ppv,
    npv = test_metrics$npv,
    auc = test_auc,
    auc_ci_lower = test_ci[1],
    auc_ci_upper = test_ci[3]
  )
  
  list(
    development = development_result,
    test = test_result,
    roc_development = roc_dev,
    roc_test = roc_test,
    cutoff = cutoff,
    development_class = dev_class,
    test_class = test_class
  )
}

### Vessel-level comparison

vessel_results <- lapply(
  names(comparison_predictors),
  function(name) {
    predictor_info <- comparison_predictors[[name]]
    evaluate_predictor(
      development_data = development_set,
      test_data = test_set,
      outcome = "label",
      predictor = predictor_info$variable,
      predictor_name = name,
      positive_if = predictor_info$positive_if
    )
  }
)

names(vessel_results) <- names(comparison_predictors)

vessel_comparison_development <- bind_rows(
  lapply(
    vessel_results,
    function(x) x$development
  )
)

vessel_comparison_test <- bind_rows(
  lapply(
    vessel_results,
    function(x) x$test
  )
)

### Patient-level comparison
# Note: Patient-level predictors are represented by the most severe vessel value.

development_patient_comparison <- development_set %>%
  group_by(patient_id) %>%
  summarise(
    patient_label = as.integer(any(label == 1)),
    predicted_prob = max(predicted_prob),
    sfr = min(sfr),
    minimal_lumen_diameter = min(minimal_lumen_diameter),
    ds_ve = max(ds_ve),
    lesion_length = max(lesion_length),
    .groups = "drop"
  )

test_patient_comparison <- test_set %>%
  group_by(patient_id) %>%
  summarise(
    patient_label = as.integer(any(label == 1)),
    predicted_prob = max(predicted_prob),
    sfr = min(sfr),
    minimal_lumen_diameter = min(minimal_lumen_diameter),
    ds_ve = max(ds_ve),
    lesion_length = max(lesion_length),
    .groups = "drop"
  )

patient_results <- lapply(
  names(comparison_predictors),
  function(name) {
    predictor_info <- comparison_predictors[[name]]
    evaluate_predictor(
      development_data = development_patient_comparison,
      test_data = test_patient_comparison,
      outcome = "patient_label",
      predictor = predictor_info$variable,
      predictor_name = name,
      positive_if = predictor_info$positive_if
    )
  }
)

names(patient_results) <- names(comparison_predictors)

patient_comparison_development <- bind_rows(
  lapply(
    patient_results,
    function(x) x$development
  )
)

patient_comparison_test <- bind_rows(
  lapply(
    patient_results,
    function(x) x$test
  )
)

### Save comparison tables

write.csv(vessel_comparison_development, "outputs/r/tables/performance_comparison_vessel_development.csv", row.names = FALSE)
write.csv(vessel_comparison_test, "outputs/r/tables/performance_comparison_vessel_test.csv", row.names = FALSE)
write.csv(patient_comparison_development, "outputs/r/tables/performance_comparison_patient_development.csv", row.names = FALSE)
write.csv(patient_comparison_test, "outputs/r/tables/performance_comparison_patient_test.csv", row.names = FALSE)

### ROC comparison figures

custom_blue <- "#013D7C"
custom_green <- "#007A49"
custom_orange <- "#D4A418"
custom_red <- "#AD2D4D"
custom_light <- "#AFC8FF"

roc_colors <- c(
  custom_blue,
  custom_green,
  custom_orange,
  custom_red,
  custom_light
)

### ROC styles

plot_roc_comparison <- function(results, roc_name, title) {
  plot(
    results[[1]][[roc_name]],
    legacy.axes = TRUE,
    main = title,
    col = roc_colors[1],
    lwd = 2
  )
  grid(
    col = "lightgrey",
    lty = "dotted"
  )
  
  for (i in 2:length(results)) {
    lines(
      results[[i]][[roc_name]],
      col = roc_colors[i],
      lwd = 2
    )
  }
  
  legend(
    "bottomright",
    legend = names(results),
    col = roc_colors,
    lwd = 2,
    cex = 0.7,
    box.lwd = 0.8
  )
}

### Vessel-level ROC

png("outputs/r/figures/roc_comparison_vessel.png",
  width = 1800,
  height = 850,
  res = 300,
  pointsize = 10,
  bg = "white"
)

layout(matrix(c(1, 2), nrow = 1))

par(
  mar = c(4.5, 4.8, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot_roc_comparison(vessel_results, "roc_development",
  paste0("Development set (n = ", nrow(development_set),")")
)

par(
  mar = c(4.5, 4.8, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot_roc_comparison(vessel_results, "roc_test",
  paste0("Test set (n = ", nrow(test_set), ")")
)

dev.off()

### Patient-level ROC

png("outputs/r/figures/roc_comparison_patient.png",
  width = 1800,
  height = 850,
  res = 300,
  pointsize = 10,
  bg = "white"
)

layout(matrix(c(1, 2), nrow = 1))

par(
  mar = c(4.5, 4.8, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot_roc_comparison(patient_results, "roc_development",
  paste0("Development set (n = ", nrow(development_patient_comparison), ")")
)

par(
  mar = c(4.5, 4.8, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot_roc_comparison(patient_results, "roc_test",
  paste0("Test set (n = ", nrow(test_patient_comparison), ")")
)

dev.off()

### Confusion matrices for final LASSO model

# Vessel-level
vessel_lasso_test_class <-
  vessel_results[["LASSO model"]]$test_class

confusion_vessel_test <- table(
  Predicted = factor(
    vessel_lasso_test_class,
    levels = c(1, 0),
    labels = c("Positive", "Negative")
  ),
  Actual = factor(
    test_set$label,
    levels = c(1, 0),
    labels = c("Positive", "Negative")
  )
)

# Patient-level
patient_lasso_test_class <-
  patient_results[["LASSO model"]]$test_class

confusion_patient_test <- table(
  Predicted = factor(
    patient_lasso_test_class,
    levels = c(1, 0),
    labels = c("Positive", "Negative")
  ),
  Actual = factor(
    test_patient_comparison$patient_label,
    levels = c(1, 0),
    labels = c("Positive", "Negative")
  )
)

format_confusion_matrix <- function(cm) {
  
  data.frame(
    Predicted = c(
      "Positive",
      "Negative"
    ),
    
    `Actual: Positive` = c(
      cm["Positive", "Positive"],
      cm["Negative", "Positive"]
    ),
    
    `Actual: Negative` = c(
      cm["Positive", "Negative"],
      cm["Negative", "Negative"]
    ),
    
    check.names = FALSE
  )
}

# Save
write.csv(format_confusion_matrix(confusion_vessel_test),"outputs/r/tables/confusion_matrix_vessel_test.csv", row.names = FALSE)
write.csv(format_confusion_matrix(confusion_patient_test), "outputs/r/tables/confusion_matrix_patient_test.csv", row.names = FALSE)