### Load packages

library(pROC)
library(dplyr)

### Define classification metrics functions

calculate_metrics <- function(actual, predicted_prob, threshold) {
  predicted_class <- ifelse(predicted_prob >= threshold, 1, 0)
  confusion_matrix <- table(
    Predicted = factor(predicted_class, levels = c(0, 1)),
    Actual = factor(actual, levels = c(0, 1))
  )
  
  TN <- confusion_matrix["0", "0"]
  FN <- confusion_matrix["0", "1"]
  FP <- confusion_matrix["1", "0"]
  TP <- confusion_matrix["1", "1"]
  
  data.frame(
    accuracy = (TP + TN) / (TP + TN + FP + FN),
    sensitivity = TP / (TP + FN),
    specificity = TN / (TN + FP),
    ppv = TP / (TP + FP),
    npv = TN / (TN + FN)
  )
}

### Vessel-level ROC analysis

roc_development <- roc(
  response = development_set$label,
  predictor = development_set$predicted_prob,
  quiet = TRUE
)

roc_test <- roc(
  response = test_set$label,
  predictor = test_set$predicted_prob,
  quiet = TRUE
)

### AUC and 95% CI

auc_development <- as.numeric(auc(roc_development))
auc_test <- as.numeric(auc(roc_test))

ci_development <- ci.auc(roc_development)
ci_test <- ci.auc(roc_test)

cat("Vessel-level AUC (development):", auc_development, "\n")
cat("95% CI:", ci_development[1], "-", ci_development[3], "\n\n")

cat("Vessel-level AUC (test):", auc_test, "\n")
cat("95% CI:", ci_test[1], "-", ci_test[3], "\n")

### Determine classification threshold on development set

youden_index <- (
  roc_development$sensitivities +
  roc_development$specificities - 1
)

optimal_index <- which.max(youden_index)
vessel_threshold <- roc_development$thresholds[optimal_index]
cat("\nOptimal threshold:", vessel_threshold, "\n")

### Vessel-level classification metrics

development_metrics <- calculate_metrics(
  actual = development_set$label,
  predicted_prob = development_set$predicted_prob,
  threshold = vessel_threshold
)

test_metrics <- calculate_metrics(
  actual = test_set$label,
  predicted_prob = test_set$predicted_prob,
  threshold = vessel_threshold
)

development_metrics$dataset <- "development"
test_metrics$dataset <- "test"

vessel_metrics <- bind_rows(
  development_metrics,
  test_metrics
) %>%
  select(
    dataset,
    everything()
  )

cat("\nVessel-level performance:\n")
print(vessel_metrics)

### Patient-level prediction
# Note: Patient-level predicted probability is represented by the
# maximum predicted probability among the patient's vessels.

development_patient_predictions <- development_set %>%
  group_by(patient_id) %>%
  summarise(
    n_vessels = n(),
    patient_label = as.integer(any(label == 1)),
    predicted_prob = max(predicted_prob),
    .groups = "drop"
  )

test_patient_predictions <- test_set %>%
  group_by(patient_id) %>%
  summarise(
    n_vessels = n(),
    patient_label = as.integer(any(label == 1)),
    predicted_prob = max(predicted_prob),
    .groups = "drop"
  )

### Patient-level ROC analysis

roc_patient_development <- roc(
  response = development_patient_predictions$patient_label,
  predictor = development_patient_predictions$predicted_prob,
  quiet = TRUE
)

roc_patient_test <- roc(
  response = test_patient_predictions$patient_label,
  predictor = test_patient_predictions$predicted_prob,
  quiet = TRUE
)

auc_patient_development <- as.numeric(auc(roc_patient_development))
auc_patient_test <- as.numeric(auc(roc_patient_test))
ci_patient_development <- ci.auc(roc_patient_development)
ci_patient_test <- ci.auc(roc_patient_test)

cat("\nPatient-level AUC (development):", auc_patient_development, "\n")
cat("Patient-level AUC (test):", auc_patient_test, "\n")

### Determine patient-level threshold on development set

patient_youden <- (
  roc_patient_development$sensitivities +
  roc_patient_development$specificities - 1
)

patient_optimal_index <- which.max(patient_youden)
patient_threshold <- roc_patient_development$thresholds[patient_optimal_index]

cat("\nPatient-level threshold:", patient_threshold, "\n")

### Patient-level classification metrics

patient_development_metrics <- calculate_metrics(
  actual = development_patient_predictions$patient_label,
  predicted_prob = development_patient_predictions$predicted_prob,
  threshold = patient_threshold
)

patient_test_metrics <- calculate_metrics(
  actual = test_patient_predictions$patient_label,
  predicted_prob = test_patient_predictions$predicted_prob,
  threshold = patient_threshold
)

patient_development_metrics$dataset <- "development"
patient_test_metrics$dataset <- "test"

patient_metrics <- bind_rows(
  patient_development_metrics,
  patient_test_metrics
) %>%
  select(
    dataset,
    everything()
  )

cat("\nPatient-level performance:\n")
print(patient_metrics)

### Add AUC nd CI to output tables

vessel_metrics$auc <- c(auc_development, auc_test)
vessel_metrics$auc_ci_lower <- c(ci_development[1], ci_test[1])
vessel_metrics$auc_ci_upper <- c(ci_development[3], ci_test[3])
patient_metrics$auc <- c(auc_patient_development, auc_patient_test)
patient_metrics$auc_ci_lower <- c(ci_patient_development[1], ci_patient_test[1])
patient_metrics$auc_ci_upper <- c(ci_patient_development[3], ci_patient_test[3])

### Save tables

write.csv(vessel_metrics, "outputs/r/tables/model_performance_vessel.csv",row.names = FALSE)
write.csv(patient_metrics, "outputs/r/tables/model_performance_patient.csv", row.names = FALSE)
write.csv(development_patient_predictions, "outputs/r/tables/patient_predictions_development.csv", row.names = FALSE)
write.csv(test_patient_predictions, "outputs/r/tables/patient_predictions_test.csv", row.names = FALSE)

### Save ROC curves

png("outputs/r/figures/roc_vessel.png",
  width = 1800,
  height = 800,
  res = 300,
  pointsize = 10,
  bg = "white"
)

layout(matrix(c(1, 2), nrow = 1))

par(
  mar = c(4.5, 4.5, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot(roc_development, legacy.axes = TRUE, main = "Development set", print.auc = TRUE)

par(
  mar = c(4.5, 4.5, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot(roc_test, legacy.axes = TRUE, main = "Test set", print.auc = TRUE)

dev.off()

png("outputs/r/figures/roc_patient.png",
  width = 1800,
  height = 800,
  res = 300,
  pointsize = 10,
  bg = "white"
)

layout(matrix(c(1, 2), nrow = 1))

par(
  mar = c(4.5, 4.5, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot(roc_patient_development, legacy.axes = TRUE, main = "Development set", print.auc = TRUE)

par(
  mar = c(4.5, 4.5, 2.5, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  pty = "s"
)

plot(roc_patient_test, legacy.axes = TRUE, main = "Test set", print.auc = TRUE)

dev.off()