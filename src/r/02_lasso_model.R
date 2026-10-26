### Load packages

library(glmnet)
library(caret)
library(dplyr)

### Define 39 features

feature_vars <- c(
  "ds_ve",
  "ds_ve_70",
  "ds_qca",
  "ds_qca_50",
  "as_qca",
  "reference_vessel_diameter",
  "reference_vessel_area",
  "minimal_lumen_diameter",
  "mla",
  "pro_diameter",
  "lesion_dia",
  "dis_diameter",
  "pro_area",
  "lesion_area",
  "dis_area",
  "pro_length",
  "lesion_length",
  "dis_length",
  "plaque_symmetry",
  "plaque_area",
  "sfr",
  "poisseuille_resistance",
  "turbulent_resistance",
  "ll_mld4",
  "inflow_angle",
  "outflow_angle",
  "age",
  "gender",
  "weight",
  "height",
  "bmi",
  "sbp",
  "dbp",
  "angina",
  "hypertension",
  "hyperlipidemia",
  "diabetes_mellitus",
  "smoking",
  "family_history_of_cad"
)

### LASSO input

x_development <- as.matrix(
  development_set %>% select(all_of(feature_vars))
)

y_development <- development_set$label

### 10-fold cross validation (still use patient-grouped)

development_patients <- patient_summary %>%
  filter(patient_id %in% development_patient_ids)

set.seed(100)
patient_folds <- createFolds(
  factor(development_patients$patient_label),
  k = 10,
  list = TRUE,
  returnTrain = FALSE
)

## Idea: first assign folds to patient rows, store the assignment results in a new column [folder_id], 
## then transfer patient assignments back to vessel rows

# Assign each patient to one fold
development_patients$fold_id <- NA_integer_
for (i in seq_along(patient_folds)) {
  development_patients$fold_id[patient_folds[[i]]] <- i
} 

# Assign every vessel the fold of its patient
fold_id <- development_patients$fold_id[
  match(
    development_set$patient_id,
    development_patients$patient_id
  )
] 

### Run cv.glmnet

set.seed(100)

lambda_grid <- 10 ^ seq(from = -3, to = 0, length.out = 100)

cv_model <- cv.glmnet(
  x = x_development,
  y = y_development,
  family = "binomial", # binary classification
  alpha = 1,           # LASSO
  lambda = lambda_grid,
  foldid = fold_id,    # 10-folded CV
  standardize = TRUE   # predictors were automatically standardized in LASSO
)

# Compare lambda.min and lambda.1se
lambda_min <- cv_model$lambda.min
lambda_1se <- cv_model$lambda.1se
cat("lambda.min:", lambda_min, "\n")
cat("lambda.1se:", lambda_1se, "\n")

### Fit LASSO coefficient paths

lasso_path <- glmnet(
  x = x_development, 
  y = y_development,
  family = "binomial",
  alpha = 1,
  lambda = lambda_grid,
  standardize = TRUE
)

### Fit best LASSO model (compare lambda.min and lambda.1se)

model_min <- glmnet(
  x = x_development,
  y = y_development,
  family = "binomial",
  alpha = 1,
  lambda = lambda_min,
  standardize = TRUE
)

model_1se <- glmnet(
  x = x_development,
  y = y_development,
  family = "binomial",
  alpha = 1,
  lambda = lambda_1se,
  standardize = TRUE
)

### Extract coefficients (both lambda.min and lambda.1se)

coef_min <- as.matrix(coef(model_min))
coef_1se <- as.matrix(coef(model_1se))

coefficient_table_min <- data.frame(
  feature = rownames(coef_min),
  coefficient = as.numeric(coef_min[, 1]),
  row.names = NULL
)

coefficient_table_1se <- data.frame(
  feature = rownames(coef_1se),
  coefficient = as.numeric(coef_1se[, 1]),
  row.names = NULL
)

selected_features_min <- coefficient_table_min %>%
  filter(feature != "(Intercept)", coefficient != 0)

selected_features_1se <- coefficient_table_1se %>%
  filter(feature != "(Intercept)", coefficient != 0)

cat("\nSelected features using lambda.min:\n")
print(selected_features_min)

cat("\nSelected features using lambda.1se:\n")
print(selected_features_1se)

### Select final model (lambda.1se)
# The 1se rule is used for the final model because it provides a more 
# parsimonious model while maintaining CV performance within one standard error of min.
final_lambda <- lambda_1se
final_model <- model_1se
final_selected_features <- selected_features_1se

### Generate development/test predictions

x_test <- as.matrix(
  test_set %>% select(all_of(feature_vars))
)

development_set$predicted_prob <- as.numeric(
  predict(final_model, newx = x_development, type = "response")
)

test_set$predicted_prob <- as.numeric(
  predict(final_model, newx = x_test, type = "response")
)

### Save outputs

dir.create("outputs/r/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/r/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  coefficient_table_min,
  "outputs/r/tables/lasso_coefficients_lambda_min.csv",
  row.names = FALSE
)

write.csv(
  coefficient_table_1se,
  "outputs/r/tables/lasso_coefficients_lambda_1se.csv",
  row.names = FALSE
)

write.csv(
  selected_features_min,
  "outputs/r/tables/lasso_selected_features_lambda_min.csv",
  row.names = FALSE
)

write.csv(
  selected_features_1se,
  "outputs/r/tables/lasso_selected_features_lambda_1se.csv",
  row.names = FALSE
)

# lasso CV model plot
png(
  "outputs/r/figures/lasso_cross_validation.png", 
  width = 1800, height = 1200, res = 300, pointsize = 10, bg = "white")
par(
  mar = c(4.5, 4.8, 1.2, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0, 
  cex.main = 1.0, 
  tcl = -0.3 
)
plot(cv_model)
dev.off()

# LASSO coefficient path plot
png(
  "outputs/r/figures/lasso_coefficient_paths.png", 
  width = 1800, height = 1200, res = 300, pointsize = 10, bg = "white")
par(
  mar = c(4.5, 4.8, 1.2, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0, 
  cex.main = 1.0, 
  tcl = -0.3 
)
plot(lasso_path, xvar = "lambda", label = TRUE)
dev.off()

# Combine 2 plots
png(
  "outputs/r/figures/lasso_combined.png",
  width = 2400,
  height = 800,
  res = 300,
  pointsize = 10,
  bg = "white"
)

layout(matrix(c(1, 2), nrow = 1))

# Left: CV plot
par(
  mar = c(4.5, 5.0, 3.2, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  tcl = -0.3,
  xpd = FALSE
)
plot(cv_model)
par(xpd = NA)

# Right: coefficient path
par(
  mar = c(4.5, 5.0, 3.2, 1.2),
  mgp = c(2.4, 0.8, 0),
  cex.axis = 0.85,
  cex.lab = 1.0,
  tcl = -0.3,
  xpd = FALSE
)
plot(
  lasso_path,
  xvar = "lambda",
  label = FALSE
)
par(xpd = NA)
dev.off()