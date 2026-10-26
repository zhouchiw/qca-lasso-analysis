### Load packages

library(readxl)
library(dplyr)
library(janitor)
library(stringr)
library(caret)

### Load data

data_path <- "data/2DQCA_total.xlsx"
vessel_data <- read_excel(data_path, sheet = "per-vessel")

### Data cleaning (simplified, mostly finished in Excel)

vessel_data <- vessel_data %>%
  clean_names() %>%
    mutate(
      center = str_trim(center),
      vessel_location = str_trim(vessel_location),
      center = recode(
        center,
        "NUH+" = "NUH",
        "NHCS+" = "NHCS"
        )
      )

# Create a unique per-patient level identifier

vessel_data <- vessel_data %>%
  mutate(
    patient_id = paste(center, patient_no, sep = "_") #e.g. NHCS_001
  )

# Dataset summary

cat("Dataset summary\n")
cat("Number of vessels:", nrow(vessel_data), "\n")
cat("Number of patients:", n_distinct(vessel_data$patient_id), "\n")
cat("Number of centers:", n_distinct(vessel_data$center), "\n")
cat("\nCenters:\n")
print(sort(unique(vessel_data$center)))
cat("\nVessel-level label distribution:\n")
print(table(vessel_data$label))

# Check missing values

missing_values <- colSums(is.na(vessel_data))

if (any(missing_values > 0)) {
  cat("\nVariables containing missing values:\n")
  print(missing_values[missing_values > 0])
} else {
  cat("\nNo missing values detected.\n")
}

### Derive patient-level outcomes

# Note: if at least one vessel of a patient is labeled positive (FFR<=0.8), then the patient-level label is positive.
patient_summary <- vessel_data %>%
  group_by(patient_id, center, patient_no) %>%
  summarise(
    n_vessels = n(),
    min_ffr = min(ffr, na.rm = TRUE),
    patient_label = as.integer(any(label == 1)),
    .groups = "drop"
  )

# validation
stopifnot(all(patient_summary$patient_label == as.integer(patient_summary$min_ffr <= 0.80)))

### Data development/test split

# Split patient rather than split vessel to prevent patient-level data leakage)
set.seed(100)
development_index <- createDataPartition(
  factor(patient_summary$patient_label), p = 0.80, list = FALSE
)

development_index <- as.vector(development_index)
development_patient_ids <- patient_summary$patient_id[development_index]
test_patient_ids <- patient_summary$patient_id[-development_index]

# Create vessel-level development and test sets

development_set <- vessel_data %>%
  filter(patient_id %in% development_patient_ids)

test_set <- vessel_data %>%
  filter(patient_id %in% test_patient_ids)

# Split summary

cat("\nDevelopment/test split\n")
cat("Development set:", nrow(development_set),"vessels from", n_distinct(development_set$patient_id),"patients\n")
cat("Test set:", nrow(test_set),"vessels from", n_distinct(test_set$patient_id),"patients\n")