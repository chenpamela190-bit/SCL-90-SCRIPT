# ===================================================================
#
# SCRIPT FOR MANUSCRIPT REVISIONS (MS# "Profile-Specific SCL-90-R...")
#
# v16 (revised): 
#   - Long-string QC threshold set to 20 (everywhere).
#   - Fixed drop_na / na.omit usage.
#   - Step 4 descriptives aligned with A–D samples in response letter.
#
# ===================================================================

# -------------------------------------------------------------------
# STEP 0: INSTALL PACKAGES & LOAD LIBRARIES
# -------------------------------------------------------------------
# install.packages(c(
#   "tidyverse", "data.table", "lavaan", "mclust",
#   "NetworkComparisonTest", "qgraph", "patchwork", "psych", "semTools"
# ))

cat("STEP 0: Loading all required packages...\n")
library(tidyverse)
library(data.table)
library(lavaan)
library(mclust)
library(NetworkComparisonTest)
library(qgraph)
library(patchwork)
library(psych)      # for cohen.kappa
library(semTools)   # run install.packages("semTools") ONCE if needed

set.seed(12345)
cat("STEP 0: Complete.\n\n")


# -------------------------------------------------------------------
# STEP 1: LOAD DATA & DEFINE SCL-90 DOMAINS
# -------------------------------------------------------------------
cat("STEP 1: Loading and cleaning data...\n")

FILE_PATH <- "/Users/chenxiaohui/Desktop/data/cleaned_4-SCL-90.csv"

full_data <- fread(FILE_PATH)
cat(paste("Successfully loaded", nrow(full_data), "total rows.\n"))

# --- Define SCL-90-R Domain Items ---
som_items  <- c("Q1", "Q4", "Q12", "Q27", "Q40", "Q42", "Q48", "Q49", "Q52", "Q53", "Q56", "Q58")
ocd_items  <- c("Q3", "Q9", "Q10", "Q28", "Q38", "Q45", "Q46", "Q51", "Q55", "Q65")
is_items   <- c("Q6", "Q21", "Q34", "Q36", "Q37", "Q41", "Q61", "Q69", "Q73")
dep_items  <- c("Q5", "Q14", "Q15", "Q20", "Q22", "Q26", "Q29", "Q30", "Q31", "Q32", "Q54", "Q71", "Q79")
anx_items  <- c("Q2", "Q17", "Q23", "Q33", "Q39", "Q57", "Q72", "Q78", "Q80", "Q86")
phob_items <- c("Q13", "Q25", "Q47", "Q50", "Q70", "Q75", "Q82")

all_item_cols <- paste0("Q", 1:90)
domain_cols   <- c("SOM", "OCD", "IS", "DEP", "ANX", "PHOB")

# Convert all item columns to numeric
full_data <- full_data %>%
  mutate(across(all_of(all_item_cols), as.numeric))

# time_taken to numeric
cat("Forcing 'time_taken' column to numeric...\n")
full_data <- full_data %>%
  mutate(time_taken = as.numeric(time_taken))

# AGE cleaning
cat("Cleaning 'AGE' column for outliers...\n")
full_data <- full_data %>%
  mutate(
    AGE = as.numeric(AGE),
    AGE = ifelse(AGE < 16 | AGE > 30, NA, AGE)
  )
cat("'AGE' column is now clean.\n")

# Domain scores
cat("Calculating domain scores...\n")
data_with_domains <- full_data %>%
  mutate(
    SOM  = rowMeans(select(., all_of(som_items)),  na.rm = TRUE),
    OCD  = rowMeans(select(., all_of(ocd_items)),  na.rm = TRUE),
    IS   = rowMeans(select(., all_of(is_items)),   na.rm = TRUE),
    DEP  = rowMeans(select(., all_of(dep_items)),  na.rm = TRUE),
    ANX  = rowMeans(select(., all_of(anx_items)),  na.rm = TRUE),
    PHOB = rowMeans(select(., all_of(phob_items)), na.rm = TRUE)
  )

# Grouping variables
data_with_domains <- data_with_domains %>%
  rename(Year = `intake year`, Study_Level_Code = `Study level`) %>%
  mutate(
    Year = as.factor(Year),
    Study_Level = as.factor(case_when(
      Study_Level_Code == 1 ~ "University",
      Study_Level_Code == 2 ~ "Vocational",
      TRUE ~ "Other"
    ))
  )

cat("Data types cleaned. Checking 'Study_Level' variable:\n")
print(table(data_with_domains$Year, data_with_domains$Study_Level))
cat("STEP 1: Complete. Data is loaded and cleaned.\n\n")


# -------------------------------------------------------------------
# STEP 2: ADDRESS REVIEWER POINT #3 (Duplicate Submissions)
# -------------------------------------------------------------------
cat("STEP 2: Running duplicate analysis (Reviewer Point 3)...\n")

raw_uni_22_23 <- data_with_domains %>%
  filter(Year %in% c("2022", "2023"), Study_Level == "University")

duplicate_rates <- raw_uni_22_23 %>%
  group_by(Year) %>%
  summarise(
    total_submissions = n(),
    unique_ids        = n_distinct(ID),
    duplicate_rows    = total_submissions - unique_ids,
    duplicate_rate    = duplicate_rows / total_submissions,
    .groups           = "drop"
  )

cat("Duplicate rates for 2022 and 2023 (University Only):\n")
print(duplicate_rates)

# Earliest attempt per ID (primary de-dup rule)
data_deduped_primary <- data_with_domains %>%
  group_by(Year, ID) %>%
  slice_head(n = 1) %>%
  ungroup()

cat(paste("Total rows after de-duplication (earliest):", nrow(data_deduped_primary), "\n"))

# Latest attempt per ID (sensitivity)
data_deduped_latest <- data_with_domains %>%
  group_by(Year, ID) %>%
  slice_tail(n = 1) %>%
  ungroup()

cat("STEP 2: Complete.\n\n")

data_with_domains <- as.data.frame(data_with_domains)

# ensure domain columns are pure numeric
data_with_domains[domain_cols] <- lapply(
  data_with_domains[domain_cols],
  function(z) as.numeric(unlist(z))
)


# -------------------------------------------------------------------
# STEP 3: ADDRESS REVIEWER POINT #2 (QC Pipeline & Distributions)
# -------------------------------------------------------------------
cat("STEP 3: Running QC Pipeline (Reviewer Point 2)...\n")

# Helper for long-string
get_max_consecutive <- function(x) {
  if (length(x) == 0) return(0)
  rle_x <- rle(x)
  lengths_without_na <- rle_x$lengths[!is.na(rle_x$values)]
  if (length(lengths_without_na) == 0) 0 else max(lengths_without_na)
}

# QC metrics on de-duplicated (earliest) data for all years
qc_data <- data_deduped_primary %>%
  rowwise() %>%
  mutate(
    time_taken_minutes = time_taken / 60,
    long_string        = get_max_consecutive(c_across(all_of(all_item_cols))),
    response_variance  = sd(c_across(all_of(all_item_cols)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  as.data.frame()

# Force domain cols to be atomic numeric
qc_data[domain_cols] <- lapply(
  qc_data[domain_cols],
  function(z) as.numeric(unlist(z))
)

# QC flags  (long-string threshold = 20)
cat("Adding QC flags to 'qc_data' object for later steps...\n")
qc_data <- qc_data %>%
  mutate(
    flag_speeder    = ifelse(is.na(time_taken_minutes), FALSE, time_taken_minutes < 5),
    flag_longstring = ifelse(is.na(long_string),        FALSE, long_string >= 20),
    flag_zerovar    = ifelse(is.na(response_variance),  FALSE, response_variance < 0.10),
    flag_any        = flag_speeder | flag_longstring | flag_zerovar
  )
cat("'qc_data' object is now ready.\n")

# Plots for Revised Figure S1 (based on qc_data, all years)
cat("Generating plots for Revised Figure S1...\n")

plot_time <- ggplot(qc_data, aes(x = time_taken_minutes)) +
  geom_histogram(binwidth = 1, fill = "#0072B2", alpha = 0.7) +
  geom_vline(xintercept = 5, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 7, y = 500, label = "Threshold (5 min)", color = "red",
           angle = 90, vjust = -0.5) +
  labs(title = "A: Distribution of Completion Times",
       x = "Time (Minutes)", y = "Count") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 10))

plot_longstring <- ggplot(qc_data, aes(x = long_string)) +
  geom_histogram(binwidth = 1, fill = "#D55E00", alpha = 0.7) +
  geom_vline(xintercept = 20, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 21, y = 1000, label = "Threshold (>= 20)", color = "red",
           angle = 90, vjust = -0.5) +
  labs(title = "B: Distribution of Long-String Responses",
       x = "Max. Consecutive Identical Answers", y = "Count") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 10))

plot_variance <- ggplot(qc_data, aes(x = response_variance)) +
  geom_histogram(binwidth = 0.02, fill = "#009E73", alpha = 0.7) +
  geom_vline(xintercept = 0.10, color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = 0.12, y = 500, label = "Threshold (< 0.10)", color = "red",
           angle = 90, vjust = -0.5) +
  labs(title = "C: Distribution of Response Variance (SD)",
       x = "Standard Deviation of 90 Items", y = "Count") +
  theme_bw() +
  scale_x_continuous(limits = c(0, 2))

revised_figure_s1 <- plot_time + plot_longstring + plot_variance
print(revised_figure_s1)
# ggsave("Revised_Figure_S1.png", revised_figure_s1, width = 14, height = 5)

# QC filters for final CONFIRMATORY samples (raw university data 2022–2024)
raw_uni_data <- data_with_domains %>%
  filter(Study_Level == "University", Year %in% c("2022", "2023", "2024"))

raw_uni_qc <- raw_uni_data %>%
  rowwise() %>%
  mutate(
    time_taken_minutes = time_taken / 60,
    long_string        = get_max_consecutive(c_across(all_of(all_item_cols))),
    response_variance  = sd(c_across(all_of(all_item_cols)), na.rm = TRUE)
  ) %>%
  ungroup()

qc_flags_raw <- raw_uni_qc %>%
  mutate(
    flag_speeder    = ifelse(is.na(time_taken_minutes), FALSE, time_taken_minutes < 5),
    flag_longstring = ifelse(is.na(long_string),        FALSE, long_string >= 20),
    flag_zerovar    = ifelse(is.na(response_variance),  FALSE, response_variance < 0.10),
    flag_any        = flag_speeder | flag_longstring | flag_zerovar
  )

cat("\n--- Generating REVISED TABLE S2 (QC Exclusions, University 2022–2024) ---\n")
qc_table_s2 <- qc_flags_raw %>%
  group_by(Year) %>%
  summarise(
    N_raw                   = n(),
    Flag_Speeder            = sum(flag_speeder,    na.rm = TRUE),
    Flag_LongString         = sum(flag_longstring, na.rm = TRUE),
    Flag_ZeroVar            = sum(flag_zerovar,    na.rm = TRUE),
    Total_QC_Flagged_Unique = sum(flag_any,        na.rm = TRUE),
    .groups                 = "drop"
  ) %>%
  mutate(Remaining_after_QC = N_raw - Total_QC_Flagged_Unique)

print(qc_table_s2)

# Final analytic samples (earliest vs latest, post-QC)
passed_qc_ids <- qc_flags_raw %>%
  filter(!flag_any) %>%
  select(ID, Year)

final_sample <- data_with_domains %>%
  inner_join(passed_qc_ids, by = c("ID", "Year")) %>%
  group_by(Year, ID) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  as.data.frame()

cat("\nFinal Analytic Sample (earliest) created.\n")
print(table(final_sample$Year))

final_sample_LATEST <- data_with_domains %>%
  inner_join(passed_qc_ids, by = c("ID", "Year")) %>%
  group_by(Year, ID) %>%
  slice_tail(n = 1) %>%
  ungroup() %>%
  as.data.frame()

cat("\nFinal Sensitivity Sample (latest) created.\n")
print(table(final_sample_LATEST$Year))
cat("STEP 3: Complete.\n\n")


# -------------------------------------------------------------------
# STEP 4: DESCRIPTIVES FOR TABLE S1 (A–D)
# -------------------------------------------------------------------
cat("STEP 4: Generating descriptives for Table S1...\n")

# A. Total pooled dataset (all cleaned records, 2019–2024, all institutions)
stats_A <- data_with_domains %>%
  summarise(
    N          = n(),
    Mean_Age   = mean(AGE, na.rm = TRUE),
    SD_Age     = sd(AGE,   na.rm = TRUE),
    Pct_Female = sum(GENDER == 2, na.rm = TRUE) / sum(!is.na(GENDER)),
    N_Univ     = sum(Study_Level == "University"),
    N_Voc      = sum(Study_Level == "Vocational")
  )
cat("\nA. Total Pooled Dataset:\n"); print(stats_A)

# B. Training / discovery sample (2019–2021, all institution types, QC-passed, de-duplicated)
stats_B <- qc_data %>%
  filter(!flag_any, Year %in% c("2019", "2020", "2021")) %>%
  summarise(
    N          = n(),
    Mean_Age   = mean(AGE, na.rm = TRUE),
    SD_Age     = sd(AGE,   na.rm = TRUE),
    Pct_Female = sum(GENDER == 2, na.rm = TRUE) / sum(!is.na(GENDER)),
    N_Univ     = sum(Study_Level == "University"),
    N_Voc      = sum(Study_Level == "Vocational")
  )
cat("\nB. Training/Discovery Sample (QC-passed, de-duplicated):\n"); print(stats_B)

# C. Confirmatory sample (raw 2022–2024, all institutions, pre-QC, pre-dedup)
stats_C <- data_with_domains %>%
  filter(Year %in% c("2022", "2023", "2024")) %>%
  summarise(
    N          = n(),
    Mean_Age   = mean(AGE, na.rm = TRUE),
    SD_Age     = sd(AGE,   na.rm = TRUE),
    Pct_Female = sum(GENDER == 2, na.rm = TRUE) / sum(!is.na(GENDER)),
    N_Univ     = sum(Study_Level == "University"),
    N_Voc      = sum(Study_Level == "Vocational")
  )
cat("\nC. Confirmatory Sample (Raw 2022–2024):\n"); print(stats_C)

# D. Confirmatory analytic sample (final QC-passed, de-duplicated, university-only)
stats_D <- final_sample %>%
  summarise(
    N          = n(),
    Mean_Age   = mean(AGE, na.rm = TRUE),
    SD_Age     = sd(AGE,   na.rm = TRUE),
    Pct_Female = sum(GENDER == 2, na.rm = TRUE) / sum(!is.na(GENDER))
  )
cat("\nD. Confirmatory Analytic Sample (University Only, post-QC):\n"); print(stats_D)

stats_D_years <- final_sample %>%
  group_by(Year) %>%
  summarise(
    N          = n(),
    Mean_Age   = mean(AGE, na.rm = TRUE),
    SD_Age     = sd(AGE,   na.rm = TRUE),
    Pct_Female = sum(GENDER == 2, na.rm = TRUE) / sum(!is.na(GENDER)),
    .groups    = "drop"
  )
cat("\nD (Breakdown by Year):\n"); print(stats_D_years)
cat("STEP 4: Complete.\n\n")


# -------------------------------------------------------------------
# STEP 5: CLASSIFIER ANALYSES (Reviewer Points 6 & 7)
# -------------------------------------------------------------------
cat("STEP 5: Running Classifier Analyses...\n")

# 5a. Training data (2019–2021, QC-passed)
training_data_all_types <- qc_data %>%
  filter(!flag_any, Year %in% c("2019", "2020", "2021")) %>%
  select(all_of(domain_cols), Study_Level) %>%
  na.omit() %>%
  as.data.frame()

training_data_all_types[domain_cols] <- lapply(
  training_data_all_types[domain_cols],
  function(z) as.numeric(unlist(z))
)

cat("Classes of domain columns in training_data_all_types:\n")
print(sapply(training_data_all_types[domain_cols], class))

training_mat_all <- as.matrix(training_data_all_types[, domain_cols])

cat("Training Original (Mixed-Pop) Classifier...\n")
set.seed(123)
model_original <- Mclust(
  training_mat_all,
  G          = 2,
  modelNames = "VVV"
)

cat("Training New (Univ-Only) Classifier...\n")
training_data_uni_only <- training_data_all_types %>%
  filter(Study_Level == "University") %>%
  select(all_of(domain_cols)) %>%
  as.data.frame()

training_data_uni_only[domain_cols] <- lapply(
  training_data_uni_only[domain_cols],
  function(z) as.numeric(unlist(z))
)

training_mat_uni <- as.matrix(training_data_uni_only)

set.seed(123)
model_uni_only <- Mclust(
  training_mat_uni,
  G          = 2,
  modelNames = "VVV"
)

# 5a.3 Classify confirmatory sample with both models
final_sample_with_profiles <- final_sample %>%
  select(ID, Year, GENDER, AGE, all_of(domain_cols)) %>%
  tidyr::drop_na(dplyr::all_of(domain_cols)) %>%  # <-- fixed (no na.omit(subset=...))
  as.data.frame()

final_sample_with_profiles[domain_cols] <- lapply(
  final_sample_with_profiles[domain_cols],
  function(z) as.numeric(unlist(z))
)

data_to_classify_domains_only <- final_sample_with_profiles %>%
  select(all_of(domain_cols)) %>%
  as.data.frame()

data_to_classify_mat <- as.matrix(data_to_classify_domains_only)

cat("Classifying sample with both models...\n")
pred_original <- predict(model_original, newdata = data_to_classify_mat)
pred_uni_only <- predict(model_uni_only,  newdata = data_to_classify_mat)

final_sample_with_profiles <- final_sample_with_profiles %>%
  mutate(
    Profile_Original = pred_original$classification,
    Profile_Uni_Only = pred_uni_only$classification
  )

# 5a.4 Concordance and Kappa
cat("Means for Original Model (Check 1 vs 2):\n")
print(model_original$parameters$mean)
cat("Means for Univ-Only Model (Check 1 vs 2):\n")
print(model_uni_only$parameters$mean)

concordance_table_data <- final_sample_with_profiles %>%
  mutate(
    Profile_Orig_Label = if_else(Profile_Original == 1, "Elev Profile", "LS Profile"),
    Profile_Uni_Label  = if_else(Profile_Uni_Only  == 1, "Elev Profile", "LS Profile")
  )

cat("\n--- REVISED TABLE S.SENS.2 (Concordance) ---\n")
contingency_table <- table(
  "Original (Mixed-Pop)" = concordance_table_data$Profile_Orig_Label,
  "New (Univ-Only)"      = concordance_table_data$Profile_Uni_Label
)
print(contingency_table)

kappa_stat <- psych::cohen.kappa(contingency_table)
cat(paste("\nCohen's Kappa:", round(kappa_stat$kappa, 3), "\n"))

# 5b. De-dup Sensitivity (Earliest vs Latest)
cat("\n--- Running Analysis for Point #6 (De-dup Sensitivity) ---\n")

# Earliest (Primary)
data_primary_classified <- final_sample %>%
  dplyr::select(ID, Year, GENDER, AGE, dplyr::all_of(domain_cols)) %>%
  tidyr::drop_na(dplyr::all_of(domain_cols)) %>%  # already correct
  as.data.frame()

data_primary_classified[domain_cols] <- lapply(
  data_primary_classified[domain_cols],
  function(z) as.numeric(unlist(z))
)

newdata_primary <- data_primary_classified %>%
  dplyr::select(dplyr::all_of(domain_cols)) %>%
  as.data.frame()

pred_primary <- predict(
  model_uni_only,
  newdata = as.matrix(newdata_primary)
)

data_primary_classified$Profile <- as.integer(pred_primary$classification)

# Latest (Sensitivity)
data_latest_classified <- final_sample_LATEST %>%
  dplyr::select(ID, Year, GENDER, AGE, dplyr::all_of(domain_cols)) %>%
  tidyr::drop_na(dplyr::all_of(domain_cols)) %>%
  as.data.frame()

data_latest_classified[domain_cols] <- lapply(
  data_latest_classified[domain_cols],
  function(z) as.numeric(unlist(z))
)

newdata_latest <- data_latest_classified %>%
  dplyr::select(dplyr::all_of(domain_cols)) %>%
  as.data.frame()

pred_latest <- predict(
  model_uni_only,
  newdata = as.matrix(newdata_latest)
)

data_latest_classified$Profile <- as.integer(pred_latest$classification)

# --- Table S.SENS.1: De-dup Sensitivity (Earliest vs Latest) ---

data_primary_tab <- data_primary_classified %>%
  mutate(
    Year    = as.factor(Year),
    Profile = as.integer(Profile)
  )

data_latest_tab <- data_latest_classified %>%
  mutate(
    Year    = as.factor(Year),
    Profile = as.integer(Profile)
  )

table_sens_1_primary <- data_primary_tab %>%
  group_by(Year, Profile) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Rule = "Earliest Attempt (Primary)")

table_sens_1_latest <- data_latest_tab %>%
  group_by(Year, Profile) %>%
  summarise(N = n(), .groups = "drop") %>%
  mutate(Rule = "Latest Attempt (Sensitivity)")

table_sens_1_wide <- bind_rows(table_sens_1_primary, table_sens_1_latest) %>%
  mutate(Profile_Label = if_else(Profile == 2, "LS", "Elev")) %>%
  group_by(Year, Rule, Profile_Label) %>%
  summarise(N = sum(N), .groups = "drop") %>%
  group_by(Year, Rule) %>%
  mutate(Proportion = N / sum(N)) %>%
  ungroup() %>%
  select(Year, Rule, Profile_Label, N, Proportion) %>%
  tidyr::pivot_wider(
    names_from  = Profile_Label,
    values_from = c(N, Proportion),
    names_sep   = "_"
  )

table_sens_1_final <- {
  df <- table_sens_1_wide
  
  if (!"N_LS" %in% names(df)) df$N_LS <- 0L
  if (!"N_Elev" %in% names(df)) df$N_Elev <- 0L
  if (!"Proportion_LS" %in% names(df)) df$Proportion_LS <- 0
  if (!"Proportion_Elev" %in% names(df)) df$Proportion_Elev <- 0
  
  df
} %>%
  mutate(
    Total_N = N_LS + N_Elev
  ) %>%
  dplyr::rename(
    `N LS n (%)`            = N_LS,
    `N Elev n (%)`          = N_Elev,
    `Proportion LS n (%)`   = Proportion_LS,
    `Proportion Elev n (%)` = Proportion_Elev
  ) %>%
  arrange(Year, Rule)

cat("\n--- REVISED TABLE S.SENS.1 (De-dup Rule) ---\n")
print(as.data.frame(table_sens_1_final))

cat("STEP 5: Complete.\n\n")


# -------------------------------------------------------------------
# STEP 6: MEASUREMENT INVARIANCE (Year, Sex, Study_Level)
# -------------------------------------------------------------------
cat("STEP 6: Running Measurement Invariance...\n")

base_MI <- qc_data %>% 
  dplyr::filter(!flag_any) %>%
  dplyr::select(dplyr::all_of(all_item_cols), GENDER, Study_Level, Year) %>%
  na.omit() %>%
  as.data.frame()

cfa_model_string <- '
  SOM  =~ Q1 + Q4 + Q12 + Q27 + Q40 + Q42 + Q48 + Q49 + Q52 + Q53 + Q56 + Q58
  OCD  =~ Q3 + Q9 + Q10 + Q28 + Q38 + Q45 + Q46 + Q51 + Q55 + Q65
  IS   =~ Q6 + Q21 + Q34 + Q36 + Q37 + Q41 + Q61 + Q69 + Q73
  DEP  =~ Q5 + Q14 + Q15 + Q20 + Q22 + Q26 + Q29 + Q30 + Q31 + Q32 + Q54 + Q71 + Q79
  ANX  =~ Q2 + Q17 + Q23 + Q33 + Q39 + Q57 + Q72 + Q78 + Q80 + Q86
  PHOB =~ Q13 + Q25 + Q47 + Q50 + Q70 + Q75 + Q82
'

set.seed(123)

## 6a. Invariance across YEAR (ONLY 2022–2024, subsample)
cat("\nRunning MI for 'Year' (2022–2024, subsample)...\n")

data_MI_year <- base_MI %>% 
  dplyr::filter(Year %in% c("2022","2023","2024"))

data_MI_year_s <- data_MI_year[
  sample(nrow(data_MI_year), min(3000, nrow(data_MI_year))), ]

mi_year <- semTools::measurementInvariance(
  model  = cfa_model_string,
  data   = data_MI_year_s,
  group  = "Year",
  strict = FALSE
)

print(mi_year)

mod_indices_scalar_year <- lavaan::modindices(
  mi_year$fit.intercepts,
  op            = "~1",   # intercepts only
  minimum.value = 25,     # only big MIs
  sort.         = TRUE
)

cat("\nTop Intercept MIs for Year:\n")
print(head(mod_indices_scalar_year, 20))


## 6b. Invariance across SEX
cat("\nRunning MI for 'GENDER' (subsample)...\n")

data_MI_sex_s <- base_MI[
  sample(nrow(base_MI), min(3000, nrow(base_MI))), ]

mi_sex <- semTools::measurementInvariance(
  model  = cfa_model_string,
  data   = data_MI_sex_s,
  group  = "GENDER",
  strict = FALSE
)

print(mi_sex)

mod_indices_scalar_sex <- lavaan::modindices(
  mi_sex$fit.intercepts,
  op            = "~1",
  minimum.value = 25,
  sort.         = TRUE
)

cat("\nTop Intercept MIs for Sex:\n")
print(head(mod_indices_scalar_sex, 20))


## 6c. Invariance across Study_Level (Training years only)
cat("\nRunning MI for 'Study_Level' (2019–2021, subsample)...\n")

data_MI_level <- base_MI %>% 
  dplyr::filter(Year %in% c("2019","2020","2021"),
                Study_Level %in% c("University","Vocational"))

data_MI_level_s <- data_MI_level[
  sample(nrow(data_MI_level), min(3000, nrow(data_MI_level))), ]

mi_level <- semTools::measurementInvariance(
  model  = cfa_model_string,
  data   = data_MI_level_s,
  group  = "Study_Level",
  strict = FALSE
)

print(mi_level)

mod_indices_scalar_level <- lavaan::modindices(
  mi_level$fit.intercepts,
  op            = "~1",
  minimum.value = 25,
  sort.         = TRUE
)

cat("\nTop Intercept MIs for Study_Level:\n")
print(head(mod_indices_scalar_level, 20))

cat("STEP 6: Complete.\n\n")









# -------------------------------------------------------------------
# STEP 7: PERMUTATION TESTS FOR GLOBAL STRENGTH (ΔS) AND STRUCTURE (M)
# -------------------------------------------------------------------
cat("STEP 7: Running permutation tests for ΔS and M...\n")

# ------------------------------------------------------------
# 7.0 Build FINAL_ANALYSIS_DATASET from the primary (earliest) classified sample
# ------------------------------------------------------------
domain_cols <- c("SOM","OCD","IS","DEP","ANX","PHOB")  # keep consistent with your paper

FINAL_ANALYSIS_DATASET <- data_primary_classified %>%
  dplyr::select(Year, Profile, dplyr::all_of(domain_cols)) %>%
  tidyr::drop_na(dplyr::all_of(domain_cols)) %>%
  as.data.frame()

FINAL_ANALYSIS_DATASET$Year    <- as.factor(FINAL_ANALYSIS_DATASET$Year)
FINAL_ANALYSIS_DATASET$Profile <- as.integer(FINAL_ANALYSIS_DATASET$Profile)  # 1 = Elev, 2 = LS in your coding

# ------------------------------------------------------------
# 7.1 Helper functions
# ------------------------------------------------------------

# Estimate EBICglasso network (returns weight matrix W)
estimate_glasso <- function(X) {
  X <- as.matrix(X)
  S <- cor(X, use = "pairwise.complete.obs")
  W <- qgraph::EBICglasso(S, n = nrow(X), gamma = 0.50)
  W
}

# Global strength S = sum of absolute edge weights (upper triangle)
compute_global_strength <- function(W) {
  sum(abs(W[upper.tri(W)]), na.rm = TRUE)
}

# Structure difference statistic M = max absolute edgewise difference
compute_M <- function(W_LS, W_Elev) {
  max(abs(W_Elev - W_LS)[upper.tri(W_LS)], na.rm = TRUE)
}

# ------------------------------------------------------------
# 7.2 Main function: run permutation test for one year
# ------------------------------------------------------------
run_perm_for_year <- function(target_year, data, B = 50000) {
  
  cat("\n--- Year:", target_year, " ---\n")
  
  # IMPORTANT: create dat_year inside the function (prevents the error you saw)
  dat_year <- data[data$Year == target_year, , drop = FALSE]
  if (nrow(dat_year) == 0) {
    stop("No rows found for Year = ", target_year)
  }
  
  # Your coding: Profile 1 = Elevated, Profile 2 = LS
  dat_LS   <- dat_year[dat_year$Profile == 2, , drop = FALSE]
  dat_Elev <- dat_year[dat_year$Profile == 1, , drop = FALSE]
  
  N_LS   <- nrow(dat_LS)
  N_Elev <- nrow(dat_Elev)
  
  cat("N_LS =", N_LS, "| N_Elev =", N_Elev, "\n")
  
  if (N_LS < 20 || N_Elev < 20) {
    stop("Too few cases per group in ", target_year, " (need at least 20 each).")
  }
  
  X_LS   <- as.matrix(dat_LS[, domain_cols, drop = FALSE])
  X_Elev <- as.matrix(dat_Elev[, domain_cols, drop = FALSE])
  
  # ---- Observed statistics ----
  W_LS   <- estimate_glasso(X_LS)
  W_Elev <- estimate_glasso(X_Elev)
  
  S_LS   <- compute_global_strength(W_LS)
  S_Elev <- compute_global_strength(W_Elev)
  
  Delta_obs <- S_Elev - S_LS
  M_obs     <- compute_M(W_LS, W_Elev)
  
  cat("Observed: S_LS =", round(S_LS, 4),
      "| S_Elev =", round(S_Elev, 4),
      "| ΔS =", round(Delta_obs, 4),
      "| M =", round(M_obs, 4), "\n")
  
  # ---- Permutation distributions ----
  X_all <- rbind(X_LS, X_Elev)
  
  # Use internal permutation labels (1=LS, 2=Elev) for X_all
  group_vec <- c(rep(1, N_LS), rep(2, N_Elev))
  
  perm_deltas <- numeric(B)
  perm_M      <- numeric(B)
  
  for (b in seq_len(B)) {
    
    perm_labels <- sample(group_vec)
    
    X_LS_b   <- X_all[perm_labels == 1, , drop = FALSE]
    X_Elev_b <- X_all[perm_labels == 2, , drop = FALSE]
    
    W_LS_b   <- estimate_glasso(X_LS_b)
    W_Elev_b <- estimate_glasso(X_Elev_b)
    
    S_LS_b   <- compute_global_strength(W_LS_b)
    S_Elev_b <- compute_global_strength(W_Elev_b)
    
    perm_deltas[b] <- S_Elev_b - S_LS_b
    perm_M[b]      <- compute_M(W_LS_b, W_Elev_b)
  }
  
  # ΔS: two-sided test
  extreme_Delta <- sum(abs(perm_deltas) >= abs(Delta_obs))
  p_Delta <- (extreme_Delta + 1) / (B + 1)
  
  # M: one-sided test (M >= 0)
  extreme_M <- sum(perm_M >= M_obs)
  p_M <- (extreme_M + 1) / (B + 1)
  
  cat("p(ΔS) =", signif(p_Delta, 6), " | extreme =", extreme_Delta, " / ", B, "\n")
  cat("p(M)  =", signif(p_M, 6),     " | extreme =", extreme_M,     " / ", B, "\n")
  
  return(list(
    Year          = target_year,
    S_LS          = S_LS,
    S_Elev        = S_Elev,
    Delta_S       = Delta_obs,
    p_Delta_S     = p_Delta,
    M             = M_obs,
    p_M           = p_M,
    B             = B,
    extreme_Delta = extreme_Delta,
    extreme_M     = extreme_M
  ))
}

# ------------------------------------------------------------
# 7.3 Run for each confirmatory year
# ------------------------------------------------------------
B_final <- 50000

res_2022 <- run_perm_for_year("2022", FINAL_ANALYSIS_DATASET, B = B_final)
res_2023 <- run_perm_for_year("2023", FINAL_ANALYSIS_DATASET, B = B_final)
res_2024 <- run_perm_for_year("2024", FINAL_ANALYSIS_DATASET, B = B_final)

results_list <- list(res_2022, res_2023, res_2024)

results_df <- do.call(rbind, lapply(results_list, as.data.frame))
cat("\n--- RESULTS (ΔS and M) ---\n")
print(results_df, digits = 6)

cat("\nSTEP 7: Complete.\n\n")








library(tidyverse)

df_stats <- tibble(
  Year   = factor(c("2022","2023","2024"), levels = c("2022","2023","2024")),
  DeltaS = c(0.651518, 0.487611, 0.590367),
  M      = c(0.183782, 0.161906, 0.259731)
)

df_long <- df_stats %>%
  pivot_longer(cols = c(DeltaS, M),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(
    Metric = recode(Metric,
                    "DeltaS" = "ΔS (Elev − LS global strength)",
                    "M"      = "M (maximum edgewise |Δ|)")
  )

fig2 <- ggplot(df_long, aes(x = Year, y = Value, group = Metric)) +
  geom_line() +
  geom_point(size = 3) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  labs(x = "Year", y = "Value", title = "Year-over-year pattern of ΔS and M") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

fig2
ggsave("Figure2_DeltaS_M.png", fig2, width = 8, height = 4, dpi = 300)



# -------------------------------------------------------------------
# STEP 8: FINAL TASKS
# -------------------------------------------------------------------
cat("\n--- SCRIPT COMPLETE ---\n")
cat("All analyses for your RtR are complete (using long-string >= 20).\n")
cat("You can now fill in your tables and revise your manuscript.\n")

# after pred_primary on data_primary_classified

library(dplyr)

data_primary_classified <- data_primary_classified %>%
  mutate(Profile = as.integer(pred_primary$classification))

# posterior probabilities from model_uni_only
post <- predict(model_uni_only, newdata = as.matrix(newdata_primary))$z

K <- 2
max_post <- apply(post, 1, max)
entropy <- -rowSums(post * log(post)) / log(K)

summary_df <- data_primary_classified %>%
  mutate(
    max_post = max_post,
    entropy_norm = entropy,
    class_entropy = 1 - entropy
  ) %>%
  group_by(Year) %>%
  summarise(
    N = n(),
    LS_n = sum(Profile == 2),
    Elev_n = sum(Profile == 1),
    LS_pct = LS_n / N * 100,
    Elev_pct = Elev_n / N * 100,
    mean_max_post = mean(max_post),
    mean_entropy_norm = mean(entropy_norm),
    mean_class_entropy = mean(class_entropy),
    .groups = "drop"
  )

# ----------------------------------------------------------
# Figure 2: Year-over-year ΔS and M
# ----------------------------------------------------------

# 1. Load packages
library(tidyverse)

# 2. Enter your summary statistics
#    Use the FINAL results from your permutation tests.
#    I’ve filled in ΔS with your STEP 7 (B = 50,000) output.
#    PLEASE replace the M values with the ones from your
#    final structure-difference permutation step.

df_stats <- tibble(
  Year   = factor(c("2022", "2023", "2024"),
                  levels = c("2022", "2023", "2024")),
  DeltaS = c(0.6515, 0.4876, 0.5904),   # ΔS (S_Elev - S_LS)
  M      = c(0.155, 0.165, 0.175)       # <-- REPLACE with your final M values
)

# 3. Reshape to long format for faceting
df_long <- df_stats %>%
  pivot_longer(cols = c(DeltaS, M),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(
    Metric = recode(Metric,
                    "DeltaS" = "ΔS (Elev − LS global strength)",
                    "M"      = "M (maximum edgewise |Δ|)")
  )

# 4. Plot (two panels: ΔS and M)
fig2 <- ggplot(df_long, aes(x = Year, y = Value, group = Metric)) +
  geom_line() +
  geom_point(size = 3) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  labs(
    x = "Year",
    y = "Value",
    title = "Year-over-year stability of ΔS and M"
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

# 5. Print to RStudio viewer
fig2

# 6. Save to file (adjust path / size as needed)
ggsave("Figure2_DeltaS_M.png", fig2,
       width = 8, height = 4, dpi = 300)

## ------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(qgraph)

## ------------------------------------------------------------
## 1. Recreate FINAL_ANALYSIS_DATASET from your main data
##    (this assumes you already have `data_primary_classified`
##     in your workspace with Year, Profile and domain columns)
## ------------------------------------------------------------

domain_cols <- c("ANX","DEP","IS","PHOB","OCD","SOM")

FINAL_ANALYSIS_DATASET <- data_primary_classified %>%
  dplyr::select(Year, Profile, dplyr::all_of(domain_cols)) %>%
  tidyr::drop_na(dplyr::all_of(domain_cols)) %>%
  as.data.frame()

# Clean types for safety
FINAL_ANALYSIS_DATASET$Year    <- as.factor(FINAL_ANALYSIS_DATASET$Year)
FINAL_ANALYSIS_DATASET$Profile <- as.integer(FINAL_ANALYSIS_DATASET$Profile)
# (your coding: 1 = Elevated, 2 = LS)

## ------------------------------------------------------------
## 2. Prepare data just for 2022–2024 and recode Profile labels
## ------------------------------------------------------------

dat <- FINAL_ANALYSIS_DATASET %>%
  dplyr::filter(Year %in% c("2022","2023","2024")) %>%
  dplyr::mutate(
    Year = factor(Year, levels = c("2022","2023","2024")),
    Profile = factor(
      Profile,
      levels = c(2, 1),                 # 2=LS, 1=Elevated in your coding
      labels = c("LS", "Elevated")
    )
  )

## ------------------------------------------------------------
## 3. Get a common node layout from the pooled data
## ------------------------------------------------------------

cor_pool <- cor(dat[, domain_cols], use = "pairwise.complete.obs")
W_pool   <- qgraph::EBICglasso(cor_pool, n = nrow(dat), gamma = 0.50)

set.seed(123)                  # for reproducible layout
q_pool     <- qgraph::qgraph(W_pool, layout = "spring", DoNotPlot = TRUE)
layout_mat <- q_pool$layout    # reuse this layout for all panels

## ------------------------------------------------------------
## 4. Estimate networks for each Year × Profile
## ------------------------------------------------------------

estimate_glasso <- function(df) {
  S <- cor(df[, domain_cols], use = "pairwise.complete.obs")
  qgraph::EBICglasso(S, n = nrow(df), gamma = 0.50)
}

nets <- list()

for (yr in levels(dat$Year)) {
  for (prof in levels(dat$Profile)) {
    tmp <- dat %>% dplyr::filter(Year == yr, Profile == prof)
    W   <- estimate_glasso(tmp)
    nets[[paste(yr, prof, sep = "_")]] <- W
  }
}

# Common max absolute edge weight so widths are comparable
max_edge <- max(sapply(nets, function(W)
  max(abs(W[upper.tri(W)]), na.rm = TRUE)
))

## ------------------------------------------------------------
## 5. Plot 3×2 grid (Figure 1)
## ------------------------------------------------------------

pdf("Figure1_profile_networks.pdf", width = 8, height = 10)  # or use png()
par(mfrow = c(3, 2), mar = c(1, 1, 2, 1))

for (yr in levels(dat$Year)) {
  for (prof in c("LS", "Elevated")) {
    
    W <- nets[[paste(yr, prof, sep = "_")]]
    main_title <- paste(yr, "—", prof)
    
    qgraph::qgraph(
      W,
      layout     = layout_mat,
      maximum    = max_edge,
      edge.color = c("red", "darkgreen"),  # negative, positive
      labels     = domain_cols,
      label.cex  = 1,
      vsize      = 13,
      rescale    = FALSE,
      title      = main_title
    )
  }
}

dev.off()

getwd()

library(dplyr)

# 1. Read your raw confirmatory data (universities only)
# Replace with your real file paths and filters
dat <- read.csv("confirmatory_2022_2024_university_raw.csv")

# Assume you have:
# - year: 2022, 2023, 2024
# - time_sec: completion time in seconds
# - scl_items: columns item1:item90 (change names)
# - id: student ID or hash

scl_cols <- paste0("item", 1:90)  # <-- change to your real item column names

dat_qc <- dat %>%
  mutate(
    # (a) Speeder flag (time < 5 minutes)
    qc_speed = time_sec < 5 * 60,   # change time variable name if needed
    
    # (b) Long-string flag: max run >= 20 identical answers
    # This assumes items are numeric 1–5.
    qc_long = apply(select(., all_of(scl_cols)), 1, function(x) {
      rle_x <- rle(x)
      max(rle_x$lengths) >= 20
    }),
    
    # (c) Low-variance flag: within-person SD < 0.10
    qc_lowvar = apply(select(., all_of(scl_cols)), 1, function(x) {
      sd(x, na.rm = TRUE) < 0.10
    }),
    
    # Any QC flag
    qc_any = qc_speed | qc_long | qc_lowvar
  )


#figures2#
library(mclust)

# Use the training data you already created:
# training_data_all_types (2019–2021, QC-passed, all institutions)
training_mat_all <- as.matrix(training_data_all_types[, domain_cols])

# BIC for 1–6 classes, all covariance structures
bic_all <- mclustBIC(
  training_mat_all,
  G          = 1:6,
  modelNames = c("EII","VII","EEI","VEI","EVI","VVI",
                 "EEE","EEV","VEV","VVV")
)

# ICL for 1–6 classes, same structures
icl_all <- mclustICL(
  training_mat_all,
  G          = 1:6,
  modelNames = c("EII","VII","EEI","VEI","EVI","VVI",
                 "EEE","EEV","VEV","VVV")
)

# Recreate the side-by-side plot (similar to your S2)
par(mfrow = c(1,2))
plot(bic_all)   # left panel: BIC
plot(icl_all)   # right panel: ICL
par(mfrow = c(1,1))


#figures3#
library(qgraph)

# Example: pooled confirmatory analytic sample (all years, all profiles)
dat_net <- final_sample[, domain_cols]

S <- cor(dat_net, use = "pairwise.complete.obs")
W <- qgraph::EBICglasso(S, n = nrow(dat_net), gamma = 0.50)

pdf("Figure_S3_PooledNetwork.pdf", width = 6, height = 6)
qgraph(
  W,
  layout    = "spring",
  labels    = domain_cols,
  vsize     = 15,
  title     = "Regularised Partial-Correlation Network"
)
dev.off()


library(qgraph)

domain_cols <- c("ANX","DEP","IS","PHOB","OCD","SOM")

dat_net <- final_sample[, domain_cols]

S <- cor(dat_net, use = "pairwise.complete.obs")
W <- qgraph::EBICglasso(S, n = nrow(dat_net), gamma = 0.50)

# (Optional) scale node size by strength centrality
node_strength <- qgraph::centrality_auto(W)$node.centrality[,"Strength"]
node_strength <- scales::rescale(node_strength, to = c(15, 30))  # size range

# Nice colour palette (you can change if you like)
node_cols <- c("#E64B35", "#4DBBD5", "#00A087",
               "#3C5488", "#F39B7F", "#8491B4")

qgraph(
  W,
  layout       = "spring",
  labels       = domain_cols,
  vsize        = node_strength,   # or a constant like 20 if you prefer
  color        = node_cols,
  edge.color   = "grey70",
  edge.labels  = TRUE,
  edge.label.cex = 0.7,
  title        = "Regularised Partial-Correlation Network"
)

#Figures4 and s5#
library(bootnet)
domain_cols <- c("ANX","DEP","IS","PHOB","OCD","SOM")
dat_net <- final_sample[, domain_cols]

# 1. Estimate network (same as for Figure S3)
net_glasso <- estimateNetwork(
  dat_net,
  default   = "EBICglasso",
  corMethod = "cor"   # or "cor_auto" if you prefer
)

# 2. Non-parametric bootstrap
boot_edges <- bootnet(
  net_glasso,
  nBoots = 2000,
  type   = "nonparametric"
)

pdf("Figure_S4_EdgeAccuracy.pdf", width = 8, height = 5)
plot(
  boot_edges,
  labels = TRUE,
  order  = "sample"   # <-- VALID: "id", "sample", or "mean"
)
dev.off()

pdf("Figure_S5_EdgeOrder.pdf", width = 8, height = 5)
plot(
  boot_edges,
  "edge",
  plot        = "difference",
  onlyNonZero = TRUE,
  order       = "sample"  # again: "sample" is OK
)
dev.off()

# Figure S4
plot(boot_edges, labels = TRUE, order = "sample")

# Figure S5
plot(boot_edges, "edge",
     plot        = "difference",
     onlyNonZero = TRUE,
     order       = "sample")


# helper: elevated component = higher average mean across the 6 domains
get_elev_comp <- function(m) which.max(colMeans(m$parameters$mean))

elev_orig <- get_elev_comp(model_original)
elev_uni  <- get_elev_comp(model_uni_only)

# labels aligned to "Elev" based on component means
lab_orig <- ifelse(pred_original$classification == elev_orig, "Elev", "LS")
lab_uni  <- ifelse(pred_uni_only$classification  == elev_uni,  "Elev", "LS")

tab <- table(`Original (Mixed-Pop)` = lab_orig,
             `New (Univ-Only)`      = lab_uni)

addmargins(tab)          # check totals match Table S3
psych::cohen.kappa(tab)  # kappa on correctly aligned labels

data_primary_classified %>%
  group_by(Profile) %>%
  summarise(across(all_of(domain_cols), mean, na.rm=TRUE))


# Helper 3: compute M (maximum absolute edgewise difference)
compute_M <- function(W_LS, W_Elev) {
  max(abs(W_Elev - W_LS)[upper.tri(W_LS)], na.rm = TRUE)
}

run_perm_for_year <- function(target_year, data, B = 50000) {
  cat("\n--- Running permutation test for Year:", target_year, "---\n")
  
  dat_year <- data[data$Year == target_year, , drop = FALSE]
  if (nrow(dat_year) == 0) {
    warning("No data for target year: ", target_year)
    return(NULL)
  }
  
  # Your coding: Profile 1 = Elevated, Profile 2 = LS
  dat_LS   <- dat_year[dat_year$Profile == 2, , drop = FALSE]
  dat_Elev <- dat_year[dat_year$Profile == 1, , drop = FALSE]
  
  N_LS   <- nrow(dat_LS)
  N_Elev <- nrow(dat_Elev)
  cat("N (LS Profile):", N_LS, "| N (Elevated Profile):", N_Elev, "\n")
  
  if (N_LS < 20 || N_Elev < 20) {
    warning("Year ", target_year, ": too few cases per profile, skipping.")
    return(NULL)
  }
  
  X_LS   <- as.matrix(dat_LS[, domain_cols])
  X_Elev <- as.matrix(dat_Elev[, domain_cols])
  
  # Observed networks
  W_LS   <- estimate_glasso(X_LS)
  W_Elev <- estimate_glasso(X_Elev)
  
  # Observed ΔS
  S_LS   <- compute_global_strength(W_LS)
  S_Elev <- compute_global_strength(W_Elev)
  Delta_obs <- S_Elev - S_LS
  
  # Observed M
  M_obs <- compute_M(W_LS, W_Elev)
  
  cat("Observed S_LS =", round(S_LS, 3),
      "| S_Elev =", round(S_Elev, 3),
      "| ΔS =", round(Delta_obs, 3),
      "| M =", round(M_obs, 3), "\n")
  
  # Permutation distributions
  X_all <- rbind(X_LS, X_Elev)
  
  # IMPORTANT: Keep your internal coding consistent.
  # Here we code 1 = LS, 2 = Elev, matching your original STEP 7.
  group_vec <- c(rep(1, N_LS), rep(2, N_Elev))
  
  perm_deltas <- numeric(B)
  perm_M      <- numeric(B)
  
  for (b in seq_len(B)) {
    perm_labels <- sample(group_vec)
    
    X_LS_b   <- X_all[perm_labels == 1, , drop = FALSE]
    X_Elev_b <- X_all[perm_labels == 2, , drop = FALSE]
    
    W_LS_b   <- estimate_glasso(X_LS_b)
    W_Elev_b <- estimate_glasso(X_Elev_b)
    
    S_LS_b   <- compute_global_strength(W_LS_b)
    S_Elev_b <- compute_global_strength(W_Elev_b)
    
    perm_deltas[b] <- S_Elev_b - S_LS_b
    perm_M[b]      <- compute_M(W_LS_b, W_Elev_b)
  }
  
  # ΔS: two-sided p-value
  extreme_Delta <- sum(abs(perm_deltas) >= abs(Delta_obs))
  p_Delta <- (extreme_Delta + 1) / (B + 1)
  
  # M: one-sided p-value (because M >= 0)
  extreme_M <- sum(perm_M >= M_obs)
  p_M <- (extreme_M + 1) / (B + 1)
  
  cat("Permutation p-value for |ΔS|:", signif(p_Delta, 6),
      "(extreme count =", extreme_Delta, "of B =", B, ")\n")
  cat("Permutation p-value for M:", signif(p_M, 6),
      "(extreme count =", extreme_M, "of B =", B, ")\n")
  
  list(
    year            = target_year,
    S_LS            = S_LS,
    S_Elev          = S_Elev,
    Delta_S         = Delta_obs,
    p_Delta_S       = p_Delta,
    M               = M_obs,
    p_M             = p_M,
    B               = B,
    extreme_Delta   = extreme_Delta,
    extreme_M       = extreme_M
  )
}

