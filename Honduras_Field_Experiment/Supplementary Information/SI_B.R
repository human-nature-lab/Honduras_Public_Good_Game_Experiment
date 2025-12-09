# =========================
# Descriptive summaries for pgg_data Tables B1, B2 and B3
# =========================

# Load minimal packages for data manipulation and reshaping
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
})

# Read in the main panel from CSV (expected to contain all variables of interest)
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# Sanity check: fail early if pgg_data was not created
stopifnot(exists("pgg_data"))

# ---- variables used in the paper table ----
# This is the ordered list of variables to appear in the descriptive tables (SI Appendix B)
vars_of_interest <- c(
  "contributing","age","gender","friends","adversaries","FI","marital_status",
  "b0100","b0200","b0600","access_routes","network_density_fr",
  "network_density_adv","network_size","finauto_q_2_temp"
)

# Keep only those that actually exist in pgg_data (robust to missing columns)
vars_of_interest <- intersect(vars_of_interest, names(pgg_data))

# ---- 1) Normalize types/codings on a *copy* (won't overwrite your original) ----
# Create a cleaned copy df with harmonized codings (binary, categorical, numeric)
df <- pgg_data %>%
  # Only keep variables of interest, in this order
  dplyr::select(all_of(vars_of_interest)) %>%
  mutate(
    # gender: normalize to binary 1 = male/man; 0 = female/woman; else try integer
    gender = case_when(
      tolower(as.character(gender)) %in% c("man","male","1","m") ~ 1L,
      tolower(as.character(gender)) %in% c("woman","female","0","f") ~ 0L,
      TRUE ~ suppressWarnings(as.integer(gender))
    ),
    # FI already 0/1 in your pipeline; coerce to integer just in case
    FI = suppressWarnings(as.integer(FI)),
    # marital_status: 1 = married/civil union, 0 = single/divorced/widowed; else integer
    marital_status = case_when(
      tolower(as.character(marital_status)) %in% c("married","civil union","1") ~ 1L,
      tolower(as.character(marital_status)) %in% c("single","divorced","widowed","0") ~ 0L,
      TRUE ~ suppressWarnings(as.integer(marital_status))
    ),
    # b0200 (indigenous status): ensure 0/1 coding
    b0200 = case_when(
      tolower(as.character(b0200)) %in% c("1","yes","yes, lenca","yes, maya chorti","other") ~ 1L,
      tolower(as.character(b0200)) %in% c("0","no") ~ 0L,
      TRUE ~ suppressWarnings(as.integer(b0200))
    ),
    # b0600 (religion): 2 = Catholic, 1 = Protestant, 0 = no religion, else integer
    b0600 = case_when(
      tolower(as.character(b0600)) %in% c("2","catholic") ~ 2L,
      tolower(as.character(b0600)) %in% c("1","protestant") ~ 1L,
      tolower(as.character(b0600)) %in% c("0","no religion") ~ 0L,
      TRUE ~ suppressWarnings(as.integer(b0600))
    ),
    # education b0100 as integer (0..13 per your coding in the paper)
    b0100 = suppressWarnings(as.integer(b0100)),
    # financial autonomy binary (0/1)
    finauto_q_2_temp = suppressWarnings(as.integer(finauto_q_2_temp)),
    # numeric coercions (safe if already numeric; otherwise try as.numeric)
    friends = suppressWarnings(as.numeric(friends)),
    adversaries = suppressWarnings(as.numeric(adversaries)),
    network_size = suppressWarnings(as.numeric(network_size)),
    network_density_fr = suppressWarnings(as.numeric(network_density_fr)),
    network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
    access_routes = suppressWarnings(as.numeric(access_routes)),
    # village_wealth_index_median_w3 = suppressWarnings(as.numeric(village_wealth_index_median_w3)),
    age = suppressWarnings(as.numeric(age)),
    contributing = suppressWarnings(as.numeric(contributing))
    #    IH = suppressWarnings(as.numeric(IH))
  )

# ---- 2) Type detection helpers ----
# Helper: check whether a vector is effectively integer-valued (within small numerical tolerance)
is_integerish <- function(x) {
  xnum <- suppressWarnings(as.numeric(x))
  all(is.na(x) | (!is.na(xnum) & abs(xnum - round(xnum)) < 1e-8))
}

# Helper: number of unique non-missing values
n_unique <- function(x) length(unique(x[!is.na(x)]))

# binary if exactly two unique non-NA values after standardization
# Works for factors/characters (two labels) or numeric/logical (two distinct codes)
is_binary <- function(x) {
  if (is.factor(x) || is.character(x)) {
    return(n_unique(trimws(as.character(x))) == 2)
  }
  if (is.numeric(x) || is.logical(x) || is.integer(x)) {
    ux <- unique(x[!is.na(x)])
    return(length(ux) == 2 || all(ux %in% c(0,1)))
  }
  FALSE
}

# numeric with a small number of integer codes -> treat as categorical (e.g., 0/1/2, Likert-like)
is_small_code_numeric <- function(x, k = 10) {
  if (!is.numeric(x) && !is.integer(x)) return(FALSE)
  if (!is_integerish(x)) return(FALSE)
  n_unique(x) <= k && !is_binary(x)
}

# Infer a semantic "type" label used to decide which summary to compute
infer_type <- function(name, x) {
  # Force specific variables to desired classes (align with your LaTeX table in the SI)
  if (name %in% c("gender","FI","marital_status","b0200","finauto_q_2_temp")) return("binary")
  if (name %in% c("b0600")) return("categorical")
  if (name %in% c("friends","adversaries","network_size")) return("count")
  if (name %in% c("contributing","age","access_routes","network_density_fr","network_density_adv")) return("continuous")
  # Fallback heuristics for everything else
  if (is_binary(x)) return("binary")
  if (is.factor(x) || is.character(x) || is_small_code_numeric(x)) return("categorical")
  if (is_integerish(x)) return("count")
  if (is.numeric(x)) return("continuous")
  "other"
}

# Map each variable in df to its inferred type
type_map <- tibble(
  variable = names(df),
  type = map2_chr(names(df), df, infer_type)
)

# ---- 3) Summaries ----

# 3a) numeric / count
# Select all variables marked as continuous or count for numeric summaries
numeric_vars <- type_map %>% filter(type %in% c("continuous","count")) %>% pull(variable)

numeric_summary <- df %>%
  dplyr::select(all_of(numeric_vars)) %>%
  # Convert to long format: one row per (variable, value)
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  mutate(value = suppressWarnings(as.numeric(value))) %>%
  group_by(variable) %>%
  # Standard numeric summary: N, mean, sd, min, max
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd   = sd(value, na.rm = TRUE),
    min  = min(value, na.rm = TRUE),
    max  = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Attach the inferred type for each variable
  left_join(type_map, by = "variable") %>%
  relocate(type, .after = variable) %>%
  # Preserve the original ordering from vars_of_interest in the output
  arrange(match(variable, vars_of_interest))

# 3b) binary (coerce to character before pivot to avoid mixed-type errors)
# Variables classified as binary
binary_vars <- type_map %>% filter(type == "binary") %>% pull(variable)

binary_summary <- df %>%
  dplyr::select(all_of(binary_vars)) %>%
  # Coerce to character to avoid issues mixing 0/1 and labels in pivot_longer
  mutate(across(everything(), ~ as.character(.))) %>%           # <-- prevent type conflict
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  # Drop missing/blank levels
  filter(!is.na(value) & value != "") %>%
  group_by(variable, value) %>%
  # Count each level of each binary variable
  summarise(n = n(), .groups = "drop_last") %>%
  # Within each variable, convert counts to proportions
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  # Order variables according to vars_of_interest, then by level frequency
  arrange(match(variable, vars_of_interest), desc(n))

# 3c) categorical
# Variables classified as categorical (more than 2 levels)
categorical_vars <- type_map %>% filter(type == "categorical") %>% pull(variable)

categorical_summary <- df %>%
  dplyr::select(all_of(categorical_vars)) %>%
  # Again, standardize to character before pivoting
  mutate(across(everything(), ~ as.character(.))) %>%           # unify types
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  # Drop missing/blank categories
  filter(!is.na(value) & value != "") %>%
  group_by(variable, value) %>%
  # Count each category for each variable
  summarise(n = n(), .groups = "drop_last") %>%
  # Convert to within-variable proportions
  mutate(pct = n / sum(n)) %>%
  ungroup() %>%
  # Order categories with the same vars_of_interest logic
  arrange(match(variable, vars_of_interest), desc(n))

# ---- 4) Print quick looks ----
# Console summary for interactive use
cat("\n== Variable types ==\n"); print(type_map, n = 100)
cat("\n== Numeric / count summary ==\n"); print(numeric_summary, n = 200)
cat("\n== Binary counts ==\n"); print(binary_summary, n = 200)
cat("\n== Categorical counts ==\n"); print(categorical_summary, n = 200)

# ---- 5) Optional: LaTeX tables (booktabs) ----
# If knitr is available, also emit LaTeX versions of the four tables (for SI / appendix)
if (requireNamespace("knitr", quietly = TRUE)) {
  cat("\n\n% LaTeX tables (booktabs)\n")
  print(knitr::kable(type_map, format = "latex", booktabs = TRUE,
                     caption = "Detected variable types in the analysis set."))
  print(knitr::kable(numeric_summary, format = "latex", booktabs = TRUE, digits = 3,
                     caption = "Summary statistics (continuous/count variables)."))
  print(knitr::kable(binary_summary, format = "latex", booktabs = TRUE, digits = 3,
                     caption = "Binary variables: counts and proportions by level."))
  print(knitr::kable(categorical_summary, format = "latex", booktabs = TRUE, digits = 3,
                     caption = "Categorical variables: counts and proportions by category."))
}

# ---- 6) Optional: write CSVs for SI/appendix ----
# Uncomment these lines if you want to export CSVs for replication package or appendix
# write.csv(type_map, "summary_types.csv", row.names = FALSE)
# write.csv(numeric_summary, "summary_numeric.csv", row.names = FALSE)
# write.csv(binary_summary, "summary_binary.csv", row.names = FALSE)
# write.csv(categorical_summary, "summary_categorical.csv", row.names = FALSE)
