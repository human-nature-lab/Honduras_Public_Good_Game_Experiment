## ------------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------------
pkgs <- c(
  "dplyr", "tidyr", "tibble",
  "fastcluster", "cluster", "factoextra",
  "RColorBrewer", "ggplot2",
  "lme4", "lmerTest", "conflicted", "modelsummary"
)

new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

## Resolve common conflicts
conflict_prefer("filter",     "dplyr")
conflict_prefer("select",     "dplyr")
conflict_prefer("lag",        "dplyr")
conflict_prefer("rename",     "dplyr")
conflict_prefer("first",      "dplyr")
conflict_prefer("summarise",  "dplyr")
conflict_prefer("count",      "dplyr")
conflict_prefer("inner_join", "dplyr")
conflict_prefer("left_join",  "dplyr")
conflicts_prefer(lmerTest::lmer)  # use lmerTest::lmer by default

## ------------------------------------------------------------------
## 1. Constants & data
## ------------------------------------------------------------------
# Game parameters (Section 2)
x <- 12L   # endowment
m <- 2L    # multiplication factor
N <- 5L    # group size

pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

## ------------------------------------------------------------------
## 2. Public–goods accounting and lags (C^t, s_i^t, u_i^t, etc.)
## ------------------------------------------------------------------

# Group totals by round: C^t and N in case of missing/variable group size
pgg_data <- pgg_data %>%
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    total_contributions = sum(contributing, na.rm = TRUE),
    N = dplyr::n()
  ) %>%
  dplyr::ungroup()

# Individual share s_i^t, payoff u_i^t, and keeping
pgg_data <- pgg_data %>%
  dplyr::mutate(
    s_i_t   = (m * total_contributions) / N,
    u_i_t   = x - contributing + s_i_t,
    keeping = x - contributing
  )

# Lagged own contribution and lagged group total
pgg_data <- pgg_data %>%
  dplyr::arrange(group, player, round_n) %>%
  dplyr::group_by(group, player) %>%
  dplyr::mutate(
    contributing_t_minus_1        = dplyr::lag(contributing, n = 1),
    total_contributions_t_minus_1 = dplyr::lag(total_contributions, n = 1)
  ) %>%
  dplyr::ungroup()

# Lagged leave‑one‑out peer mean c_peer^{t-1}
pgg_data <- pgg_data %>%
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    average_contribution_others_t_minus_1 =
      (total_contributions_t_minus_1 - contributing_t_minus_1) / (N - 1)
  ) %>%
  dplyr::ungroup()

# First round has no lag: set peer lag to 0 (or drop those obs later)
pgg_data$average_contribution_others_t_minus_1[
  is.na(pgg_data$average_contribution_others_t_minus_1)
] <- 0

# Total welfare W^t = sum of squared contributions (used elsewhere)
pgg_data <- pgg_data %>%
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(W_t = sum(contributing^2)) %>%
  dplyr::ungroup()

## ------------------------------------------------------------------
## 3. Basic type conversions (no scaling / no centering)
## ------------------------------------------------------------------

# Make sure numeric vars are numeric
pgg_data$friends             <- as.numeric(pgg_data$friends)
pgg_data$adversaries         <- as.numeric(pgg_data$adversaries)
pgg_data$age                 <- as.numeric(pgg_data$age)
pgg_data$network_density_fr  <- as.numeric(pgg_data$network_density_fr)
pgg_data$network_density_adv <- as.numeric(pgg_data$network_density_adv)
pgg_data$network_size        <- as.numeric(pgg_data$network_size)
pgg_data$b0100               <- as.numeric(pgg_data$b0100)
pgg_data$access_routes       <- as.numeric(pgg_data$access_routes)

# Factors / dummies
pgg_data$gender         <- factor(pgg_data$gender)
pgg_data$marital_status <- factor(pgg_data$marital_status)

# Religion: b0600 (2 = Catholic as reference, following Table 1)
pgg_data$b0600 <- stats::relevel(
  factor(pgg_data$b0600),
  ref = "2"   # Catholic
)

## ------------------------------------------------------------------
## 4. Build analysis dataframe for §5.2 (Eq. 10)
## ------------------------------------------------------------------
## c_own_lag   = own contribution at t-1
## c_group_lag = lagged LOO peer mean in Lempiras (0–12)
## round_n     = numeric round counter
## round_f     = factor for village×round random intercepts

dat <- pgg_data %>%
  dplyr::arrange(village_code, group, player, round_n) %>%
  dplyr::mutate(
    c_own_lag   = contributing_t_minus_1,
    c_group_lag = average_contribution_others_t_minus_1,
    round_n     = as.numeric(round_n),
    round_f     = factor(round_n),
    player      = factor(player),
    group       = factor(group),
    village_code   = factor(village_code),
    gender         = factor(gender),
    marital_status = factor(marital_status),
    b0600          = stats::relevel(factor(b0600), ref = "2") # Catholic
  ) %>%
  # Drop rows with missing lags (i.e., round 1 and any incomplete histories)
  dplyr::filter(!is.na(c_own_lag), !is.na(c_group_lag))

## Quick sanity checks on lags
stopifnot(all(dat$c_own_lag   >= 0 & dat$c_own_lag   <= x, na.rm = TRUE))
stopifnot(all(dat$c_group_lag >= 0 & dat$c_group_lag <= x, na.rm = TRUE))

## ------------------------------------------------------------------
## 5. Baseline random‑slope learning model (Eq. 10)
## ------------------------------------------------------------------
## c_igt = β0 + β1 c_own_{ig,t-1} + β2 c_peer_{ig,t-1} + β3 r_t + X'γ
##       + δ_g + α_v + κ_{v,t} + a_i
##       + (by‑player random slopes on round, own lag, peer lag)

mod_rs_learn_full <- lmer(
  contributing ~
    # AR(1) with peer effects + time trend
    c_own_lag + c_group_lag + round_n+
    # Covariates (as in Table 1)
    age + gender + friends + adversaries + FI + marital_status +
    network_density_fr + network_density_adv + network_size +
    b0100 + b0600 + b0200 + access_routes +
    # Random intercepts
    (1 | village_code) +               # α_v
    (1 | group) +                      # δ_g
    (1 | village_code:round_n) +       # κ_{v,t} village×round shocks
    (1 | player) +                     # a_i
    # By‑player uncorrelated random slopes (diagonal var‑cov)
    (0 + round_n    || player) +
    (0 + c_own_lag  || player) +
    (0 + c_group_lag|| player),
  data = dat,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl   = list(maxfun = 2e5)
  )
)

summary(mod_rs_learn_full)

## ----------------------------------------------------------
## 5.1. Map internal variable names -> readable labels
## ----------------------------------------------------------
## Names here match the typical lmer output (e.g., gender1, b06000).

coef_map <- c(
  # Intercept & AR terms
  "(Intercept)"  = "Intercept (baseline contribution, L)",
  "c_own_lag"    = "Own contribution (t−1)",
  "c_group_lag"  = "Group mean contribution (t−1)",
  "round_n"      = "Round (time trend)",
  
  # Individual covariates
  "age"              = "Age (years)",
  "gender1"          = "Male (vs. female)",
  "friends"          = "Number of named friends",
  "adversaries"      = "Number of named adversaries",
  "FI"               = "Food insufficiencies",
  "marital_status1"  = "Married / cohabiting",
  
  # Education / identity / religion
  "b0100"   = "Education (0–13 years)",
  "b0200"   = "Indigenous",
  "b06000"  = "Relative not religious",
  "b06001"  = "Relative Protestant",
  
  # Village / network covariates
  "network_density_fr"  = "Friendship network density (village)",
  "network_density_adv" = "Adversarial network density (village)",
  "network_size"        = "Village network size",
  "access_routes"       = "Road access (1–4)"
)

## ----------------------------------------------------------
## 5.2. Produce regression table for mod_rs_learn_full only
## ----------------------------------------------------------

modelsummary(
  mod_rs_learn_full,
  coef_map  = coef_map,           # apply new labels
  statistic = "({std.error})",    # show SEs in parentheses
  stars     = TRUE,               # significance stars
  gof_omit  = "IC|Log.Lik|RMSE",  # optional: hide some GOF rows
  title     = "Random-effects ARX model with renamed covariates",
  output    = "markdown"          # use "latex" or "html" for paper
)


## ------------------------------------------------------------------
## 6. OPTIONAL: Heterogeneity in slopes (gender, friends, FI)
## ------------------------------------------------------------------
## This corresponds to the interaction specification you sketched:
## allowing average slopes on own/peer lags to vary with gender,
## number of friends, and food insecurity (FI).

#mod_rs_het_full <- lmer(
#  contributing ~
#    # main AR terms
#    c_own_lag + c_group_lag +
#    # heterogeneity in average slopes
#    c_own_lag:gender + c_group_lag:gender +
#    c_own_lag:friends + c_group_lag:friends +
#    c_own_lag:FI      + c_group_lag:FI +
#    # time trend
#    round_n +
#    # main covariates (same set as above)
#    age + gender + friends + adversaries + FI + marital_status +
#    network_density_fr + network_density_adv + network_size +
#    b0100 + b0600 + b0200 + access_routes +
#    # random effects structure identical to baseline model
#    (1 | village_code) +
#    (1 | group) +
#    (1 | village_code:round_f) +
#    (1 | player) +
#    (0 + round_n    || player) +
#    (0 + c_own_lag  || player) +
#    (0 + c_group_lag|| player),
#  data = dat,
#  control = lmerControl(
#    optimizer = "bobyqa",
#    optCtrl   = list(maxfun = 2e5)
#  )
#)

#summary(mod_rs_het_full)

## ------------------------------------------------------------------
## 7. OPTIONAL: Centered version for interpretability
## ------------------------------------------------------------------
## If you want the intercept and random slopes to be interpreted at
## “average” round and average lagged contributions, uncomment below.

# mu_round <- mean(dat$round_n,     na.rm = TRUE)  # ~5.5 for 1..10
# mu_own   <- mean(dat$c_own_lag,   na.rm = TRUE)
# mu_peer  <- mean(dat$c_group_lag, na.rm = TRUE)
#
# dat_c <- dat %>%
#   mutate(
#     round_c       = round_n    - mu_round,
#     c_own_lag_c   = c_own_lag  - mu_own,
#     c_group_lag_c = c_group_lag- mu_peer
#   )
#
# mod_rs_learn_centered <- lmer(
#   contributing ~
#     c_own_lag_c + c_group_lag_c + round_c +
#     age + gender + friends + adversaries + FI + marital_status +
#     network_density_fr + network_density_adv + network_size +
#     b0100 + b0600 + b0200 + access_routes +
#     (1 | village_code) +
#     (1 | group) +
#     (1 | village_code:round_f) +
#     (1 | player) +
#     (0 + round_c       || player) +
#     (0 + c_own_lag_c   || player) +
#     (0 + c_group_lag_c || player),
#   data    = dat_c,
#   control = lmerControl(
#     optimizer = "bobyqa",
#     optCtrl   = list(maxfun = 2e5)
#   )
# )
#
# summary(mod_rs_learn_centered)
