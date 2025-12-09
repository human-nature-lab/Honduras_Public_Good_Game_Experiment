###############################################################################
## 0) PACKAGES ---------------------------------------------------------------
###############################################################################
# Define the set of packages required for the full replication:
# - dplyr, tidyr, tibble, stringr, purrr: data wrangling
# - ggplot2: plotting
# - lme4, broom, broom.mixed, splines: mixed models & spline terms
# - ggeffects, marginaleffects, emmeans: marginal effects and predictions
# - yardstick, pROC, DescTools, rsample: model performance tools
# - rlang: programming helpers
# - geepack, sandwich, lmtest: GEE and cluster-robust inference
# - survival, survminer: Cox models and KM plots
# - depmixS4: hidden Markov models
pkgs <- c(
  "dplyr","tidyr","tibble","ggplot2","stringr","purrr",
  "lme4","broom","broom.mixed","splines","ggeffects","marginaleffects","emmeans",
  "yardstick","pROC","DescTools","rsample","rlang",
  "geepack","survival","survminer","depmixS4"
)

# Install any missing packages (with dependencies), then load all of them
new <- setdiff(pkgs, rownames(installed.packages()))
if(length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Fix random seed for reproducibility of any stochastic routines (e.g. HMM fitting)
set.seed(123)

# Main experimental panel at the individual × round level
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)


###############################################################################
## 1) DATA PREP --------------------------------------------------------------
###############################################################################
# Threshold c*: contribution level that defines the "High" state (H) at baseline.
# Here c* is set to the mean round-1 contribution over all players.
c_star <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::summarise(m = mean(contributing, na.rm = TRUE)) %>%
  dplyr::pull(m)

# Baseline covariates: one row per player (round 1 only).
# - Coerce IDs and factors to the intended types.
# - Construct z-scores for continuous covariates to ease interpretation.
baseline <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::distinct(player, .keep_all = TRUE) %>%   # one baseline row per player
  dplyr::transmute(
    player        = as.character(player),
    village_code  = as.character(village_code),
    gender        = factor(gender),
    marital_status= factor(marital_status),
    b0600         = factor(b0600),     # religion (re-leveled below)
    access_routes = factor(access_routes),
    age           = suppressWarnings(as.numeric(age)),
    friends       = suppressWarnings(as.numeric(friends)),
    adversaries   = suppressWarnings(as.numeric(adversaries)),
    network_density_fr  = suppressWarnings(as.numeric(network_density_fr)),
    network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
    network_size        = suppressWarnings(as.numeric(network_size)),
    b0100              = suppressWarnings(as.numeric(b0100)),  # education (years)
    b0200              = as.numeric(b0200),  # indigenous indicator (as numeric here)
    FI                 = suppressWarnings(as.numeric(FI))      # financial autonomy index
  ) %>%
  # Standardize all continuous covariates used later in the models.
  dplyr::mutate(dplyr::across(
    c(age,friends,adversaries,network_density_fr,network_density_adv,network_size,b0100,FI),
    ~ as.numeric(scale(.x)),
    .names = "{.col}_z"
  )) %>%
  # Set religion reference category to "2" (Catholic) if that level exists.
  # All other levels are interpreted relative to Catholic.
  dplyr::mutate(b0600 = stats::relevel(b0600, ref = "2"))


# Construct the panel used in the transition analyses:
# - One row per player × round with the binary state s_t (H vs L)
# - Group-level leave-one-out (LOO) peer mean contributions, lagged by one round
# - Spell length (duration already spent in current regime at t-1)
panel <- pgg_data %>%
  dplyr::filter(is.finite(contributing)) %>%
  dplyr::mutate(
    player       = as.character(player),
    village_code = as.character(village_code),
    group        = as.character(group),
    # State indicator: s_t = 1 if contribution ≥ c*, 0 otherwise
    s            = as.integer(contributing >= c_star),
    b0600        = as.factor(b0600),
    gender       = as.factor(gender)
  ) %>%
  # Within each group-round, compute leave-one-out peer mean contribution at time t
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    G = dplyr::n(),
    peer_mean_c_t = ifelse(
      G > 1,
      (sum(contributing, na.rm = TRUE) - contributing) / (G - 1),
      NA_real_
    )
  ) %>% dplyr::ungroup() %>%
  # Lag the peer mean by group to obtain peer_mean_c_{t-1}
  dplyr::arrange(group, round_n) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(peer_mean_c_tm1 = dplyr::lag(peer_mean_c_t)) %>%
  dplyr::ungroup() %>%
  # Within-player lags and duration (spell) indices
  dplyr::arrange(player, round_n) %>%
  dplyr::group_by(player) %>%
  dplyr::mutate(
    # Previous-period state s_{t-1}
    s_lag   = dplyr::lag(s),
    # Time index (round number)
    r_t     = as.integer(round_n),
    # Identify “runs” of consecutive same-state periods
    run_id  = cumsum(dplyr::if_else(
      is.na(dplyr::lag(s)), TRUE, s != dplyr::lag(s)
    )),
    # Position within the run (1, 2, 3, …)
    spell_index = dplyr::row_number(),
    # Spell length at t-1 (duration already spent in current regime)
    spell_prev  = dplyr::lag(spell_index)
  ) %>%
  dplyr::ungroup() %>%
  # Merge baseline covariates (one row per player) into the full panel
  dplyr::left_join(
    baseline %>% dplyr::select(
      player, village_code,
      gender, marital_status, b0600, b0200, access_routes,
      ends_with("_z")
    ),
    by = c("player", "village_code")
  ) %>%
  # Scale peer mean by the per-round endowment (assumed to be 12 units)
  dplyr::mutate(peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12)

# Construct “at risk” data sets for the two transition types:
# - drop_df: players currently High (s_{t-1}=1), at risk of dropping H → L
# - rise_df: players currently Low (s_{t-1}=0), at risk of rising L → H
drop_df <- panel %>%      # At risk of H->L transitions
  dplyr::filter(!is.na(s_lag), s_lag == 1) %>%
  dplyr::mutate(y_drop = as.integer(s == 0)) %>%  # event indicator: drop if s_t=0
  dplyr::filter(
    is.finite(peer_mean_c_tm1_scaled),
    is.finite(spell_prev)
  )

rise_df <- panel %>%      # At risk of L->H transitions
  dplyr::filter(!is.na(s_lag), s_lag == 0) %>%
  dplyr::mutate(y_rise = as.integer(s == 1)) %>%  # event indicator: rise if s_t=1
  dplyr::filter(
    is.finite(peer_mean_c_tm1_scaled),
    is.finite(spell_prev)
  )

# Ensure religion factor in the at-risk sets uses Catholic ("2") as reference
drop_df$b0600.x <- relevel(drop_df$b0600.x, ref = "2")   # Catholic
rise_df$b0600.x <- relevel(rise_df$b0600.x, ref = "2")


###############################################################################
## 2) TRANSITION COUNTS & AT-RISK SETS ---------------------------------------
###############################################################################
# Compute empirical 2×2 transition counts (from s_{t-1} to s_t) by round.
# This is purely descriptive: it shows how many L→L, L→H, H→L, H→H transitions occur.
trans_by_round <- panel %>%
  dplyr::filter(!is.na(s_lag)) %>%
  dplyr::count(round_n, s_lag, s, name = "n") %>%
  tidyr::complete(round_n, s_lag = 0:1, s = 0:1, fill = list(n = 0)) %>%
  dplyr::arrange(round_n, s_lag, s)

# Convert transition counts for each round into an explicit 2×2 matrix
# (rows: from L/H; columns: to L/H).
trans_list <- trans_by_round %>%
  dplyr::group_by(round_n) %>%
  dplyr::group_split() %>%
  purrr::map(~{
    M <- matrix(0, nrow = 2, ncol = 2,
                dimnames = list(from = c("L","H"), to = c("L","H")))
    for(i in seq_len(nrow(.x))) M[.x$s_lag[i] + 1, .x$s[i] + 1] <- .x$n[i]
    M
  })
names(trans_list) <- sort(unique(trans_by_round$round_n))

# Transition counts aggregated over all rounds
trans_total <- panel %>%
  dplyr::filter(!is.na(s_lag)) %>%
  dplyr::count(s_lag, s, name = "n") %>%
  tidyr::complete(s_lag = 0:1, s = 0:1, fill = list(n = 0)) %>%
  dplyr::arrange(s_lag, s)

cat("\n--- Transition counts by round: list of 2x2 matrices (from rows, to cols) ---\n")
print(trans_list[[1]])   # Example: transitions between round 1 and 2
cat("\n--- Total transitions over all rounds (from s_{t-1} to s_t) ---\n")
print(trans_total)


###############################################################################
## 3) TRANSITION-SPECIFIC GLMMs (Drop & Rise) + DURATION SPLINES -------------
###############################################################################
# Define the common right-hand side (RHS) of the transition models.
# It includes:
# - peer_mean_c_tm1_scaled: peer signal at t-1
# - r_t: round index (time trend)
# - ns(spell_prev, df=3): flexible effect of prior duration in current regime
# - baseline covariates: gender, marital status, religion, ethnicity, access, and
#   standardized individual/ network characteristics.
rhs <- paste(c(
  "peer_mean_c_tm1_scaled", "r_t", "ns(spell_prev, df=3)",  # dynamics & duration
  "gender.x","marital_status.x","b0600.x","b0200.x","access_routes.x",
  "age_z","friends_z","adversaries_z","network_density_fr_z",
  "network_density_adv_z","network_size_z","b0100_z","FI_z"
), collapse = " + ")

# Drop model: conditional hazard of leaving the High state (H→L),
# estimated only on observations with s_{t-1} = 1.
form_drop <- as.formula(paste0("y_drop ~ ", rhs, " + (1|village_code)"))
mod_drop  <- glmer(
  form_drop,
  data   = drop_df,
  family = binomial,
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)

# Rise model: conditional hazard of leaving the Low state (L→H),
# estimated only on observations with s_{t-1} = 0.
form_rise <- as.formula(paste0("y_rise ~ ", rhs, " + (1|village_code)"))
mod_rise  <- glmer(
  form_rise,
  data   = rise_df,
  family = binomial,
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)

cat("\n--- DROP model (HL vs HH) summary ---\n");  print(summary(mod_drop))
cat("\n--- RISE model (LH vs LL) summary ---\n");  print(summary(mod_rise))

# Helper function: extract odds ratios (fixed effects), random-effect variances,
# and information criteria from a fitted GLMM.
or_table <- function(fit){
  fx <- broom.mixed::tidy(
    fit,
    effects      = "fixed",
    conf.int     = TRUE,
    exponentiate = TRUE
  )
  re <- broom.mixed::tidy(fit, effects = "ran_pars")
  list(
    OR = fx %>% dplyr::select(term, estimate, conf.low, conf.high, p.value),
    RE = re,
    IC = c(AIC = AIC(fit), BIC = BIC(fit))
  )
}
drop_out <- or_table(mod_drop)
rise_out <- or_table(mod_rise)

cat("\n--- DROP ORs ---\n"); print(drop_out$OR)
cat("\n--- RISE ORs ---\n"); print(rise_out$OR)
cat("\n--- Random-effects (variances) ---\n"); print(drop_out$RE); print(rise_out$RE)
cat("\n--- AIC/BIC ---\n"); print(drop_out$IC); print(rise_out$IC)


###############################################################################
## 4) FIGURE 1: MARGINAL EFFECTS ---------------------------------------------
###############################################################################
###############################################################################
##  FIGURE 1: MARGINAL EFFECTS (Drop vs friends × gender; Rise vs religion)  ##
###############################################################################
# This section constructs marginal effect plots corresponding to Figure 1:
# - Panel A: Drop hazard as a function of friends_z, by gender
# - Panel B: Rise hazard by religion

# Ensure plotting packages used here are available
if (!requireNamespace("ggeffects", quietly = TRUE)) install.packages("ggeffects")
if (!requireNamespace("ggplot2",  quietly = TRUE)) install.packages("ggplot2")
library(ggeffects)
library(ggplot2)

# Sensible defaults for conditioning in the marginal effects:
# For each model, we hold r_t, spell_prev, and the peer signal at representative
# values, while other covariates are set by ggeffects (typically at means or
# baseline factor levels).

# Conditioning values for the DROP model (H->L)
cond_drop <- list(
  r_t        = median(drop_df$r_t, na.rm = TRUE),
  spell_prev = median(drop_df$spell_prev, na.rm = TRUE),
  peer_mean_c_tm1_scaled = mean(drop_df$peer_mean_c_tm1_scaled, na.rm = TRUE)
  # Remaining covariates are set to reference or mean values by ggeffects
)

# Conditioning values for the RISE model (L->H)
cond_rise <- list(
  r_t        = median(rise_df$r_t, na.rm = TRUE),
  spell_prev = median(rise_df$spell_prev, na.rm = TRUE),
  peer_mean_c_tm1_scaled = mean(rise_df$peer_mean_c_tm1_scaled, na.rm = TRUE)
)

# ------------------------------ PLOT A: DROP ----------------------------------
# Predicted Pr(H->L) over a range of friends_z, for each gender.
# friends_z is evaluated on a symmetric grid [-2, 2] in z-score units.
me_drop_friends <- ggpredict(
  model     = mod_drop,
  terms     = c("friends_z [-2:2 by=0.2]", "gender.x"),  # slope × grouping factor
  type      = "fixed",                                   # fixed effects only
  condition = cond_drop
)

# If gender.x is coded as 0/1 (numeric levels), relabel them for readability.
if ("group" %in% names(me_drop_friends)) {
  me_drop_friends$group <- factor(
    me_drop_friends$group,
    levels = levels(drop_df$gender.x),
    labels = c("Female","Male")[seq_len(nlevels(drop_df$gender.x))]
  )
}

p_drop <- ggplot(me_drop_friends,
                 aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = .15, colour = NA) +
  geom_line(size = 1.1) +
  labs(
    x      = "Friends (z-scored)",
    y      = "Pr(drop H→L)",
    colour = "Gender",
    fill   = "Gender",
    title  = "Drop hazard vs. friends",
    subtitle = "Fixed-effects predictions; 95% CI"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------ PLOT B: RISE ----------------------------------
# Predicted Pr(L->H) as a function of religion (b0600.x), evaluated at
# the conditioning values defined above. ggpredict creates one prediction
# per factor level.
me_rise_relig <- ggpredict(
  model     = mod_rise,
  terms     = "b0600.x",     # religion as a factor
  type      = "fixed",
  condition = cond_rise
)

# Relabel simple numeric religion codes (0/1/2) if coding is:
# 0 = No religion, 1 = Protestant, 2 = Catholic (baseline).
pretty_labs <- c("0" = "No religion", "1" = "Protestant", "2" = "Catholic")
if (all(me_rise_relig$x %in% names(pretty_labs))) {
  me_rise_relig$x <- factor(
    pretty_labs[as.character(me_rise_relig$x)],
    levels = pretty_labs[intersect(c("0","1","2"), names(pretty_labs))]
  )
}

p_rise <- ggplot(me_rise_relig, aes(x = x, y = predicted)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = .08) +
  labs(
    x      = "Religion (baseline: Catholic)",
    y      = "Pr(rise L→H)",
    title  = "Rise hazard by religion",
    subtitle = "Fixed-effects predictions; 95% CI"
  ) +
  theme_minimal(base_size = 12)

# Display the two separate panels; they can be combined later if desired
p_drop
p_rise

# Example two-panel layout (commented out):
# if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
# library(patchwork)
# p_drop + p_rise + plot_annotation(title = "Figure 1. Marginal effects")


######## 5. Figure 2 (duration hazard curves)

# ==============================================
# Figure 2: Duration hazard curves (Drop & Rise)
# ==============================================
# This section plots the discrete hazard (transition probability) as a function
# of spell length (time already spent in the current regime), holding other
# covariates at representative values.

library(dplyr)
library(ggplot2)
library(ggeffects)   # ggpredict handles ns() spline terms
library(forcats)

# 0) Small helpers -----------------------------------------------------------
# majority_level(): return the most frequent level of a factor (excluding NA).
majority_level <- function(f) {
  tab <- table(f, useNA = "no")
  if (length(tab) == 0) return(levels(f)[1])
  names(tab)[which.max(tab)]
}

# build_condition(): construct a list of conditioning values for ggpredict,
# setting numeric variables to their sample means and factors to their majority
# level in the training frame.
build_condition <- function(df, num_vars, fac_vars) {
  cond <- list()
  for (v in num_vars) if (v %in% names(df)) cond[[v]] <- mean(df[[v]], na.rm = TRUE)
  for (v in fac_vars) if (v %in% names(df)) cond[[v]] <- majority_level(df[[v]])
  cond
}

# 1) DROP duration curve  ----------------------------------------------------
# Variables used in mod_drop; we pass them to build_condition() to construct
# a representative covariate profile for duration curves.
num_vars_drop <- c("peer_mean_c_tm1_scaled","r_t",
                   "age_z","friends_z","adversaries_z",
                   "network_density_fr_z","network_density_adv_z",
                   "network_size_z","b0100_z","FI_z")
fac_vars_drop <- c("gender.x","marital_status.x","b0600.x","b0200.x","access_routes.x")

cond_drop <- build_condition(drop_df, num_vars_drop, fac_vars_drop)

# ggeffects can expand ns(spell_prev, 3) if we ask for "[all]" to loop over the
# distinct observed values of spell_prev.
me_drop_dur <- ggpredict(
  model     = mod_drop,
  terms     = c("spell_prev [all]"),  # discrete durations observed in the data
  condition = cond_drop,
  type      = "fixed"                 # fixed effects only
)

p_drop_dur <- plot(me_drop_dur) +
  labs(
    x     = "Time already spent High (t−1)",
    y     = "Pr(drop H→L)",
    title = "Discrete hazard of drop by spell length"
  ) +
  theme_minimal()

# 2) RISE duration curve  ----------------------------------------------------
num_vars_rise <- c("peer_mean_c_tm1_scaled","r_t",
                   "age_z","friends_z","adversaries_z",
                   "network_density_fr_z","network_density_adv_z",
                   "network_size_z","b0100_z","FI_z")
fac_vars_rise <- c("gender.x","marital_status.x","b0600.x","b0200.x","access_routes.x")

cond_rise <- build_condition(rise_df, num_vars_rise, fac_vars_rise)

me_rise_dur <- ggpredict(
  model     = mod_rise,
  terms     = c("spell_prev [all]"),
  condition = cond_rise,
  type      = "fixed"
)

p_rise_dur <- plot(me_rise_dur) +
  labs(
    x     = "Time already spent Low (t−1)",
    y     = "Pr(rise L→H)",
    title = "Discrete hazard of rise by spell length"
  ) +
  theme_minimal()

# Show the two duration plots (Drop and Rise)
p_drop_dur
p_rise_dur


# ----- Fixed-effects predictions aligned to a new data frame -----------------
# pred_fixed_glmer_aligned():
# Given a fitted glmer object and a new data frame, compute fixed-effects
# predictions (no random effects) with careful alignment of:
# - factor levels
# - design matrix columns
# This is important for robust out-of-sample prediction and for calibration.
pred_fixed_glmer_aligned <- function(mod, df, type = c("response","link")) {
  type <- match.arg(type)
  
  # Extract the fixed-effects formula (remove random effects from the full formula)
  f_fix <- lme4::nobars(formula(mod))
  
  # Training model frame and its contrast specification
  mf_train <- model.frame(mod)
  contr    <- attr(mf_train, "contrasts")
  
  # Ensure that factor variables in df share the same levels as in training data
  fac_train <- names(which(vapply(mf_train, is.factor, logical(1))))
  for (v in fac_train) {
    if (v %in% names(df)) {
      df[[v]] <- factor(df[[v]], levels = levels(mf_train[[v]]))
    }
  }
  
  # Construct a model frame allowing NA rows, then the corresponding design matrix
  mf_new <- model.frame(f_fix, df, na.action = na.pass)
  Xnew   <- model.matrix(f_fix, data = mf_new, contrasts.arg = contr)
  
  # Align columns of Xnew with the fixed-effect coefficient vector beta
  beta <- fixef(mod)
  
  # If the model has columns that are missing in Xnew, append zero-columns for them
  missing_cols <- setdiff(names(beta), colnames(Xnew))
  if (length(missing_cols)) {
    Xnew <- cbind(
      Xnew,
      `colnames<-`(
        matrix(0, nrow = nrow(Xnew), ncol = length(missing_cols)),
        missing_cols
      )
    )
  }
  Xnew <- Xnew[, names(beta), drop = FALSE]
  
  # Predict only for rows with complete covariates in Xnew
  ok  <- stats::complete.cases(Xnew)
  eta <- rep(NA_real_, nrow(df))
  if (any(ok)) eta[ok] <- as.numeric(Xnew[ok, , drop = FALSE] %*% beta)
  
  if (type == "link") return(eta)
  plogis(eta)  # convert log-odds to probabilities
}

# ---- Utility: coerce types in new data to match a reference vector ----------
# This is used below to ensure factor/ numeric types are consistent before
# calling pred_fixed_glmer_aligned().
coerce_like <- function(x_new, x_ref){
  if (is.factor(x_ref)) factor(x_new, levels = levels(x_ref))
  else if (is.numeric(x_ref)) suppressWarnings(as.numeric(x_new))
  else if (is.integer(x_ref)) suppressWarnings(as.integer(x_new))
  else x_new
}

# Ensure factor columns in drop_df conform to the factor structure the model
# was actually fit on (from model.frame(mod_drop)).
mf_train_drop <- model.frame(mod_drop)
for (v in c("gender.x","marital_status.x","b0600.x","b0200.x","access_routes.x")) {
  if (v %in% names(drop_df) && v %in% names(mf_train_drop)) {
    drop_df[[v]] <- coerce_like(drop_df[[v]], mf_train_drop[[v]])
  }
}

# Same for rise_df and mod_rise
mf_train_rise <- model.frame(mod_rise)
for (v in c("gender.x","marital_status.x","b0600.x","b0200.x","access_routes.x")) {
  if (v %in% names(rise_df) && v %in% names(mf_train_rise)) {
    rise_df[[v]] <- coerce_like(rise_df[[v]], mf_train_rise[[v]])
  }
}

# Fixed-effects predicted probabilities for all rows in the at-risk sets
drop_df$.phat <- pred_fixed_glmer_aligned(mod_drop, drop_df, type = "response")
rise_df$.phat <- pred_fixed_glmer_aligned(mod_rise, rise_df, type  = "response")



###############################################################################
## 6) MODEL FIT & CALIBRATION (AUC, Brier, Calibration curves) ---------------
###############################################################################
# calc_metrics_safe():
# - computes fixed-effects predictions,
# - keeps rows where both y and p are defined,
# - returns AUC, Brier score, and calibration intercept & slope.
calc_metrics_safe <- function(df, fit, yvar){
  # 1) fixed-effects predictions aligned to the model (no random effects)
  p <- as.numeric(pred_fixed_glmer_aligned(fit, df, type = "response"))
  
  # 2) keep only rows where both y and p are defined
  y  <- df[[yvar]]
  ok <- is.finite(p) & !is.na(y)
  y  <- y[ok]; p <- p[ok]
  
  # 3) performance metrics
  auc   <- as.numeric(pROC::auc(y, p))
  brier <- mean((y - p)^2)
  
  # 4) calibration on the logit scale (intercept and slope from y ~ logit(p))
  eps   <- 1e-8
  lp    <- qlogis(pmin(pmax(p, eps), 1 - eps))
  cal_fit <- glm(y ~ lp, family = binomial())
  cal_int <- unname(coef(cal_fit)[1])
  cal_slo <- unname(coef(cal_fit)[2])
  
  list(
    auc          = auc,
    brier        = brier,
    cal_intercept= cal_int,
    cal_slope    = cal_slo,
    pred         = p,
    y            = y
  )
}

# Evaluate fit for Drop and Rise models
drop_fit <- calc_metrics_safe(drop_df, mod_drop, "y_drop")
rise_fit <- calc_metrics_safe(rise_df,  mod_rise,  "y_rise")

cat("\n--- DROP fit ---\n",
    "AUC =", round(drop_fit$auc, 3),
    "| Brier =", round(drop_fit$brier, 4),
    "| Cal. intercept =", round(drop_fit$cal_intercept, 3),
    "| Cal. slope =", round(drop_fit$cal_slope, 3), "\n")

cat("\n--- RISE fit ---\n",
    "AUC =", round(rise_fit$auc, 3),
    "| Brier =", round(rise_fit$brier, 4),
    "| Cal. intercept =", round(rise_fit$cal_intercept, 3),
    "| Cal. slope =", round(rise_fit$cal_slope, 3), "\n")

# Decile-based calibration plots:
# - Sort observations by predicted probability,
# - Split into deciles,
# - Compare mean predicted vs. mean observed outcome per decile.
plot_cal <- function(y, p, title, ylab){
  tibble(y, p) |>
    dplyr::mutate(dec = dplyr::ntile(p, 10)) |>
    dplyr::group_by(dec) |>
    dplyr::summarise(pred = mean(p), obs = mean(y), .groups = "drop") |>
    ggplot(aes(pred, obs)) +
    geom_point(size = 2) + geom_line() +
    geom_abline(linetype = 2, colour = "grey50") +
    labs(x = "Predicted (bin mean)", y = ylab, title = title) +
    theme_minimal()
}

print(plot_cal(drop_fit$y, drop_fit$pred,
               "Calibration: H→L (Drop)", "Observed drop rate"))
print(plot_cal(rise_fit$y, rise_fit$pred,
               "Calibration: L→H (Rise)", "Observed rise rate"))



###############################################################################
## 7) TRANSITION SYSTEM CHECK: INTEGRATED MARKOV LOGIT -----------------------
###############################################################################
# This section fits a single integrated mixed-effects logit in which:
# - s_t is regressed on s_{t-1}, peer signals, duration, and covariates.
# - All main covariates are interacted with s_{t-1}, so the model nests both
#   Drop (H→L) and Rise (L→H) hazards in one system.
# - Random intercepts are included at village and player levels.

# Fit / reuse integrated Markov logit if not already present in the workspace
if (!exists("mod_markov")) {
  markov_form <- as.formula(paste0(
    "s ~ s_lag + ",
    "peer_mean_c_tm1_scaled*s_lag + r_t*s_lag + ns(spell_prev,3)*s_lag + ",
    "gender.x*s_lag + marital_status.x*s_lag + b0600.x*s_lag + b0200.x*s_lag + access_routes.x*s_lag + ",
    "age_z*s_lag + friends_z*s_lag + adversaries_z*s_lag + ",
    "network_density_fr_z*s_lag + network_density_adv_z*s_lag + network_size_z*s_lag + ",
    "b0100_z*s_lag + FI_z*s_lag + (1|village_code) + (1|player)"
  ))
  markov_df <- panel %>% dplyr::filter(
    !is.na(s_lag),
    is.finite(peer_mean_c_tm1_scaled),
    is.finite(spell_prev)
  )
  mod_markov <- lme4::glmer(
    markov_form,
    data   = markov_df,
    family = binomial,
    nAGQ   = 0,
    control = lme4::glmerControl(
      optimizer   = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl     = list(maxfun = 2e4)
    )
  )
  
  summary(mod_markov)
}

# Model frame used in the integrated Markov logit
mf_markov <- model.frame(mod_markov)

# In-sample predicted probabilities (fixed effects only)
mf_markov$.phat <- as.numeric(
  predict(mod_markov, type = "response", re.form = NA)
)

# Row-mean transition probabilities by previous state s_{t-1}
P_hat <- mf_markov |>
  dplyr::group_by(s_lag) |>
  dplyr::summarise(
    P_q1 = mean(.phat),   # Pr(s_t=1 | s_{t-1})
    P_q0 = 1 - P_q1,
    .groups = "drop"
  ) |>
  dplyr::arrange(s_lag)

# Construct a 2×2 transition matrix P_mat from row means (L/H as states)
P_mat <- matrix(NA_real_, 2, 2, dimnames = list(from = c("L","H"), to = c("L","H")))
P_mat[1,2] <- P_hat$P_q1[P_hat$s_lag == 0]; P_mat[1,1] <- 1 - P_mat[1,2]
P_mat[2,2] <- P_hat$P_q1[P_hat$s_lag == 1]; P_mat[2,1] <- 1 - P_mat[2,2]
cat("\n--- Mean transition matrix (Markov logit) ---\n"); print(round(P_mat, 3))

# Forward iteration of the Markov chain to round 10:
# Start from round-1 share of High types; iteratively apply P_mat for 9 transitions.
pi1 <- panel |>
  dplyr::filter(round_n == 1) |>
  dplyr::summarise(pH = mean(s == 1)) |>
  dplyr::pull()
pi <- c(L = 1 - pi1, H = pi1)
for (k in 1:9) pi <- as.numeric(pi %*% P_mat)
cat("Implied round-10 High share:", round(pi["H"], 3), "\n")

# Compare with the observed share of High types in round 10
obs_H10 <- panel |>
  dplyr::filter(round_n == 10) |>
  dplyr::summarise(pH = mean(s == 1)) |>
  dplyr::pull()
cat("Observed round-10 High share:", round(obs_H10, 3), "\n")

# Ensure factor types in markov_df match the model’s training frame.
# This is useful if one recomputes predictions from markov_df using
# pred_fixed_glmer_aligned().
mf_train_mk <- model.frame(mod_markov)
fac_vars_mk <- names(which(vapply(mf_train_mk, is.factor, logical(1))))
for (v in intersect(fac_vars_mk, names(markov_df))) {
  markov_df[[v]] <- factor(markov_df[[v]], levels = levels(mf_train_mk[[v]]))
}

# Aligned fixed-effects predictions for the integrated Markov model
markov_df$.phat <- pred_fixed_glmer_aligned(mod_markov, markov_df, type = "response")

# Recompute row-mean transition probabilities based on markov_df
P_hat <- markov_df |>
  dplyr::group_by(s_lag) |>
  dplyr::summarise(
    P_q1 = mean(.phat, na.rm = TRUE),
    P_q0 = 1 - P_q1,
    .groups = "drop"
  ) |>
  dplyr::arrange(s_lag)

# Build P_mat and iterate as above (structure already set up in P_mat)

# Heatmap of the transition matrix with probabilities labeled in each cell
P_ci <- tibble(
  from = rep(c("L","H"), each = 2),
  to   = rep(c("L","H"), times = 2),
  p    = c(P_mat[1,1], P_mat[1,2], P_mat[2,1], P_mat[2,2])
)
ggplot(P_ci, aes(to, from, fill = p)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", p)),
            colour = "white", fontface = "bold") +
  scale_fill_gradient(low = "#9ecae1", high = "#08519c") +
  labs(
    title = "Mean transition matrix (Markov logit)",
    x     = "to",
    y     = "from"
  ) +
  theme_minimal()


###############################################################################
## 8) TARGETED HETEROGENEITY (INTERACTIONS) ----------------------------------
###############################################################################
# In this section, we augment the Drop and Rise models with selected interaction
# terms to capture specific heterogeneous effects.

# Drop model: add interaction terms
# - friends_z × gender.x  (heterogeneous peer effect by gender)
# - network_density_adv_z × b0100_z (network of adversaries × education)
form_drop_int <- update(
  form_drop,
  . ~ . + friends_z:gender.x + network_density_adv_z:b0100_z
)
mod_drop_int <- glmer(
  form_drop_int,
  data   = drop_df,
  family = binomial,
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)
summary(mod_drop_int)

# Rise model: heterogeneous effect of financial autonomy (FI_z) by "No religion".
# We create an indicator no_relig = 1 if religion code corresponds to "No religion".
# (Here we assume "0" encodes No religion; adjust if coding differs.)
rise_df <- rise_df %>%
  dplyr::mutate(no_relig = as.integer(as.character(b0600.x) == "0"))
form_rise_int <- update(
  form_rise,
  . ~ . + no_relig:FI_z
)
mod_rise_int <- glmer(
  form_rise_int,
  data   = rise_df,
  family = binomial,
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)
summary(mod_rise_int)

# Figure: simple effects for Drop by gender (friends_z × gender.x)
plot_drop_int <- ggeffects::ggpredict(
  mod_drop_int,
  terms = c("friends_z[-2:2 by=0.25]","gender.x")
)
ggplot(plot_drop_int, aes(x, predicted, colour = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = .12, colour = NA) +
  labs(
    x      = "Friends (z)",
    y      = "P(H→L | at risk)",
    colour = "Gender",
    fill   = "Gender",
    title  = "Drop (H→L): friends × gender"
  ) +
  theme_minimal()

# Figure: simple effects for Rise by FI_z, comparing No religion vs others
plot_rise_int <- ggeffects::ggpredict(
  mod_rise_int,
  terms = c("FI_z[-2:2 by=0.25]","no_relig")
)
ggplot(plot_rise_int, aes(x, predicted, colour = group)) +
  geom_line(size = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = .12, colour = NA) +
  scale_color_discrete(labels = c("Catholic/other","No religion")) +
  scale_fill_discrete(labels = c("Catholic/other","No religion")) +
  labs(
    x      = "Financial autonomy (z)",
    y      = "P(L→H | at risk)",
    colour = "Religion",
    fill   = "Religion",
    title  = "Rise (L→H): no-religion × financial autonomy"
  ) +
  theme_minimal()


###############################################################################
## 9) ROBUSTNESS: THRESHOLD SENSITIVITY --------------------------------------
###############################################################################
# Robustness check: re-estimate Drop and Rise models for alternative thresholds
# c* used to define High/Low states. Here we consider:
# - mean (baseline),
# - median,
# - mean ± 0.5 SD (evaluated on round-1 contributions).

# Thresholds computed directly from round-1 contributions
th_mean   <- c_star
th_median <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::summarise(m = median(contributing, na.rm = TRUE)) %>%
  dplyr::pull(m)
th_sd     <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::summarise(s = sd(contributing, na.rm = TRUE)) %>%
  dplyr::pull(s)
th_grid   <- c(
  mean      = th_mean,
  median    = th_median,
  mean_m05  = th_mean - 0.5 * th_sd,
  mean_p05  = th_mean + 0.5 * th_sd
)

# First version of fit_at_threshold(): fits Drop/Rise GLMMs at a numeric threshold
# and returns selected odds ratios and CIs. It is later superseded by a more
# general version that also handles named methods (mean, median, etc.).
fit_at_threshold <- function(th){
  pan <- pgg_data %>%
    dplyr::filter(is.finite(contributing)) %>%
    dplyr::mutate(
      player       = as.character(player),
      village_code = as.character(village_code),
      group        = as.character(group),
      s            = as.integer(contributing >= th)
    ) %>%
    # LOO peer mean at time t
    dplyr::group_by(group, round_n) %>%
    dplyr::mutate(
      G = dplyr::n(),
      peer_mean_c_t = ifelse(
        G > 1,
        (sum(contributing, na.rm = TRUE) - contributing) / (G - 1),
        NA_real_
      )
    ) %>%
    dplyr::ungroup() %>%
    # Lag peer mean by group
    dplyr::arrange(group, round_n) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(peer_mean_c_tm1 = dplyr::lag(peer_mean_c_t)) %>%
    dplyr::ungroup() %>%
    # Compute own lag and spell_prev
    dplyr::arrange(player, round_n) %>%
    dplyr::group_by(player) %>%
    dplyr::mutate(
      s_lag       = dplyr::lag(s),
      run_id      = cumsum(dplyr::if_else(
        is.na(dplyr::lag(s)), TRUE, s != dplyr::lag(s)
      )),
      spell_index = dplyr::row_number(),
      spell_prev  = dplyr::lag(spell_index)
    ) %>%
    dplyr::ungroup() %>%
    # Merge baseline covariates and scale peer mean by endowment
    dplyr::left_join(
      baseline %>% dplyr::select(
        player, village_code,
        gender, marital_status, b0600, b0200, access_routes, ends_with("_z")
      ),
      by = c("player","village_code")
    ) %>%
    dplyr::mutate(peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12)
  
  # At-risk sets under this alternative threshold
  drop_df2 <- pan %>%
    dplyr::filter(
      !is.na(s_lag), s_lag == 1,
      is.finite(peer_mean_c_tm1_scaled), is.finite(spell_prev)
    ) %>%
    dplyr::mutate(y_drop = as.integer(s == 0))
  rise_df2 <- pan %>%
    dplyr::filter(
      !is.na(s_lag), s_lag == 0,
      is.finite(peer_mean_c_tm1_scaled), is.finite(spell_prev)
    ) %>%
    dplyr::mutate(y_rise = as.integer(s == 1))
  
  # Fit Drop model with the same form_drop formula (defined earlier)
  f_drop <- glmer(
    form_drop,
    data   = drop_df2,
    family = binomial,
    nAGQ   = 0,
    control = glmerControl(
      optimizer   = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl     = list(maxfun = 2e4)
    )
  )
  summary(f_drop)
  
  # Fit Rise model with the same form_rise formula
  f_rise <- glmer(
    form_rise,
    data   = rise_df2,
    family = binomial,
    nAGQ   = 0,
    control = glmerControl(
      optimizer   = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl     = list(maxfun = 2e4)
    )
  )
  
  summary(f_rise)
  
  # Extract ORs (and CIs) for selected terms from the two models
  or_drop <- broom.mixed::tidy(
    f_drop, effects = "fixed", conf.int = TRUE, exponentiate = TRUE
  )
  or_rise <- broom.mixed::tidy(
    f_rise, effects = "fixed", conf.int = TRUE, exponentiate = TRUE
  )
  
  tibble(
    threshold = th,
    drop_friends_OR = or_drop$estimate[or_drop$term == "friends_z"],
    drop_friends_lo = or_drop$conf.low[or_drop$term == "friends_z"],
    drop_friends_hi = or_drop$conf.high[or_drop$term == "friends_z"],
    # The no-religion coefficient name must match the coding of b0600;
    # here "b06000" is used as a placeholder and may need adjustment.
    rise_norelig_OR = or_rise$estimate[or_rise$term == "b06000"],
    rise_norelig_lo = or_rise$conf.low[or_rise$term == "b06000"],
    rise_norelig_hi = or_rise$conf.high[or_rise$term == "b06000"]
  )
}

# Load a minimal subset of packages for the robustness routine (no printing)
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(lme4)
  library(splines)
})

# --- helper: compute a numeric threshold from round-1 contributions ----------
# threshold_from(): maps a label ("mean", "median", "mean_m0.5sd", "mean_p0.5sd")
# to a numeric c* computed from the round-1 contribution distribution.
threshold_from <- function(x, method = c("mean","median","mean_m0.5sd","mean_p0.5sd")) {
  method <- match.arg(method)
  m  <- mean(x, na.rm = TRUE)
  sd <- stats::sd(x, na.rm = TRUE)
  switch(
    method,
    mean        = m,
    median      = stats::median(x, na.rm = TRUE),
    mean_m0.5sd = m - 0.5 * sd,
    mean_p0.5sd = m + 0.5 * sd
  )
}

# --- helper: coerce new (x) to same type/levels as reference (ref) -----------
# (Re-defined here; same idea as above but kept local to this section)
coerce_like <- function(x, ref) {
  if (is.factor(ref)) factor(x, levels = levels(ref))
  else if (is.numeric(ref)) suppressWarnings(as.numeric(x))
  else if (is.integer(ref)) suppressWarnings(as.integer(x))
  else x
}

# --- MAIN ROBUSTNESS FUNCTION -------------------------------------------------
# Second (main) definition of fit_at_threshold(): accepts either a numeric
# threshold or a label indicating how to compute c* from round-1 contributions.
fit_at_threshold <- function(th) {
  # 1) Numeric threshold c_star: either directly given or computed from a label
  c_star <- if (is.numeric(th)) {
    th
  } else {
    r1 <- pgg_data %>%
      filter(round_n == 1, is.finite(contributing)) %>%
      pull(contributing)
    threshold_from(r1, method = th)
  }
  
  # 2) Baseline covariates with names/types matching the working models
  baseline <- pgg_data %>%
    filter(round_n == 1) %>%
    distinct(player, village_code, .keep_all = TRUE) %>%
    transmute(
      player       = as.character(player),
      village_code = as.character(village_code),
      gender.x         = factor(gender),
      marital_status.x = factor(marital_status),
      b0600.x          = factor(b0600),   # religion
      b0200.x          = factor(b0200),   # ethnicity (factor here)
      access_routes.x  = factor(access_routes),
      age_z                 = as.numeric(scale(as.numeric(age))),
      friends_z             = as.numeric(scale(as.numeric(friends))),
      adversaries_z         = as.numeric(scale(as.numeric(adversaries))),
      network_density_fr_z  = as.numeric(scale(as.numeric(network_density_fr))),
      network_density_adv_z = as.numeric(scale(as.numeric(network_density_adv))),
      network_size_z        = as.numeric(scale(as.numeric(network_size))),
      b0100_z               = as.numeric(scale(as.numeric(b0100))),
      FI_z                  = as.numeric(scale(as.numeric(FI)))
    )
  
  # 3) Panel with lagged state, lagged peer mean, and spell_prev under this c*
  panel <- pgg_data %>%
    filter(is.finite(contributing)) %>%
    mutate(
      player       = as.character(player),
      village_code = as.character(village_code),
      group        = as.character(group),
      s            = as.integer(contributing >= c_star)
    ) %>%
    # LOO peer mean at time t
    group_by(group, round_n) %>%
    mutate(
      G = n(),
      peer_mean_c_t = ifelse(
        G > 1,
        (sum(contributing, na.rm = TRUE) - contributing) / (G - 1),
        NA_real_
      )
    ) %>%
    ungroup() %>%
    arrange(group, round_n) %>%
    group_by(group) %>%
    mutate(peer_mean_c_tm1 = dplyr::lag(peer_mean_c_t)) %>%
    ungroup() %>%
    # own lag, time index, and duration up to t-1
    arrange(player, round_n) %>%
    group_by(player) %>%
    mutate(
      s_lag       = dplyr::lag(s),
      r_t         = as.integer(round_n),
      run_id      = cumsum(dplyr::if_else(
        is.na(dplyr::lag(s)), TRUE, s != dplyr::lag(s)
      )),
      spell_index = dplyr::row_number(),
      spell_prev  = dplyr::lag(spell_index)
    ) %>%
    ungroup() %>%
    mutate(peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12) %>%
    left_join(baseline, by = c("player","village_code"))
  
  # 4) At-risk sets for Drop (s_lag==1) and Rise (s_lag==0)
  drop_df <- panel %>%
    filter(s_lag == 1) %>%
    mutate(y_drop = as.integer(s == 0))
  
  rise_df <- panel %>%
    filter(s_lag == 0) %>%
    mutate(y_rise = as.integer(s == 1))
  
  # 5) Ensure factor columns in the at-risk sets are well-formed
  fac_vars <- c("gender.x","marital_status.x","b0600.x","b0200.x","access_routes.x")
  drop_df[fac_vars] <- lapply(drop_df[fac_vars], function(x) factor(x))
  rise_df[fac_vars] <- lapply(rise_df[fac_vars], function(x) factor(x))
  
  # 6) Fit GLMMs at this threshold with the same fixed effects as the working
  # specifications (but defined locally here).
  form_drop <- y_drop ~
    peer_mean_c_tm1_scaled + r_t + ns(spell_prev, 3) +
    gender.x + marital_status.x + b0600.x + b0200.x + access_routes.x +
    age_z + friends_z + adversaries_z + network_density_fr_z + network_density_adv_z +
    network_size_z + b0100_z + FI_z +
    (1 | village_code)
  
  form_rise <- y_rise ~
    peer_mean_c_tm1_scaled + r_t + ns(spell_prev, 3) +
    gender.x + marital_status.x + b0600.x + b0200.x + access_routes.x +
    age_z + friends_z + adversaries_z + network_density_fr_z + network_density_adv_z +
    network_size_z + b0100_z + FI_z +
    (1 | village_code)
  
  mod_drop <- glmer(
    form_drop,
    data   = drop_df,
    family = binomial,
    nAGQ   = 0,
    control = glmerControl(
      optimizer   = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl     = list(maxfun = 2e4)
    )
  )
  
  summary(mod_drop)
  
  mod_rise <- glmer(
    form_rise,
    data   = rise_df,
    family = binomial,
    nAGQ   = 0,
    control = glmerControl(
      optimizer   = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl     = list(maxfun = 2e4)
    )
  )
  
  summary(mod_rise)
  
  # 7) Extract key odds ratios and random-effect SDs plus AIC/BIC
  pick_or <- function(mod, term) {
    b  <- suppressWarnings(coef(summary(mod))[term, "Estimate"])
    se <- suppressWarnings(coef(summary(mod))[term, "Std. Error"])
    if (is.na(b) || is.na(se)) return(c(OR = NA_real_, lo = NA_real_, hi = NA_real_))
    c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
  }
  
  out <- tibble(
    threshold     = if (is.numeric(th)) th else th,
    # Drop model fit and selected ORs
    drop_AIC      = AIC(mod_drop),
    drop_BIC      = BIC(mod_drop),
    drop_RE_SD    = as.numeric(sqrt(VarCorr(mod_drop)$village_code[1])),
    drop_OR_friends = pick_or(mod_drop, "friends_z")["OR"],
    drop_OR_gender  = pick_or(mod_drop, "gender.x1")["OR"],
    # Rise model fit and selected OR (e.g., for no religion)
    rise_AIC      = AIC(mod_rise),
    rise_BIC      = BIC(mod_rise),
    rise_RE_SD    = as.numeric(sqrt(VarCorr(mod_rise)$village_code[1])),
    rise_OR_norel  = {
      # If b0600.x has multiple levels, obtain one of the non-baseline levels
      # as a representative "no religion" coefficient.
      levs <- levels(rise_df$b0600.x)
      coef_name <- paste0("b0600.x",
                          levs[levs != levels(rise_df$b0600.x)[1]][1])
      pick_or(mod_rise, coef_name)["OR"]
    }
  )
  
  out
}

# ----------------- Run the robustness grid -----------------------------------
# th_grid here is a set of *methods* that threshold_from() can handle.
th_grid <- c("mean","median","mean_m0.5sd","mean_p0.5sd")
robust_tbl <- purrr::map_df(th_grid, fit_at_threshold)

robust_tbl

###############################################################################

###############################################################################
## 10) ALTERNATIVE LINKS & GEE  (robust, GEE-ready) ---------------------------
###############################################################################
# This section explores alternative link functions and population-averaged
# specifications:
# - Probit GLMMs
# - GLMs with player-clustered robust standard errors (sandwich estimator)

# --- Probit GLMMs with the same specification as the logit models -----------
mod_drop_probit <- glmer(
  update(form_drop, . ~ .),
  data   = drop_df,
  family = binomial(link = "probit"),
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)
summary(mod_drop_probit)

mod_rise_probit <- glmer(
  update(form_rise, . ~ .),
  data   = rise_df,
  family = binomial(link = "probit"),
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)
summary(mod_rise_probit)

# As a population-averaged robustness check we also estimate logistic GLMs with
# player-clustered (CR) standard errors (sandwich estimator). The resulting odds
# ratios are compared to those from the mixed-effects logit.

library(forcats)

# Simpler duration specification for GEE/GLM: spell_prev and its square,
# instead of a spline, to avoid potential issues with spline-based design matrices.
form_drop_gee <- y_drop ~
  peer_mean_c_tm1_scaled + r_t + spell_prev + I(spell_prev^2) +
  gender.x + marital_status.x + b0600.x + b0200.x + access_routes.x +
  age_z + friends_z + adversaries_z +
  network_density_fr_z + network_density_adv_z + network_size_z +
  b0100_z + FI_z

form_rise_gee <- y_rise ~
  peer_mean_c_tm1_scaled + r_t + spell_prev + I(spell_prev^2) +
  gender.x + marital_status.x + b0600.x + b0200.x + access_routes.x +
  age_z + friends_z + adversaries_z +
  network_density_fr_z + network_density_adv_z + network_size_z +
  b0100_z + FI_z


# ===== Robust alternative to GEE: GLM + cluster-robust SEs ====================
library(sandwich)
library(lmtest)

# GLM for Drop with the GEE-style specification; cluster-robust SEs by player
glm_drop <- glm(form_drop_gee, data = drop_df, family = binomial)
V_drop   <- sandwich::vcovCL(glm_drop, cluster = ~ player)  # clustered on player
ct_drop  <- lmtest::coeftest(glm_drop, vcov = V_drop)

# GLM for Rise with cluster-robust SEs
glm_rise <- glm(form_rise_gee, data = rise_df, family = binomial)
V_rise   <- sandwich::vcovCL(glm_rise, cluster = ~ player)
ct_rise  <- lmtest::coeftest(glm_rise, vcov = V_rise)

# Odds ratios with robust 95% CI from cluster-robust GLMs
or_table <- function(ct) {
  est <- ct[, "Estimate"]; se <- ct[, "Std. Error"]
  OR  <- exp(est)
  lo  <- exp(est - 1.96 * se)
  hi  <- exp(est + 1.96 * se)
  data.frame(
    term = rownames(ct),
    OR, lo, hi,
    z = ct[, "z value"],
    p = ct[, "Pr(>|z|)"]
  )
}
or_drop_cr <- or_table(ct_drop)
or_rise_cr <- or_table(ct_rise)

cat("\n--- GLM + cluster-robust SEs (DROP): first 10 rows ---\n")
print(head(or_drop_cr, 10))
cat("\n--- GLM + cluster-robust SEs (RISE): first 10 rows ---\n")
print(head(or_rise_cr, 10))



###############################################################################
## 11) EVENT-HISTORY (Cox) WITH VILLAGE FRAILTY + KM BY FRIENDS TERCILES -----
###############################################################################
# This section recasts the transition problem in a continuous-time event-history
# framework:
# - time_to_event() constructs spell-level time-to-first-event (drop/rise)
# - Cox proportional hazards models with village-level frailty (random effects)
# - Kaplan-Meier curves by terciles of friends_z for Drop

# Helper to compute time-to-first event for a given starting state:
# - start_state: initial regime (1=High, 0=Low)
# - target_state: state triggering the event (0 or 1)
time_to_event <- function(df, start_state = 1, target_state = 0){
  # Order by player and round; mark when the player is at risk and when event occurs
  dat <- df %>%
    dplyr::arrange(player, round_n) %>%
    dplyr::group_by(player) %>%
    dplyr::mutate(
      start_ok = first(s) == start_state,
      event    = dplyr::if_else(
        !is.na(dplyr::lag(s)) & dplyr::lag(s) == start_state & s == target_state,
        1L, 0L
      ),
      at_risk  = dplyr::lag(s) == start_state
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(at_risk))
  
  # Within each player, define at-risk spells and count time steps within them
  dat <- dat %>%
    dplyr::group_by(player) %>%
    dplyr::mutate(
      at_risk_id = cumsum(dplyr::if_else(
        is.na(dplyr::lag(at_risk)) | at_risk < dplyr::lag(at_risk), 1L, 0L
      ))
    ) %>%
    dplyr::group_by(player, at_risk_id) %>%
    dplyr::mutate(t = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # For each player, extract time to the first event (if any) and event indicator
  fe <- dat %>%
    dplyr::group_by(player) %>%
    dplyr::summarise(
      time  = ifelse(any(event == 1), min(t[event == 1]), max(t)),
      event = as.integer(any(event == 1)),
      village_code = first(village_code),
      friends_z    = first(friends_z),
      .groups      = "drop"
    )
  fe
}

# Spell-level time-to-first-event data for drops and rises
cox_drop <- time_to_event(panel, start_state = 1, target_state = 0)  # H->L
cox_rise <- time_to_event(panel, start_state = 0, target_state = 1)  # L->H

# Cox proportional hazards with village-level frailty (gamma) for Drop:
# baseline hazard of time to first drop, conditioning on friends_z, gender, religion
cox_drop_fit <- survival::coxph(
  Surv(time, event) ~ friends_z.x + gender.x + b0600.x + frailty(village_code),
  data = panel %>%
    dplyr::left_join(cox_drop, by = c("player","village_code")) %>%
    dplyr::distinct(player, .keep_all = TRUE)
)
summary(cox_drop_fit)

# Cox model for Rise: time to first rise, with FI_z and religion as covariates
cox_rise_fit <- survival::coxph(
  Surv(time, event) ~ b0600.x + FI_z + frailty(village_code),
  data = panel %>%
    dplyr::left_join(cox_rise, by = c("player","village_code")) %>%
    dplyr::distinct(player, .keep_all = TRUE)
)
summary(cox_rise_fit)

# Kaplan-Meier survival curves for time to first drop (H→L) by terciles of friends_z
km_df <- panel %>%
  dplyr::distinct(player, friends_z) %>%
  dplyr::mutate(
    friends_terc = cut(
      friends_z,
      breaks = quantile(friends_z, probs = c(0, .33, .66, 1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Low","Med","High")
    )
  ) %>%
  dplyr::left_join(cox_drop, by = "player")
fit_km <- survival::survfit(Surv(time, event) ~ friends_terc, data = km_df)
survminer::ggsurvplot(
  fit_km,
  data       = km_df,
  risk.table = TRUE,
  ggtheme    = theme_minimal(),
  title      = "Time to first drop (H→L) by friends terciles"
)


###############################################################################
## 12) HMM CROSS-CHECK (2-state, no covariates) ------------------------------
###############################################################################
# Hidden Markov model (HMM) cross-check:
# - Contributions are z-scored within rounds.
# - A 2-state Gaussian HMM is fit with depmixS4.
# - The Viterbi path is used to estimate a transition matrix.
# - The evolution of the high-emission state share by round is plotted.

library(dplyr)
library(tidyr)
library(ggplot2)
library(depmixS4)

# 1) Prepare panel for HMM: contribution standardized within each round
pan_hmm <- pgg_data %>%
  dplyr::select(player, round_n, contributing) %>%
  dplyr::filter(is.finite(contributing)) %>%
  dplyr::group_by(round_n) %>%
  dplyr::mutate(y_z = as.numeric(scale(contributing))) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(player, round_n) %>%
  dplyr::mutate(player = as.character(player))

# Number of observations per player (required by depmixS4)
ntimes <- as.integer(table(pan_hmm$player))

# 2) Fit a 2-state Gaussian HMM to standardized contributions
m_hmm   <- depmixS4::depmix(
  response = y_z ~ 1,
  data     = pan_hmm,
  nstates  = 2,
  family   = gaussian(),
  ntimes   = ntimes
)
fit_hmm <- depmixS4::fit(
  m_hmm,
  emcontrol = depmixS4::em.control(maxit = 300),
  verbose   = FALSE
)

# 3) Viterbi sequence of most likely hidden states
post <- depmixS4::posterior(fit_hmm, type = "viterbi")
pan_hmm$hmm_state <- as.integer(post$state)

# 4) Transition matrix from the Viterbi path (empirical Markov chain over states)
get_P_from_viterbi <- function(df, nstates = 2, laplace = 1e-6) {
  # Build player-wise transitions (from t to t+1)
  T <- df %>%
    dplyr::group_by(player) %>%
    dplyr::arrange(round_n, .by_group = TRUE) %>%
    dplyr::transmute(
      from = hmm_state,
      to   = dplyr::lead(hmm_state)
    ) %>%
    dplyr::filter(!is.na(to)) %>%
    dplyr::ungroup()
  
  # Force factor levels 1:nstates
  T$from <- factor(T$from, levels = 1:nstates)
  T$to   <- factor(T$to,   levels = 1:nstates)
  
  # nstates × nstates count matrix
  M <- with(T, table(from, to))
  M <- M + laplace                      # Laplace smoothing to avoid zero rows
  P <- sweep(M, 1, rowSums(M), "/")     # row-normalized transition matrix
  dimnames(P) <- list(
    from = paste0("S", 1:nstates),
    to   = paste0("S", 1:nstates)
  )
  P
}

P_vit <- get_P_from_viterbi(pan_hmm, 2)
cat("\n--- Viterbi-based transition matrix ---\n")
print(round(P_vit, 3))

# 5) Parametric transition matrix from the fitted HMM (if available)
get_P_param <- function(fm) {
  n <- length(fm@transition)
  P <- matrix(NA_real_, n, n)
  for (i in seq_len(n)) {
    # Predict transition probabilities for state i
    P[i, ] <- as.numeric(predict(fm@transition[[i]], type = "response"))
  }
  dimnames(P) <- list(
    from = paste0("S", 1:n),
    to   = paste0("S", 1:n)
  )
  P
}
P_par <- try(get_P_param(fit_hmm), silent = TRUE)
if (inherits(P_par, "try-error") || any(!is.finite(P_par))) {
  message("Parametric P not available; falling back to Viterbi-based P.")
  P_par <- P_vit
}
cat("\n--- HMM transition matrix (parametric if available) ---\n")
print(round(P_par, 3))

# 6) Share in state 2 (high-emission state) by round
hmm_share <- pan_hmm %>%
  dplyr::group_by(round_n) %>%
  dplyr::summarise(
    p_state2 = mean(hmm_state == 2),
    .groups  = "drop"
  )

ggplot(hmm_share, aes(round_n, p_state2)) +
  geom_line(size = 1.1) + geom_point() +
  labs(
    x     = "Round",
    y     = "Share in state 2",
    title = "HMM: share in high-emission state by round"
  ) +
  theme_minimal()


###############################################################################
## 13) PLACEBO LEAD TEST (peer_{t+1}) ----------------------------------------
###############################################################################
# Placebo test for dynamic peer effects:
# - Construct peer_mean_c_{t+1} (lead)
# - Regress s_t on s_{t-1}, peer_{t-1}, peer_{t+1}, controls, and random effects.
# - If the temporal ordering is correctly specified, the lead term should have
#   no predictive power once peer_{t-1} and controls are included.

# Build peer mean at t and lead it to obtain peer_mean_c_{t+1}
panel_lead <- panel %>%
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    peer_mean_c_t = ifelse(
      dplyr::n() > 1,
      (sum(contributing, na.rm = TRUE) - contributing) / (dplyr::n() - 1),
      NA_real_
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, round_n) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(peer_mean_c_tp1 = dplyr::lead(peer_mean_c_t)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(peer_mean_c_tp1_scaled = peer_mean_c_tp1 / 12) %>%
  dplyr::filter(
    !is.na(s_lag),
    is.finite(peer_mean_c_tm1_scaled),
    is.finite(peer_mean_c_tp1_scaled)
  )

# Placebo specification: s_t on s_{t-1} + peer_{t-1} + peer_{t+1} + controls,
# with random intercepts at village and player levels.
placebo_form <- as.formula(paste0(
  "s ~ s_lag + peer_mean_c_tm1_scaled + peer_mean_c_tp1_scaled + r_t + ",
  "gender.x + age_z + friends_z + adversaries_z + b0100_z + FI_z + ",
  "(1|village_code) + (1|player)"
))
mod_placebo <- glmer(
  placebo_form,
  data   = panel_lead,
  family = binomial,
  nAGQ   = 0,
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)
summary(mod_placebo)$coefficients["peer_mean_c_tp1_scaled",]
cat("\n--- Placebo lead peer effect (should be ~0) ---\n")
print(
  broom.mixed::tidy(mod_placebo, effects = "fixed") %>%
    dplyr::filter(term == "peer_mean_c_tp1_scaled")
)
