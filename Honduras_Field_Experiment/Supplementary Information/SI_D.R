## ==============================================================
##  Dynamics of switching and learning during sessions
##  (A) volatility        (B) direction among switchers
##  (C) learning in levels (D) robustness & counterfactuals
## ==============================================================

## ---- packages -------------------------------------------------
## Load main packages used throughout the switching / learning analysis.
pkgs <- c(
  "dplyr","tidyr","ggplot2","broom","broom.mixed","splines",
  "lme4","pROC","rsample","yardstick","ggeffects","purrr"
)

## Install any missing packages
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)

## Attach packages
invisible(lapply(pkgs, library, character.only = TRUE))

## Optional / faster FE + sequence analytics
## These are used for robustness and visualization but not strictly required.
if (!requireNamespace("fixest", quietly=TRUE))    try(install.packages("fixest"))
if (!requireNamespace("TraMineR", quietly=TRUE))  try(install.packages("TraMineR"))
if (!requireNamespace("performance", quietly=TRUE)) try(install.packages("performance"))

## Reproducible random numbers (for CV splits, k-means, etc.)
set.seed(123)

## Load the dataset used throughout the IV analysis (player × round panel).
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)


## ---- helpers you used earlier (safe re-define if missing) ----
## pred_fixed_glmer_aligned(): get predictions from a glmer model using
## only the fixed effects, carefully aligning factor levels and columns.
## This is used for model diagnostics, CV, and the Markov transition system.
if (!exists("pred_fixed_glmer_aligned")) {
  pred_fixed_glmer_aligned <- function(mod, df, type = c("response","link")) {
    type <- match.arg(type)
    
    ## Remove random effects from model formula to obtain fixed-effects part
    f_fix   <- lme4::nobars(formula(mod))
    
    ## Recover training model frame and contrasts to reuse factor encodings
    mf_tr   <- model.frame(mod)
    contr   <- attr(mf_tr,"contrasts")
    
    ## Ensure that factors in new data have the same levels as in training
    fac_tr  <- names(which(vapply(mf_tr, is.factor, logical(1))))
    for (v in fac_tr) if (v %in% names(df)) df[[v]] <- factor(df[[v]], levels = levels(mf_tr[[v]]))
    
    ## Build model matrix for new data with aligned contrasts
    mf_new  <- model.frame(f_fix, df, na.action = na.pass)
    Xnew    <- model.matrix(f_fix, data = mf_new, contrasts.arg = contr)
    
    ## Extract fixed effects and ensure Xnew has the same columns
    beta    <- fixef(mod)
    miss    <- setdiff(names(beta), colnames(Xnew))
    if (length(miss)) {
      ## Add zero-columns for any missing dummy variables
      Xnew <- cbind(
        Xnew,
        `colnames<-`(matrix(0, nrow = nrow(Xnew), ncol = length(miss)), miss)
      )
    }
    ## Column order must match the order of coefficients
    Xnew <- Xnew[, names(beta), drop = FALSE]
    
    ## Compute linear predictor for rows with complete covariates
    ok  <- stats::complete.cases(Xnew)
    eta <- rep(NA_real_, nrow(df))
    if (any(ok)) eta[ok] <- as.numeric(Xnew[ok, , drop = FALSE] %*% beta)
    
    ## Return on link or response scale
    if (type == "link") return(eta)
    plogis(eta)
  }
}

## calc_metrics_safe(): convenient summary of binary prediction quality:
## AUC, Brier score, and logistic calibration (intercept, slope).
if (!exists("calc_metrics_safe")) {
  calc_metrics_safe <- function(df, fit, yvar){
    ## Fixed-effects-only predictions for the binary outcome
    p  <- as.numeric(pred_fixed_glmer_aligned(fit, df, type = "response"))
    y  <- df[[yvar]]
    
    ## Keep only finite predictions and non-missing outcomes
    ok <- is.finite(p) & !is.na(y)
    y <- y[ok]; p <- p[ok]
    
    ## AUC and Brier score
    auc   <- as.numeric(pROC::auc(y, p))
    brier <- mean((y - p)^2)
    
    ## Calibration (regression of y on log-odds of p)
    eps   <- 1e-8
    lp    <- qlogis(pmin(pmax(p, eps), 1-eps))
    cal   <- glm(y ~ lp, family=binomial())
    
    list(
      auc           = auc,
      brier         = brier,
      cal_intercept = coef(cal)[1],
      cal_slope     = coef(cal)[2],
      y             = y,
      pred          = p
    )
  }
}

## or_table(): tidy table of odds ratios and confidence intervals
## for fixed effects from lme4::glmer models.
or_table <- function(fit){
  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) |>
    dplyr::select(term, estimate, conf.low, conf.high, p.value)
}


## ---- 0. Data slices we need ----------------------------------
## Compute c* = average round-1 contribution (used as H/L cut-off).
c_star <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::summarise(m = mean(contributing, na.rm = TRUE)) %>%
  dplyr::pull(m)

## Baseline covariates: one row per player (round 1 only).
## - Coerce IDs and factors to the intended types.
## - Construct z-scores for continuous covariates to ease interpretation.
baseline <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::distinct(player, .keep_all = TRUE) %>%  # one baseline row per player
  dplyr::transmute(
    player        = as.character(player),
    village_code  = as.character(village_code),
    gender        = factor(gender),
    marital_status= factor(marital_status),
    b0600         = factor(b0600),    # religion (re-leveled below)
    access_routes = factor(access_routes),
    age           = suppressWarnings(as.numeric(age)),
    friends       = suppressWarnings(as.numeric(friends)),
    adversaries   = suppressWarnings(as.numeric(adversaries)),
    network_density_fr  = suppressWarnings(as.numeric(network_density_fr)),
    network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
    network_size        = suppressWarnings(as.numeric(network_size)),
    b0100              = suppressWarnings(as.numeric(b0100)),  # education (years)
    b0200              = as.numeric(b0200),   # indigenous indicator (numeric dummy)
    FI                 = suppressWarnings(as.numeric(FI))      # financial autonomy index
  ) %>%
  ## Standardize continuous covariates used later in the models.
  dplyr::mutate(dplyr::across(
    c(
      age, friends, adversaries,
      network_density_fr, network_density_adv, network_size,
      b0100, FI
    ),
    ~ as.numeric(scale(.x)),
    .names = "{.col}_z"
  )) %>%
  ## Set reference religion to "2" (Catholic) if that level exists.
  ## Other levels are interpreted relative to Catholic.
  dplyr::mutate(b0600 = stats::relevel(b0600, ref = "2"))

## Construct the panel used in the transition analyses:
## - One row per player × round with state s_t: High (H) vs Low (L)
## - Group-level leave-one-out peer mean contributions at t and t-1
## - Spell length: time already spent in current regime at t-1
panel <- pgg_data %>%
  dplyr::filter(is.finite(contributing)) %>%                      # keep valid contributions
  dplyr::mutate(
    player       = as.character(player),
    village_code = as.character(village_code),
    group        = as.character(group),
    
    ## State indicator: s_t = 1 if contribution ≥ c*, 0 otherwise.
    s            = as.integer(contributing >= c_star),
    
    ## Ensure religion and gender are factors for later models.
    b0600        = as.factor(b0600),
    gender       = as.factor(gender)
  ) %>%
  ## Within each group-round, compute leave-one-out peer mean contribution at time t
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    G = dplyr::n(),
    peer_mean_c_t = ifelse(
      G > 1,
      (sum(contributing, na.rm = TRUE) - contributing) / (G - 1),  # LOO mean
      NA_real_
    )
  ) %>%
  dplyr::ungroup() %>%
  ## Lag the peer mean at the group level: peer_mean_c_{t-1}
  dplyr::arrange(group, round_n) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(peer_mean_c_tm1 = dplyr::lag(peer_mean_c_t)) %>%
  dplyr::ungroup() %>%
  ## Within-player lags and duration (spell) indices
  dplyr::arrange(player, round_n) %>%
  dplyr::group_by(player) %>%
  dplyr::mutate(
    ## Previous-period state s_{t-1}
    s_lag   = dplyr::lag(s),
    
    ## Time index (round number)
    r_t     = as.integer(round_n),
    
    ## Identify “runs” of consecutive same-state periods:
    ## increment run_id every time the state changes or at the first observation.
    run_id  = cumsum(dplyr::if_else(
      is.na(dplyr::lag(s)), TRUE, s != dplyr::lag(s)
    )),
    
    ## Position within the run (1, 2, 3, …) – here equal to cumulative length.
    spell_index = dplyr::row_number(),
    
    ## Spell length at t-1 (duration already spent in current regime)
    spell_prev  = dplyr::lag(spell_index)
  ) %>%
  dplyr::ungroup() %>%
  ## Merge baseline covariates (one row per player) into the full panel.
  dplyr::left_join(
    baseline %>% dplyr::select(
      player, village_code,
      gender, marital_status, b0600, b0200, access_routes,
      ends_with("_z")
    ),
    by = c("player", "village_code")
  ) %>%
  ## Scale peer mean by the per-round endowment (assumed to be 12 units)
  dplyr::mutate(peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12)

## Potential transitions (for rounds t ≥ 2; i.e., where s_lag is observed).
switch_df <- panel %>%
  dplyr::arrange(player, round_n) %>%
  dplyr::filter(!is.na(s_lag)) %>%
  dplyr::mutate(
    ## Binary indicator for any switch in this round (H↔L).
    y_switch = as.integer(s != s_lag),
    
    ## Direction of switch: HL (1) vs LH (0); non-switches set to NA.
    dir_hl   = dplyr::case_when(
      s_lag==1 & s==0 ~ 1L,   # High → Low
      s_lag==0 & s==1 ~ 0L,   # Low → High
      TRUE ~ NA_integer_
    ),
    
    ## Covariates mirroring the GLMMs used in the main text
    peer_mean_c_tm1_scaled = peer_mean_c_tm1_scaled,
    r_t = as.integer(round_n)
  ) %>%
  ## Keep only observations with well-defined scaled peer mean
  dplyr::filter(is.finite(peer_mean_c_tm1_scaled))


## ==============================================================
## A1. Volatility: probability of ANY switch this round
## ==============================================================

## Logit mixed model for probability of any switch in round t:
## - Outcome: y_switch (1 if state changed, 0 otherwise)
## - Key predictors: peer mean (scaled), round index, spell length (spline)
## - Controls: baseline socio-demographics and network covariates
## - Random intercept: village_code (to capture village-level heterogeneity)
form_any <- y_switch ~ peer_mean_c_tm1_scaled + r_t + ns(spell_prev,3) +
  gender.x + marital_status.x + b0600.x + b0200.x + access_routes.x +
  age_z + friends_z + adversaries_z +
  network_density_fr_z + network_density_adv_z + network_size_z +
  b0100_z + FI_z + (1|village_code)

mod_any  <- lme4::glmer(
  form_any, data = switch_df, family = binomial,
  nAGQ = 0,
  control = lme4::glmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)

summary(mod_any)

## Odds-ratio table for fixed effects
any_or   <- or_table(mod_any)

## In-sample predictive performance using fixed-effects predictions only
any_fit  <- calc_metrics_safe(switch_df, mod_any, "y_switch")

cat(
  "\n[ANY switch] AUC =", round(any_fit$auc,3),
  "| Brier =", round(any_fit$brier,4),
  "| Cal(interc,slope) =",
  paste(round(any_fit$cal_intercept,3), round(any_fit$cal_slope,3), sep="/"), "\n"
)

## Figure: distribution of # switches per player
## Counts how many times each player changes High/Low state over the session.
fig_switch_hist <- switch_df %>%
  dplyr::group_by(player) %>% dplyr::summarise(n_switch = sum(y_switch), .groups="drop") %>%
  ggplot(aes(n_switch)) +
  geom_histogram(binwidth=1, boundary=-0.5) +
  labs(
    x = "# switches per player (t=2..10)",
    y = "Count",
    title = "Distribution of switches"
  ) +
  theme_minimal()

## Figure: stacked state-distribution by round (HH/HL/LH/LL)
## For each t and state pair (s_{t-1}, s_t), compute the share of HH, HL, LH, LL.
trans_stack <- switch_df %>%
  dplyr::mutate(
    pair = dplyr::case_when(
      s_lag==1 & s==1 ~ "HH",
      s_lag==1 & s==0 ~ "HL",
      s_lag==0 & s==1 ~ "LH",
      TRUE            ~ "LL"
    )
  ) %>%
  dplyr::count(round_n, pair, name="n") %>%
  dplyr::group_by(round_n) %>%
  dplyr::mutate(p = n/sum(n)) %>%
  dplyr::ungroup()

fig_state_stack <- ggplot(trans_stack, aes(round_n, p, fill = pair)) +
  geom_area(position="fill", colour=NA) +
  scale_y_continuous(labels=scales::percent) +
  labs(
    x="Round", y="Share", fill="Transition",
    title="State distribution by round"
  ) +
  theme_minimal()


## ==============================================================
## A2. Direction among switchers: HL (=1) vs LH (=0)
## ==============================================================

## Restrict to true switchers with defined direction.
dir_df <- switch_df %>% dplyr::filter(!is.na(dir_hl))

## Logit mixed model for switch direction:
## - Outcome: dir_hl = 1 if switch is High→Low, 0 if Low→High.
## - Same covariates as in the volatility model.
form_dir <- dir_hl ~ peer_mean_c_tm1_scaled + r_t + ns(spell_prev,3) +
  gender.x + marital_status.x + b0600.x + b0200.x + access_routes.x +
  age_z + friends_z + adversaries_z +
  network_density_fr_z + network_density_adv_z + network_size_z +
  b0100_z + FI_z + (1|village_code)

mod_dir  <- lme4::glmer(
  form_dir, data = dir_df, family = binomial,
  nAGQ = 0,
  control = lme4::glmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)

summary(mod_dir)

## In-sample predictive fit for direction model (HL vs LH)
dir_or   <- or_table(mod_dir)
dir_fit  <- calc_metrics_safe(dir_df, mod_dir, "dir_hl")

cat(
  "\n[HL vs LH] AUC =", round(dir_fit$auc,3),
  "| Brier =", round(dir_fit$brier,4), "\n"
)

## Figure: heatmap of HL vs LH counts by round
dir_heat <- dir_df %>%
  dplyr::mutate(direction = ifelse(dir_hl==1, "HL (drop)", "LH (rise)")) %>%
  dplyr::count(round_n, direction, name="n")

fig_dir_heat <- ggplot(dir_heat, aes(round_n, direction, fill = n)) +
  geom_tile() +
  geom_text(aes(label=n), colour="white", fontface="bold") +
  scale_fill_gradient(low="#c6dbef", high="#08519c") +
  labs(
    x="Round", y="",
    title="Timing of switches by direction"
  ) +
  theme_minimal()


## ==============================================================
## A3. Learning in levels: updating model
## ==============================================================

## --- (3) Learning in levels: player FE + village×round FE, clustered at group ---

library(fixest)

## Build the continuous level panel (contributions, not just High/Low)
## - c_tm1: own lagged contribution
## - peer_mean_c_tm1: leave-one-out group mean at t-1
## - Factorize IDs for FE estimation
level_df <- panel %>%
  dplyr::arrange(player, round_n) %>%
  dplyr::group_by(player) %>%
  dplyr::mutate(c_tm1 = dplyr::lag(contributing)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(is.finite(c_tm1), is.finite(peer_mean_c_tm1)) %>%
  dplyr::mutate(
    player       = factor(player),
    village_code = factor(village_code),
    round_n      = factor(round_n)
  )

## OPTION A (explicit): create a village×round FE and include it.
## This captures village shocks that vary by round (e.g. local norms updated each round).
level_df <- level_df %>%
  dplyr::mutate(vill_round = interaction(village_code, round_n, drop = TRUE))

## Fixed-effects model:
## - Outcome: contributing (in L)
## - Covariates: peer mean at t-1, own contribution at t-1
## - FE: player and village×round
## - Clustered s.e.: group (matching experimental grouping).
level_fe <- feols(
  contributing ~ peer_mean_c_tm1 + c_tm1 | player + vill_round,
  cluster = ~ group, data = level_df
)

## OPTION B (one‑liner): use FE interaction operator ^ in the FE slot.
## Equivalent to vill_round, but specified directly in the formula.
level_fe_alt <- feols(
  contributing ~ peer_mean_c_tm1 + c_tm1 | player + village_code^round_n,
  cluster = ~ group, data = level_df
)

## Tidy table of fixed-effects coefficients and confidence intervals.
summ_fe <- broom::tidy(level_fe, conf.int = TRUE)
print(level_fe)


## ==============================================================
## A4. Cross-level interactions (single-stage)
## ==============================================================

## Mixed model for updating in levels with cross-level interactions:
## - Random slopes for own and peer lag by player.
## - Village×round random intercept.
## - Interactions between lagged contributions and baseline covariates.
form_xlev <- contributing ~ c_tm1 + peer_mean_c_tm1 +
  gender.x + b0100_z + FI_z +
  c_tm1:gender.x + peer_mean_c_tm1:gender.x +
  c_tm1:b0100_z   + peer_mean_c_tm1:b0100_z +
  c_tm1:FI_z      + peer_mean_c_tm1:FI_z +
  (c_tm1 + peer_mean_c_tm1 | player) + (1 | village_code:round_n)

mod_xlev <- lme4::lmer(
  form_xlev, data = level_df,
  control = lme4::lmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)

summary(mod_xlev)

## Extract only interaction terms (cross-level interactions).
xlev_tab <- broom.mixed::tidy(mod_xlev, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(grepl(":", term))  # show interactions only


## ==============================================================
## B5. Time-varying peer response
## ==============================================================

## Mixed model allowing the peer effect to vary over rounds:
## - Interaction peer_mean_c_tm1 * r_t.
## - Random slopes in peer_mean_c_tm1 and own-lag c_tm1 by player.
## - Village×round random intercept.
mod_timevary <- lme4::lmer(
  contributing ~ peer_mean_c_tm1*r_t + c_tm1 +
    (peer_mean_c_tm1 + c_tm1 | player) + (1|village_code:round_n),
  data = level_df,
  control = lme4::lmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)

summary(mod_timevary)

## Peer slope over rounds (marginal effect of peer mean on contribution).
me_peer_round <- ggeffects::ggpredict(
  mod_timevary,
  terms = c("r_t [2:10]","peer_mean_c_tm1")
)

fig_peer_round <- ggplot(me_peer_round, aes(x, predicted)) +
  geom_line(size=1.1) +
  labs(
    x="Round",
    y="Marginal effect of peer mean",
    title="Peer response over rounds (marginal)"
  ) +
  theme_minimal()


## ==============================================================
## B6. Nonlinearity checks (splines)
## ==============================================================

## Mixed model with spline terms for own and peer lag:
## - ns(peer_mean_c_tm1,3) and ns(c_tm1,3) allow flexible nonlinearities.
## - Random slopes by player, village×round random intercept.
mod_spline <- lme4::lmer(
  contributing ~ ns(peer_mean_c_tm1,3) + ns(c_tm1,3) +
    (peer_mean_c_tm1 + c_tm1 | player) + (1|village_code:round_n),
  data = level_df,
  control = lme4::lmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)

summary(mod_spline)

## Marginal predicted contribution as a function of peer mean (nonlinear).
me_nl_group <- ggeffects::ggpredict(mod_spline, terms = "peer_mean_c_tm1 [all]")

fig_nonlinear <- ggplot(me_nl_group, aes(x, predicted)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.15) +
  geom_line(size=1.1) +
  labs(
    x="LOO peer mean (t-1)",
    y="Predicted contribution",
    title="Nonlinear peer response (spline)"
  ) +
  theme_minimal()


## ==============================================================
## B7. Lead–lag falsification in the level model
## ==============================================================

## Construct lead and lag peer mean variables to test for spurious lead effects.
panel_lead_lvl <- panel %>%
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    peer_mean_c_t = ifelse(
      dplyr::n()>1,
      (sum(contributing, na.rm=TRUE)-contributing)/(dplyr::n()-1),
      NA_real_
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, round_n) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(peer_mean_c_tp1 = dplyr::lead(peer_mean_c_t)) %>%  # lead(t+1)
  dplyr::ungroup() %>%
  dplyr::arrange(player, round_n) %>%
  dplyr::group_by(player) %>%
  dplyr::mutate(c_tm1 = dplyr::lag(contributing)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    is.finite(c_tm1),
    is.finite(peer_mean_c_tm1),
    is.finite(peer_mean_c_tp1)
  )

## Lead-lag model:
## - Include both peer_mean_c_tm1 (lag) and peer_mean_c_tp1 (lead).
## - Significant lead effect would suggest reverse causality / omitted dynamics.
mod_lead <- lme4::lmer(
  contributing ~ peer_mean_c_tm1 + peer_mean_c_tp1 + c_tm1 +
    (peer_mean_c_tm1 + c_tm1 | player) + (1|village_code:round_n),
  data = panel_lead_lvl,
  control = lme4::lmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)

summary(mod_lead )

## Extract coefficient on the lead peer mean for reporting.
lead_coef <- broom.mixed::tidy(mod_lead, effects="fixed") %>%
  dplyr::filter(term=="peer_mean_c_tp1")
print(lead_coef)


## ==============================================================
## C8. Sequence analytics (TraMineR) – optional visual glue
## ==============================================================

## Sequence visualizations of High/Low paths, if TraMineR is available.
if (requireNamespace("TraMineR", quietly=TRUE)) {
  library(TraMineR)
  
  ## Build wide dataset with one row per player and one column per round (H/L).
  seq_wide <- panel %>%
    dplyr::select(player, round_n, s) %>%
    dplyr::mutate(state = ifelse(s==1,"H","L")) %>%
    dplyr::select(-s) %>%
    tidyr::pivot_wider(
      names_from = round_n, values_from = state,
      values_fill = "L", names_prefix = "r"
    ) %>%
    dplyr::arrange(player) %>% dplyr::distinct(player, .keep_all = TRUE)
  
  alph <- c("L","H"); labs <- c("Low","High")
  
  ## Define sequence object for TraMineR
  seqobj <- seqdef(seq_wide[,-1], alphabet = alph, states = alph, labels = labs)
  
  ## Figure: state distribution across rounds
  TraMineR::seqdplot(seqobj, with.legend = "right",
                     main = "State distribution by round (TraMineR)")
  
  ## Figure: transition rates (printed to console)
  TraMineR::seqtransn(seqobj)               # overall transition matrix
  TraMineR::seqtrate(seqobj, with.missing=FALSE) # per-round transition rates
}


## ==============================================================
## D9. Out-of-sample prediction (10-fold CV)
## ==============================================================

## --- Fixed-effects predictions aligned to a glmer fit (no RE) ---
## (Redefinition here is intentional for the CV block; same logic as above.)
pred_fixed_glmer_aligned <- function(mod, df, type = c("response","link")) {
  type  <- match.arg(type)
  
  ## Fixed-effects portion of the fitted glmer model
  f_fix <- lme4::nobars(formula(mod))
  
  ## Training model frame used to extract factor levels / contrasts
  mf_tr <- model.frame(mod)
  contr <- attr(mf_tr, "contrasts")
  
  ## Align factor levels in new data with those used in training
  fac_tr <- names(which(vapply(mf_tr, is.factor, logical(1))))
  for (v in fac_tr) if (v %in% names(df)) df[[v]] <- factor(df[[v]], levels = levels(mf_tr[[v]]))
  
  ## Build model matrix in new data with aligned contrasts
  mf_new <- model.frame(f_fix, df, na.action = na.pass)
  Xnew   <- model.matrix(f_fix, data = mf_new, contrasts.arg = contr)
  beta   <- fixef(mod)
  
  ## Add zero-columns for missing factor dummies if needed
  miss <- setdiff(names(beta), colnames(Xnew))
  if (length(miss)) {
    Xnew <- cbind(
      Xnew,
      `colnames<-`(matrix(0, nrow = nrow(Xnew), ncol = length(miss)), miss)
    )
  }
  Xnew <- Xnew[, names(beta), drop = FALSE]
  
  ## Compute fixed-effects-only linear predictor
  ok  <- stats::complete.cases(Xnew)
  eta <- rep(NA_real_, nrow(df))
  if (any(ok)) eta[ok] <- as.numeric(Xnew[ok, , drop = FALSE] %*% beta)
  
  if (type == "link") return(eta)
  plogis(eta)
}

## (i) Any-switch model: 10-fold CV for AUC/Brier via fixed-effects predictions.
set.seed(123)
vfolds <- rsample::vfold_cv(switch_df, v = 10)

cv_bin_metrics <- purrr::map_df(vfolds$splits, function(sp){
  ## Split into training and test folds
  train <- rsample::analysis(sp)
  test  <- rsample::assessment(sp)
  
  ## Refit the any-switch model on the training fold
  fit <- lme4::glmer(
    form_any, data = train, family = binomial, nAGQ = 0,
    control = lme4::glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 2e4)
    )
  )
  
  ## Fixed-effects-only predictions on the test fold
  p <- pred_fixed_glmer_aligned(fit, test, "response")
  y <- test$y_switch
  ok <- is.finite(p) & !is.na(y)
  
  tibble(
    auc   = as.numeric(pROC::auc(y[ok], p[ok])),
    brier = mean((y[ok] - p[ok])^2)
  )
})

## Aggregate CV metrics across folds.
cv_any <- cv_bin_metrics %>% summarise(AUC = mean(auc), Brier = mean(brier))
print(cv_any)


## ==============================================================
## D10. Counterfactual round-10 High share (transition system)
## ==============================================================

## Integrated Markov logit for transitions:
## - Outcome: s_t (High/Low)
## - Predictors: s_{t-1} and its interactions with covariates, plus random effects.
## - This gives round-specific transition probabilities used in Markov simulations.
if (!exists("mod_markov")) {
  markov_form <- as.formula(paste0(
    "s ~ s_lag + ",
    "peer_mean_c_tm1_scaled*s_lag + r_t*s_lag + ns(spell_prev,3)*s_lag + ",
    "gender.x*s_lag + marital_status.x*s_lag + b0600.x*s_lag + b0200.x*s_lag + access_routes.x*s_lag + ",
    "age_z*s_lag + friends_z*s_lag + adversaries_z*s_lag + ",
    "network_density_fr_z*s_lag + network_density_adv_z*s_lag + network_size_z*s_lag + ",
    "b0100_z*s_lag + FI_z*s_lag + (1|village_code) + (1|player)"
  ))
  
  markov_df <- panel %>%
    dplyr::filter(
      !is.na(s_lag),
      is.finite(peer_mean_c_tm1_scaled),
      is.finite(spell_prev)
    )
  
  mod_markov <- lme4::glmer(
    markov_form, data=markov_df, family=binomial,
    nAGQ=0,
    control = lme4::glmerControl(
      optimizer="nloptwrap",
      calc.derivs=FALSE,
      optCtrl=list(maxfun=2e4)
    )
  )
  
  summary(mod_markov)
}

## Helper: construct 2×2 transition matrix P from a fitted Markov logit
## by averaging predicted P(s_t=1 | s_{t-1}) for s_{t-1}=0 and s_{t-1}=1.
rowmean_P <- function(fit, df){
  p <- pred_fixed_glmer_aligned(fit, df, "response")
  out <- df %>% dplyr::mutate(.phat = p) %>%
    dplyr::group_by(s_lag) %>%
    dplyr::summarise(P_q1 = mean(.phat, na.rm=TRUE), .groups="drop") %>%
    dplyr::arrange(s_lag)
  
  P <- matrix(NA_real_, 2, 2, dimnames = list(from=c("L","H"), to=c("L","H")))
  ## For s_lag=0, probability of moving to H is P[1,2]; staying Low is 1 - this.
  P[1,2] <- out$P_q1[out$s_lag==0]; P[1,1] <- 1 - P[1,2]
  ## For s_lag=1, probability of staying High is P[2,2]; dropping to Low is 1 - this.
  P[2,2] <- out$P_q1[out$s_lag==1]; P[2,1] <- 1 - P[2,2]
  P
}

## Baseline transition matrix
P_base <- rowmean_P(mod_markov, markov_df)

## Counterfactual 1: +1 SD in friends_z for H-at-risk (players currently High).
cf1 <- markov_df
sd_f <- stats::sd(cf1$friends_z, na.rm=TRUE)
cf1$friends_z[cf1$s_lag==1] <- cf1$friends_z[cf1$s_lag==1] + sd_f
P_cf1 <- rowmean_P(mod_markov, cf1)

## Counterfactual 2: switch Catholic -> no religion for L-at-risk (players currently Low)
cf2 <- markov_df
if (is.factor(cf2$b0600.x)) {
  lv <- levels(cf2$b0600.x)
  if ("0" %in% lv) cf2$b0600.x[cf2$s_lag==0] <- factor("0", levels=lv)
}
P_cf2 <- rowmean_P(mod_markov, cf2)

## Counterfactual 3: +1 SD peer mean in early rounds (t <= 5)
cf3 <- markov_df
sd_m <- stats::sd(cf3$peer_mean_c_tm1_scaled, na.rm=TRUE)
cf3$peer_mean_c_tm1_scaled[cf3$r_t <= 5] <- cf3$peer_mean_c_tm1_scaled[cf3$r_t <= 5] + sd_m
P_cf3 <- rowmean_P(mod_markov, cf3)

## Forward iterate from round 1 to round 10 to get expected share High.
## pi1 = observed share High in round 1.
pi1 <- panel %>% dplyr::filter(round_n==1) %>%
  dplyr::summarise(pH = mean(s==1)) %>% dplyr::pull()

iterate10 <- function(P){
  ## Initial distribution (L,H) at round 1
  pi <- c(L=1-pi1, H=pi1)
  ## Apply transition matrix 9 times to reach round 10
  for(k in 1:9) pi <- as.numeric(pi %*% P)
  names(pi) <- c("L","H")
  pi
}

## High-share at round 10, baseline and counterfactuals.
s_base <- iterate10(P_base)["H"]
s_cf1  <- iterate10(P_cf1)["H"]
s_cf2  <- iterate10(P_cf2)["H"]
s_cf3  <- iterate10(P_cf3)["H"]

round10_tbl <- tibble(
  scenario = c("Baseline","Friends +1 SD (H at risk)","No‑religion for L at risk","Peer +1 SD (early rounds)"),
  pH_round10 = c(s_base, s_cf1, s_cf2, s_cf3)
)
print(round10_tbl)

## (Optional) bar plot of round‑10 High share under each scenario
fig_counterf <- ggplot(round10_tbl, aes(reorder(scenario, pH_round10), pH_round10)) +
  geom_col() + coord_flip() +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  labs(
    x="", y="Pr(High) at round 10",
    title="Counterfactual round‑10 High share"
  ) +
  theme_minimal()


## ==============================================================
## Save or print figures (uncomment to save)
## ==============================================================
## ggsave("fig_switch_hist.png", fig_switch_hist, width=6, height=4, dpi=300)
## ggsave("fig_state_stack.png", fig_state_stack, width=6, height=4, dpi=300)
## ggsave("fig_dir_heat.png",   fig_dir_heat,   width=6, height=3.8, dpi=300)
## ggsave("fig_slopes_density.png", fig_slopes_density, width=6, height=3.8, dpi=300)
## ggsave("fig_partial_group.png",   fig_partial_group, width=6, height=4, dpi=300)
## ggsave("fig_peer_round.png",      fig_peer_round,    width=6, height=4, dpi=300)
## ggsave("fig_nonlinear.png",       fig_nonlinear,     width=6, height=4, dpi=300)
## ggsave("fig_counterf.png",        fig_counterf,      width=6, height=4, dpi=300)


## ==============================================================
## Peer‑response heterogeneity by social position (friends / adversaries)
## --------------------------------------------------------------
## Question: Do players with more friends (or adversaries) respond
## differently to the peer signal?
## ==============================================================

## --- peer response × network position in the level model ---
mod_peer_x_net <- lme4::lmer(
  contributing ~ peer_mean_c_tm1*friends_z + peer_mean_c_tm1*adversaries_z + c_tm1 +
    (peer_mean_c_tm1 + c_tm1 | player) + (1 | village_code:round_n),
  data = level_df,
  control = lme4::lmerControl(
    optimizer = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl = list(maxfun = 2e4)
  )
)
summary(mod_peer_x_net)

## Marginal effect of peer mean at low/med/high friends and adversaries
library(ggeffects)

me_peer_friends <- ggpredict(
  mod_peer_x_net,
  terms = c("peer_mean_c_tm1 [-2:2 by=0.2]", "friends_z [-1.5,0,1.5]")
)

p_peer_friends <- ggplot(me_peer_friends, aes(x, predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "LOO peer mean (t-1)",
    y = "Predicted contribution",
    colour = "Friends (z)", fill = "Friends (z)",
    title = "Peer response by social embeddedness (friends)"
  ) +
  theme_minimal()

me_peer_advers <- ggpredict(
  mod_peer_x_net,
  terms = c("peer_mean_c_tm1 [-2:2 by=0.2]", "adversaries_z [-1.5,0,1.5]")
)

p_peer_advers <- ggplot(me_peer_advers, aes(x, predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, colour = NA) +
  geom_line(size = 1.1) +
  labs(
    x = "LOO peer mean (t-1)",
    y = "Predicted contribution",
    colour = "Adversaries (z)", fill = "Adversaries (z)",
    title = "Peer response by social tension (adversaries)"
  ) +
  theme_minimal()

p_peer_friends; p_peer_advers


## ==============================================================
## H2) State‑dependent updating (does last state change imitation?)
## --------------------------------------------------------------
## Question: Is peer‑responsiveness different after being High vs Low at t-1?
## ==============================================================

## Recompute last state (s_lag) directly in the level panel.
## Assumes: c_star defined; level_df has player, round_n (factor),
## contributing, c_tm1, peer_mean_c_tm1.
level_df2 <- level_df %>%
  dplyr::mutate(
    player        = as.character(player),
    round_n_int   = as.integer(as.character(round_n)),  # ensure numeric order
    village_code  = as.character(village_code),
    
    ## New village×round factor used as group random effect
    vill_round    = interaction(village_code, round_n_int, drop = TRUE),
    
    ## Current state s_now from continuous contribution using c_star cutoff
    s_now         = as.integer(contributing >= c_star)      # current state H(1)/L(0)
  ) %>%
  dplyr::arrange(player, round_n_int) %>%
  dplyr::group_by(player) %>%
  dplyr::mutate(s_lag = dplyr::lag(s_now)) %>%             # last state
  dplyr::ungroup() %>%
  dplyr::mutate(
    ## Factor label for last state: Low vs High at t-1
    s_lag_f = factor(s_lag, levels = c(0,1), labels = c("L(t-1)","H(t-1)"))
  ) %>%
  dplyr::filter(!is.na(s_lag_f))      # drop round-1 rows (no lag)

## State-dependent peer/own response model:
## - Interactions peer_mean_c_tm1*s_lag_f and c_tm1*s_lag_f.
## - Random slopes by player, village×round random intercept.
mod_state_dep <- lme4::lmer(
  contributing ~ peer_mean_c_tm1*s_lag_f + c_tm1*s_lag_f +
    (peer_mean_c_tm1 + c_tm1 | player) + (1 | vill_round),
  data = level_df2,
  control = lme4::lmerControl(
    optimizer="nloptwrap", calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)
summary(mod_state_dep)

## Plot: peer → contribution, by last state (L(t-1) vs H(t-1)).
me_state_peer <- ggeffects::ggpredict(
  mod_state_dep,
  terms = c("peer_mean_c_tm1 [-2:2 by=0.2]","s_lag_f")
)

p_state_peer <- ggplot(me_state_peer, aes(x, predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.15, colour=NA) +
  geom_line(size=1.1) +
  labs(
    x="LOO peer mean (t-1)",
    y="Predicted contribution",
    colour="Last state", fill="Last state",
    title="State-dependent peer response"
  ) +
  theme_minimal()
p_state_peer


## ==============================================================
## H3) “Learner types’’ from random‑slope clustering
## --------------------------------------------------------------
## Question: Are there distinct behavioral types based on (own‑lag vs peer‑lag) slopes?
## ==============================================================

## Use the cross-level random-slopes model (mod_xlev):
##   - Random slopes for c_tm1 and peer_mean_c_tm1 by player.
coef_mat <- coef(mod_xlev)$player[, c("c_tm1","peer_mean_c_tm1")]
colnames(coef_mat) <- c("slope_own","slope_peer")

## K-means clustering on standardized slopes to identify learner types.
set.seed(123)
km <- kmeans(scale(coef_mat), centers = 3, nstart = 50)

learner_types <- tibble(
  player = rownames(coef(mod_xlev)$player),
  slope_own  = coef_mat[,1],
  slope_peer = coef_mat[,2],
  type = factor(km$cluster, labels = c("Peer‑driven","Own‑sticky","Low‑responsive"))
)

## Scatter plot of learner types in (own-slope, peer-slope) space.
p_types <- ggplot(learner_types, aes(slope_own, slope_peer, colour = type)) +
  geom_point(alpha = .7) +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  labs(
    x = "Own‑lag slope", y = "Peer‑lag slope", colour = "Type",
    title = "Learner types from random‑slope clustering"
  ) +
  theme_minimal()
p_types

## Profile the types by selected baseline covariates.
baseline_players <- panel %>%
  dplyr::arrange(player, round_n) %>%
  dplyr::group_by(player) %>%
  dplyr::summarise(
    gender = dplyr::first(gender.x),
    b0100_z = dplyr::first(b0100_z),
    FI_z = dplyr::first(FI_z),
    friends_z = dplyr::first(friends_z),
    adversaries_z = dplyr::first(adversaries_z),
    .groups = "drop"
  )

type_profile <- learner_types %>%
  dplyr::left_join(baseline_players, by = "player") %>%
  dplyr::group_by(type) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_own  = mean(slope_own),
    mean_peer = mean(slope_peer),
    male_share = mean(as.numeric(as.character(gender))==1, na.rm=TRUE),
    educ = mean(b0100_z, na.rm=TRUE),
    FI   = mean(FI_z, na.rm=TRUE),
    friends = mean(friends_z, na.rm=TRUE),
    advers = mean(adversaries_z, na.rm=TRUE),
    .groups = "drop"
  )
print(type_profile)


## ==============================================================
## H4) Heterogeneity in switch direction (HL vs LH) by gender & friends
## --------------------------------------------------------------
## Question: Among switchers, do gender and friend ties tilt the direction?
## ==============================================================

## Extend direction model with interactions:
## - gender.x * friends_z
## - gender.x * peer_mean_c_tm1_scaled
## - friends_z * peer_mean_c_tm1_scaled
form_dir_het <- update(
  form_dir,
  . ~ . + gender.x*friends_z + gender.x*peer_mean_c_tm1_scaled + friends_z*peer_mean_c_tm1_scaled
)

mod_dir_het <- lme4::glmer(
  form_dir_het, data = dir_df, family = binomial,
  nAGQ = 0,
  control = lme4::glmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)
summary(mod_dir_het)

## Predicted Pr(HL) over peer mean, by gender.
me_dir_gender <- ggpredict(
  mod_dir_het,
  terms = c("peer_mean_c_tm1_scaled [-2:2 by=0.2]","gender.x")
)

p_dir_gender <- ggplot(me_dir_gender, aes(x, predicted, colour = group, fill = group)) +
  geom_ribbon(aes(ymin=conf.low, ymax=conf.high), alpha=.15, colour=NA) +
  geom_line(size=1.1) +
  labs(
    x="Peer mean (t-1, scaled)",
    y="Pr(HL | switch)",
    colour="Gender", fill="Gender",
    title="Direction among switchers: HL vs LH by gender × peer signal"
  ) +
  theme_minimal()
p_dir_gender


## ==============================================================
## H5) Subgroup counterfactuals (who benefits most?)
## --------------------------------------------------------------
## Question: Under the same counterfactuals, how does the round‑10 High
## share change by subgroup (e.g., gender or friend terciles)?
## ==============================================================

## -----------------------------------------------------------------
## Helper: fixed-effects-only predictions for lme4::glmer objects
## (redeclared here for convenience in the subgroup Markov analysis)
## -----------------------------------------------------------------
if (!exists("pred_fixed_glmer_aligned")) {
  pred_fixed_glmer_aligned <- function(mod, df, type = c("response","link")) {
    type <- match.arg(type)
    
    f_fix  <- lme4::nobars(stats::formula(mod))
    mf_tr  <- stats::model.frame(mod)
    contr  <- attr(mf_tr, "contrasts")
    
    ## Align factor levels
    fac_tr <- names(which(vapply(mf_tr, is.factor, logical(1))))
    for (v in fac_tr) if (v %in% names(df)) df[[v]] <- factor(df[[v]], levels = levels(mf_tr[[v]]))
    
    ## Build model matrix with same contrasts and order
    mf_new <- stats::model.frame(f_fix, df, na.action = stats::na.pass)
    Xnew   <- stats::model.matrix(f_fix, data = mf_new, contrasts.arg = contr)
    beta   <- lme4::fixef(mod)
    
    miss   <- setdiff(names(beta), colnames(Xnew))
    if (length(miss)) {
      Xnew <- cbind(
        Xnew,
        `colnames<-`(matrix(0, nrow = nrow(Xnew), ncol = length(miss)), miss)
      )
    }
    Xnew <- Xnew[, names(beta), drop = FALSE]
    
    ok   <- stats::complete.cases(Xnew)
    eta  <- rep(NA_real_, nrow(df))
    if (any(ok)) eta[ok] <- as.numeric(Xnew[ok, , drop = FALSE] %*% beta)
    
    if (type == "link") return(eta)
    plogis(eta)
  }
}

## -----------------------------------------------------------------
## Row-mean transition matrix from a fitted integrated Markov logit,
## for a given subset of rows (subgroup selection via idx).
## -----------------------------------------------------------------
rowmean_P_by <- function(fit, df, idx) {
  df_sub <- df[idx, , drop = FALSE]
  p <- pred_fixed_glmer_aligned(fit, df_sub, "response")
  
  out <- df_sub %>%
    dplyr::mutate(.phat = p) %>%
    dplyr::group_by(s_lag) %>%
    dplyr::summarise(P_q1 = mean(.phat, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(s_lag)
  
  P <- matrix(NA_real_, 2, 2, dimnames = list(from = c("L","H"), to = c("L","H")))
  P[1,2] <- out$P_q1[out$s_lag == 0]; P[1,1] <- 1 - P[1,2]
  P[2,2] <- out$P_q1[out$s_lag == 1]; P[2,1] <- 1 - P[2,2]
  P
}

## Forward iteration to round 10 given initial High share pi1
iterate10 <- function(P, pi1) {
  pi <- c(L = 1 - pi1, H = pi1)
  for (k in 1:9) pi <- as.numeric(pi %*% P)
  unname(pi[2])
}

## -----------------------------------------------------------------
## Build ‘friends’ terciles safely (no length mismatch, no name clash)
## -----------------------------------------------------------------
friends_baseline <- panel %>%
  dplyr::mutate(player = as.character(player)) %>%
  dplyr::group_by(player) %>%
  dplyr::summarise(friends_z_base = dplyr::first(friends_z), .groups = "drop")

## Join baseline friends_z to markov_df and compute terciles on the joined column
## Use ntile() -> robust when quantiles tie, includes lowest by design.
markov_df2 <- markov_df %>%
  dplyr::mutate(player = as.character(player)) %>%
  dplyr::left_join(friends_baseline, by = "player") %>%
  dplyr::mutate(
    fr_terc_id = dplyr::ntile(friends_z_base, 3),
    fr_terc    = factor(fr_terc_id, levels = 1:3, labels = c("Low","Med","High"))
  )

## (Optional sanity check)
## table(is.na(markov_df2$friends_z_base))

## -----------------------------------------------------------------
## Round‑10 High share by gender (robust coding)
## -----------------------------------------------------------------
male_levels <- c("1","male","Male","M","m")

## Construct a robust male indicator from gender.x (handles factor/character/numeric).
is_male <- if (is.factor(markov_df2$gender.x) || is.character(markov_df2$gender.x)) {
  tolower(as.character(markov_df2$gender.x)) %in% tolower(male_levels)
} else {
  as.numeric(markov_df2$gender.x) == 1
}

## Transition matrices by gender
P_male   <- rowmean_P_by(mod_markov, markov_df2, is_male)
P_female <- rowmean_P_by(mod_markov, markov_df2, !is_male)

## Initial High share (same pi1 as above)
pi1 <- panel %>% dplyr::filter(round_n == 1) %>%
  dplyr::summarise(pH = mean(s == 1)) %>% dplyr::pull()

tbl_gender <- tibble::tibble(
  subgroup    = c("Male","Female"),
  pH_round10  = c(iterate10(P_male, pi1), iterate10(P_female, pi1))
)
print(tbl_gender)

## -----------------------------------------------------------------
## Round‑10 High share by baseline friends tercile
## -----------------------------------------------------------------
res_terc <- dplyr::bind_rows(lapply(levels(markov_df2$fr_terc), function(L) {
  idx <- which(markov_df2$fr_terc == L)
  if (!length(idx))
    return(tibble::tibble(subgroup = paste("Friends", L), pH_round10 = NA_real_))
  P <- rowmean_P_by(mod_markov, markov_df2, idx)
  tibble::tibble(subgroup = paste("Friends", L), pH_round10 = iterate10(P, pi1))
}))
print(res_terc)


## ==============================================================
## Additional robustness: alternative random-effects structures
## ==============================================================

## Option 1: keep only village RE and round FE (factor round_n_int).
mod_state_dep2 <- lme4::lmer(
  contributing ~ peer_mean_c_tm1*s_lag_f + c_tm1*s_lag_f +
    (peer_mean_c_tm1 + c_tm1 | player) + (1 | village_code) + factor(round_n_int),
  data = level_df2,
  control = lme4::lmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)
summary(mod_state_dep2)

## Option 2: separate REs for village and round (instead of their interaction)
mod_state_dep3 <- lme4::lmer(
  contributing ~ peer_mean_c_tm1*s_lag_f + c_tm1*s_lag_f +
    (peer_mean_c_tm1 + c_tm1 | player) + (1 | village_code) + (1 | round_n_int),
  data = level_df2,
  control = lme4::lmerControl(
    optimizer="nloptwrap",
    calc.derivs=FALSE,
    optCtrl=list(maxfun=2e4)
  )
)
summary(mod_state_dep3)


## ==============================================================
## Alternative specifications for switching and levels
## (explicit village/group/player RE structures)
## ==============================================================

library(lme4); library(splines); library(dplyr)

## --------------------------------------------------------------
## B.1 Any switch (HL or LH) with village/group/player RE
## --------------------------------------------------------------
sw_df <- panel %>%
  arrange(player, round_n) %>%
  mutate(
    ## y_sw = 1 if current state differs from previous round,
    ##         0 otherwise (including first round, treated as 0).
    y_sw  = as.integer(s != lag(s)),
    y_sw  = ifelse(is.na(y_sw), 0L, y_sw),
    
    ## Ensure factors for random-effects grouping variables
    village_code = factor(village_code),
    group        = factor(group),
    player       = factor(player)
  )

mod_anyswitch <- glmer(
  y_sw ~ peer_mean_c_tm1_scaled + r_t + ns(spell_prev, 3) +
    gender.x + age_z + friends_z + adversaries_z +
    network_density_fr_z + network_density_adv_z + network_size_z +
    b0100_z + FI_z + marital_status.x + b0200.x + b0600.x + access_routes.x +
    (1 | village_code) + (1 | group) + (1 | player),
  data = sw_df %>% filter(!is.na(y_sw), is.finite(peer_mean_c_tm1_scaled)),
  family = binomial,
  control = glmerControl(
    optimizer = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl = list(maxfun = 2e4)
  )
)
summary(mod_anyswitch)

## --------------------------------------------------------------
## B.2 Direction among switchers (HL=1 vs LH=0) with village/group/player RE
## --------------------------------------------------------------
dir_df <- panel %>%
  arrange(player, round_n) %>%
  mutate(
    y_sw  = as.integer(s != lag(s)),
    ## y_dir = 1 for HL (drop), 0 for LH (rise), NA for non-switchers.
    y_dir = ifelse(y_sw == 1L, as.integer(lag(s) == 1L & s == 0L), NA_integer_),
    village_code = factor(village_code),
    group        = factor(group),
    player       = factor(player)
  ) %>%
  filter(y_sw == 1L, is.finite(peer_mean_c_tm1_scaled))

mod_dir <- glmer(
  y_dir ~ peer_mean_c_tm1_scaled + r_t + ns(spell_prev, 3) +
    gender.x + age_z + friends_z + adversaries_z +
    network_density_fr_z + network_density_adv_z + network_size_z +
    b0100_z + FI_z + marital_status.x + b0200.x + b0600.x + access_routes.x +
    (1 | village_code) + (1 | group) + (1 | player),
  data = dir_df,
  family = binomial,
  control = glmerControl(
    optimizer = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl = list(maxfun = 2e4)
  )
)
summary(mod_dir)

## Reporting tip (kept as a code comment for clarity):
## Odds ratios are per one-unit change in the covariate. For variables measured
## in units of the endowment (12 Lempiras), the implied OR per 1 L is OR^(1/12).


## --------------------------------------------------------------
## B.3 Levels (updating) models — FE in lempiras + LMM with RE
## --------------------------------------------------------------

## Fixed‑effects model (player FE and village×round FE; cluster s.e. at group).
## If you use fixest (recommended for speed with many FEs):
##   install.packages("fixest")
library(fixest)

lev_df <- panel %>%
  filter(
    is.finite(contributing),
    is.finite(peer_mean_c_tm1),
    is.finite(lag(contributing))
  ) %>%
  mutate(
    c_own_lag   = lag(contributing),        # own contribution at t-1 (in L)
    c_group_lag = peer_mean_c_tm1,          # peer mean at t-1 (in L)
    vill_round  = interaction(village_code, round_n, drop = TRUE)
  )

fe_levels <- feols(
  contributing ~ c_group_lag + c_own_lag | player + vill_round,
  cluster = ~ group,
  data = lev_df
)
summary(fe_levels)

## Random‑effects LMM with player random slopes and village/group/village×round RE:
library(lme4)

lmm_levels <- lmer(
  contributing ~ c_group_lag + c_own_lag +
    (1 | village_code) + (1 | group) + (1 | village_code:round_n) +
    (1 | player) +
    (0 + c_group_lag || player) + (0 + c_own_lag || player),
  data = lev_df,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
summary(lmm_levels)
