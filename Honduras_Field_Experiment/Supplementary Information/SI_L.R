# ---------------------------------------------------------------------------
# Variance decomposition across village / group / player
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr); library(lme4)})

# Main experimental panel at the individual × round level
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

GAME_GROUP_VAR <- "group"
CONTRIB_VAR    <- "contributing"
ROUND_VAR      <- "round_n"
PLAYER_VAR     <- "player"
VILLAGE_VAR    <- "village_code"

FRIENDS_VAR    <- "friends"
AUTONOMY_VAR   <- "finauto_q_2_temp"
GENDER_VAR     <- "gender"

MAX_CONTRIB    <- 12
ALPHA_PARAM    <- 0.5   # exponent in the altruism utility term c^α

## ---------------------------------------------------------------------------
## Helper functions (panel / game configuration)
## ---------------------------------------------------------------------------

## Simple column-existence check used throughout.
check_cols <- function(df, vars) {
  miss <- setdiff(vars, names(df))
  if(length(miss)) stop("Missing: ", paste(miss, collapse=", "))
  invisible(TRUE)
}

## Infer group size N as the modal per-round group size.
infer_N <- function(df) {
  dt <- as.data.table(df)
  dt[, .N, by=c(GAME_GROUP_VAR, ROUND_VAR)
  ][, .N, by=N][order(-N)][1,N]
}

## Add leave-one-out peer mean and lags, using data.table for speed.
add_peer_means <- function(df) {
  check_cols(df, c(GAME_GROUP_VAR,ROUND_VAR,PLAYER_VAR,CONTRIB_VAR))
  dt <- as.data.table(df)
  
  ## Group-level totals and sizes
  dt[, grp_total := sum(get(CONTRIB_VAR)), by=c(GAME_GROUP_VAR,ROUND_VAR)]
  dt[, grp_n     := .N,                        by=c(GAME_GROUP_VAR,ROUND_VAR)]
  
  ## LOO peer mean at time t
  dt[, peer_mean_t := fifelse(
    grp_n>1,
    (grp_total - get(CONTRIB_VAR))/pmax(grp_n-1,1),
    NA_real_) ]
  
  ## Order within player and add lagged own contribution and lagged peer mean
  setorderv(dt, c(PLAYER_VAR,ROUND_VAR))
  dt[, `:=`(
    self_lag      = shift(get(CONTRIB_VAR),1L,type="lag"),
    peer_mean_tm1 = shift(peer_mean_t,1L,type="lag")
  ), by=PLAYER_VAR]
  
  as_tibble(dt)
}

library(dplyr)

## ------------------------------------------------------------------
## Helper: check that a data frame has the required columns
## ------------------------------------------------------------------
.check_cols <- function(df, vars) {
  vars <- unname(vars)
  missing <- setdiff(vars, names(df))
  if (length(missing) > 0) {
    stop("`df` is missing required columns: ",
         paste(missing, collapse = ", "),
         call. = FALSE)
  }
  invisible(TRUE)
}

## ------------------------------------------------------------------
## Link generic variable names used in downstream modules (finite-horizon,
## bifurcation, threshold-sensitivity, etc.) to your COL_* config.
## (only set them if they don't already exist)
## ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
  library(tidyr)
  library(rlang)
  library(data.table)
  library(purrr)     # <-- add this
})

library(depmixS4)   # add this near the top of your script

# ---------------------------------------------------------------------------
# Variance decomposition: village / group / player levels
# ---------------------------------------------------------------------------
variance_decomp_levels <- function(df, nested = TRUE){
  stopifnot(all(c("village_code","group","player","contributing") %in% names(df)))
  
  dat <- df |>
    transmute(y = contributing,
              village = factor(village_code),
              group   = factor(group),
              player  = factor(player))
  
  if (nested) {
    # Explicit nesting: village : group : player
    dat <- dat |>
      mutate(group_v  = interaction(village, group,  drop = TRUE),
             player_v = interaction(village, group, player, drop = TRUE))
    m <- lmer(y ~ 1 + (1|village) + (1|group_v) + (1|player_v),
              data = dat, REML = TRUE)
  } else {
    # Crossed (use this if IDs are globally unique across villages)
    m <- lmer(y ~ 1 + (1|village) + (1|group) + (1|player),
              data = dat, REML = TRUE)
  }
  
  vc <- as.data.frame(VarCorr(m)) |>
    mutate(share = vcov / sum(vcov))
  
  out_clean <- vc |>
    dplyr::select(grp, vcov, sdcor, share) |>
    mutate(share_pct = round(100*share, 1))
  
  list(model = m,
       varcomp = vc,
       varcomp_clean = out_clean,
       singular = isSingular(m))
}

# ---------------------------------------------------------------------------
# Welfare accounting (observed, full cooperation, structural prediction, policy)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({library(dplyr)})

welfare_accounting <- function(df,
                               b = 2, kappa = 1, c_max = 12,
                               struct_params = NULL, policy_m = 0){
  
  stopifnot(all(c("group","round_n","player","contributing") %in% names(df)))
  
  # Detect N from the data (modal per-group size across rounds)
  N <- df |>
    count(group, round_n, name = "n") |>
    count(n, sort = TRUE) |>
    slice(1) |>
    pull(n)
  
  # Observed per-player payoff
  obs <- df |>
    group_by(group, round_n) |>
    mutate(grp_total = sum(contributing, na.rm = TRUE)) |>
    ungroup() |>
    mutate(material = (b/N)*grp_total - kappa*contributing) |>
    summarise(payoff_per_player = mean(material, na.rm = TRUE), .groups = "drop") |>
    summarise(payoff_per_player = mean(payoff_per_player, na.rm = TRUE)) |>
    mutate(scenario = "Observed", m = 0)
  
  # Full cooperation (everyone contributes c_max)
  full <- tibble(
    scenario = "FullCoop",
    payoff_per_player = b*c_max - kappa*c_max,
    m = 0
  )
  
  # Structural prediction (if structural params provided)
  pred <- NULL
  if (!is.null(struct_params)) {
    sp <- struct_params
    
    if ("d_i"   %in% names(sp) && !"d_hat"   %in% names(sp)) sp$d_hat   <- sp$d_i
    if ("phi_i" %in% names(sp) && !"phi_hat" %in% names(sp)) sp$phi_hat <- sp$phi_i
    stopifnot(all(c("player","d_hat","phi_hat") %in% names(sp)))
    
    # Compute implied c^H = [ (0.5 d) / ( (kappa - b/N) + 2*phi ) ]^2 truncated to [0, c_max]
    Delta <- (kappa - b/N)
    if (Delta <= 0) {
      warning("kappa - b/N <= 0; structural c^H undefined. Skipping structural scenario.")
    } else {
      sp <- sp |>
        transmute(
          player,
          cH = pmin(
            c_max,
            pmax(
              0,
              ((0.5 * d_hat) / pmax(Delta + 2 * phi_hat, 1e-8))^2
            )
          )
        )
      
      gmap <- df |>
        filter(round_n == 1) |>
        distinct(player, group)
      
      sp <- sp |>
        left_join(gmap, by = "player") |>
        group_by(group) |>
        mutate(grp_mean_cH = mean(cH, na.rm = TRUE)) |>
        ungroup() |>
        mutate(pay = (b/N)*(cH + (N-1)*grp_mean_cH) - kappa*cH)
      
      pred <- tibble(
        scenario = "StructPred",
        payoff_per_player = mean(sp$pay, na.rm = TRUE),
        m = 0
      )
    }
  }
  
  # Policy-m: per-player payoff if everyone contributes c_max under a subsidy m
  kappa_p <- max(kappa - policy_m, 0)
  pol <- tibble(
    scenario = "Policy_m",
    payoff_per_player = b*c_max - kappa_p*c_max,
    m = policy_m
  )
  
  out <- bind_rows(obs, full, pred, pol) |>
    dplyr::select(scenario, payoff_per_player, m)
  
  return(out)
}

# --- Run welfare accounting -----------------------------------------------

# Pick whichever structural object you have available:
struct_guess <-
  if (exists("pgg_param_player")) pgg_param_player else
    if (exists("results") && "struct_params" %in% names(results)) results$struct_params else
      NULL

pgg_data$player <- as.character(pgg_data$player)

ad3 <- welfare_accounting(pgg_data,
                          b = 2, kappa = 1, c_max = 12,
                          struct_params = struct_guess,
                          policy_m = 0.5)

# ---------------------------------------------------------------------------
# Peer norm models with friend / adversary moderation
# ---------------------------------------------------------------------------

# (requires edge lists with columns from,to)
suppressPackageStartupMessages({
  library(dplyr)
  library(lme4)
})

# -- Peer mean (group excluding self), lagged --
if (!exists("add_peer_means")) {
  add_peer_means <- function(df){
    stopifnot(all(c("group","round_n","player","contributing") %in% names(df)))
    df %>%
      group_by(group, round_n) %>%
      mutate(g_sum = sum(contributing, na.rm=TRUE),
             g_n   = n(),
             peer_mean_t = ifelse(g_n>1,
                                  (g_sum - contributing)/pmax(g_n-1,1),
                                  NA_real_)) %>%
      ungroup() %>%
      arrange(player, round_n) %>%
      group_by(player) %>%
      mutate(peer_mean_tm1 = dplyr::lag(peer_mean_t)) %>%
      ungroup() %>%
      dplyr::select(-g_sum, -g_n)
  }
}

# -- Tiny model comparison helper (AIC/BIC/R2/ RMSE) --
compare_models_simple <- function(models, data, y = "contributing"){
  get_r2 <- function(m){
    r2m <- NA_real_; r2c <- NA_real_
    if (requireNamespace("MuMIn", quietly = TRUE)) {
      r <- try(MuMIn::r.squaredGLMM(m), silent = TRUE)
      if (!inherits(r, "try-error")) {
        r2m <- as.numeric(r[1, "R2m"]); r2c <- as.numeric(r[1, "R2c"])
      }
    }
    c(r2m, r2c)
  }
  
  tibble(
    model = names(models),
    AIC   = vapply(models, AIC, numeric(1)),
    BIC   = vapply(models, BIC, numeric(1)),
    R2m   = vapply(models, function(m) get_r2(m)[1], numeric(1)),
    R2c   = vapply(models, function(m) get_r2(m)[2], numeric(1)),
    RMSE  = vapply(models, function(m){
      pr <- as.numeric(predict(m))  # include random intercepts
      sqrt(mean((data[[y]] - pr)^2, na.rm = TRUE))
    }, numeric(1))
  ) %>% arrange(AIC)
}

# -- Norm model comparison: peer norm with friend/adversary moderation (counts fallback) --
norm_fit_comparison_counts <- function(df){
  need <- c("player","round_n","contributing","friends","adversaries")
  stopifnot(all(need %in% names(df)))
  
  dd <- add_peer_means(df) %>%
    filter(!is.na(peer_mean_tm1)) %>%
    mutate(
      # center peer norm, z-score moderators
      pm = as.numeric(scale(peer_mean_tm1, center = TRUE, scale = FALSE)),
      fz = as.numeric(scale(friends)),
      az = as.numeric(scale(adversaries))
    )
  
  mG   <- lmer(contributing ~ pm                 + (1|player), data = dd, REML = FALSE)
  mGF  <- lmer(contributing ~ pm*fz              + (1|player), data = dd, REML = FALSE)
  mGA  <- lmer(contributing ~ pm*az              + (1|player), data = dd, REML = FALSE)
  mGFA <- lmer(contributing ~ pm*fz + pm*az      + (1|player), data = dd, REML = FALSE)
  
  cmp <- compare_models_simple(list(mG=mG, mGF=mGF, mGA=mGA, mGFA=mGFA), data = dd)
  list(models = list(mG=mG, mGF=mGF, mGA=mGA, mGFA=mGFA), compare = cmp, data = dd)
}

# =============================================================================
# Honduras Repeated Public‑Goods Games – Full Game‑Theoretic Extensions Toolkit
# -----------------------------------------------------------------------------
# =============================================================================

# =============================================================================
# 00 CONFIGURATION -------------------------------------------------------------
# =============================================================================
# Edit below to match your dataset.

## Game parameters
B_PARAM             <- 2       # public return per 1 Lempira (b)
KAPPA_PARAM         <- 1       # private cost per 1 Lempira (κ)
N_GROUP_DEFAULT     <- 5       # group size (5-person sessions)
ALPHA               <- 0.5     # altruism exponent (√·)
C_MAX               <- 12      # max contribution
K_GLOBAL            <- 1.0     # scale in norm penalty

## Column mappings for pgg_data
COL_PLAYER          <- "player"
COL_ROUND           <- "round_n"
COL_CONTRIB         <- "contributing"
COL_GROUP           <- "group"
COL_VILLAGE         <- "village_code"
COL_GENDER          <- "gender"
COL_AUTON           <- "finauto_q_2_temp"
COL_FRIENDS         <- "friends"
COL_ADVS            <- "adversaries"
#COL_IH              <- "IH"
COL_AGE             <- "age"
COL_EDU             <- "b0100"
COL_RELIG           <- "b0600"
COL_MARITAL         <- "marital_status"
COL_FI              <- "FI"
COL_ACCESS          <- "access_routes"
#COL_VILL_WEALTH     <- "village_wealth_index_median_w3"

## Optional adjacency matrices (rownames = player IDs)
FRIEND_ADJ          <- NULL
ADVERSARY_ADJ       <- NULL

## Estimation options
MIN_ROUNDS_PER_PLAYER <- 6
PARAM_BOUNDS_DI       <- c(0, 10)
PARAM_BOUNDS_PHI      <- c(0, 10)
NLS_START_DI          <- 0.5
NLS_START_PHI         <- 0.5

## Simulation options
SIM_SEED            <- 123
SIM_STEPS_MORAN     <- 5000
SIM_BURNIN          <- 1000
SIM_REPS            <- 200

# =============================================================================
# 09 HIDDEN MARKOV MODEL (2 states)
# =============================================================================

message("[C] Fitting 2‑state HMM … (this may take a minute)")

set.seed(SIM_SEED)
hmm_df <- pgg_data %>%
  dplyr::select(
    id   = .data[[COL_PLAYER]],
    time = .data[[COL_ROUND]],
    y    = .data[[COL_CONTRIB]]
  )

nt     <- as.numeric(table(hmm_df$id))

hmm_mod <- depmixS4::depmix(
  y ~ 1,
  data    = hmm_df,
  nstates = 2,
  family  = gaussian(),
  ntimes  = nt
)
hmm_fit <- suppressMessages(fit(hmm_mod, verbose = FALSE))
assign("pgg_hmm_fit", hmm_fit, envir = .GlobalEnv)

post <- depmixS4::posterior(hmm_fit)
hmm_df$state <- factor(post$state)
means     <- tapply(hmm_df$y, hmm_df$state, mean)
hi_state  <- names(means)[which.max(means)]

hmm_states <- hmm_df %>%
  group_by(id) %>%
  summarise(
    start    = dplyr::first(state),
    end      = dplyr::last(state),
    hi_share = mean(state == hi_state),
    .groups  = "drop"
  )

# ---------------------------------------------------------------------------
# Structural parameter construction (d_i, phi_i) from observed paths
# ---------------------------------------------------------------------------

slope_ols <- function(y, x){
  ok <- is.finite(y) & is.finite(x)
  if(sum(ok) < 3) return(NA_real_)
  vx <- stats::var(x[ok])
  if(!is.finite(vx) || vx <= 0) return(NA_real_)
  stats::cov(x[ok], y[ok]) / vx
}

build_param_df <- function(df, b = 2, kappa = 1, N = 5, c_max = 12){
  stopifnot(all(c("player","group","round_n","contributing") %in% names(df)))
  
  dt <- df %>%
    dplyr::mutate(player = as.character(player)) %>%
    dplyr::group_by(group, round_n) %>%
    dplyr::mutate(
      gsum = sum(contributing, na.rm = TRUE),
      gN   = dplyr::n(),
      peer_mean_t = dplyr::if_else(gN > 1, (gsum - contributing) / pmax(gN - 1, 1), NA_real_)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(player, round_n) %>%
    dplyr::group_by(player) %>%
    dplyr::mutate(peer_mean_tm1 = dplyr::lag(peer_mean_t)) %>%
    dplyr::ungroup() %>%
    dplyr::select(player, contributing, peer_mean_tm1)
  
  per <- dt %>%
    dplyr::group_by(player) %>%
    dplyr::summarise(
      n_used    = sum(is.finite(contributing)),
      c_bar     = mean(contributing, na.rm = TRUE),
      beta_peer = slope_ols(contributing, peer_mean_tm1),
      rss       = sum((contributing - mean(contributing, na.rm = TRUE))^2, na.rm = TRUE),
      .groups   = "drop"
    )
  
  Delta <- kappa - b / N
  
  per %>%
    dplyr::mutate(
      c_bar   = pmin(c_max, pmax(0, c_bar)),
      phi_hat = pmax(beta_peer, 0),
      d_hat   = ifelse(
        is.finite(Delta) && Delta > 0,
        2 * (Delta + 2 * phi_hat) * sqrt(pmax(c_bar, 0)),
        NA_real_
      ),
      d_i     = d_hat,
      phi_i   = phi_hat
    ) %>%
    dplyr::select(player, d_i, phi_i, d_hat, phi_hat, rss, n_used)
}

param_df <- if(exists("pgg_param_player") && is.data.frame(pgg_param_player)) {
  pgg_param_player %>%
    dplyr::mutate(player = as.character(player)) %>%
    { if("d_i" %in% names(.) && !"d_hat" %in% names(.)) dplyr::mutate(., d_hat = d_i) else . } %>%
    { if("d_hat" %in% names(.) && !"d_i" %in% names(.)) dplyr::mutate(., d_i = d_hat) else . } %>%
    { if("phi_i" %in% names(.) && !"phi_hat" %in% names(.)) dplyr::mutate(., phi_hat = phi_i) else . } %>%
    { if("phi_hat" %in% names(.) && !"phi_i" %in% names(.)) dplyr::mutate(., phi_i = phi_hat) else . }
} else {
  build_param_df(pgg_data, b = B_PARAM, kappa = KAPPA_PARAM, N = N_GROUP_DEFAULT, c_max = C_MAX)
}

hmm_states <- hmm_states %>% dplyr::mutate(id = as.character(id))
param_df   <- param_df   %>% dplyr::mutate(player = as.character(player))

param_df <- dplyr::left_join(param_df, hmm_states, by = c("player" = "id"))

# =============================================================================
#  SAVE KEY OBJECTS TO GLOBAL ENVIRONMENT -------------------------------------
# =============================================================================

assign("pgg_param_player",      param_df,       envir = .GlobalEnv)

## ============================================================================
## Finite-horizon unravel check vs observed decline
## ============================================================================

finite_horizon_check <- function(df, thresh_high=NULL){
  .check_cols(df, c(GAME_GROUP_VAR, ROUND_VAR, CONTRIB_VAR))
  if(is.null(thresh_high)){
    thresh_high <- df %>% filter(!!sym(ROUND_VAR)==1) %>%
      summarise(th = mean(!!sym(CONTRIB_VAR), na.rm=TRUE)) %>% pull(th)
  }
  g_init <- df %>% filter(!!sym(ROUND_VAR)==1) %>%
    group_by(!!sym(GAME_GROUP_VAR)) %>%
    summarise(init_share = mean(!!sym(CONTRIB_VAR)>thresh_high, na.rm=TRUE), .groups="drop")
  
  g_path <- df %>%
    group_by(!!sym(GAME_GROUP_VAR), !!sym(ROUND_VAR)) %>%
    summarise(g_mean = mean(!!sym(CONTRIB_VAR), na.rm=TRUE), .groups="drop") %>%
    left_join(g_init, by=GAME_GROUP_VAR) %>%
    mutate(q = cut(init_share, breaks=quantile(init_share, probs=seq(0,1,0.25), na.rm=TRUE),
                   include.lowest=TRUE, labels=paste0("Q",1:4)))
  
  p <- ggplot(g_path, aes(x=.data[[ROUND_VAR]], y=g_mean, group=.data[[GAME_GROUP_VAR]], colour=q)) +
    geom_line(alpha=0.2) +
    stat_summary(aes(group=q), fun=mean, geom="line", linewidth=1.1) +
    theme_minimal() +
    labs(x="Round", y="Group mean contribution", colour="Init High quartile",
         title="Finite-horizon unravel vs observed paths")
  list(data=g_path, plot=p, thresh_high=thresh_high)
}

## ============================================================================
## Bootstrap the Low-path share of the bifurcation
## ============================================================================

bootstrap_bifurcation <- function(df, B=200, seed=123){
  set.seed(seed)
  .check_cols(df, c(VILLAGE_VAR, PLAYER_VAR, ROUND_VAR, CONTRIB_VAR))
  vil_ids <- unique(df[[VILLAGE_VAR]])
  out <- vector("list", B)
  for(b in seq_len(B)){
    samp_v <- sample(vil_ids, length(vil_ids), replace = TRUE)
    d_b <- purrr::map_dfr(samp_v, function(v) {
      df %>%
        filter(!!sym(VILLAGE_VAR) == v)
    }) %>%
      group_by(player = !!sym(PLAYER_VAR), round = !!sym(ROUND_VAR)) %>%
      summarise(c = mean(!!sym(CONTRIB_VAR), na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = round, values_from = c, names_prefix = "r") %>%
      drop_na()
    mat <- as.matrix(d_b[,-1])
    cl  <- cutree(hclust(dist(scale(mat)), method="ward.D2"), k=2)
    mu  <- rowMeans(mat); dfc <- tibble(player=d_b$player, cl=factor(cl), mu)
    labs <- dfc %>% group_by(cl) %>% summarise(mu=mean(mu), .groups="drop") %>%
      arrange(desc(mu)) %>% pull(cl)
    out[[b]] <- tibble(iter=b, low_share = mean(dfc$cl==labs[2]))
  }
  bind_rows(out)
}

## ============================================================================
## Threshold-sensitivity of logit predictors across cutoffs
## ============================================================================

threshold_sensitivity_grid <- function(df, thresh_grid=c(5,6,7), add_tipping=TRUE){
  .check_cols(df, c(PLAYER_VAR, ROUND_VAR, CONTRIB_VAR, FRIENDS_VAR, AUTONOMY_VAR, GENDER_VAR))
  if(add_tipping){
    tip <- df %>% filter(!!sym(ROUND_VAR)==1) %>%
      summarise(th = mean(!!sym(CONTRIB_VAR), na.rm=TRUE)) %>% pull(th)
    thresh_grid <- unique(c(thresh_grid, tip))
  }
  purrr::map_dfr(thresh_grid, function(th){
    ydat <- df %>% filter(!!sym(ROUND_VAR)==10) %>%
      mutate(high10 = !!sym(CONTRIB_VAR) > th) %>%
      dplyr::select(player = !!sym(PLAYER_VAR), high10)
    base <- df %>% filter(!!sym(ROUND_VAR)==1) %>%
      distinct(player=!!sym(PLAYER_VAR),
               friends=!!sym(FRIENDS_VAR),
               autonomy=!!sym(AUTONOMY_VAR),
               gender=!!sym(GENDER_VAR))
    dat <- left_join(base, ydat, by="player") %>% filter(!is.na(high10))
    broom::tidy(glm(high10 ~ friends + autonomy + gender, data=dat, family=binomial)) %>%
      mutate(thresh = th)
  })
}

## ============================================================================
## Conceptual flow diagram (DiagrammeR)
## ============================================================================

conceptual_flow_diagram <- function(){
  
  DiagrammeR::grViz("
   digraph flow {
     graph [rankdir=LR]; node [shape=box, style=rounded];
     Trait   [label='Observed Traits\\n(friends, autonomy, gender, humility)'];
     Param   [label='Behavioral Params\\n(d_i, phi_i, omega_i)'];
     Dyn     [label='Myopic Update & Imitation\\n(partial-adjust)'];
     Basin   [label='High vs Low Basin\\n(bifurcation)'];
     Village [label='Village Outcomes\\n(mean contrib, stability)'];
     Policy  [label='Policy Levers\\n(subsidy, visibility, seeding)'];
     Trait -> Param -> Dyn -> Basin -> Village -> Policy;
     Policy -> Param [style=dashed, label='affect incentives'];
   }")
}

# --- Run main analytical modules --------------------------------------------

ad1 <- variance_decomp_levels(pgg_data, nested = TRUE)

# Tidy table (percent shares)
ad1$varcomp_clean

# If singular due to ~0 group variance and you want a leaner spec:
# m_trim <- lme4::lmer(contributing ~ 1 + (1|village) + (1|player), data = pgg_data, REML = TRUE)
# as.data.frame(VarCorr(m_trim)) |> mutate(share = vcov / sum(vcov))

ad3

# ----- RUN norm-model comparison --------------------------------------------
rm(list = c("ad4"))        # clean up any half-created object from the error
ad4 <- norm_fit_comparison_counts(pgg_data)
ad4$compare                # AIC/BIC/R2/ RMSE table
summary(ad4$models$mGFA)   # see interactions: pm:fz (friends), pm:az (adversaries)

# Finite-horizon check
ad5 <- finite_horizon_check(pgg_data)

ad5$plot

# Bootstrap bifurcation share and threshold-sensitivity grid
ad6 <- bootstrap_bifurcation(pgg_data); hist(ad6$low_share, breaks=20)
ad7 <- threshold_sensitivity_grid(pgg_data)

# Conceptual flow diagram
conceptual_flow_diagram()

# ---------------------------------------------------------------------------
# Reporting tables and plots
# ---------------------------------------------------------------------------

# Common deps
suppressPackageStartupMessages({
  library(dplyr); library(knitr); library(broom); library(broom.mixed)
})

## Variance components table
ad1_tab <- if ("varcomp_clean" %in% names(ad1)) {
  base <- ad1$varcomp_clean %>%
    transmute(Level = grp, Variance = vcov, `Std.Dev` = sdcor, `Share (%)` = round(share_pct, 1))
  total <- tibble(Level = "Total",
                  Variance = sum(base$Variance),
                  `Std.Dev` = sqrt(sum(base$Variance)),
                  `Share (%)` = 100)
  bind_rows(base, total)
} else {
  ad1$varcomp %>%
    mutate(share_pct = 100 * vcov / sum(vcov)) %>%
    transmute(Level = grp, Variance = vcov, `Std.Dev` = sdcor, `Share (%)` = round(share_pct, 1))
}
kable(ad1_tab, digits = 3, caption = "Variance components and shares")

## Welfare accounting table
kable(ad3 %>% mutate(across(payoff_per_player, round, 2)),
      caption = "Per-player welfare across scenarios")

## Norm-fit comparison (model selection table)
kable(ad4$compare, digits = 3, caption = "Peer vs friend vs adversary norm models")

# (Optional) fixed-effects table for the richest spec if you want it:
# fe <- broom::tidy(ad4$models$mGFA, effects="fixed")
# kable(fe, digits=3, caption="Fixed effects (peer × friends/adversaries)")

suppressPackageStartupMessages({
  library(ggplot2); library(scales); library(DiagrammeR)
})

## Finite-horizon check plot
print(ad5$plot)

## Bootstrap of bifurcation share
ggplot(ad6, aes(low_share)) +
  geom_histogram(bins = 20) +
  labs(x = "Share in Low-path cluster (bootstrap draw)",
       y = "Count",
       title = "Bootstrap distribution of low-path share") +
  theme_minimal()

## Threshold sensitivity (effects vs threshold)
ggplot(ad7, aes(x = thresh, y = estimate, colour = term)) +
  geom_line() +
  facet_wrap(~term, scales = "free_y") +
  labs(x = "Threshold (Lempiras)", y = "Coefficient",
       title = "Sensitivity of end-state regressions to threshold") +
  theme_minimal()

## Conceptual flow diagram (DiagrammeR)
conceptual_flow_diagram()  # renders DiagrammeR graph
