## Produces:
##   – Tipping‑point drift plot (Δc vs c)
##   – Village critical‑mass logistic curve
##   – HMM transition matrix heatmap
##   – Distributions of behavioural slopes (β_peer, β_self)
##   – Structural altruism histogram (d_i)
##   – Early‑warning ROC curve and associated logit coefficients

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(mgcv)
  library(broom); library(pROC); library(lme4)
})

pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

################################################################################
## HONDURAS PUBLIC-GOODS EXTENSIONS: WHAT THIS SCRIPT DOES
## ----------------------------------------------------------------------------- 
## This script layers game-theoretic structure onto the cleaned long panel
## `pgg_data` (10-round public-goods game, about 2,600 players in 134 villages).
##
## MAIN MODULES
## ------------
##   A1  Structural back-out of player-level parameters:
##       - d_hat : altruism weight on own contribution
##
##   A3  Behavioural humility slopes from a mixed model:
##       contribution_it ~ peer_mean_{t-1} + own_{t-1} with random slopes by i.
##
##   A4  Empirical tipping point c* from a smooth drift function
##       Δc = E[c_{t+1}-c_t | c_t], estimated via a GAM in c_t.
##
##   B   Village-level critical mass:
##       P(village ends high) as a logistic function of initial high share,
##       plus friend-density and autonomy rate.
##
##   C   2-state Hidden-Markov model (HMM):
##       latent Low vs High contribution states with (optionally) covariate-
##       dependent transition probabilities (friends, autonomy, gender).
##       Emissions are Gaussian on the (z-scored) contribution scale.
##
## INPUTS
## ------
## * pgg_data          : long panel in player × round format
## * Stage-game params : b_param, kappa_param (here fixed at b = 2, κ = 1)
## * Column mappings   : group, round, contrib, player, village, etc.
##
## MAIN CALL
## ---------
## results <- run_extensions(
##   pgg_data,
##   b_param     = 2,
##   kappa_param = 1
## )
##
## OBJECTS RETURNED
## ----------------
## results$struct_params : player-level d_hat, phi_hat, fit info, predicted cH
## results$behav_hum     : player-level β_peer (imitation) and β_self (habit)
## results$tipping_root  : empirical tipping point c* (Δc = 0)
## results$village_cm    : list(model=glm, crit=critical mass, data=village df)
## results$hmm_fit       : fitted depmixS4 2-state HMM
## results$hmm_post      : per-row posterior state probabilities and Viterbi path
##
################################################################################

pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

################################################################################
#  Honduras Repeated Public‑Goods Games – Extensions (A, B, C)  – FINAL
################################################################################

## ---------------------------------------------------------------------------
## 0 · PACKAGES
## ---------------------------------------------------------------------------
## pkgs lists all packages used by the extension machinery.
## If some are missing, they are installed; then all are loaded.
## If the 'conflicted' package is already loaded, we explicitly prefer
## data.table’s := operator to avoid ambiguity.

pkgs <- c(
  "dplyr","tidyr","purrr","tibble","rlang",
  "data.table",
  "minpack.lm",
  "lme4","broom","broom.mixed",
  "mgcv",
  "depmixS4",
  "ggplot2"
)
new <- setdiff(pkgs, rownames(installed.packages()))
if(length(new)) install.packages(new, dependencies=TRUE)
invisible(lapply(pkgs, library, character.only=TRUE))
if("conflicted" %in% loadedNamespaces())
  conflicts_prefer(data.table::`:=`)

## ---------------------------------------------------------------------------
## 1 · COLUMN MAP
## ---------------------------------------------------------------------------
## String constants so that the rest of the code can be written in terms of
## GAME_GROUP_VAR, CONTRIB_VAR, etc., instead of hard-coding column names.

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
## 2 · HELPERS
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

## ---------------------------------------------------------------------------
## 3 · A1 STRUCTURAL BACK‑OUT
## ---------------------------------------------------------------------------
## estimate_structural():
##   For each player i, fit (d_i, φ_i) so that the first-order condition from
##   the simple structural model is approximately satisfied in each round.
##   This uses nonlinear least squares (minpack.lm::nls.lm) in (d_i, φ_i),
##   subject to non-negativity constraints. It also computes the implied
##   high-state contribution c^H_i = [(0.5 d_i) / (κ - b/N)]^2 truncated
##   to [0, MAX_CONTRIB].

estimate_structural <- function(df,b_param,kappa_param,eps=1e-2,min_r=3){
  N <- infer_N(df)
  df <- add_peer_means(df)
  
  ## Keep only rows with a defined peer mean at t−1.
  keep <- df %>%
    dplyr::select(player=all_of(PLAYER_VAR),
           c=all_of(CONTRIB_VAR),
           norm=peer_mean_tm1) %>%
    filter(!is.na(norm)) %>%
    mutate(c_safe=pmax(c,eps))
  
  ## Constant term in the FOC: b/N − κ
  const <- (b_param/N)-kappa_param
  
  ## Per-player NLS fit of (d_i, φ_i).
  fit_one <- function(pid){
    sub <- keep[keep$player==pid,]
    if(nrow(sub)<min_r)
      return(tibble(player=pid,d_hat=NA,phi_hat=NA,n_obs=nrow(sub),
                    conv_code=NA,rss=NA))
    
    res_fun <- function(par){
      d_i<-par[1]; phi_i<-par[2]
      const + d_i*ALPHA_PARAM*sub$c_safe^(ALPHA_PARAM-1) -
        phi_i*(sub$c - sub$norm)
    }
    
    ## Initial value for φ_i from a simple regression; truncated at 0.
    phi0<-tryCatch(-coef(lm(I(const)~I(sub$c-sub$norm)))[2],
                   error=function(e)0.1)
    
    fit<-tryCatch(
      minpack.lm::nls.lm(c(0.5,max(phi0,0)),fn=res_fun,
                         lower=c(0,0),upper=c(Inf,Inf),
                         control=nls.lm.control(maxiter=100)),
      error=function(e)NULL)
    
    if(is.null(fit))
      return(tibble(player=pid,d_hat=NA,phi_hat=NA,n_obs=nrow(sub),
                    conv_code=1,rss=NA))
    
    tibble(player=pid,
           d_hat=fit$par[1],
           phi_hat=fit$par[2],
           n_obs=nrow(sub),
           conv_code=fit$info,
           rss=sum(res_fun(fit$par)^2))
  }
  
  ## Implied high fixed point c^H under the simple closed-form rule.
  denom<-(kappa_param - b_param/N)
  map_dfr(unique(keep$player),fit_one) %>%
    mutate(cH_pred=ifelse(denom>0,
                          pmin(MAX_CONTRIB,
                               pmax(0,((0.5*d_hat)/denom)^2)),
                          NA))
}

## ---------------------------------------------------------------------------
## 4 · A3 BEHAVIOURAL HUMILITY
## ---------------------------------------------------------------------------
## behav_humility():
##   Mixed model:
##     contributing_it ~ peer_mean_{t-1} + own_{t-1}
##   with random slopes in peer/self by player. Returns per-player estimates
##   β_peer and β_self used as empirical “behavioural humility” indices.

behav_humility <- function(df){
  df<-add_peer_means(df)%>%
    filter(!is.na(peer_mean_tm1),!is.na(self_lag))%>%
    mutate(peer_c=scale(peer_mean_tm1,TRUE,FALSE),
           self_c=scale(self_lag,TRUE,FALSE))
  
  form<-reformulate(c("peer_c","self_c",
                      paste0("(peer_c+self_c|",PLAYER_VAR,")")),
                    response=CONTRIB_VAR)
  
  mod<-lme4::lmer(form,data=df,REML=FALSE)
  
  broom.mixed::tidy(mod,effects="ran_vals")%>%
    filter(term%in%c("peer_c","self_c"))%>%
    dplyr::select(level,term,estimate)%>%
    pivot_wider(names_from=term,values_from=estimate)%>%
    rename(player=level,beta_peer=peer_c,beta_self=self_c)
}

## ---------------------------------------------------------------------------
## 5 · B VILLAGE CRITICAL MASS
## ---------------------------------------------------------------------------
## village_crit_mass():
##   Village-level logistic regression:
##     final_high_v (ends above threshold) ~ init_share_v + friend density + autonomy rate
##   Returns fitted model, estimated “critical share” where P=0.5, and the
##   underlying village-level dataset.

village_crit_mass<-function(df,thresh){
  init<-df%>%
    filter(.data[[ROUND_VAR]]==1)%>%
    mutate(high1=.data[[CONTRIB_VAR]]>thresh)%>%
    group_by(.data[[VILLAGE_VAR]])%>%
    summarise(init_share=mean(high1),
              fr_dens=mean(.data[[FRIENDS_VAR]],na.rm=TRUE),
              auton_rate=mean(.data[[AUTONOMY_VAR]],na.rm=TRUE),
              .groups="drop")
  
  final<-df%>%
    filter(.data[[ROUND_VAR]]==10)%>%
    group_by(.data[[VILLAGE_VAR]])%>%
    summarise(final_mean=mean(.data[[CONTRIB_VAR]],na.rm=TRUE),
              .groups="drop")
  
  vil<-left_join(init,final,by=VILLAGE_VAR)%>%
    filter(!is.na(final_mean))%>%
    mutate(final_high=final_mean>thresh)
  
  mod<-glm(final_high~init_share+fr_dens+auton_rate,
           data=vil,family=binomial)
  
  ## Critical mass s^crit solves P(final_high | init_share = s) = 0.5.
  root<-tryCatch(
    uniroot(function(x){
      nd<-with(vil,
               data.frame(init_share=x,
                          fr_dens=mean(fr_dens,na.rm=TRUE),
                          auton_rate=mean(auton_rate,na.rm=TRUE)))
      predict(mod,nd,type="response")-0.5
    },c(0,1))$root,
    error=function(e)NA_real_)
  
  list(model=mod,crit=root,data=vil)
}

## ---------------------------------------------------------------------------
## 6 · C HIDDEN‑MARKOV MODEL
## ---------------------------------------------------------------------------
## prep_hmm():
##   Reshape panel into per-player time series with covariates.
## fit_hmm():
##   Fit a 2-state HMM with Gaussian emissions (on z-scored contributions).
##   Transitions depend on friends, autonomy, gender; if that fails,
##   fall back to a homogeneous transition model.

prep_hmm<-function(df){
  df%>%
    dplyr::select(player=all_of(PLAYER_VAR),
           round=all_of(ROUND_VAR),
           contrib=all_of(CONTRIB_VAR),
           friends=all_of(FRIENDS_VAR),
           autonomy=all_of(AUTONOMY_VAR),
           gender=all_of(GENDER_VAR))%>%
    arrange(player,round)%>%
    group_by(player)%>%
    mutate(ts=row_number())%>%
    ungroup()
}

fit_hmm<-function(df){
  # clean NAs: we require complete sequences in all covariates
  covs<-c("contrib","friends","autonomy","gender")
  df2<-df%>%
    group_by(player)%>%
    filter(!any(is.na(across(all_of(covs)))))%>%
    ungroup()
  if(nrow(df2)==0) stop("No complete sequences after dropping NAs.")
  
  df2$contrib_z<-scale(df2$contrib)
  nts<-as.numeric(table(df2$player))
  
  make_mod<-function(trans){
    depmixS4::depmix(
      response=list(contrib_z~1),
      data=df2,
      nstates=2,
      family=list(gaussian()),
      transition=as.formula(trans),
      ntimes=nts)
  }
  
  ## Try covariate-dependent transitions; fall back to homogeneous if needed.
  mod<-tryCatch(make_mod("~ friends+autonomy+gender"),
                error=function(e)NULL)
  if(is.null(mod)){
    message("Falling back to homogeneous transitions.")
    mod<-make_mod("~1")
  }
  
  set.seed(123)
  fit<-depmixS4::fit(mod,verbose=FALSE)
  list(model=fit, df=df2)
}

## ---------------------------------------------------------------------------
## 7 · MASTER DRIVER
## ---------------------------------------------------------------------------
## run_extensions():
##   One-stop function that:
##     1) checks columns and infers group size N
##     2) adds peer-means/lagged variables
##     3) runs structural back-out 
##     4) runs behavioural humility mixed model 
##     5) estimates empirical tipping point c* 
##     6) runs village critical mass logistic 
##     7) fits the HMM and builds posterior state series 
##   Returns all objects in a named list.

run_extensions<-function(panel,b_param,kappa_param){
  message("\n[Extensions] Starting …")
  check_cols(panel,c(GAME_GROUP_VAR,CONTRIB_VAR,ROUND_VAR,PLAYER_VAR))
  N<-infer_N(panel); message("Detected group size N=",N)
  
  panel<-add_peer_means(panel)
  
  message("A1 · structural back‑out …")
  pars<-estimate_structural(panel,b_param,kappa_param)
  
  message("A3 · behavioural humility …")
  hum<-behav_humility(panel)
  
  message("A4 · empirical tipping …")
  tip_df<-panel%>%
    mutate(c=.data[[CONTRIB_VAR]])%>%
    arrange(.data[[PLAYER_VAR]],.data[[ROUND_VAR]])%>%
    group_by(.data[[PLAYER_VAR]])%>%
    mutate(c_lead=lead(c),delta=c_lead-c)%>%
    ungroup()%>%
    filter(!is.na(delta))
  tip_mod<-mgcv::gam(delta~s(c,k=5),data=tip_df)
  
  ## Empirical tipping root c* where predicted drift crosses zero.
  root<-tryCatch(
    uniroot(function(x)
      predict(tip_mod,data.frame(c=x),type="response"),
      c(0,MAX_CONTRIB))$root,
    error=function(e)mean(panel[[CONTRIB_VAR]],na.rm=TRUE))
  message("  Tipping root ≈ ",round(root,2)," Lps")
  
  message("B · village critical mass …")
  vil<-village_crit_mass(panel,root)
  
  message("C · hidden‑Markov model …")
  hmm_full<-prep_hmm(panel)
  hmm_res <- fit_hmm(hmm_full)
  hmm_fit<-hmm_res$model
  hmm_df <- hmm_res$df
  post <- depmixS4::posterior(hmm_fit)
  hmm_post<-cbind(hmm_df,post)
  
  message("[Extensions] Done.\n")
  list(
    struct_params=pars,
    behav_hum    =hum,
    tipping_root =root,
    village_cm   =vil,
    hmm_fit      =hmm_fit,
    hmm_post     =hmm_post
  )
}

## Run extensions once (if not already in memory) and store tipping point.
if (!exists("results")) {
  results <- run_extensions(pgg_data, b_param = 2, kappa_param = 1)
}
pgg_tipping <- results$tipping_root

# ---------- Helper: peer means/lag for later figures ----------
add_peer_means <- function(df){
  df %>%
    group_by(group, round_n) %>%
    mutate(g_sum = sum(contributing), g_n = n(),
           peer_mean_t = ifelse(g_n>1, (g_sum - contributing)/pmax(g_n-1,1), NA_real_)) %>%
    ungroup() %>%
    arrange(player, round_n) %>%
    group_by(player) %>%
    mutate(self_lag = dplyr::lag(contributing),
           peer_mean_tm1 = dplyr::lag(peer_mean_t)) %>%
    ungroup() %>%
    dplyr::select(-g_sum, -g_n)
}

# ================================ Fig 2 ======================================
# Drift function Δc vs c with GAM fit and empirical tipping point c*
# Used for the tipping-point figure: E[Δc | c] and zero-crossing at c*.

tip_df <- pgg_data %>%
  arrange(player, round_n) %>%
  group_by(player) %>%
  mutate(c = contributing, c_lead = dplyr::lead(c), delta = c_lead - c) %>%
  ungroup() %>%
  filter(!is.na(delta))

gam_fit <- mgcv::gam(delta ~ s(c, k = 5), data = tip_df)
grid <- data.frame(c = seq(0, 12, length.out = 400))
pr   <- predict(gam_fit, newdata = grid, se.fit = TRUE)
df_g <- transform(grid, fit = pr$fit, lo = pr$fit - 1.96*pr$se.fit, hi = pr$fit + 1.96*pr$se.fit)

# Helper to compute approximate zero-crossing of a smooth curve.
zero_cross <- function(x, y){
  idx <- which(diff(sign(y)) != 0)
  if(!length(idx)) return(NA_real_)
  i <- idx[1]
  x[i] - y[i]*(x[i+1]-x[i])/(y[i+1]-y[i])
}

c_star    <- pgg_tipping
c_star_lo <- zero_cross(df_g$c, df_g$hi)  # lower bound from upper band
c_star_hi <- zero_cross(df_g$c, df_g$lo)  # upper bound from lower band

cat(sprintf("\n[Fig 2] Empirical tipping c* ≈ %.2f L (approx CI: [%.2f, %.2f])\n",
            c_star, c_star_lo, c_star_hi))

fig2_tipping <- ggplot(df_g, aes(c, fit)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(size = 1) +
  geom_vline(xintercept = c_star, linetype = 2) +
  labs(title = "Drift Δc vs c with empirical tipping point c*",
       x = "Current contribution c", y = "Predicted drift E[Δc | c]") +
  theme_minimal()
# Interpretation: for c below c*, drift is positive (tends upward); above c*,
# drift is negative (tends downward).

# ================================ Fig 4 ======================================
# Village critical mass logistic:
# P(village ends above tipping point) as a function of initial high share.

vil_init <- pgg_data %>%
  filter(round_n == 1) %>%
  mutate(high1 = contributing > c_star) %>%
  group_by(village_code) %>% summarise(init_share = mean(high1), .groups="drop")

vil_final <- pgg_data %>%
  filter(round_n == 10) %>%
  group_by(village_code) %>%
  summarise(final_mean = mean(contributing), .groups="drop") %>%
  mutate(final_high = final_mean > c_star)

vil_df <- left_join(vil_init, vil_final, by="village_code")
logit_m <- glm(final_high ~ init_share, data = vil_df, family = binomial)

nd <- data.frame(init_share = seq(0, 1, length.out = 501))
nd$prob <- predict(logit_m, nd, type = "response")

# 50–50 crossing and its delta-method SE
co <- coef(logit_m)               # intercept and slope
V  <- vcov(logit_m)
scrit <- -co[1]/co[2]
g <- c(-1/co[2], co[1]/(co[2]^2)) # gradient wrt (intercept,slope)
scrit_se <- sqrt(t(g) %*% V %*% g)
scrit_lo <- scrit - 1.96*scrit_se
scrit_hi <- scrit + 1.96*scrit_se

cat(sprintf("[Fig 4] Critical mass s^crit ≈ %.2f (95%% CI [%.2f, %.2f])\n", scrit, scrit_lo, scrit_hi))
cat("Logit coefficient for init_share:\n"); print(summary(logit_m)$coef["init_share", , drop=FALSE])

fig4_cm <- ggplot(vil_df, aes(init_share, as.numeric(final_high))) +
  geom_jitter(height = .05, width = 0, alpha = .5) +
  geom_line(data = nd, aes(init_share, prob), linewidth = 1) +
  geom_vline(xintercept = scrit, linetype = 2) +
  scale_y_continuous("P(village finishes high)", limits = c(0,1)) +
  labs(title = "Village critical mass",
       x = "Initial share high (round 1)") +
  theme_minimal()
# Interpretation: villages starting above s^crit are much more likely
# to end in the high-contribution basin.

# ================================ Fig 3 ======================================
# HMM diagnostics: emission means by state, empirical transition matrix,
# and share of observations in the high-contribution state.

hmm_post <- results$hmm_post

# If only z-scores are present, use them as 'contrib' for summaries.
if (!"contrib" %in% names(hmm_post) && "contrib_z" %in% names(hmm_post)) {
  hmm_post$contrib <- hmm_post$contrib_z
}

emit_means <- hmm_post %>%
  group_by(state) %>%
  summarise(emission_mean = mean(contrib, na.rm = TRUE), .groups = "drop") %>%
  arrange(emission_mean)
print(emit_means)

# Empirical transition probabilities from Viterbi path.
trans_df <- hmm_post %>%
  arrange(player, round) %>%
  group_by(player) %>%
  mutate(state_next = dplyr::lead(state)) %>%
  ungroup() %>%
  filter(!is.na(state), !is.na(state_next)) %>%
  count(state, state_next, name = "n") %>%
  group_by(state) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

# Transition matrix and “stickiness” (diagonal probabilities).
Pmat <- trans_df %>%
  dplyr::select(state, state_next, p) %>%
  tidyr::pivot_wider(names_from = state_next, values_from = p, values_fill = 0) %>%
  tibble::column_to_rownames("state") %>% as.matrix()
cat("\n[Fig 3] Transition matrix (rows=current, cols=next):\n"); print(round(Pmat, 3))
cat(sprintf("Diagonal stickiness: %.3f and %.3f\n", Pmat[1,1], Pmat[2,2]))

# Share of rounds in the higher-mean state.
hi_state <- emit_means$state[which.max(emit_means$emission_mean)]
share_hi <- hmm_post %>%
  mutate(is_hi = as.integer(state == hi_state)) %>%
  summarise(share_high = mean(is_hi, na.rm = TRUE)) %>% pull(share_high)
cat(sprintf("Share of observations in High state: %.3f\n", share_hi))

fig3_hmm <- ggplot(trans_df %>% mutate(state = paste0("From ", state),
                                       state_next = paste0("To ", state_next)),
                   aes(state_next, state, fill = p)) +
  geom_tile() + geom_text(aes(label = sprintf("%.2f", p))) +
  scale_fill_continuous(name = "Probability") +
  labs(title = "Empirical transition matrix (Viterbi states)") +
  theme_minimal() + theme(axis.title = element_blank())
# Interpretation: diagonals show persistence of each state; off-diagonals show
# Low→High and High→Low switching rates.

# ================================ Fig 5 ======================================
# Behavioural humility: empirical distributions of β_peer and β_self.

behav <- results$behav_hum %>% rename(beta_peer = beta_peer, beta_self = beta_self)

bp <- stats::quantile(behav$beta_peer, c(.1,.5,.9), na.rm=TRUE)
bs <- stats::quantile(behav$beta_self, c(.1,.5,.9), na.rm=TRUE)
share_anti  <- mean(behav$beta_peer < 0, na.rm=TRUE)
share_strong<- mean(behav$beta_peer > bp["90%"], na.rm=TRUE)

cat(sprintf("\n[Fig 5] beta_peer: median=%.3f (P10=%.3f, P90=%.3f); anti-imitators=%.1f%%; strong copiers=%.1f%%\n",
            bp["50%"], bp["10%"], bp["90%"], 100*share_anti, 100*share_strong))
cat(sprintf("[Fig 5] beta_self: median=%.3f (P10=%.3f, P90=%.3f)\n",
            bs["50%"], bs["10%"], bs["90%"]))

fig5_peer <- ggplot(behav, aes(beta_peer)) +
  geom_histogram(bins = 40, alpha = .85) +
  labs(title = expression("Distribution of " * beta[peer]),
       x = expression(beta[peer]), y = "Count") +
  theme_minimal()

fig5_self <- ggplot(behav, aes(beta_self)) +
  geom_histogram(bins = 40, alpha = .85) +
  labs(title = expression("Distribution of " * beta[self]),
       x = expression(beta[self]), y = "Count") +
  theme_minimal()

fig5 <- fig5_peer / fig5_self
# Interpretation: β_self tends to be larger in magnitude; β_peer has substantial
# mass near zero with both “anti-imitators” (β_peer < 0) and strong copiers.

# ================================ Fig 6 ======================================
# Structural altruism d_i: histogram and summary shares for
#   (i) players whose c^H hits the contribution cap
#   (ii) players with near-zero φ_i.

struct <- results$struct_params
share_at_cap <- mean(struct$cH_pred >= 12 - 1e-9, na.rm=TRUE)
phi_near0    <- mean(struct$phi_hat <= 0.1, na.rm=TRUE)  # threshold for “near 0”

cat(sprintf("\n[Fig 6] Share with c^H at contribution cap (12 L): %.1f%%\n", 100*share_at_cap))
cat(sprintf("[Fig 6] Share with φ_i ≤ 0.1 (near-zero norm pull): %.1f%%\n", 100*phi_near0))
cat("[Fig 6] d_i summary:\n"); print(summary(struct$d_hat))

fig6_d <- ggplot(struct, aes(d_hat)) +
  geom_histogram(bins = 40, alpha = .85) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.6,
           label = sprintf("c^H at cap: %.0f%%\nφ near 0: %.0f%%", 100*share_at_cap, 100*phi_near0)) +
  labs(title = expression("Structural altruism " * d[i] * " (histogram)"),
       x = expression(d[i]), y = "Count") +
  theme_minimal()

# ================================ Table 2 / ROC ===============================
# Early-warning prediction of ending in the high state:
#   – Logistic regression using early-round behaviour and traits
#   – Train/test split with out-of-sample ROC and AUC
#   – Coefficients exported as odds ratios with confidence intervals.

set.seed(123)
early_rounds <- 1:3; final_round <- 10; thresh_high <- c_star

## Early-round summary features at the player level.
feat_early <- pgg_data %>%
  filter(round_n %in% early_rounds) %>%
  group_by(player) %>%
  summarise(
    c_mean_e = mean(contributing, na.rm = TRUE),
    c_sd_e   = sd  (contributing, na.rm = TRUE),
    c_r1     = contributing[round_n == min(early_rounds)][1],
    c_r3     = contributing[round_n == max(early_rounds)][1],
    c_slope  = (c_r3 - c_r1)/(max(early_rounds)-min(early_rounds)),
    .groups = "drop"
  )

## Baseline covariates from round 1 and final high/low outcome.
base1 <- pgg_data %>%
  filter(round_n == 1) %>%
  distinct(player, village_code, friends, finauto_q_2_temp, gender, age)

ydat <- pgg_data %>%
  filter(round_n == final_round) %>%
  transmute(player, final_high = contributing > thresh_high)

## Assemble analysis dataset and split into training/testing folds.
dat <- feat_early %>%
  left_join(base1, by = "player") %>%
  left_join(ydat,  by = "player") %>%
  filter(!is.na(final_high))

idx   <- sample.int(nrow(dat), floor(0.7*nrow(dat)))
train <- dat[idx,]; test <- dat[-idx,]

ew_glm <- glm(final_high ~ c_mean_e + c_sd_e + c_slope + friends + finauto_q_2_temp + gender + age ,
              data = train, family = binomial)

test$pred <- predict(ew_glm, newdata = test, type = "response")

roc_obj <- pROC::roc(test$final_high, test$pred)
auc_val <- as.numeric(pROC::auc(roc_obj))
cat(sprintf("\n[Table 2] Early-warning AUC (holdout): %.3f\n", auc_val))

# Optimal classification threshold by Youden’s J statistic (for summary Sens/Spec).
opt <- pROC::coords(roc_obj, "best", best.method="youden",
                    ret=c("threshold","sensitivity","specificity"))
cat(sprintf("[Table 2] Optimal threshold=%.3f, Sens=%.2f, Spec=%.2f\n",
            opt["threshold"], opt["sensitivity"], opt["specificity"]))

table2_coefs <- broom::tidy(ew_glm, conf.int = TRUE, exponentiate = TRUE) %>%
  mutate(term = dplyr::recode(term,
                              "(Intercept)" = "Intercept",
                              "c_mean_e"    = "Early mean",
                              "c_sd_e"      = "Early SD",
                              "c_slope"     = "Early slope",
                              "friends"     = "Friends",
                              "finauto_q_2_temp" = "Financial autonomy",
                              "gender"      = "Gender (1=male)",
                              "age"         = "Age"))
# "IH"          = "Intellectual humility"))
print(table2_coefs)

## ROC curve for the early-warning classifier.
figROC <- ggplot(data.frame(tpr = roc_obj$sensitivities,
                            fpr = 1 - roc_obj$specificities),
                 aes(fpr, tpr)) +
  geom_line() + geom_abline(linetype=2) +
  labs(title = sprintf("ROC (AUC = %.3f)", auc_val),
       x = "False positive rate", y = "True positive rate") +
  theme_minimal()
# Interpretation: plots the trade-off between Sensitivity and Specificity;
# AUC summarizes overall early-warning performance.

# ============================== SAVE/PRINT (optional) =========================

print(fig2_tipping); print(fig3_hmm); print(fig4_cm); print(fig5); print(fig6_d); print(figROC)

ggsave("fig2_tipping.pdf", fig2_tipping, width=6.5, height=4)
ggsave("fig3_hmm.pdf",     fig3_hmm,     width=5.5, height=4.5)
ggsave("fig4_cm.pdf",      fig4_cm,      width=6.5, height=4)
ggsave("fig5_behav.pdf",   fig5,         width=6.5, height=6.5)
ggsave("fig6_d.pdf",       fig6_d,       width=6.5, height=4)
ggsave("figROC.pdf",       figROC,       width=5.5, height=4.5)
write.csv(table2_coefs, "table2_early_warning_coeffs.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Additional setup for alternative structural estimation (per-player FOC fits).
# -----------------------------------------------------------------------------
# 1) Build per-person panel with leave-one-out peer mean in later rounds (2–10).
pan <- pgg_data |>
  arrange(group, player, round_n) |>
  group_by(group, round_n) |>
  mutate(G = n(),
         peer_mean_t = ifelse(G > 1, (sum(contributing) - contributing) / (G - 1), NA_real_)) |>
  ungroup() |>
  group_by(player) |>
  mutate(peer_mean_tm1 = lag(peer_mean_t),
         c_tm1         = lag(contributing)) |>
  ungroup() |>
  filter(round_n >= 2, is.finite(peer_mean_tm1), is.finite(contributing))

# 2) Guard against log/FOC problems from zeros in contributions.
pan <- pan |> mutate(c_eff = pmax(contributing, 1e-3))

# 3) For each player i, one can, if desired, solve the minimum-distance FOC for
#    (d_i, phi_i) using (c_eff, peer_mean_tm1); residuals can be winsorized
#    before summarising. The implementation of that step is not included here,
#    since the structural back-out above already provides d_hat, phi_hat.
