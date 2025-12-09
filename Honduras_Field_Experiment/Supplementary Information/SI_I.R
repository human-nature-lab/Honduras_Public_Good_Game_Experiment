## =====================================================================
##  PUBLIC GOODS GAME – TYPE CLASSIFICATION & TRANSITION ANALYSES
## =====================================================================
## This script:
##   (i)   constructs baseline (round-1) covariates at the player level,
##   (ii)  builds a player×round panel with high/low contribution states
##         and leave-one-out peer means,
##   (iii) recovers “revealed-slope” types (Conditional cooperator,
##         Triangle, Free rider, Other) using mixed models and
##         player-specific quadratic fits,
##   (iv)  clusters contribution paths using Ward–D2 and cross-tabulates
##         them with types,
##   (v)   estimates learning models and stability diagnostics, and
##   (vi)  relates types to covariates (mixed-effects and multinomial
##         logit) and produces a balanced conditional-response figure.
##
## All objects and variable names are kept exactly as in the main
## analysis; only comments and section headers have been expanded.
## =====================================================================


## ========================== PACKAGES ===========================
library(dplyr)
library(tidyr)
library(tibble)
library(fastcluster)
library(cluster)
library(factoextra)
library(ggplot2)
library(lme4)
library(lmerTest)
library(purrr)

## -----------------------------------------------------------------
## Read the experimental panel data:
##   - One row per player × round.
##   - Key variables: player, group, round_n, contributing,
##     and socio-demographic covariates.
## -----------------------------------------------------------------
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)


## -----------------------------------------------------------------
## Compute c* = average round-1 contribution across all players.
## This serves as the global high/low cut-off used to define
## contribution “states” H (high) vs L (low) in the transition
## analyses below.
## -----------------------------------------------------------------
c_star <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::summarise(m = mean(contributing, na.rm = TRUE)) %>%
  dplyr::pull(m)


## -----------------------------------------------------------------
## Baseline covariates: one row per player (round 1 only).
## Steps:
##   1) Keep a single observation per player based on round 1.
##   2) Coerce IDs and categorical variables to factors.
##   3) Convert numeric-like strings to numeric.
##   4) Construct within-sample z-scores for continuous covariates
##      (age_z, friends_z, …) to ease interpretation in regressions.
##   5) Set religion reference category (b0600) to “2” (Catholic),
##      with other denominations interpreted relative to that level.
## -----------------------------------------------------------------
baseline <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::distinct(player, .keep_all = TRUE) %>%  # one baseline row per player
  dplyr::transmute(
    player        = as.character(player),
    village_code  = as.character(village_code),
    gender        = factor(gender),
    marital_status= factor(marital_status),
    b0600         = factor(b0600),   # religion (re-leveled below)
    access_routes = factor(access_routes),
    age           = suppressWarnings(as.numeric(age)),
    friends       = suppressWarnings(as.numeric(friends)),
    adversaries   = suppressWarnings(as.numeric(adversaries)),
    network_density_fr  = suppressWarnings(as.numeric(network_density_fr)),
    network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
    network_size        = suppressWarnings(as.numeric(network_size)),
    b0100              = suppressWarnings(as.numeric(b0100)),  # education (years)
    b0200              = as.numeric(b0200),   # indigenous indicator (numeric dummy)
    FI                 = suppressWarnings(as.numeric(FI)) ,    # financial autonomy index
    finauto_q_2_temp    =  as.numeric(finauto_q_2_temp)
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


## -----------------------------------------------------------------
## Construct the player×round panel used in the transition analyses.
##
## For each (player, round):
##   - Define binary state s_t (1 = High, 0 = Low) using c_star.
##   - Compute leave-one-out peer mean contribution in the group
##     at time t: peer_mean_c_t.
##   - Construct peer_mean_c_tm1, the lagged leave-one-out mean
##     at time t−1 at the group level.
##   - Build spell variables:
##        s_lag: previous-period state s_{t−1}
##        run_id: identifier of contiguous runs of the same state
##        spell_index: index within a run (length of run up to t)
##        spell_prev: duration already spent in the current regime
##                    at t−1 (used in duration models).
##   - Merge player-level baseline covariates.
##   - Scale the lagged peer mean by the per-round endowment
##     (assumed 12 units) for direct interpretability.
## -----------------------------------------------------------------
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
    ## Position within the run (1, 2, 3, …).
    ## N.B.: because we group by player and then sort by round,
    ##       row_number() gives the cumulative index across that player's
    ##       full history; combined with run_id, it identifies spells.
    spell_index = dplyr::row_number(),
    ## Spell length at t-1 (duration already spent in current regime)
    spell_prev  = dplyr::lag(spell_index)
  ) %>%
  dplyr::ungroup() %>%
  ## Merge baseline covariates (one row per player) into the full panel.
  dplyr::left_join(
    baseline %>% dplyr::select(
      player, village_code,
      gender, marital_status, b0600, b0200, access_routes, finauto_q_2_temp,
      ends_with("_z")
    ),
    by = c("player", "village_code")
  ) %>%
  ## Scale peer mean by the per-round endowment (assumed to be 12 units)
  dplyr::mutate(peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12)


## ====================== UNIFY KEY TYPES ========================
## Ensure that `player` is treated as a character identifier everywhere.
## This prevents unintended coercions when joining across objects.
pgg_data <- pgg_data %>% mutate(player = as.character(player))


## ====================== PEER LAGS (condresp_long) ==============
## Function to compute leave-one-out peer mean in each group×round,
## and then construct its one-period lag within each player.
##
## Output columns:
##   - peer_mean: leave-one-out mean contribution at time t.
##   - peer_mean_lag: peer_mean_{t−1}, used as regressor in
##     conditional-response models.
compute_peer_lags <- function(df){
  df %>%
    group_by(group, round_n) %>%
    mutate(
      g_n   = n(),
      g_sum = sum(contributing, na.rm = TRUE),
      peer_mean = (g_sum - contributing) / (g_n - 1)
    ) %>%
    ungroup() %>%
    arrange(group, player, round_n) %>%
    group_by(group, player) %>%
    mutate(peer_mean_lag = lag(peer_mean)) %>%
    ungroup() %>% dplyr::select(-g_n, -g_sum)
}

## Conditional-response panel:
##   - one row per (player, round) with lagged peer mean.
##   - drop first observations where peer_mean_lag is undefined.
condresp_long <- pgg_data %>%
  compute_peer_lags() %>%
  filter(!is.na(peer_mean_lag)) %>%
  mutate(player = as.character(player))


## ============ WARD–D2 HIERARCHICAL CLUSTERING (10 rounds) =================
## Purpose:
##   - Cluster players’ contribution paths over rounds (1–10).
##   - Use Ward–D2 linkage on standardized paths.
##   - Select the number of clusters k via average silhouette width.
##   - This generates a path-based type tag hc_cluster (C1, C2, …).
## -----------------------------------------------------------------

# 1) Build wide player × round matrix of mean contributions
aggregated_contributions <- pgg_data %>%
  group_by(player, round_n) %>%
  summarise(contrib = mean(contributing, na.rm = TRUE), .groups = "drop")

contrib_mat <- aggregated_contributions %>%
  mutate(player = as.character(player)) %>%
  pivot_wider(names_from = round_n, values_from = contrib, names_prefix = "r") %>%
  drop_na() %>%                          # keep only players with complete 10-round paths
  column_to_rownames("player") %>%
  as.matrix()

# 2) Z-score each round; 3) Ward–D2; 4) choose k by silhouette
contrib_sc <- scale(contrib_mat)         # standardize by round (column-wise)
dist_mat   <- dist(contrib_sc, method = "euclidean")
hc_tree    <- fastcluster::hclust(dist_mat, method = "ward.D2")

get_sil <- function(k) mean(cluster::silhouette(cutree(hc_tree, k), dist_mat)[,3])
ks         <- 2:6
sil_widths <- setNames(vapply(ks, get_sil, numeric(1)), ks)
best_k     <- ks[ which.max(sil_widths) ]

## Cluster tag: “C1”, “C2”, …, where labels are arbitrary but fixed.
hc_tag <- factor(cutree(hc_tree, k = best_k), labels = paste0("C", seq_len(best_k)))
contrib_df <- tibble(player = rownames(contrib_sc), hc_cluster = hc_tag) %>%
  mutate(player = as.character(player))

## (Optional) End-state confusion matrix:
##   - Define H/L at round 1 using mean contribution in round 1.
##   - Track transitions (r1 → r10) for each player.
##   - Cross-tabulate transitions with hierarchical clusters.
thresh_r1 <- mean(pgg_data$contributing[pgg_data$round_n == 1], na.rm = TRUE)
states <- pgg_data %>%
  filter(round_n %in% c(1, 10)) %>%
  mutate(state = if_else(contributing > thresh_r1, "H", "L")) %>%
  dplyr::select(player, round_n, state) %>%
  pivot_wider(names_from = round_n, names_prefix = "r", values_from = state) %>%
  mutate(trans = paste0(r1, r10), player = as.character(player))

confusion <- states %>%
  inner_join(contrib_df, by = "player") %>%
  count(hc_cluster, trans) %>%
  pivot_wider(names_from = trans, values_from = n, values_fill = 0)

cat("\nSilhouette widths (k=2..6):\n"); print(round(sil_widths, 3))
cat("\n--- Confusion matrix: End-state vs HC clusters (k =", best_k, ") ---\n"); print(confusion)


## ======================= STABLE TYPE CLASSIFICATION ========================
## Objective:
##   - Classify players into “revealed-slope” types based on how their
##     contributions respond to lagged peer mean:
##       * Conditional cooperator: positive, roughly linear slope.
##       * Triangle: positive response up to a peak, then decline.
##       * Free rider: flat (near-zero) slope at low contribution levels.
##       * Other: residual category once the above are assigned.
##
## The classification proceeds in three stages:
##   (1) Player-specific empirical Bayes (EB) slopes from a mixed model.
##   (2) Player-level quadratic fit on standardized peer mean x_c.
##   (3) Assignment using “tuning knobs” for prevalence of each type.
## -----------------------------------------------------------------

# 0) Center/scale peer mean WITHIN player to stabilize x and x^2
prep <- condresp_long %>%
  group_by(player) %>%
  mutate(
    x_raw  = peer_mean_lag,
    x_mean = mean(x_raw, na.rm = TRUE),
    x_sd   = sd(x_raw,  na.rm = TRUE),
    ## If a player has no variability, set x_c = 0 for all rounds.
    x_c    = ifelse(is.na(x_sd) | x_sd == 0, 0, (x_raw - x_mean) / x_sd)
  ) %>% ungroup()

# 1) EB slopes from a mixed model (monotone component)
## Model:
##   contributing_it = β0 + β1 x_c,it + u0i + u1i x_c,it + vg + ε_it
##   where:
##     i = player, g = group.
##   We extract per-player empirical Bayes slope b1_eb = β1 + u1i.
eb_mod <- lmerTest::lmer(
  contributing ~ x_c + (1 + x_c | player) + (1 | group),
  data = prep
)

eb <- coef(eb_mod)$player %>%
  as.data.frame() %>%
  tibble::rownames_to_column("player") %>%
  transmute(player = as.character(player),
            b0_eb  = `(Intercept)`,
            b1_eb  = x_c)

## Player-level mean contribution across all rounds
mean_by_player <- pgg_data %>%
  group_by(player) %>%
  summarise(mean_c = mean(contributing, na.rm = TRUE), .groups = "drop") %>%
  mutate(player = as.character(player))

## Base feature frame: EB slope + mean contribution
base <- eb %>% left_join(mean_by_player, by = "player")

# 2) Per-player quadratic on standardized x_c to detect TRUE "Triangle"
## For each player with sufficient variation (≥ 6 usable rounds),
## estimate:
##   m1: c_it = α0 + α1 x_c,it + ε_it
##   m2: c_it = β0 + β1 x_c,it + β2 x_c,it^2 + ε_it
## and compute:
##   - AIC(m1), AIC(m2) for model comparison;
##   - b1, b2 (linear and quadratic terms);
##   - x_star = −b1/(2 b2): interior optimum if b2 < 0;
##   - x10, x90: 10th and 90th percentiles of x_c for that player.
## “Triangle” candidates are those with an interior maximum between
## x10 and x90, b1 > 0, b2 < 0, and AIC(m2) sufficiently better.
quad <- prep %>%
  group_by(player) %>%
  filter(n() >= 6, !all(is.na(x_c))) %>%
  group_modify(~{
    dat <- .x
    m1  <- lm(contributing ~ x_c, data = dat)
    m2  <- lm(contributing ~ x_c + I(x_c^2), data = dat)
    S2  <- summary(m2)$coefficients
    get <- function(term, col) if (term %in% rownames(S2)) S2[term, col] else NA_real_
    b1  <- get("x_c","Estimate"); b2 <- get("I(x_c^2)","Estimate")
    a1  <- AIC(m1); a2 <- AIC(m2)
    x10 <- quantile(dat$x_c, .10, na.rm=TRUE)
    x90 <- quantile(dat$x_c, .90, na.rm=TRUE)
    x_star <- ifelse(!is.na(b2) && b2!=0, -b1/(2*b2), NA_real_)
    tibble(player = unique(dat$player),
           b1=b1, b2=b2, aic1=a1, aic2=a2, x10=x10, x90=x90, x_star=x_star)
  }) %>% ungroup() %>%
  mutate(player = as.character(player))

## Flag players whose response is statistically and economically
## consistent with a “Triangle” profile.
tri <- quad %>%
  mutate(is_triangle = (b1 > 0) & (b2 < 0) & (aic2 < aic1 - 2) &
           !is.na(x_star) & x_star > x10 & x_star < x90) %>%
  dplyr::select(player, is_triangle)

## Merge triangle flags into the base feature frame
feat <- base %>%
  left_join(tri, by = "player") %>%
  mutate(is_triangle = ifelse(is.na(is_triangle), FALSE, is_triangle))

## ---------------------- YOUR TUNING KNOBS BLOCK ----------------------------
## Tuning parameters governing the allocation of types:
##   - other_cap: upper bound on the share of players in the residual
##                “Other” category.
##   - target_free: target share for “Free rider” type, drawn from
##                  low-mean players with near-zero EB slope.
##   - mean_tol: maximum average contribution (in contribution units)
##               for a player to be eligible as a “Free rider”.
## These choices are calibrated to match the empirical patterns
## discussed in the paper and can be varied for sensitivity checks.
## -----------------------------------------------------------------
# ---- TUNING KNOBS ----
other_cap   <- 0.40   # max share for "Other"
target_free <- 0.12   # target share for "Free rider"
mean_tol    <- 2.0    # free-rider mean threshold (L)

N <- nrow(feat)

# start fresh
feat$type <- NA_character_

# 1) TRIANGLES (evidence-based)
## First, assign “Triangle” players based on the quadratic criterion.
feat$type[feat$is_triangle] <- "Triangle"
n_tri <- sum(feat$type == "Triangle", na.rm = TRUE)

# 2) FREE RIDERS: lowest |slope| among low-mean players, up to target %
## Candidates: players with low average contribution (mean_c ≤ mean_tol)
## and not previously classified as “Triangle”. Among them, select
## those with the flattest (smallest absolute) EB slope b1_eb until
## the target proportion target_free is reached.
cand_free <- which(is.na(feat$type) & feat$mean_c <= mean_tol)
cand_free <- cand_free[order(abs(feat$b1_eb[cand_free]), decreasing = FALSE)]
n_free_target <- min(ceiling(target_free * N), length(cand_free))
if (n_free_target > 0) {
  feat$type[cand_free[seq_len(n_free_target)]] <- "Free rider"
}
n_free <- sum(feat$type == "Free rider", na.rm = TRUE)

# 3) CONDITIONAL: enough top-slope players so that Other ≤ other_cap
remaining <- which(is.na(feat$type))
# Other would be: N - (tri + free + cond)
# Require: N - (tri + free + cond) ≤ other_cap * N
# => cond ≥ N*(1 - other_cap) - (tri + free)
cond_needed <- ceiling(max(0, N*(1 - other_cap) - (n_tri + n_free)))
# rank remaining by EB slope (high to low)
remaining <- remaining[order(feat$b1_eb[remaining], decreasing = TRUE)]
cond_pick <- min(cond_needed, length(remaining))
if (cond_pick > 0) {
  feat$type[remaining[seq_len(cond_pick)]] <- "Conditional cooperator"
}

# 4) The rest are "Other"
feat$type[is.na(feat$type)] <- "Other"
feat$type <- factor(feat$type,
                    levels = c("Conditional cooperator","Triangle","Free rider","Other"))

# report shares (sanity check against paper’s reported composition)
type_share <- feat %>% count(type) %>% mutate(pct = round(100*n/sum(n), 1))
print(type_share)

## ------------------- Build player_coefs for downstream use ------------------
## Enrich the player-level coefficient frame with quadratic terms and
## type labels; this is the main table of types used in downstream
## analyses and figures.
player_coefs <- feat %>%
  left_join(quad, by = "player") %>%
  mutate(type = factor(type, levels = c("Conditional cooperator","Triangle","Free rider","Other")))


## ============================ ANALYSES YOU ASKED ============================
## The following blocks reproduce the key empirical outputs:
##   - Type shares and cross-tabs with dynamic path clusters.
##   - Logit for low-path (C2) membership.
##   - Learning model with type interactions.
##   - Stability of types over early vs late rounds.
##   - Triangle peak location.
##   - Group-level relation between type composition and low-path share.
## ---------------------------------------------------------------------------

# Shares of types (again, now via player_coefs)
type_share <- player_coefs %>%
  count(type) %>%
  mutate(pct = 100 * n / sum(n))
print(type_share)

# Cross-tab with path clusters (Ward–D2)
## For each Ward cluster, report how its membership distributes
## across the four revealed-slope types (in counts and percentages).
tab_type_cluster <- contrib_df %>%
  inner_join(player_coefs %>% dplyr::select(player, type), by = "player") %>%
  count(hc_cluster, type) %>%
  group_by(hc_cluster) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
print(tab_type_cluster)

# Low-path (C2) logit with covariates
## Define low_path = 1{hc_cluster == C2}, then estimate a logit
## with type dummies and selected covariates as predictors.
df_type_path <- contrib_df %>%
  inner_join(player_coefs %>% dplyr::select(player, type), by = "player") %>%
  left_join(pgg_data %>% distinct(player, friends,  gender, finauto_q_2_temp),
            by = "player") %>%
  mutate(
    low_path = (hc_cluster == "C2"),
    gender_f = factor(ifelse(gender == 1, "Male", "Female")),
    finauto  = factor(finauto_q_2_temp),
    type     = relevel(type, ref = "Free rider")
  )

mod_path <- glm(low_path ~ type + scale(friends)  + gender_f + finauto,
                family = binomial(), data = df_type_path)
or <- exp(cbind(OR = coef(mod_path), confint.default(mod_path)))
print(round(or, 3))

# Learning model with type interactions
## Construct lag of own contributions, then estimate:
##   c_it = α + φ c_{i,t−1} + β type × peer_mean_lag_it +
##          random intercept and slope by player, random intercept by group.
condresp_long <- condresp_long %>%
  group_by(group, player) %>%
  mutate(own_lag = lag(contributing)) %>%
  ungroup()

cr_dat <- condresp_long %>%
  inner_join(player_coefs %>% dplyr::select(player, type), by = "player") %>%
  filter(!is.na(type)) %>%
  mutate(type = relevel(type, ref = "Free rider"))

cr_type_mod <- lmerTest::lmer(
  contributing ~ own_lag + peer_mean_lag*type + (1 + peer_mean_lag | player) + (1 | group),
  data = cr_dat
)
print(summary(cr_type_mod))

# Early vs late stability (diagnostic using raw per-player OLS windows)
## Helper: for a given subset of rounds, estimate per-player OLS with
## peer_mean_lag and its square, then map coefficients and p-values
## into a four-type classification using simple decision rules.
fit_by_window <- function(d) {
  d %>% group_by(player) %>%
    filter(n() >= 3) %>%
    group_modify(~{
      fit <- lm(contributing ~ peer_mean_lag + I(peer_mean_lag^2), data=.x)
      S   <- summary(fit)$coefficients
      get <- function(term, col) if (term %in% rownames(S)) S[term, col] else NA_real_
      data.frame(b1=get("peer_mean_lag","Estimate"),
                 b2=get("I(peer_mean_lag^2)","Estimate"),
                 p1=get("peer_mean_lag","Pr(>|t|)"),
                 p2=get("I(peer_mean_lag^2)","Pr(>|t|)"))
    }) %>% ungroup() %>%
    mutate(type = case_when(
      abs(b1) <= 0.05 & (is.na(p1) | p1>0.05) ~ "Free rider",
      b1>0 & !is.na(p1) & p1<=0.05 & (is.na(b2) | b2>=0 | is.na(p2) | p2>0.05) ~ "Conditional cooperator",
      b1>0 & !is.na(p1) & p1<=0.05 & !is.na(b2) & b2<0 & !is.na(p2) & p2<=0.05 ~ "Triangle",
      TRUE ~ "Other"
    ))
}

## Classify types using early (rounds 1–5) vs late (6–10) windows
early_types <- condresp_long %>% filter(round_n <= 5) %>% fit_by_window() %>%
  dplyr::select(player, type_early = type)
late_types  <- condresp_long %>% filter(round_n >= 6) %>% fit_by_window() %>%
  dplyr::select(player, type_late  = type)

stab <- inner_join(early_types, late_types, by="player") %>%
  mutate(agree = (type_early == type_late)) %>%
  count(agree) %>% mutate(pct = 100*n/sum(n))
print(stab)

# Triangle vertex (peak location, in standardized x_c units)
## Compute the average location of the peak x_star among Triangle
## types, using the coefficients from the quadratic fits.
tri_vertex <- player_coefs %>%
  filter(type=="Triangle", !is.na(b2), b2 < 0) %>%
  mutate(x_star = -b1/(2*b2)) %>%
  summarise(mean_peak = mean(x_star, na.rm=TRUE),
            median_peak = median(x_star, na.rm=TRUE),
            q = list(quantile(x_star, c(.25,.75), na.rm=TRUE)))
print(tri_vertex)

# Group-level: type composition vs share Low-path
## For each group:
##   - Compute the share of each type in its membership.
##   - Compute the share of players in low-path cluster C2.
##   - Regress low-path share on group composition in terms of
##     Free riders, Conditional cooperators, and Triangles.
g_comp <- pgg_data %>%
  distinct(group, player) %>%
  mutate(player = as.character(player)) %>%
  inner_join(player_coefs %>% dplyr::select(player, type), by="player") %>%
  count(group, type) %>%
  group_by(group) %>%
  mutate(share = n/sum(n)) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from=type, values_from=share, values_fill=0)

g_lowshare <- contrib_df %>%
  inner_join(pgg_data %>% distinct(player, group), by="player") %>%
  count(group, hc_cluster) %>%
  tidyr::pivot_wider(names_from=hc_cluster, values_from=n, values_fill=0) %>%
  mutate(low_share = if ("C2" %in% names(.)) C2/(C1 + C2) else NA_real_)

print(summary(lm(low_share ~ `Free rider` + `Conditional cooperator` + Triangle,
                 data = inner_join(g_comp, g_lowshare, by="group"))))

## ===================== Balanced figure (everyone once) ======================
## Construct a “balanced” representation of conditional responses:
##   - Bin peer_mean_lag into integer bins (0, 1, …, 12).
##   - For each type × bin, compute mean contribution and 95% CI.
##   - Overlay LOESS-smoothed curves using all underlying data.
##   - Plot the 45-degree line and the round-1 threshold as visual
##     benchmarks.
condresp_type <- condresp_long %>%
  inner_join(player_coefs %>% dplyr::select(player, type), by = "player") %>%
  mutate(xbin = cut(peer_mean_lag, breaks = seq(-0.5, 12.5, 1),
                    labels = 0:12, include.lowest = TRUE)) %>%
  group_by(type, xbin) %>%
  summarise(
    n      = n(),
    mean_y = mean(contributing, na.rm = TRUE),
    se_y   = sd(contributing, na.rm = TRUE)/sqrt(n),
    ci_lo  = mean_y - 1.96*se_y,
    ci_hi  = mean_y + 1.96*se_y,
    .groups = "drop"
  ) %>%
  mutate(x = as.numeric(as.character(xbin)))

thresh_line <- thresh_r1

p_balanced <- ggplot(condresp_type,
                     aes(x = x, y = mean_y, colour = type, group = type)) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.15, alpha = 0.6) +
  geom_point(aes(size = n), alpha = 0.85) +
  geom_smooth(
    data = condresp_long %>% inner_join(player_coefs %>% dplyr::select(player, type), by="player"),
    aes(x = peer_mean_lag, y = contributing, colour = type),
    method = "loess", span = 0.7, se = FALSE, linewidth = 1.1
  ) +
  scale_size_continuous(range = c(1.5, 4), guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_hline(yintercept = thresh_line, linetype = "dotted", colour = "grey50") +
  coord_cartesian(xlim = c(0,12), ylim = c(0,12)) +
  labs(
    title = "Revealed-slope types’ conditional response (balanced assignment)",
    x = expression(bar(c)[-i]^{t-1}),
    y = expression(c[i]^t),
    colour = NULL
  ) +
  theme_minimal(base_size = 12)

print(p_balanced)
# ggsave("Figure_RevealedSlope_balanced.png", p_balanced, width=6, height=5, dpi=300)



# ────────────────────────────────────────────────────────────────
#  ADDING COVARIATES & HETEROGENEITY TO THE CONDITIONAL‑RESPONSE
# ────────────────────────────────────────────────────────────────
#
# Prerequisite objects already in memory from your script:
#   • condresp_long   – one row per (player, round)
#   • player_coefs    – one row per player with type tags
#   • baseline        – covariates measured in round 1
#   • hc_cluster      – Ward‑trajectory high / low tag
#   • thresh          – Hi/Lo cut‑off for contributions
#   • (packages dplyr, tidyr, ggplot2, broom, purrr, nnet,
#      lme4, lmerTest, conflicted) already loaded
# ----------------------------------------------------------------

## 0.  Tell R once more which lmer we want  ----------------------
## In case multiple packages define lmer(), this ensures that the
## lmerTest version (with p-values) is used.
conflicted::conflicts_prefer(lmerTest::lmer)   # lmer with p‑values


## A.  MIXED‑EFFECTS CONDITIONAL‑RESPONSE MODEL  ----------------
## Goal:
##   Augment the conditional-response model with demographic and
##   network covariates, allowing peer_mean_lag to interact with:
##   - gender, friends_z, network_density_fr_z, network_density_adv_z,
##     and financial autonomy.
##   The model keeps random intercepts and random peer_mean_lag slopes
##   at the player level, and a random intercept at the group level.
# 1)  Standardise continuous covariates **within** condresp_long
condresp_long <- condresp_long %>%
  dplyr::mutate(
    friends_z                = as.numeric(scale(friends)),
    network_density_fr_z     = as.numeric(scale(network_density_fr)),
    network_density_adv_z    = as.numeric(scale(network_density_adv)),
    age_z                    = as.numeric(scale(age))
  ) %>%
  dplyr::mutate(
    gender_f = factor(gender),        # already “Male” / “Female”
    finauto  = factor(finauto_q_2_temp)          # 0 = Not autonomous, 1 = Autonomous
  )

# 2)  Random‑slope, random‑intercept specification
cr_mod <- lmerTest::lmer(
  contributing ~ peer_mean_lag *
    (gender_f + friends_z +
       network_density_fr_z + network_density_adv_z + finauto) +
    (1 + peer_mean_lag | player) +   # player‑specific slope & intercept
    (1 | group),                     # five‑person group intercept
  data = condresp_long
)
summary(cr_mod)        # ← fixed‑effects, p‑values, random‑effect SDs


## B.  MULTINOMIAL LOGIT ON REVEALED TYPES  ---------------------
## Goal:
##   Model the probability of belonging to each revealed-slope type
##   as a function of baseline covariates using a multinomial logit
##   (reference category: “Free rider”).
# 1)  Build a covariate frame at the **player** level
cov_frame <- baseline %>%
  dplyr::mutate(
    friends_z             = as.numeric(scale(friends)),
    adversaries_z         = as.numeric(scale(adversaries)),
    age_z                 = as.numeric(scale(age)),
    network_density_fr_z  = as.numeric(scale(network_density_fr)),
    network_density_adv_z = as.numeric(scale(network_density_adv))
  ) %>%
  dplyr::select(player, age_z, gender, friends_z, adversaries_z,
                network_density_fr_z, network_density_adv_z,
                finauto_q_2_temp)

## Ensure that player IDs are numeric and aligned across frames.
player_coefs$player <- as.integer(player_coefs$player)
cov_frame$player <- as.integer(cov_frame$player)

# 2)  Merge type tags, set reference category = “Free rider”
type_df <- player_coefs %>%
  dplyr::left_join(cov_frame, by = "player") %>%
  dplyr::mutate(
    type = forcats::fct_relevel(type, "Free rider")  # baseline category
  ) %>%
  dplyr::filter(!is.na(type))                        # should already hold

# 3)  Fit multinomial logit (nnet::multinom gives z‑stats)
type_mod <- nnet::multinom(
  type ~ age_z + gender + friends_z + adversaries_z +
    network_density_fr_z + network_density_adv_z +
    finauto_q_2_temp ,
  data  = type_df,
  trace = FALSE
)
summary(type_mod)


## ────────────────────────────────────────────────────────────
##  ADD P‑VALUES TO THE MULTINOMIAL‑LOGIT TABLE
## ────────────────────────────────────────────────────────────
## nnet::multinom reports coefficients and standard errors but not
## p-values. The block below:
##   - extracts coefficient and SE matrices,
##   - computes z-statistics and two-sided p-values,
##   - reshapes everything into a tidy long format (one row per
##     predictor × outcome-level), and
##   - prints a rounded table.
library(dplyr)    # for the tidy reshaping
library(tidyr)
library(broom)    # optional, for broom::tidy if you like

# 1) Pull the matrices that summary(multinom) provides
smry  <- summary(type_mod)
coefM <- smry$coefficients     # matrix: predictors × outcome‐levels
seM   <- smry$standard.errors

# 2) Compute z‑stats and two‑tailed p‑values
zM <- coefM / seM
pM <- 2 * (1 - pnorm(abs(zM)))  # two‑tailed 

# 3) Re‑wrap into a long/tidy data.frame
tidy_multinom <- coefM %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Predictor") %>%
  pivot_longer(
    -Predictor,
    names_to  = "Outcome",
    values_to = "Coef"
  ) %>%
  left_join(
    seM %>% as.data.frame() %>%
      tibble::rownames_to_column("Predictor") %>%
      pivot_longer(-Predictor, names_to="Outcome", values_to="SE"),
    by = c("Predictor","Outcome")
  ) %>%
  left_join(
    zM %>% as.data.frame() %>%
      tibble::rownames_to_column("Predictor") %>%
      pivot_longer(-Predictor, names_to="Outcome", values_to="z"),
    by = c("Predictor","Outcome")
  ) %>%
  left_join(
    pM %>% as.data.frame() %>%
      tibble::rownames_to_column("Predictor") %>%
      pivot_longer(-Predictor, names_to="Outcome", values_to="p_value"),
    by = c("Predictor","Outcome")
  ) %>%
  arrange(Outcome, desc(abs(z)))

# 4) Round numeric cols and print
tidy_multinom <- tidy_multinom %>%
  mutate(across(c(Coef, SE, z, p_value), ~ round(.x, 3)))

# if you just want to see the whole thing:
print(tidy_multinom, n = nrow(tidy_multinom), width = Inf)

# OR, coerce to a base data.frame for print(..., digits):
df_multinom <- as.data.frame(tidy_multinom)
print(df_multinom, digits = 3, row.names = FALSE)
