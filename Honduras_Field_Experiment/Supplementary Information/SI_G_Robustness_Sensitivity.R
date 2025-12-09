###############################################################################
# ROBUSTNESS CHECKS – supplementary sensitivity analyses for LL/LH/HL/HH
# transition labels and multinomial models (Section 5; SI G.1–G.3).
###############################################################################
# This block is self‑contained: it reloads the data and reconstructs the
# necessary derived objects used in the robustness exercises below.

pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

conflicted::conflicts_prefer(nnet::multinom)

pgg_data$player <- as.character(pgg_data$player)

thresh <- 6

###############################################################################
# 0. PACKAGES ----------------------------------------------------------------
###############################################################################
# Ensure that all packages required for the robustness checks are installed
# and attached (beyond those used in the main analysis script).
pkgs <- c(
  "dplyr", "tidyr", "ggplot2", "nnet", "lme4",
  "tableone", "marginaleffects", "glmnet", "randomForest",
  "rpart", "rpart.plot"
)
new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
lapply(pkgs, library, character.only = TRUE)

###############################################################################
# 1. END-STATE LABELS (LL LH HL HH) + COVARIATES -----------------------------
###############################################################################
# Rebuild High/Low states using the round‑1 mean as the threshold, then
# construct end‑state transitions LL/LH/HL/HH and merge in baseline covariates.
thresh <- mean(
  pgg_data$contributing[pgg_data$round_n == 1],
  na.rm = TRUE
)

states <- pgg_data %>%
  filter(round_n %in% c(1, 10)) %>%
  mutate(state = if_else(contributing > thresh, "H", "L")) %>%
  dplyr::select(player, village_code, round_n, state) %>%
  pivot_wider(
    names_from = round_n,
    names_prefix = "r",
    values_from = state
  ) %>%
  mutate(trans = paste0(r1, r10))

baseline <- pgg_data %>%
  filter(round_n == 1) %>%
  distinct(player, .keep_all = TRUE) %>%
  dplyr::select(
    player, village_code, age, gender, friends, adversaries, FI,
    marital_status, network_density_fr, network_density_adv, network_size,
    b0100, b0600, b0200, access_routes, finauto_q_2_temp
  )

df <- states %>%
  left_join(baseline, by = c("player", "village_code")) %>%
  mutate(
    gender = factor(gender),
    marital_status = factor(marital_status),
    b0600 = relevel(factor(b0600), ref = "2")
  )

###############################################################################
# 2. STANDARDISE CONTINUOUS PREDICTORS ---------------------------------------
###############################################################################
# Standardise continuous covariates (z‑scores) so coefficients are comparable
# across variables and the multinomial logit is numerically well‑behaved.
cont_vars <- c(
  "age", "friends", "adversaries",
  "network_density_fr", "network_density_adv", "network_size",
  "b0100"
)

df <- df %>%
  mutate(across(all_of(cont_vars), ~ as.numeric(scale(.x)), .names = "{.col}_z"))

###############################################################################
# 3. COUNTS -------------------------------------------------------------------
###############################################################################
# Report the empirical distribution of end‑state transitions: LL, LH, HL, HH.
cat("\n--- Counts of LL / LH / HL / HH ---\n")
print(states %>% count(trans))

###############################################################################
# 4. MULTINOMIAL LOGIT (baseline = LL) ---------------------------------------
###############################################################################
# Set up the multinomial logit specification for transitions, using LL as
# the reference category and including standardised covariates plus factors.
df$trans <- factor(df$trans, levels = c("LL", "LH", "HL", "HH"))

multi_formula <- as.formula(
  paste(
    "trans ~",
    paste0(cont_vars, "_z", collapse = " + "),
    "+ gender + marital_status + b0600"
  )
)

# Classify each player as H/L in rounds 1 and 10, then define transition type
thresh_r1 <- mean(pgg_data$contributing[pgg_data$round_n == 1], na.rm = TRUE)
states <- pgg_data %>%
  filter(round_n %in% c(1, 10)) %>%
  mutate(state = if_else(contributing > thresh_r1, "H", "L")) %>%
  dplyr::select(player, round_n, state) %>%
  pivot_wider(names_from = round_n, names_prefix = "r", values_from = state) %>%
  mutate(trans = paste0(r1, r10), player = as.character(player))

# ---------------------------------------------------------------------------
# 0. House‑keeping: enforce a single, consistent type for the player ID
#    across all data frames used in joins or model fitting.
# ---------------------------------------------------------------------------
pgg_data$player <- as.character(pgg_data$player)
states$player         <- as.character(states$player)
#df$player             <- as.character(df$player)

# ---------------------------------------------------------------------------
# 5.1 Plateau of the “High” share over the 10 rounds
# ---------------------------------------------------------------------------
# Compute, and (optionally) plot, the fraction of decisions classified as High
# in each round, to verify that the aggregate H/L composition stabilises well
# before the final period.
prop_high <- pgg_data %>%
  mutate(state = if_else(contributing > thresh, "H", "L")) %>%
  group_by(round_n) %>%
  summarise(prop_H = mean(state == "H"), .groups = "drop")

ggplot(prop_high, aes(round_n, prop_H)) +
  geom_line() + geom_point() +
  scale_x_continuous(breaks = 1:10) +
  labs(title = "Proportion of High contributors by round",
       subtitle = sprintf("Threshold = %.2f (mean of round 1)", thresh),
       x = "Round", y = "Proportion High") +
  theme_minimal()

# ---------------------------------------------------------------------------
# 5.2 Threshold‑choice sensitivity
#     • Compare the original round‑1‑mean cut‑off with (i) the global mean
#       across all observations and (ii) the theoretical cut‑off at 6.
#     • Check how often transition labels change and whether multinomial
#       coefficients preserve their sign under alternative thresholds.
# ---------------------------------------------------------------------------
mk_trans <- function(thr){
  pgg_data %>%
    filter(round_n %in% c(1, 10)) %>%
    mutate(state = if_else(contributing > thr, "H", "L")) %>%
    dplyr::select(player, round_n, state) %>%
    pivot_wider(names_from = round_n, names_prefix = "r", values_from = state) %>%
    mutate(trans  = paste0(r1, r10),
           player = as.character(player))
}

trans_global <- mk_trans(mean(pgg_data$contributing, na.rm = TRUE))  # global mean
trans_theory <- mk_trans(6)                                          # threshold = 6

cat("\n--- Confusion: original vs global-mean cut-off ---\n")
print(table(original = states$trans,
            global   = trans_global$trans))

cat("\n--- Confusion: original vs theoretical (6-point) cut-off ---\n")
print(table(original = states$trans,
            theory   = trans_theory$trans))

# ---------- Sign agreement of multinomial logit coefficients ----------------
# Refit the multinomial model under alternative thresholds and compare the
# sign pattern of coefficients with the baseline specification.
get_coef <- function(df_trans){
  tmp <- df %>%
    dplyr::select(-trans) %>%                           # drop old transition label
    left_join(df_trans %>% dplyr::select(player, trans), by = "player") %>%
    mutate(trans = factor(trans, levels = c("LL","LH","HL","HH")))
  coef(multinom(multi_formula, data = tmp, trace = FALSE))
}

sign_same <- sign(get_coef(states)) == sign(get_coef(trans_global))
cat("\nNumber of coefficients with identical sign (orig vs global): ",
    sum(sign_same), " / ", length(sign_same), "\n", sep = "")

# ---------------------------------------------------------------------------
# 5.3 Multiple‑flip statistics (how often players cross the H/L threshold)
# ---------------------------------------------------------------------------
# Count, for each player, how many times they switch between High and Low
# across the ten rounds, then summarise the prevalence of single‑ and
# multiple‑flip behaviour in the sample.
flip_stats <- pgg_data %>%
  mutate(state = if_else(contributing > thresh, "H", "L")) %>%
  arrange(player, round_n) %>%
  group_by(player) %>%
  mutate(flip = state != lag(state, default = first(state))) %>%
  summarise(n_flips = sum(flip), .groups = "drop") %>%
  summarise(
    total_players = n(),
    flipped_once  = sum(n_flips == 1),
    flipped_multi = sum(n_flips  > 1),
    prop_multi    = flipped_multi / total_players
  )

cat("\n--- Multiple-flip statistics (threshold = round-1 mean) ---\n")
print(flip_stats) 

# ---------------------------------------------------------------------------
# ROBUSTNESS 5.4 – End‑state based on the last TWO rounds (9 & 10)
# ---------------------------------------------------------------------------
# 1. Helper: classify each player as High/Low in round 1 and in the “last
#    period” defined over rounds 9 and 10 (High if high in at least one of
#    the two rounds, or in both rounds when both = TRUE).
# ---------------------------------------------------------------------------
mk_two_round_trans <- function(thr, both = FALSE){
  last_state <- pgg_data %>%                   # get r9 & r10 states
    filter(round_n %in% c(9,10)) %>%
    mutate(state = if_else(contributing > thr, "H", "L")) %>%
    group_by(player, village_code) %>%
    summarise(end_high = if (both) all(state == "H") else any(state == "H"),
              .groups = "drop") %>%
    mutate(rLast = if_else(end_high, "H", "L"))
  
  pgg_data %>%                                # get r1 state
    filter(round_n == 1) %>%
    mutate(r1 = if_else(contributing > thr, "H", "L")) %>%
    dplyr::select(player, village_code, r1) %>%
    left_join(last_state, by = c("player","village_code")) %>%
    mutate(trans = paste0(r1, rLast))             # LL, LH, HL, HH
}

# 2. Construct the new two‑round transition labels and compare them to the
#    original r10‑based labels stored in `states`.
# ---------------------------------------------------------------------------
states_two <- mk_two_round_trans(thresh)             # lenient version
cat("\n--- Confusion: original (r10) vs two‑round endpoint (r9|r10) ---\n")
print(table(original = states$trans, twoRound = states_two$trans))

# 3. Re‑run the multinomial logit with the two‑round endpoint labels
#    to check whether the direction and size of covariate effects are stable.
# ---------------------------------------------------------------------------
df_two <- df %>%
  dplyr::select(-trans) %>%
  left_join(states_two %>% dplyr::select(player, trans), by = "player") %>%
  mutate(trans = factor(trans, levels = c("LL","LH","HL","HH")))

multi_two <- multinom(multi_formula, data = df_two, trace = FALSE)
cat("\n--- Multinomial Logit (endpoint = r9 | r10 High, LL baseline) ---\n")
print(summary(multi_two)) 

###############################################################################
# END OF ROBUSTNESS BLOCK
###############################################################################

###############################################################################
## ROBUSTNESS CHECK — INTERPRETATION OF THE NUMERICAL OUTPUT ABOVE
###############################################################################

## ── 5.1 Plateau of “High” share ────────────────────────────────────────────
## • The prop_high series tracks, by round, the fraction of decisions above
##   the round‑1 mean threshold. The line plot (not shown here) declines
##   modestly in early rounds and then flattens, with only small movements
##   over rounds 4–10.
## • This stabilisation of the aggregate High vs. Low mix well before the
##   final period justifies treating round 10 as an empirical “end‑state”.

## ── 5.2 Threshold‑choice sensitivity ───────────────────────────────────────
## • Using the round‑1 mean as the cut‑off yields 1,055 HH, 472 HL, 195 LH and
##   869 LL trajectories.
## • Re‑defining High/Low with either the global mean of all contributions or
##   the theoretical 6‑point cut‑off produces exactly the same LL/LH/HL/HH
##   labels: both confusion matrices are purely diagonal.
## • The multinomial logit refitted under the global‑mean cut‑off preserves
##   the sign of every coefficient (36 / 36). Thus, both the classification
##   and the qualitative pattern of covariate effects are invariant to these
##   reasonable alternative thresholds.

## ── 5.3 Multiple‑flip statistics ───────────────────────────────────────────
##   total players     : 2 591
##   flipped exactly 1 :   266  (≈10 %)
##   flipped > 1 times : 1 110  (≈43 %)
##   share multi‑flip  :  0.428
## • A substantial share of subjects move back and forth across the High/Low
##   line at least twice, indicating active experimentation at the individual
##   level.
## • Combined with the aggregate plateau in 5.1, this supports the view that
##   these flips reflect short‑run “search” behaviour around a stable long‑run
##   mix of High and Low contributors, rather than a third enduring regime.

## ── 5.4 Endpoint based on the last TWO rounds (9 & 10) ─────────────────────
## • When the endpoint is defined as “High in at least one of rounds 9 or 10”
##   instead of “High in round 10 only”, HH and LH cases are fully stable.
##   Re‑classification occurs only among HL and LL trajectories (HL → HH and
##   LL → LH), i.e. for players close to the threshold in late rounds.
## • The multinomial logit with these two‑round endpoints (multi_two) yields
##   coefficient patterns closely aligned with the baseline model: male gender
##   and friend degree remain positively associated with High outcomes, and
##   no key covariate reverses sign.

## ── Bottom‑line robustness take‑away ───────────────────────────────────────
## • Endpoint labels (LL / LH / HL / HH) are highly stable under plausible
##   alternative High/Low cut‑offs and under a two‑round endpoint definition.
## • The direction (and most significance patterns) of the multinomial‑logit
##   coefficients is preserved across these perturbations.
## • Many players explore both sides of the threshold, but by round 10 almost
##   all have settled into one of two persistent paths – a High‑contribution
##   regime and a Low‑contribution regime – consistent with the bifurcation
##   narrative developed in the main text.
###############################################################################
