#### PACKAGES ------------------------------------------------------------------
# Load and attach all required packages.
# If a package is missing, install it first (with dependencies).
pkgs <- c("dplyr","tidyr","stringr","AER","sandwich","lmtest")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))

# Fix random seed for reproducibility of any random operations later.
set.seed(123)

# Load the dataset used throughout the IV analysis.
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)


#### 1) BUILD CLEAN FRAME + LAGGED LOO PEERS -----------------------------------
# Goal of this block:
#   • Harmonize coding of key covariates (gender, marital status, religion, etc.).
#   • Construct group-level leave-one-out (LOO) contemporaneous peer means.
#   • Construct lagged peer means peer_mean_tm1 at the group level.
#
# The resulting object `dat` is the working panel used for the IV construction.

dat <- pgg_data %>%
  mutate(
    # Binary sex indicator: 1 = male, 0 = female (all other/unrecognized → NA)
    male = case_when(
      is.numeric(gender)   ~ as.integer(gender == 1),
      is.character(gender) ~ as.integer(tolower(gender) %in% c("1","man","male","m")),
      is.factor(gender)    ~ as.integer(tolower(as.character(gender)) %in% c("1","man","male","m")),
      TRUE ~ NA_integer_
    ),
    # Marital status: 1 = married / civil union, 0 = otherwise
    married = case_when(
      marital_status %in% c(1,"1","Married","Civil Union","married","civil union") ~ 1L,
      marital_status %in% c(0,"0","Single","Separated","Widowed","Divorced","Other") ~ 0L,
      TRUE ~ NA_integer_
    ),
    # Religion dummies based on b0600 coding in the paper:
    #   rel_protestant = 1{Protestant}, rel_none = 1{no religion};
    #   Catholics form the omitted category.
    rel_protestant = as.integer(b0600 == 1),
    rel_none       = as.integer(b0600 == 0),
    
    # Coerce various survey covariates to numeric, suppressing warnings if coercion fails.
    FI          = suppressWarnings(as.numeric(FI)),
    b0200       = suppressWarnings(as.numeric(b0200)),
    b0100       = suppressWarnings(as.numeric(b0100)),
    age         = suppressWarnings(as.numeric(age)),
    friends     = suppressWarnings(as.numeric(friends)),
    adversaries = suppressWarnings(as.numeric(adversaries)),
    
    # Game indexing variables
    round_n      = as.integer(round_n),
    player       = as.factor(player),
    village_code = as.factor(village_code),
    # Village-by-round factor (for later FE / demeaning).
    vill_round   = interaction(village_code, round_n, drop = TRUE)
  ) %>%
  # In each group×round, compute LOO peer means.
  group_by(group, round_n) %>%
  mutate(
    grp_n       = dplyr::n(),
    # contemporaneous leave-one-out peer mean contribution
    # (sum of others' contributions divided by group size minus one)
    peer_mean_t = ifelse(grp_n > 1,
                         (sum(contributing, na.rm = TRUE) - contributing) / (grp_n - 1),
                         NA_real_)
  ) %>% ungroup() %>%
  # Sort by group and round so that lags are well-defined.
  arrange(group, round_n) %>%
  group_by(group) %>%
  # Lag the contemporaneous peer mean by one round within group.
  mutate(peer_mean_tm1 = dplyr::lag(peer_mean_t)) %>%
  ungroup()


#### 2) TIME-VARYING SHIFT–SHARE INSTRUMENT ------------------------------------
# Goal of this block:
#   • Construct a time-varying shift–share IV Z_shiftshare for the lagged peer mean.
#   • Shares are based on leave-one-out group composition in baseline traits.
#   • Shifters are round-specific cross-sectional means of the same traits.
#
# base_z: the list of baseline covariates used in the shift–share construction.

base_z <- c("male","age","friends","adversaries","FI","married",
            "b0100","rel_protestant","rel_none","b0200")

# gmem: one row per (group, player), with the relevant baseline covariates.
gmem <- dat %>% dplyr::select(group, player, dplyr::all_of(base_z)) %>% distinct()

# gsum: for each group, sum of each baseline covariate and the group size G.
gsum <- gmem %>% group_by(group) %>%
  summarise(across(all_of(base_z), ~ sum(.x, na.rm = TRUE), .names = "sum_{.col}"),
            G = dplyr::n(), .groups = "drop")

# gshare: leave-one-out group shares for each player and each baseline trait.
#   share_{k,ig} = (sum_k_g - x_{i,k}) / (G - 1)
gshare <- gmem %>%
  tidyr::pivot_longer(cols = all_of(base_z), names_to = "var", values_to = "val") %>%
  left_join(
    gsum %>% tidyr::pivot_longer(cols = dplyr::starts_with("sum_"),
                                 names_to = "sumvar", values_to = "sumval") %>%
      mutate(var = sub("^sum_","",sumvar)) %>%
      dplyr::select(group, var, sumval, G),
    by = c("group","var")
  ) %>%
  mutate(share = (sumval - val) / pmax(G - 1, 1),
         share_name = paste0("share_", var)) %>%
  dplyr::select(group, player, share_name, share) %>%
  tidyr::pivot_wider(names_from = share_name, values_from = share)

# shifts: time-varying “shifters” = round-specific cross-village means of base_z.
# For each round, compute µ_t (mean of each baseline covariate), then lag by 1 round.
shifts <- dat %>% group_by(round_n) %>%
  summarise(across(all_of(base_z), ~ mean(.x, na.rm = TRUE), .names = "mu_{.col}"),
            .groups = "drop") %>%
  arrange(round_n) %>%
  mutate(across(dplyr::starts_with("mu_"), ~ dplyr::lag(.)))

# Bring shares and shifters together for each player×round observation.
ss <- dat %>%
  left_join(gshare, by = c("group","player")) %>%
  left_join(shifts, by = "round_n")

# Construct the scalar shift–share instrument Z_shiftshare as:
#   Z_shiftshare_igt = Σ_k share_{k,ig} × µ_{k,t-1}
share_vars <- grep("^share_", names(ss), value = TRUE)
mu_vars    <- sub("^share_", "mu_", share_vars)
ss$Z_shiftshare <- rowSums(ss[, share_vars] * ss[, mu_vars], na.rm = TRUE)

# Final IV analysis frame:
#   • outcome: contributing
#   • endogenous regressor: peer_mean_tm1
#   • instrument: Z_shiftshare
#   • controls: male, age, friends, adversaries, FI, married, b0100, rel_protestant, rel_none, b0200
iv_df <- ss %>%
  dplyr::select(group, player, vill_round,
         contributing, peer_mean_tm1, Z_shiftshare,
         male, age, friends, adversaries, FI, married,
         b0100, rel_protestant, rel_none, b0200) %>%
  mutate(peer_mean_tm1 = as.numeric(peer_mean_tm1),
         Z_shiftshare  = as.numeric(Z_shiftshare)) %>%
  # Keep only observations with complete information for the main variables.
  filter(complete.cases(peer_mean_tm1, Z_shiftshare, contributing, group, player, vill_round))


#### 3) TWO-WAY DEMEANING (player FE + village×round FE) -----------------------
# Goal:
#   • Implement the “within” transformation to absorb individual FE and village×round FE
#     via explicit two-way demeaning (Frisch–Waugh–Lovell logic).
#   • tw_demean(x, f1, f2) = x - mean(x|f1) - mean(x|f2) + overall mean(x)

tw_demean <- function(x, f1, f2){
  mu <- mean(x, na.rm = TRUE)                                        # overall mean
  m1 <- ave(x, f1, FUN = function(z) mean(z, na.rm = TRUE))          # mean by first factor (e.g. player)
  m2 <- ave(x, f2, FUN = function(z) mean(z, na.rm = TRUE))          # mean by second factor (e.g. village×round)
  x - m1 - m2 + mu                                                   # residual after removing both FE
}

# Two-way demean the outcome, endogenous regressor, and instrument.
y_w <- tw_demean(iv_df$contributing,  iv_df$player, iv_df$vill_round)
e_w <- tw_demean(iv_df$peer_mean_tm1, iv_df$player, iv_df$vill_round)
z_w <- tw_demean(iv_df$Z_shiftshare,  iv_df$player, iv_df$vill_round)

# List of time-varying controls that may survive the fixed effects transformation.
ctrl_all <- c("male","age","friends","adversaries","FI","married",
              "b0100","rel_protestant","rel_none","b0200")

# Two-way demeaning of each control; drop controls that become constant or degenerate
# after demeaning (to avoid singularity in the regression).
ctrl_w_df <- lapply(ctrl_all, function(vn){
  xv <- iv_df[[vn]]; if (is.null(xv)) return(NULL)
  res <- tw_demean(as.numeric(xv), iv_df$player, iv_df$vill_round)
  if (all(!is.finite(res)) || sd(res[is.finite(res)], na.rm = TRUE) == 0) return(NULL)
  res
})
names(ctrl_w_df) <- ctrl_all
ctrl_w_df <- ctrl_w_df[!sapply(ctrl_w_df, is.null)]

# Assemble the working regression data frame.
#   y = demeaned contribution
#   e = demeaned lagged peer mean (endogenous regressor)
#   Z = demeaned shift–share instrument
#   group = group ID (for clustering)
reg <- data.frame(y = y_w, e = e_w, Z = z_w, group = iv_df$group)
if (length(ctrl_w_df)) reg <- cbind(reg, as.data.frame(ctrl_w_df))
reg <- reg[complete.cases(reg[, c("y","e","Z")]), , drop = FALSE]


#### 4) 2SLS WITH CLUSTERED SEs (AER::ivreg) -----------------------------------
# Goal:
#   • Estimate the causal effect of lagged peer mean (e) on contributions (y),
#     controlling for time-varying covariates, after partialling out player and
#     village×round FE via demeaning.
#   • Use shift–share Z as the instrument and cluster-robust variance by group.

# Construct 2SLS formula:
#   • If there are no controls: y ~ e | Z
#   • With controls X:        y ~ e + X | Z + X
if (ncol(reg) == 4){
  f_iv <- y ~ e | Z
  Xnames <- character(0)
} else {
  Xnames <- setdiff(names(reg), c("y","e","Z","group"))
  f_iv <- as.formula(
    paste0("y ~ e + ", paste(Xnames, collapse = " + "),
           " | Z + ", paste(Xnames, collapse = " + "))
  )
}

# 2SLS estimation using ivreg.
m_iv  <- AER::ivreg(f_iv, data = reg)

# Cluster-robust (HC1) variance–covariance matrix, clustering by group.
vcovG <- sandwich::vcovCL(m_iv, cluster = reg$group, type = "HC1")

cat("\n=== 2SLS (player & village×round FE via demeaning) — cluster-robust coefs ===\n")
print(lmtest::coeftest(m_iv, vcov = vcovG))

cat("\n=== Diagnostics (Weak-IV F, Wu–Hausman, Sargan if over-ID) ===\n")
print(summary(m_iv, diagnostics = TRUE, vcov. = vcovG))


#### 5) CLUSTER-ROBUST FIRST STAGE ---------------------------------------------
# Goal:
#   • Report first-stage regression e ~ Z + X (on demeaned variables),
#     cluster-robust SEs, and a Kleibergen–Paap-type F-test on Z.

# First-stage specification (on the same demeaned reg data).
fs_form <- as.formula(
  paste0("e ~ Z", if (length(Xnames)) paste0(" + ", paste(Xnames, collapse = " + ")) else "")
)
fs_lm <- lm(fs_form, data = reg)

# Cluster-robust VCOV for first-stage.
vcov_fs <- sandwich::vcovCL(fs_lm, cluster = reg$group, type = "HC1")

cat("\n=== First-stage (cluster-robust) ===\n")
print(lmtest::coeftest(fs_lm, vcov = vcov_fs))

# KP-style F-test on the instrument(s), if car is available.
if (requireNamespace("car", quietly = TRUE)) {
  cat("\nKP-style F on instrument(s):\n")
  print(car::linearHypothesis(fs_lm, c("Z = 0"), vcov. = vcov_fs, test = "F"))
}


#### 6) PLACEBO: EARLIEST ROUND WITH NON-ZERO INSTRUMENT VARIATION --------------
# Goal:
#   • Identify the earliest round with non-degenerate Z_shiftshare variation.
#   • Run a cross-sectional placebo regression of contributing on Z_shiftshare
#     in that round (demeaned by village), clustered by village.
#   • This checks for spurious cross-sectional relationships at early rounds.

# Within-round standard deviation of Z_shiftshare, by round_n.
ss$Z_var <- ave(ss$Z_shiftshare, ss$round_n, FUN = function(x) sd(x, na.rm = TRUE))

# Earliest round with finite Z_shiftshare and positive within-round variation.
r_placebo <- min(ss$round_n[is.finite(ss$Z_shiftshare) & ss$Z_var > 0], na.rm = TRUE)

# Take one observation per player in the placebo round with complete data.
placebo <- ss %>%
  dplyr::filter(round_n == r_placebo) %>%
  dplyr::distinct(player, .keep_all = TRUE) %>%
  dplyr::filter(complete.cases(contributing, Z_shiftshare, village_code))

# Within-village fixed effects by demeaning contributing and Z_shiftshare,
# then regress without intercept, clustering by village_code.
y_tilde <- with(placebo, contributing  - ave(contributing,  village_code, FUN = mean))
z_tilde <- with(placebo, Z_shiftshare - ave(Z_shiftshare, village_code, FUN = mean))
m_pl <- lm(y_tilde ~ z_tilde - 1, data = placebo)
vc_pl <- sandwich::vcovCL(m_pl, cluster = placebo$village_code, type = "HC1")

cat("\n=== Placebo (round ", r_placebo, ", demeaned by village) ===\n", sep = "")
print(lmtest::coeftest(m_pl, vcov = vc_pl))


#### 7) OVER-ID VARIANT: TWO SHIFT–SHARE IVs -----------------------------------
# Goal:
#   • Construct two distinct shift–share instruments Z1 and Z2 using disjoint
#     sets of baseline traits, yielding an over-identified IV system.
#   • Re-run the FE–demeaned 2SLS with both IVs and report over-ID tests.

# Split the baseline traits into two disjoint groups for Z1 and Z2.
traits1 <- c("male","married","b0100","rel_protestant","rel_none","b0200")
traits2 <- setdiff(base_z, traits1)

# Build Z1, Z2 *inside ss* using the same share × shifter structure, but
# with separate trait subsets.
S1 <- as.matrix(ss[, paste0("share_", traits1), drop = FALSE])
M1 <- as.matrix(ss[, paste0("mu_"   , traits1), drop = FALSE])
S2 <- as.matrix(ss[, paste0("share_", traits2), drop = FALSE])
M2 <- as.matrix(ss[, paste0("mu_"   , traits2), drop = FALSE])
ss$Z1 <- rowSums(S1 * M1, na.rm = TRUE)
ss$Z2 <- rowSums(S2 * M2, na.rm = TRUE)

# Join Z1 and Z2 into iv_df so that rows align exactly; keep only finite Z1, Z2.
iv_df2 <- iv_df %>%
  dplyr::left_join(ss %>% dplyr::select(group, player, vill_round, Z1, Z2),
                   by = c("group","player","vill_round")) %>%
  dplyr::filter(is.finite(Z1), is.finite(Z2))

# Apply the same two-way demeaning transformation to outcome, regressor,
# and both instruments under the same FE structure.
y2  <- tw_demean(iv_df2$contributing,  iv_df2$player, iv_df2$vill_round)
e2  <- tw_demean(iv_df2$peer_mean_tm1, iv_df2$player, iv_df2$vill_round)
Z1w <- tw_demean(iv_df2$Z1,           iv_df2$player, iv_df2$vill_round)
Z2w <- tw_demean(iv_df2$Z2,           iv_df2$player, iv_df2$vill_round)

# Two-way demean the controls again for this subsample.
ctrl_w2 <- lapply(ctrl_all, function(vn){
  xv <- iv_df2[[vn]]; if (is.null(xv)) return(NULL)
  res <- tw_demean(as.numeric(xv), iv_df2$player, iv_df2$vill_round)
  if (all(!is.finite(res)) || sd(res[is.finite(res)], na.rm = TRUE) == 0) return(NULL)
  res
})
names(ctrl_w2) <- ctrl_all
ctrl_w2 <- ctrl_w2[!sapply(ctrl_w2, is.null)]

# Regression frame for the two-IV specification.
reg2 <- data.frame(y = y2, e = e2, Z1 = Z1w, Z2 = Z2w, group = iv_df2$group)
if (length(ctrl_w2)) reg2 <- cbind(reg2, as.data.frame(ctrl_w2))
reg2 <- reg2[complete.cases(reg2), , drop = FALSE]

# Over-identified 2SLS:
#   • If no controls:  y ~ e | Z1 + Z2
#   • With controls X: y ~ e + X | Z1 + Z2 + X
if (ncol(reg2) == 5) {
  f_iv2 <- y ~ e | Z1 + Z2
  Xnames2 <- character(0)
} else {
  Xnames2 <- setdiff(names(reg2), c("y","e","Z1","Z2","group"))
  f_iv2 <- as.formula(
    paste0("y ~ e + ", paste(Xnames2, collapse = " + "),
           " | Z1 + Z2 + ", paste(Xnames2, collapse = " + "))
  )
}

# Estimate two-IV 2SLS and cluster-robust VCOV.
m_iv2 <- AER::ivreg(f_iv2, data = reg2)
vc2   <- sandwich::vcovCL(m_iv2, cluster = reg2$group, type = "HC1")

cat("\n=== 2SLS with TWO shift–share IVs (over-ID) — cluster-robust ===\n")
print(lmtest::coeftest(m_iv2, vcov = vc2))
cat("\nDiagnostics (Weak-IV, Wu–Hausman, Sargan):\n")
print(summary(m_iv2, diagnostics = TRUE, vcov. = vc2))

# Cluster-robust first stage for the two-IV specification, with both instruments.
fs2_form <- as.formula(
  paste0("e ~ Z1 + Z2", if (length(Xnames2)) paste0(" + ", paste(Xnames2, collapse = " + ")) else "")
)
fs2 <- lm(fs2_form, data = reg2)
vcfs2 <- sandwich::vcovCL(fs2, cluster = reg2$group, type = "HC1")

cat("\n=== First stage (two IVs, cluster-robust) ===\n")
print(lmtest::coeftest(fs2, vcov = vcfs2))
if (requireNamespace("car", quietly = TRUE)) {
  cat("\nKP-style joint F on Z1,Z2:\n")
  print(car::linearHypothesis(fs2, c("Z1 = 0", "Z2 = 0"), vcov. = vcfs2, test = "F"))
}


#### 8) OPTIONAL FIRST-STAGE PLOTS (demeaned) -----------------------------------
# Goal:
#   • Provide diagnostic plots for the first-stage relationships after demeaning.
#   • (a) e vs Z in the one-IV case.
#   • (b) e vs fitted values from Z1, Z2 in the two-IV case.

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  # One-IV first stage: scatter + fitted line of e on Z (both demeaned).
  ggplot(reg, aes(x = Z, y = e)) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = "lm", linewidth = 0.6, se = TRUE) +
    labs(title = "First stage (demeaned): e ~ Z",
         x = "Z (demeaned)",
         y = "peer_mean_tm1 (demeaned)") +
    theme_minimal()
  
  # Two-IV first stage: scatter + fitted line of e on predicted values from Z1 and Z2.
  ehat <- fitted(lm(e ~ Z1 + Z2, data = reg2))
  ggplot(data.frame(e = reg2$e, ehat = ehat), aes(x = ehat, y = e)) +
    geom_point(alpha = 0.15) +
    geom_smooth(method = "lm", linewidth = 0.6, se = TRUE) +
    labs(title = "First stage fit (demeaned): e ~ Z1 + Z2",
         x = "Fitted from (Z1,Z2)",
         y = "e (demeaned)") +
    theme_minimal()
}
