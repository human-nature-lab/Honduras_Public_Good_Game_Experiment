## ================== STRONGER IVs + DIAGNOSTICS ==========================
## Load required packages, install any that are missing, and set seed
pkgs <- c("dplyr","tidyr","AER","sandwich","lmtest","glmnet","Matrix")
new  <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))
set.seed(123)

## Import main panel dataset
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

## Helper: two-way demeaning (player FE + village×round FE)
## Removes player and village×round means and adds back grand mean
tw_demean <- function(x, f1, f2){
  mu <- mean(x, na.rm = TRUE)
  m1 <- ave(x, f1, FUN = function(z) mean(z, na.rm = TRUE))
  m2 <- ave(x, f2, FUN = function(z) mean(z, na.rm = TRUE))
  x - m1 - m2 + mu
}

library(dplyr)
library(tidyr)

## -----------------------------------------------------------------------
## Construct leave-one-out peer mean in period t and deeper lags
## -----------------------------------------------------------------------

## First pass: leave-one-out peer mean (within group and round)
dat_iv <- pgg_data %>%
  arrange(group, round_n) %>%
  group_by(group) %>%
  mutate(
    # Leave-one-out peer mean at time t (within group)
    peer_mean_t = ifelse(
      n() > 1,
      (sum(contributing, na.rm = TRUE) - contributing) / (n() - 1),
      NA_real_
    )
  ) %>%
  ungroup()

## Second pass: create lags t-2 and t-3 of peer_mean_t within group
dat_iv <- dat_iv %>%
  arrange(group, round_n) %>%
  group_by(group) %>%
  mutate(
    peer_mean_tm2 = dplyr::lag(peer_mean_t, 2),
    peer_mean_tm3 = dplyr::lag(peer_mean_t, 3)
  ) %>%
  ungroup()

## Quick check that peer_mean_t and deeper lags are present
head(dat_iv)

## -----------------------------------------------------------------------
## 0) Ensure predetermined traits exist and are numeric/binary as needed
## -----------------------------------------------------------------------
dat_iv <- dat_iv %>%
  mutate(
    ## Binary male indicator from various gender codings
    male = case_when(
      is.numeric(gender) ~ as.integer(gender == 1),
      is.character(gender) ~ as.integer(tolower(gender) %in% c("1","male","man","m")),
      is.factor(gender) ~ as.integer(tolower(as.character(gender)) %in% c("1","male","man","m")),
      TRUE ~ NA_integer_
    ),
    ## Binary married indicator from marital_status
    married = case_when(
      is.numeric(marital_status) ~ as.integer(marital_status == 1),
      is.character(marital_status) ~ as.integer(tolower(marital_status) %in% c("1","married","civil union","union")),
      is.factor(marital_status) ~ as.integer(tolower(as.character(marital_status)) %in% c("1","married","civil union","union")),
      TRUE ~ NA_integer_
    ),
    ## Religion: Protestant vs other
    rel_protestant = case_when(
      is.numeric(b0600) ~ as.integer(b0600 == 1),
      is.character(b0600) ~ as.integer(tolower(b0600) %in% c("1","protestant")),
      is.factor(b0600) ~ as.integer(tolower(as.character(b0600)) %in% c("1","protestant")),
      TRUE ~ NA_integer_
    ),
    ## Religion: none vs other
    rel_none = case_when(
      is.numeric(b0600) ~ as.integer(b0600 == 0),
      is.character(b0600) ~ as.integer(tolower(b0600) %in% c("0","none","no religion","nonreligious")),
      is.factor(b0600) ~ as.integer(tolower(as.character(b0600)) %in% c("0","none","no religion","nonreligious")),
      TRUE ~ NA_integer_
    ),
    ## Convert remaining traits to numeric (robust to character storage)
    age = suppressWarnings(as.numeric(age)),
    friends = suppressWarnings(as.numeric(friends)),
    adversaries = suppressWarnings(as.numeric(adversaries)),
    FI = suppressWarnings(as.numeric(FI)),
    b0100 = suppressWarnings(as.numeric(b0100)),
    b0200 = suppressWarnings(as.numeric(b0200))
  )

## List of baseline traits used in shift–share constructions
base_z <- c("male","age","friends","adversaries","FI","married","b0100",
            "rel_protestant","rel_none","b0200")

## -----------------------------------------------------------------------
## 1) Rebuild peer means: contemporaneous LOO and deeper lags
## -----------------------------------------------------------------------
dat_iv <- dat_iv %>%
  arrange(group, round_n) %>%
  group_by(group, round_n) %>%
  mutate(
    grp_n = dplyr::n(),
    ## Leave-one-out peer mean c̄_{-i,gt} with protection against n=1
    peer_mean_t = ifelse(
      grp_n > 1,
      (sum(contributing, na.rm = TRUE) - contributing) / pmax(grp_n - 1, 1),
      NA_real_
    )
  ) %>% ungroup() %>%
  arrange(group, round_n) %>%
  group_by(group) %>%
  mutate(
    ## Lags of peer mean within group
    peer_mean_tm1 = dplyr::lag(peer_mean_t, 1),
    peer_mean_tm2 = dplyr::lag(peer_mean_t, 2),
    peer_mean_tm3 = dplyr::lag(peer_mean_t, 3)
  ) %>% ungroup()

## -----------------------------------------------------------------------
## 2) Build group LOO shares s^(k) (time-invariant within group)
## -----------------------------------------------------------------------

## One row per (group, player) with a single value for each baseline trait
gmem <- dat_iv %>%
  group_by(group, player) %>%
  summarise(
    across(
      all_of(base_z),
      ~ { v <- na.omit(.); if (length(v)==0) NA_real_ else as.numeric(v[1]) }
    ),
    .groups = "drop"
  )

## Group totals and sizes for each trait
gsum <- gmem %>%
  group_by(group) %>%
  summarise(
    across(all_of(base_z), ~ sum(.x, na.rm = TRUE), .names = "sum_{.col}"),
    G = dplyr::n(),   # group size
    .groups = "drop"
  )

## Leave-one-out group means (shares) for each trait
gshare <- gmem %>%
  pivot_longer(all_of(base_z), names_to = "var", values_to = "val") %>%
  left_join(
    gsum %>% pivot_longer(starts_with("sum_"), names_to = "sumvar", values_to = "sumval") %>%
      mutate(var = sub("^sum_","",sumvar)) %>%
      dplyr::select(group, var, sumval, G),
    by = c("group","var")
  ) %>%
  mutate(
    share = (sumval - val) / pmax(G - 1, 1),
    share_name = paste0("share_", var)
  ) %>%
  dplyr::select(group, player, share_name, share) %>%
  pivot_wider(names_from = share_name, values_from = share)

## -----------------------------------------------------------------------
## 3) Build outside-village shifters μ_{−v,t-1} and lag them
## -----------------------------------------------------------------------

## Totals across all villages by round
totals <- dat_iv %>%
  group_by(round_n) %>%
  summarise(
    across(all_of(base_z), ~ sum(.x, na.rm = TRUE), .names = "sum_{.col}"),
    across(all_of(base_z), ~ sum(is.finite(.x)), .names = "n_{.col}"),
    .groups = "drop"
  )

## Totals within village by round
vill_tot <- dat_iv %>%
  group_by(round_n, village_code) %>%
  summarise(
    across(all_of(base_z), ~ sum(.x, na.rm = TRUE), .names = "sumv_{.col}"),
    across(all_of(base_z), ~ sum(is.finite(.x)), .names = "nv_{.col}"),
    .groups = "drop"
  ) %>%
  left_join(totals, by = "round_n")

## Compute leave-one-village averages for each trait
for (v in base_z) {
  vill_tot[[paste0("mux_", v)]] <-
    (vill_tot[[paste0("sum_", v)]] - vill_tot[[paste0("sumv_", v)]]) /
    pmax(vill_tot[[paste0("n_", v)]] - vill_tot[[paste0("nv_", v)]], 1)
}

## Lag these shifters by one round within village
vill_shifts_lag <- vill_tot %>%
  arrange(village_code, round_n) %>%
  group_by(village_code) %>%
  mutate(across(starts_with("mux_"), ~ dplyr::lag(.x))) %>%
  ungroup() %>%
  dplyr::select(village_code, round_n, starts_with("mux_"))

## -----------------------------------------------------------------------
## 4) Join shares & shifters; construct shift–share instrument Z_LOV
## -----------------------------------------------------------------------
dat_iv <- dat_iv %>%
  left_join(gshare, by = c("group","player")) %>%
  left_join(vill_shifts_lag, by = c("village_code","round_n"))

share_mat <- as.matrix(dat_iv[, paste0("share_", base_z), drop = FALSE])
mux_mat   <- as.matrix(dat_iv[, paste0("mux_",   base_z), drop = FALSE])
dat_iv$Z_LOV <- rowSums(share_mat * mux_mat, na.rm = TRUE)

## -----------------------------------------------------------------------
## 5) Build estimation frame with FE keys
## -----------------------------------------------------------------------
iv_base <- dat_iv %>%
  transmute(
    group,
    player = as.factor(player),
    vill_round = interaction(village_code, round_n, drop = TRUE),
    y  = contributing,      # outcome: contribution in Lempiras
    e  = peer_mean_tm1,     # endogenous regressor: lagged peer mean
    Z_tm2 = peer_mean_tm2,  # internal deeper-lag instrument
    Z_tm3 = peer_mean_tm3,  # third lag (used in CF–IV)
    Z_LOV = Z_LOV           # external shift–share instrument
  ) %>%
  filter(complete.cases(y, e, Z_tm2, Z_LOV, group, player, vill_round))

## 3) Two-way demeaning (player and village×round FE)
y_w  <- tw_demean(iv_base$y,  iv_base$player, iv_base$vill_round)
e_w  <- tw_demean(iv_base$e,  iv_base$player, iv_base$vill_round)
z2_w <- tw_demean(iv_base$Z_tm2, iv_base$player, iv_base$vill_round)
z3_w <- if(all(is.na(iv_base$Z_tm3))) NULL else tw_demean(iv_base$Z_tm3, iv_base$player, iv_base$vill_round)
zL_w <- tw_demean(iv_base$Z_LOV,  iv_base$player, iv_base$vill_round)

## Estimation frame for IV regressions on demeaned variables
reg_tm2 <- data.frame(y = y_w, e = e_w, Z2 = z2_w, ZL = zL_w, group = iv_base$group)
reg_tm2 <- reg_tm2[complete.cases(reg_tm2), , drop = FALSE]

## -----------------------------------------------------------------------
## 4) 2SLS (just-ID using Z_{t-2})
## -----------------------------------------------------------------------
m_tm2   <- AER::ivreg(y ~ e | Z2, data = reg_tm2)
vc_tm2  <- sandwich::vcovCL(m_tm2, cluster = reg_tm2$group, type = "HC1")
cat("\n=== 2SLS: instrument = peer_mean_{t-2} (demeaned) ===\n")
print(lmtest::coeftest(m_tm2, vcov = vc_tm2))

## Cluster-robust first stage & KP-like F (approximate, via car::linearHypothesis)
fs_tm2  <- lm(e ~ Z2, data = reg_tm2)
vc_fs2  <- sandwich::vcovCL(fs_tm2, cluster = reg_tm2$group, type = "HC1")
cat("\nFirst stage (cluster-robust):\n"); print(lmtest::coeftest(fs_tm2, vcov = vc_fs2))
if (requireNamespace("car", quietly = TRUE)) {
  cat("KP-style F on Z2:\n"); print(car::linearHypothesis(fs_tm2, "Z2 = 0", vcov. = vc_fs2, test = "F"))
}

## -----------------------------------------------------------------------
## 5) Over-ID spec: Z_{t-2} + LOV shift–share
## -----------------------------------------------------------------------
m_over  <- AER::ivreg(y ~ e | Z2 + ZL, data = reg_tm2)
vc_over <- sandwich::vcovCL(m_over, cluster = reg_tm2$group, type = "HC1")
cat("\n=== 2SLS: instruments = {Z_{t-2}, Z_LOV} (demeaned) ===\n")
print(lmtest::coeftest(m_over, vcov = vc_over))

cat("\nDiagnostics (AER; uses non-clustered formulas internally for some tests):\n")
print(summary(m_over, diagnostics = TRUE, vcov. = vc_over))

## First-stage regression for over-identified specification
fs_over <- lm(e ~ Z2 + ZL, data = reg_tm2)
vc_fso  <- sandwich::vcovCL(fs_over, cluster = reg_tm2$group, type = "HC1")
cat("\nFirst stage (two IVs, cluster-robust):\n"); print(lmtest::coeftest(fs_over, vcov = vc_fso))
if (requireNamespace("car", quietly = TRUE)) {
  cat("KP-style joint F on {Z2, ZL}:\n"); print(car::linearHypothesis(fs_over, c("Z2 = 0","ZL = 0"), vcov. = vc_fso, test = "F"))
}

## -----------------------------------------------------------------------
## 6) Weak-IV robust AR CI for the over-ID spec
##     (heteroskedastic-robust, non-clustered, using ivmodel)
## -----------------------------------------------------------------------
if (!requireNamespace("ivmodel", quietly = TRUE)) {
  install.packages("ivmodel")
}
library(ivmodel)

## Build Y, D, Z from the demeaned over-ID frame
Y <- reg_tm2$y                    # outcome (two-way demeaned)
D <- reg_tm2$e                    # endogenous regressor
Z <- as.matrix(reg_tm2[, c("Z2","ZL")])  # instruments (demeaned Z_{t-2}, Z_LOV)

## Fit IV model object for ivmodel
iv_obj <- ivmodel::ivmodel(Y = Y, D = D, Z = Z)

## Anderson–Rubin 95% CI (heteroskedastic-robust, non-clustered)
ar_out <- ivmodel::AR.test(iv_obj, alpha = 0.05)

cat("\n=== Heteroskedastic-robust Anderson–Rubin 95% CI (non-clustered) ===\n")
print(ar_out$ci)
cat(ar_out$ci.info, "\n")

## -----------------------------------------------------------------------
## 7) Cross-fitted “optimal IV” (CF–IV) for e = peer_mean_{t-1}
##     Use Z_LOV, Z_{t-2}, Z_{t-3} as predictors for e
## -----------------------------------------------------------------------
reg_cf <- reg_tm2 %>%
  mutate(Z3 = if(!is.null(z3_w)) z3_w else NA_real_)

## Matrix of exogenous predictors for CF–IV (lasso/ridge)
exog_mat <- as.matrix(reg_cf[, c("Z2","ZL","Z3")])
exog_mat[,"Z3"][!is.finite(exog_mat[,"Z3"])] <- 0  # allow missing Z3

K <- 5  # number of folds based on village×round grouping
fold_ids <- as.integer(as.factor(iv_base$vill_round[complete.cases(iv_base$y,iv_base$e,iv_base$Z_tm2)]))
## Map many vill_rounds into K folds
u_fold <- unique(fold_ids); map <- setNames(sample(rep(1:K, length.out = length(u_fold))), u_fold)
fold_k  <- map[as.character(fold_ids)]

## Cross-fitted predicted e
Z_cf <- rep(NA_real_, nrow(reg_cf))
for (k in 1:K) {
  trn <- which(fold_k != k); tst <- which(fold_k == k)
  ## Lasso; if it fails, fallback to ridge
  fit <- tryCatch(
    glmnet::cv.glmnet(
      x = exog_mat[trn, , drop = FALSE],
      y = reg_cf$e[trn], alpha = 1, nfolds = 5,
      standardize = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- glmnet::cv.glmnet(
      x = exog_mat[trn, , drop = FALSE],
      y = reg_cf$e[trn], alpha = 0, nfolds = 5,
      standardize = TRUE
    )
  }
  Z_cf[tst] <- as.numeric(predict(fit, newx = exog_mat[tst, , drop = FALSE], s = "lambda.min"))
}

## CF–IV 2SLS and first-stage diagnostics
if (sd(Z_cf, na.rm = TRUE) > 0) {
  reg_cf2 <- reg_cf %>% mutate(Zcf = Z_cf) %>% filter(is.finite(Zcf))
  m_cf  <- AER::ivreg(y ~ e | Zcf, data = reg_cf2)
  vc_cf <- sandwich::vcovCL(m_cf, cluster = reg_cf2$group, type = "HC1")
  cat("\n=== 2SLS: instrument = Cross-fitted predicted e (CF–IV) ===\n")
  print(lmtest::coeftest(m_cf, vcov = vc_cf))
  fs_cf <- lm(e ~ Zcf, data = reg_cf2)
  vc_fs_cf <- sandwich::vcovCL(fs_cf, cluster = reg_cf2$group, type = "HC1")
  cat("\nFirst stage (CF–IV, cluster-robust):\n"); print(lmtest::coeftest(fs_cf, vcov = vc_fs_cf))
  if (requireNamespace("car", quietly = TRUE)) {
    cat("KP-style F on Zcf:\n"); print(car::linearHypothesis(fs_cf, "Zcf = 0", vcov. = vc_fs_cf, test = "F"))
  }
} else {
  cat("\nCF–IV had near-zero variance; skip. (Try K=3 or include additional predetermined features.)\n")
}
## ================== END: STRONGER IVs ==================================

## -----------------------------------------------------------------------
## Baseline balance of instruments on covariates (levels and FE-demeaned)
## -----------------------------------------------------------------------

## Build covariate frame aligned with IV sample
Xbal <- dat_iv %>%
  transmute(
    player = as.factor(player),
    vill_round = interaction(village_code, round_n, drop = TRUE),
    group,
    male, age, friends, adversaries, FI, married, b0100, rel_protestant, rel_none, b0200
  )

bal_df <- iv_base %>%
  dplyr::select(player, vill_round, group, Z2 = Z_tm2, ZL = Z_LOV) %>%
  left_join(Xbal, by = c("player","vill_round","group")) %>%
  filter(complete.cases(.))

## Demean instruments by player and vill_round for FE-style balance check
demean_by <- function(v) tw_demean(v, bal_df$player, bal_df$vill_round)
Z2w <- demean_by(bal_df$Z2); ZLw <- demean_by(bal_df$ZL)

## Design matrix for baseline covariates
Xmat <- model.matrix(
  ~ male + age + friends + adversaries + FI + married +
    b0100 + rel_protestant + rel_none + b0200,
  data = bal_df
)[, -1]

## Balance test for Z2 (just-identified instrument)
bal_Z <- lm(Z2w ~ Xmat)
vcZ   <- sandwich::vcovCL(bal_Z, cluster = bal_df$group, type = "HC1")

## Coefficient names are prefixed with "Xmat"
tested <- paste0("Xmat", colnames(Xmat))
if (requireNamespace("car", quietly = TRUE)) {
  print(car::linearHypothesis(bal_Z, tested, vcov. = vcZ, test = "F"))
}

## Balance test for Z_LOV (shift–share instrument)
bal_L <- lm(ZLw ~ Xmat)
vcL   <- sandwich::vcovCL(bal_L, cluster = bal_df$group, type = "HC1")

tested_L <- paste0("Xmat", colnames(Xmat))  # prefix with "Xmat"
if (requireNamespace("car", quietly = TRUE)) {
  print(car::linearHypothesis(bal_L, tested_L, vcov. = vcL, test = "F"))
}

## -----------------------------------------------------------------------
## Permutation test for first-stage strength of Z_{t-2}
## -----------------------------------------------------------------------
B <- 500L
perm_F <- numeric(B)
for (b in 1:B) {
  ## Permute Z2 within vill_round cells to break structural link
  Zperm <- ave(
    reg_tm2$Z2,
    iv_base$vill_round[complete.cases(iv_base$y,iv_base$e,iv_base$Z_tm2)],
    FUN = function(z) sample(z, length(z))
  )
  fs_b  <- lm(e ~ Zperm, data = transform(reg_tm2, Zperm = Zperm))
  vc_b  <- sandwich::vcovCL(fs_b, cluster = reg_tm2$group, type = "HC1")
  if (requireNamespace("car", quietly = TRUE)) {
    perm_F[b] <- car::linearHypothesis(fs_b, "Zperm = 0", vcov.=vc_b, test="F")$`F`[2]
  } else perm_F[b] <- NA_real_
}
cat("\nPermutation p-value for first-stage (Z_{t-2}):\n")
print(mean(
  !is.na(perm_F) &
    perm_F >= car::linearHypothesis(fs_tm2, "Z2=0", vcov.=vc_fs2, test="F")$`F`[2]
))

## -----------------------------------------------------------------------
## ========== RCT SEEDING PLAN (BLOCKED BY VILLAGE), ANALYSIS FOR t=2 ==========
## -----------------------------------------------------------------------
pkgs2 <- c("randomizr","estimatr")
new2  <- setdiff(pkgs2, rownames(installed.packages()))
if (length(new2)) install.packages(new2, dependencies = TRUE)
invisible(lapply(pkgs2, library, character.only = TRUE))

## Group-level frame for round 1 (one row per group)
g1 <- pgg_data %>%
  filter(round_n == 1) %>%
  group_by(village_code, group) %>%
  summarise(n = n(), .groups = "drop")

## Block-randomize seeding at the group level within village
set.seed(123)
g1$seed_assign <- randomizr::block_ra(blocks = g1$village_code, prob = 0.5)

## Merge assignment back; compute peer means in r1 and outcomes in r2
panel12 <- pgg_data %>%
  dplyr::select(village_code, group, player, round_n, contributing) %>%
  left_join(g1, by = c("village_code","group")) %>%
  group_by(group, round_n) %>%
  mutate(peer_mean_t = (sum(contributing) - contributing) / (n - 1)) %>%
  ungroup() %>%
  group_by(group, player) %>%
  summarise(
    peer_mean_r1 = peer_mean_t[round_n == 1][1],
    y_r2         = contributing[round_n == 2][1],
    seed_assign  = seed_assign[round_n == 1][1],
    village_code = village_code[round_n == 1][1]
  ) %>% ungroup() %>%
  filter(!is.na(y_r2), !is.na(peer_mean_r1), !is.na(seed_assign))

## ITT: seed assignment effect on round-2 contribution (clustered by group)
itt <- estimatr::lm_robust(y_r2 ~ seed_assign, clusters = group, data = panel12)
cat("\n=== ITT: Seed assignment effect on round-2 contribution ===\n")
print(summary(itt))

## First stage: assignment → peer_mean_r1
fs  <- estimatr::lm_robust(peer_mean_r1 ~ seed_assign, clusters = group, data = panel12)
cat("\n=== First stage: assignment -> peer_mean_r1 ===\n")
print(summary(fs))

## IV: effect of early peer exposure on later contribution (cluster-robust)
iv  <- estimatr::iv_robust(y_r2 ~ peer_mean_r1 | seed_assign, clusters = group, data = panel12)
cat("\n=== 2SLS: effect of early peer exposure (IV with seed assignment) ===\n")
print(summary(iv))

## -----------------------------------------------------------------------
## FE-demeaned balance: Z_{t-2}, Z_LOV vs baseline covariates
## -----------------------------------------------------------------------
## Demean each baseline covariate by (player, vill_round)
demean_col <- function(v) tw_demean(v, bal_df$player, bal_df$vill_round)

Xw_df <- bal_df %>%
  transmute(
    male_w         = demean_col(male),
    age_w          = demean_col(age),
    friends_w      = demean_col(friends),
    adversaries_w  = demean_col(adversaries),
    FI_w           = demean_col(FI),
    married_w      = demean_col(married),
    b0100_w        = demean_col(b0100),
    rel_prot_w     = demean_col(rel_protestant),
    rel_none_w     = demean_col(rel_none),
    b0200_w        = demean_col(b0200),
    group          = group
  )

## Join demeaned X's with demeaned instruments
bal_FE <- cbind.data.frame(Z2w = Z2w, ZLw = ZLw, Xw_df)

## Joint balance test for Z2w ~ demeaned X
balZ_fe <- lm(Z2w ~ . - ZLw - group, data = bal_FE)
vcZ_fe  <- sandwich::vcovCL(balZ_fe, cluster = bal_FE$group, type = "HC1")
if (requireNamespace("car", quietly = TRUE)) {
  tested <- setdiff(names(coef(balZ_fe)), "(Intercept)")
  cat("\nBalance (FE‑demeaned): does Z_{t-2} load on baseline covariates?\n")
  print(car::linearHypothesis(balZ_fe, tested, vcov. = vcZ_fe, test = "F"))
}

## Joint balance test for ZLw ~ demeaned X
balL_fe <- lm(ZLw ~ . - Z2w - group, data = bal_FE)
vcL_fe  <- sandwich::vcovCL(balL_fe, cluster = bal_FE$group, type = "HC1")
if (requireNamespace("car", quietly = TRUE)) {
  tested <- setdiff(names(coef(balL_fe)), "(Intercept)")
  cat("\nBalance (FE‑demeaned): does Z_LOV load on baseline covariates?\n")
  print(car::linearHypothesis(balL_fe, tested, vcov. = vcL_fe, test = "F"))
}

## -----------------------------------------------------------------------
## Anderson–Rubin CI via ivmodel wrapper (same interface as before)
## -----------------------------------------------------------------------
if (!requireNamespace("ivmodel", quietly = TRUE)) {
  install.packages("ivmodel")
}
library(ivmodel)

## Wrapper so ar_ci_cluster(reg_tm2, ...) returns AR CI bounds
ar_ci_cluster <- function(reg_df, alpha = 0.05,
                          y = "y", e = "e",
                          Zs = c("Z2","ZL"),
                          cluster = "group") {
  
  ## Build Y, D, Z from reg_df
  stopifnot(all(c(y, e, Zs) %in% names(reg_df)))
  
  Y <- reg_df[[y]]
  D <- reg_df[[e]]
  Z <- as.matrix(reg_df[, Zs])
  
  iv_obj <- ivmodel::ivmodel(Y = Y, D = D, Z = Z)
  ar_out <- ivmodel::AR.test(iv_obj, alpha = alpha)
  
  c(
    AR_CI_low  = ar_out$ci[1, "lower"],
    AR_CI_high = ar_out$ci[1, "upper"]
  )
}

cat("\n=== Anderson–Rubin 95% CI (heteroskedastic-robust, non-clustered) for over-ID spec ===\n")
print(ar_ci_cluster(reg_tm2, alpha = 0.05, y = "y", e = "e", Zs = c("Z2","ZL"), cluster = "group"))

## -----------------------------------------------------------------------
## Associational regression: outcome vs baseline traits (FE residualized)
## -----------------------------------------------------------------------

## Build a clean covariate frame at the observation level
covars <- c("male","age","friends","adversaries","FI","married","b0100",
            "rel_protestant","rel_none","b0200")

## Residualize outcome and covariates by village×round FE
resid_by <- function(v, key) v - ave(v, key, FUN=function(z) mean(z, na.rm=TRUE))

## 0) Rebuild iv_base2 if missing (expects dat_iv to exist)
if (!exists("iv_base2")) {
  if (!exists("dat_iv")) stop("dat_iv not found. Re-create dat_iv first, then rerun this block.")
  iv_base2 <- dat_iv %>%
    transmute(
      group,
      player      = as.factor(player),
      vill_round  = interaction(village_code, round_n, drop = TRUE),
      round_n,
      y  = contributing,
      e  = peer_mean_tm1,
      Z2 = peer_mean_tm2,
      ZL = Z_LOV
    ) %>%
    filter(complete.cases(y, e, Z2, ZL, group, player, vill_round, round_n))
}

covars <- c("male","age","friends","adversaries","FI",
            "married","b0100","rel_protestant","rel_none","b0200")

## Collapse player-level covariates for associational regressions
covars_player <- dat_iv %>%
  dplyr::group_by(player) %>%
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(covars),
      ~ { v <- stats::na.omit(.x); if (length(v) == 0) NA_real_ else v[1] }
    ),
    .groups = "drop"
  )

iv_base2$player <- as.numeric(iv_base2$player)
covars_player$player <- as.numeric(covars_player$player)

df_asso <- iv_base2 %>%
  dplyr::mutate(vr = vill_round) %>%
  dplyr::left_join(covars_player, by = "player") %>%
  dplyr::mutate(
    y_r = resid_by(y, vr),
    dplyr::across(dplyr::all_of(covars), ~ resid_by(.x, vr), .names = "{.col}_r")
  ) %>%
  dplyr::filter(complete.cases(y_r, dplyr::across(dplyr::ends_with("_r"))))

## Regress residualized y on residualized covariates; cluster by group
f_asso <- as.formula(paste("y_r ~", paste(paste0(covars, "_r"), collapse=" + ")))
m_asso <- lm(f_asso, data = df_asso)
vc_asso <- sandwich::vcovCL(m_asso, cluster = df_asso$group, type = "HC1")
print(lmtest::coeftest(m_asso, vcov = vc_asso))

## -----------------------------------------------------------------------
## Example trait-level IV (placeholder; requires real exogenous instruments)
## -----------------------------------------------------------------------
iv_map <- list(
  # "b0100"   = c("policy_offer"),         # for schooling (hypothetical)
  # "friends" = c("preassigned_seating"),  # for network exposure (hypothetical)
  # "FI"      = c("lottery_win")          # for financial index (hypothetical)
)

trait_iv <- function(df, covar, iv_names){
  stopifnot(covar %in% names(df))
  stopifnot(all(iv_names %in% names(df)))
  ## Two-way demeaning to absorb player and village×round shocks
  df2 <- df %>%
    dplyr::mutate(
      y_w = tw_demean(y,  player, vill_round),
      x_w = tw_demean(.data[[covar]], player, vill_round)
    )
  ## If trait is time-invariant within player, cannot estimate FE effect
  if (sd(df2$x_w, na.rm=TRUE) == 0) {
    message(sprintf("'%s' has no within-player variation; cannot estimate its level effect with player FE.", covar))
    return(NULL)
  }
  ## Demean instruments by the same FE structure
  for (nm in iv_names) df2[[paste0(nm,"_w")]] <- tw_demean(df2[[nm]], df2$player, df2$vill_round)
  
  keep <- c("y_w","x_w", paste0(iv_names,"_w"), "group")
  df2 <- df2[stats::complete.cases(df2[, keep]), keep, drop=FALSE]
  
  if (nrow(df2) == 0) { message("No complete rows for trait IV."); return(NULL) }
  
  rhs_iv <- paste(paste0(iv_names,"_w"), collapse = " + ")
  fml <- stats::as.formula(paste0("y_w ~ x_w | ", rhs_iv))
  fit <- AER::ivreg(fml, data=df2)
  vc  <- sandwich::vcovCL(fit, cluster = df2$group, type="HC1")
  print(lmtest::coeftest(fit, vcov = vc))
  invisible(fit)
}

## Run trait-level IVs if instruments are specified in iv_map
if (length(iv_map)) {
  for (nm in names(iv_map)) {
    cat("\n=== Trait-level IV for", nm, "===\n")
    trait_iv(iv_base2, nm, iv_map[[nm]])
  }
} else {
  message("No trait instruments provided. Supply iv_map with real exogenous instruments to identify trait effects.")
}

## -----------------------------------------------------------------------
## Moderator analysis: IV with interaction terms (local-centered moderator)
## -----------------------------------------------------------------------

library(dplyr); library(AER); library(sandwich); library(lmtest)

## Helper for z-scoring with NA- and zero-variance handling
scale_na <- function(x){
  x <- as.numeric(x); mu <- mean(x, na.rm=TRUE); sdv <- stats::sd(x, na.rm=TRUE)
  if (!is.finite(sdv) || sdv == 0) x - mu else (x - mu)/sdv
}

## 0) Rebuild iv_base2 if missing (expects dat_iv to exist)
if (!exists("iv_base2")) {
  if (!exists("dat_iv")) stop("dat_iv not found. Re-create dat_iv first, then rerun this block.")
  iv_base2 <- dat_iv %>%
    transmute(
      group,
      player      = as.factor(player),
      vill_round  = interaction(village_code, round_n, drop = TRUE),
      round_n,
      y  = contributing,
      e  = peer_mean_tm1,
      Z2 = peer_mean_tm2,
      ZL = Z_LOV
    ) %>%
    filter(complete.cases(y, e, Z2, ZL, group, player, vill_round, round_n))
}

## 1) Rebuild ivw from iv_base2 so player is present in demeaned frame
ivw <- iv_base2 %>%
  mutate(
    y_w  = tw_demean(y,  player, vill_round),
    e_w  = tw_demean(e,  player, vill_round),
    Z2_w = tw_demean(Z2, player, vill_round),
    ZL_w = tw_demean(ZL, player, vill_round)
  )

## Player-level moderator construction from pgg_data
mods3 <- pgg_data %>%
  group_by(player) %>%
  summarise(
    player_key = as.character(first(player)),
    
    ## male from gender
    male = case_when(
      is.numeric(first(gender))   ~ as.integer(first(gender) == 1),
      is.character(first(gender)) ~ as.integer(tolower(first(gender)) %in% c("1","male","man","m")),
      is.factor(first(gender))    ~ as.integer(tolower(as.character(first(gender))) %in% c("1","male","man","m")),
      TRUE ~ NA_integer_
    ),
    
    ## baseline continuous traits
    educ_raw        = suppressWarnings(as.numeric(first(b0100))),
    friends_raw     = suppressWarnings(as.numeric(first(friends))),
    adversaries_raw = suppressWarnings(as.numeric(first(adversaries))),
    age_raw         = suppressWarnings(as.numeric(first(age))),
    FI_raw          = suppressWarnings(as.numeric(first(FI))),
    b0100_raw       = suppressWarnings(as.numeric(first(b0100))),
    b0200_raw       = suppressWarnings(as.numeric(first(b0200))),
    b0600_raw       = suppressWarnings(as.numeric(first(b0600))),
    
    ## married from marital_status
    married_raw = case_when(
      is.numeric(first(marital_status)) ~ as.integer(first(marital_status) == 1),
      is.character(first(marital_status)) ~ as.integer(tolower(first(marital_status)) %in%
                                                         c("1","married","civil union","union")),
      is.factor(first(marital_status)) ~ as.integer(tolower(as.character(first(marital_status))) %in%
                                                      c("1","married","civil union","union")),
      TRUE ~ NA_integer_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    ## z-score and scale moderators
    educ_z        = scale_na(educ_raw),
    friends_z     = scale_na(friends_raw),
    adversaries_z = scale_na(adversaries_raw),
    age_z         = scale_na(age_raw),
    FI_z          = FI_raw,
    b0100_z       = scale_na(b0100_raw),
    b0200_z       = b0200_raw,
    b0600_z       = b0600_raw,
    married_z     = married_raw
  ) %>%
  dplyr::select(
    player_key, male,
    educ_z, friends_z, adversaries_z, age_z, FI_z, b0200_z, b0600_z, married_z
  )

## Join moderators to ivw and coalesce in case of duplicates
ivw <- ivw %>%
  dplyr::mutate(player_key = as.character(player)) %>%
  dplyr::left_join(mods3, by = "player_key", suffix = c("", ".mods")) %>%
  dplyr::mutate(
    male      = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("male",      "male.mods")))),
    educ_z    = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("educ_z",    "educ_z.mods")))),
    friends_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("friends_z", "friends_z.mods")))),
    adversaries_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("adversaries_z", "adversaries_z.mods")))),
    age_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("age_z", "age_z.mods")))),
    FI_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("FI_z", "FI_z.mods")))),
    b0200_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("b0200_z", "b0200_z.mods")))),
    b0600_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("b0600_z", "b0600_z.mods")))),
    married_z = dplyr::coalesce(!!!dplyr::select(., dplyr::any_of(c("married_z", "married_z.mods")))),
  ) %>%
  dplyr::select(-dplyr::ends_with(".mods"), -player_key)

## Check availability of each moderator in ivw
print(sapply(ivw[, c("male","educ_z","friends_z","adversaries_z","age_z","FI_z",
                     "b0200_z","b0600_z","married_z")],
             function(x) sum(is.finite(x))))

## IV with instrumented interaction (local-centered moderator)
iv_mod <- function(data, Mname, use_overid=TRUE){
  if (!(Mname %in% names(data))) return(NULL)
  
  ## Center moderator within village×round to avoid collinearity with FE
  Mloc <- with(data, get(Mname) -
                 ave(get(Mname), vill_round, FUN=function(z) mean(z, na.rm=TRUE)))
  
  df <- data.frame(
    y  = data$y_w, e = data$e_w, Z2 = data$Z2_w, ZL = data$ZL_w,
    Mloc = as.numeric(Mloc), grp = data$group
  )
  df <- df[is.finite(rowSums(df[, c("y","e","Z2","ZL","Mloc")])), , drop=FALSE]
  if (!nrow(df)) return(NULL)
  
  ## Interaction terms and corresponding instruments
  df$eM  <- df$e  * df$Mloc
  df$Z2M <- df$Z2 * df$Mloc
  df$ZLM <- df$ZL * df$Mloc
  if (use_overid && sd(df$ZLM, na.rm=TRUE) == 0) use_overid <- FALSE
  
  fml <- if (use_overid) y ~ e + eM | Z2 + Z2M + ZL + ZLM else y ~ e + eM | Z2 + Z2M
  fit <- tryCatch(AER::ivreg(fml, data=df), error=function(e) NULL)
  if (is.null(fit)) return(NULL)
  
  vc  <- sandwich::vcovCL(fit, cluster = df$grp, type="HC1")
  tab <- lmtest::coeftest(fit, vcov=vc)
  
  take <- function(rn, cn) if (!is.null(tab) && rn %in% rownames(tab)) unname(tab[rn, cn]) else NA_real_
  
  data.frame(
    moderator = Mname,
    beta_peer = take("e",  "Estimate"),
    se_peer   = take("e",  "Std. Error"),
    beta_int  = take("eM", "Estimate"),
    se_int    = take("eM", "Std. Error"),
    p_int     = take("eM", "Pr(>|t|)"),
    n         = nrow(df),
    row.names = NULL
  )
}

## Run IV moderator regressions for dplyr::selected moderators
all_mods <- c("male","educ_z","friends_z","adversaries_z","age_z",
              "FI_z","b0200_z","b0600_z","married_z")
res_list <- lapply(all_mods, function(m) try(iv_mod(ivw, m), silent = TRUE))
mod_results <- dplyr::bind_rows(
  Filter(function(x) is.data.frame(x) && nrow(x) > 0, res_list)
)
print(mod_results)
