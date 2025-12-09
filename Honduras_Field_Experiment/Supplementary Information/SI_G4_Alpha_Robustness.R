# =============================================================================
# Master helper and robustness script for the public-good game analysis
# -----------------------------------------------------------------------------
# This script assumes that:
#   - The main panel data are stored in `pgg_data` (one row per player-round).
#   - Structural routines (e.g. `run_extensions()`) and main results objects
#     (e.g. `results`) may already exist in the workspace, but the script is
#     robust to their absence and uses sensible fallbacks.
#   - The script generates tables and diagnostics for:
#       * Ward-cluster vs end-state confusion,
#       * HMM / Markov transitions and High-state occupation,
#       * Behavioral humility distributions,
#       * Structural parameters and proxies,
#       * Early-warning model,
#       * Dynamic state logit models,
#       * α–robustness of structural primitives and policy predictions (Sα.1–Sα.4).
# =============================================================================

# =============================================================================
# Helpers (put once, early in your tables script)
# =============================================================================

# Convenience wrapper for package management.
# - Input: character vector of package names.
# - For each package p:
#     * If p is not installed, install it with dependencies.
#     * Then load it (quietly) into the current session.
# - This is called once near the top of the script so that all downstream
#   code can assume the packages are available.
safe_lib <- function(pkgs) {
  for (p in pkgs) {
    if (!suppressWarnings(require(p, character.only = TRUE))) {
      install.packages(p, dependencies = TRUE)
      library(p, character.only = TRUE)
    }
  }
}

# Core packages used throughout:
# - dplyr, tidyr, tibble : data manipulation and reshaping.
# - gt, janitor          : table creation, formatting, and totals.
# - broom, broom.mixed   : tidying regression and mixed-model output.
# - pROC                 : ROC curve construction and AUC calculations.
# - lme4                 : mixed-effects (multilevel) models.
safe_lib(c("dplyr","tidyr","tibble","gt","janitor",
           "broom","broom.mixed","pROC","lme4"))

# Main panel data from the public-good game.
# Expected core columns (used below) include:
#   - player        : individual identifier.
#   - group         : group identifier in the game.
#   - round_n       : period index.
#   - contributing  : contribution in the round (0–endowment).
#   - village_code  : village identifier (for higher-level clustering).
# Additional covariates are used when available (demographics, networks, etc.).
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# Optional packages used only if available:
# - depmixS4     : hidden Markov models (HMM) for continuous contributions.
# - modelsummary : compact and flexible regression table formatting.
if (requireNamespace("depmixS4", quietly = TRUE)) library(depmixS4)
if (requireNamespace("modelsummary", quietly = TRUE)) library(modelsummary)

# Light-weight HTML saver for gt tables.
# - save_gt(gt_tbl, stem, outdir):
#     * Writes an HTML file named "<stem>.html" in `outdir`.
#     * Returns the filename invisibly and prints a message.
# - This helper keeps table-writing syntax short in later blocks.
if (!exists("save_gt")) {
  save_gt <- function(gt_tbl, stem, outdir = getwd(), ...) {
    fn <- file.path(outdir, paste0(stem, ".html"))
    gt::gtsave(gt_tbl, filename = fn, ...)
    message("Saved: ", fn)
    invisible(fn)
  }
}

# ---- Utilities fallbacks --------------------------------

# Compute the contribution threshold c* as the mean contribution in round 1.
# - c* is used:
#     * to define High/Low contribution states (e.g. s_t = 1 if c_t >= c*),
#     * to construct early-warning outcomes (High in the final round),
#     * as a fallback for Markov-chain state definitions when HMMs are not used.
get_c_star <- function(df) {
  df %>%
    dplyr::filter(round_n == 1) %>%
    dplyr::summarise(c_star = mean(contributing, na.rm = TRUE)) %>%
    dplyr::pull(c_star)
}

# Ensure global Markov transition matrix Pmat (2x2 Low/High) and share_hi exist.
# Priority order for constructing Pmat:
#   1) Use HMM posterior in results$hmm_post (preferred: 2-state HMM on contrib).
#   2) If that is missing but depmixS4 is installed, fit a new 2-state Gaussian HMM.
#   3) Otherwise, use a simple threshold-based 2-state Markov chain (High/Low via c*).
#
# Output (stored globally):
#   - Pmat    : 2x2 matrix with rows = current state, cols = next state.
#               States are ordered as c("Low","High").
#   - share_hi: share of time (across all player-rounds) spent in the High state.
ensure_Pmat_sharehi <- function() {
  # If Pmat already exists in the environment and is a matrix, we reuse it.
  # If share_hi is missing, it is set to NA (so that downstream code can run).
  if (exists("Pmat", inherits = TRUE) && is.matrix(Pmat)) {
    if (!exists("share_hi", inherits = TRUE)) share_hi <<- NA_real_
    return(invisible(NULL))
  }
  
  # ----- Case 1: use existing HMM posterior from results$hmm_post ----------
  # This assumes:
  #   - An HMM has already been estimated elsewhere in the analysis.
  #   - Its posterior (including state assignments and contributions) is stored
  #     as `results$hmm_post`.
  if (exists("results", inherits = TRUE) && is.list(results) && "hmm_post" %in% names(results)) {
    hmm_post <- results$hmm_post
    
    # Align contribution column name to 'contrib' if needed.
    # Some earlier code might have stored standardized contributions as contrib_z.
    if (!("contrib" %in% names(hmm_post)) && "contrib_z" %in% names(hmm_post)) {
      hmm_post$contrib <- hmm_post$contrib_z
    }
    # Align round column name to 'round' if needed.
    if (!("round" %in% names(hmm_post)) && "round_n" %in% names(hmm_post)) {
      hmm_post$round <- hmm_post$round_n
    }
    
    # Emission means by state to identify which latent state corresponds to
    # higher contributions (the “High” state in the two-state Markov chain).
    emit <- hmm_post %>%
      dplyr::group_by(state) %>%
      dplyr::summarise(mu = mean(contrib, na.rm = TRUE), .groups = "drop")
    hi_state <- emit$state[which.max(emit$mu)]
    
    # Empirical state transition probabilities from the posterior state sequence.
    # For each player, we compute transitions state_t → state_{t+1}.
    trans <- hmm_post %>%
      dplyr::arrange(player, round) %>%
      dplyr::group_by(player) %>%
      dplyr::mutate(state_next = dplyr::lead(state)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(state_next)) %>%
      dplyr::count(state, state_next, name = "n") %>%
      dplyr::group_by(state) %>% dplyr::mutate(p = n / sum(n)) %>% dplyr::ungroup()
    
    # Transition matrix (rows: current, cols: next), ordered by emission means.
    # We first create a state-by-state matrix in original labels, then re-order
    # rows and columns so that "Low" comes first and "High" second.
    P <- trans %>%
      dplyr::select(state, state_next, p) %>%
      tidyr::pivot_wider(names_from = state_next, values_from = p, values_fill = 0) %>%
      tibble::column_to_rownames("state") %>% as.matrix()
    ord <- emit$state[order(emit$mu)]  # low-emission state first, high-emission second
    P <- P[as.character(ord), as.character(ord), drop = FALSE]
    colnames(P) <- rownames(P) <- c("Low","High")
    
    Pmat <<- P
    # share_hi is the fraction of posterior states classified as the high-emission state.
    share_hi <<- mean(hmm_post$state == hi_state, na.rm = TRUE)
    return(invisible(NULL))
  }
  
  # ----- Case 2: fit a 2-state Gaussian HMM using depmixS4 -------------------
  # Used when no posterior state sequence is supplied in `results`, but the
  # depmixS4 package is installed. The HMM is estimated on standardized
  # contributions, treating each player as a separate sequence.
  if (requireNamespace("depmixS4", quietly = TRUE)) {
    pan <- pgg_data %>%
      dplyr::select(player, round_n, contributing) %>%
      dplyr::filter(is.finite(contributing)) %>%
      dplyr::mutate(player = as.character(player)) %>%
      dplyr::arrange(player, round_n) %>%
      # Standardize contributions within each round (across players).
      dplyr::group_by(round_n) %>%
      dplyr::mutate(y_z = as.numeric(scale(contributing))) %>%
      dplyr::ungroup()
    
    # ntimes: number of observations per player, required by depmixS4 for
    # multi-sequence HMMs (panel structure).
    ntimes <- as.integer(table(pan$player))
    
    # Two-state Gaussian HMM on standardized contributions.
    m_hmm  <- depmixS4::depmix(response = y_z ~ 1, data = pan, nstates = 2,
                               family = gaussian(), ntimes = ntimes)
    fit_hmm <- depmixS4::fit(m_hmm, emcontrol = depmixS4::em.control(maxit = 300), verbose = FALSE)
    post <- depmixS4::posterior(fit_hmm)
    pan$state <- post$state
    
    # Identify which state is Low vs High based on mean standardized contribution.
    emit <- pan %>%
      dplyr::group_by(state) %>%
      dplyr::summarise(mu = mean(y_z, na.rm = TRUE), .groups = "drop")
    hi_state <- emit$state[which.max(emit$mu)]
    
    # Transition matrix from latent state sequences (empirical transitions).
    trans <- pan %>%
      dplyr::arrange(player, round_n) %>%
      dplyr::group_by(player) %>%
      dplyr::mutate(state_next = dplyr::lead(state)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(!is.na(state_next)) %>%
      dplyr::count(state, state_next, name = "n") %>%
      dplyr::group_by(state) %>% dplyr::mutate(p = n / sum(n)) %>% dplyr::ungroup()
    
    P <- trans %>%
      dplyr::select(state, state_next, p) %>%
      tidyr::pivot_wider(names_from = state_next, values_from = p, values_fill = 0) %>%
      tibble::column_to_rownames("state") %>% as.matrix()
    ord <- emit$state[order(emit$mu)]
    P <- P[as.character(ord), as.character(ord), drop = FALSE]
    colnames(P) <- rownames(P) <- c("Low","High")
    
    Pmat <<- P
    share_hi <<- mean(pan$state == hi_state, na.rm = TRUE)
    return(invisible(NULL))
  }
  
  # ----- Case 3: threshold-based two-state Markov chain -----------------------
  # Fallback when no HMM machinery is available.
  # States are defined deterministically by whether contributions exceed c*
  # (the mean contribution in round 1). We then compute empirical transitions.
  c_star <- get_c_star(pgg_data)
  pan <- pgg_data %>%
    dplyr::mutate(player = as.character(player),
                  st = ifelse(contributing > c_star, "High","Low")) %>%
    dplyr::arrange(player, round_n)
  
  emp <- pan %>%
    dplyr::group_by(player) %>%
    dplyr::mutate(st_next = dplyr::lead(st)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(st_next))
  
  # Transition counts Low/High → Low/High from the empirical High/Low states.
  tab <- xtabs(~ st + st_next, data = emp)
  lev <- c("Low","High")
  full <- matrix(0, nrow = 2, ncol = 2, dimnames = list(st = lev, st_next = lev))
  if (length(dimnames(tab)) == 2) full[rownames(tab), colnames(tab)] <- tab
  
  # Row-normalize to obtain transition probabilities.
  P <- sweep(full, 1, rowSums(full), FUN = "/")
  Pmat <<- P
  # share_hi: unconditional fraction of High states in the full panel.
  share_hi <<- mean(pan$st == "High", na.rm = TRUE)
  invisible(NULL)
}

# Safe scalar null-coalesce: returns 'a' if it is a finite scalar; otherwise 'b'.
# - Used when extracting fixed effects from lmer objects where some coefficients
#   may be absent or not identifiable. This avoids errors in downstream code.
`%||%` <- function(a, b) if (!is.null(a) && length(a) == 1 && is.finite(a)) a else b

# Build per-player imitation (β_peer) and habit (β_self) slopes.
# Sources, in order of preference:
#   - Preferred: an existing lmer object `mod_rs_learn_full` with random slopes
#                on own and group lagged contributions (main heterogeneity model).
#   - Otherwise: fit a minimal learning model internally (random intercept +
#                uncorrelated random slopes for own and group lags).
#
# Returns:
#   - A tibble with columns: player, beta_peer, beta_self
#     where:
#       * beta_self captures habit / persistence on own past contributions,
#       * beta_peer captures imitation / responsiveness to peers.
ensure_behav_hum <- function() {
  # Case 1: use pre-computed behavior summary in results$behav_hum.
  # This path assumes the main analysis has already produced the per-player
  # β_peer and β_self estimates and stored them as `results$behav_hum`.
  if (exists("results") && is.list(results) && "behav_hum" %in% names(results)) {
    return(dplyr::rename(results$behav_hum, beta_peer = beta_peer, beta_self = beta_self))
  }
  
  # Case 2: use an already-fitted random-slope learning model.
  # The model `mod_rs_learn_full` is expected to be an lmer with random slopes
  # on lagged own and group contributions at the player level.
  if (exists("mod_rs_learn_full")) {
    rn <- lme4::ranef(mod_rs_learn_full)$player
    cn <- colnames(rn)
    # Identify the relevant random-slope columns (first matching column).
    own_col  <- cn[grepl("^c_own_lag",   cn)][1]
    peer_col <- cn[grepl("^c_group_lag", cn)][1]
    fx <- lme4::fixef(mod_rs_learn_full)
    fx_own  <- unname(fx[grep("^c_own_lag",   names(fx))[1]])
    fx_peer <- unname(fx[grep("^c_group_lag", names(fx))[1]])
    
    # β_self and β_peer are fixed effect + player-specific random deviation.
    # If for some reason the relevant column is missing, set the deviation to 0.
    return(tibble::tibble(
      player     = rownames(rn),
      beta_self  = (fx_own  %||% 0) + if (!is.na(own_col))  rn[[own_col]]  else 0,
      beta_peer  = (fx_peer %||% 0) + if (!is.na(peer_col)) rn[[peer_col]] else 0
    ))
  } else {
    # Case 3: minimal learning model fitted within this script.
    # The model regresses contributions on:
    #   - c_own_lag   : own lagged contribution,
    #   - c_group_lag : lagged mean of peers in the group (excluding self),
    # with player-level random intercepts and (uncorrelated) random slopes.
    
    # Build panel with own and peer lags.
    dat <- pgg_data %>%
      dplyr::arrange(group, player, round_n) %>%
      dplyr::group_by(group, round_n) %>%
      dplyr::mutate(
        G = dplyr::n(),
        # Peer mean at time t excluding the focal player.
        peer_mean_t = ifelse(G > 1, (sum(contributing) - contributing)/pmax(G-1,1), NA_real_)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(c_group_lag = dplyr::lag(peer_mean_t)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(player) %>%
      dplyr::mutate(c_own_lag   = dplyr::lag(contributing)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is.finite(c_own_lag), is.finite(c_group_lag))
    
    # Random-intercept and random-slope lmer: own and group lags.
    # Random effects structure:
    #   (1 + c_own_lag + c_group_lag || player)
    # uses independent random effects for intercept and slopes at the player level.
    mod <- lme4::lmer(contributing ~ c_own_lag + c_group_lag +
                        (1 + c_own_lag + c_group_lag || player),
                      data = dat,
                      control = lme4::lmerControl(optimizer = "bobyqa",
                                                  optCtrl   = list(maxfun = 2e5)))
    rn <- lme4::ranef(mod)$player
    fx <- lme4::fixef(mod)
    
    # Extract per-player slopes: fixed effect + random deviations.
    return(tibble::tibble(
      player    = rownames(rn),
      beta_self = as.numeric(fx["c_own_lag"])    + rn[["c_own_lag"]],
      beta_peer = as.numeric(fx["c_group_lag"])  + rn[["c_group_lag"]]
    ))
  }
}

# =============================================================================
# Confusion matrix: Ward clusters × end states (LL/LH/HL/HH)
# Place: after the Ward clustering block that produced `confusion`; before ARI.
# =============================================================================

# This block summarizes how hierarchical clusters (hc_cluster) relate to
# four end-state patterns (HH, HL, LH, LL), as defined in the main analysis:
#   - A table of row percentages (cluster-by-end-state distribution),
#   - A counts-only confusion matrix with row/column totals.
if (exists("confusion")) {
  
  # Order of end states in tables (matches notation in the main analysis).
  end_order <- c("HH","HL","LH","LL")
  
  # Long format: one row per (cluster, end_state), including row totals and
  # row-wise percentages within each cluster.
  conf_long <- confusion %>%
    tidyr::pivot_longer(cols = -hc_cluster, names_to = "end_state", values_to = "n") %>%
    dplyr::mutate(end_state = factor(end_state, levels = end_order)) %>%
    dplyr::group_by(hc_cluster) %>%
    dplyr::mutate(
      row_total = sum(n, na.rm = TRUE),
      row_pct   = dplyr::if_else(row_total > 0, n / row_total, NA_real_)
    ) %>%
    dplyr::ungroup()
  
  # Wide format: counts and row percentages by (HC cluster × end state).
  # Columns:
  #   - n_HH, n_HL, ...        : counts,
  #   - row_pct_HH, row_pct_HL : cluster-specific shares.
  conf_wide <- conf_long %>%
    dplyr::select(hc_cluster, end_state, n, row_pct) %>%
    tidyr::pivot_wider(names_from = end_state,
                       values_from = c(n, row_pct),
                       values_fill = 0)
  
  # Table with row percentages by cluster (main confusion summary).
  gt_conf <- conf_wide %>%
    gt::gt() %>%
    gt::fmt_percent(columns = dplyr::starts_with("row_pct_"), decimals = 1) %>%
    gt::tab_header(title = "Confusion: Ward clusters vs. end-state transitions")
  
  save_gt(gt_conf, "Table_B2_confusion_ward_vs_endstate")
  
  # Counts-only confusion matrix with row/column totals added by janitor.
  counts_only <- confusion %>%
    dplyr::select(hc_cluster, dplyr::all_of(end_order)) %>%
    janitor::adorn_totals(where = c("row","col")) %>%
    as.data.frame()
  
  gt_counts <- counts_only %>%
    gt::gt() %>%
    gt::tab_header(title = "Confusion Matrix (counts) — Ward clusters × {HH, HL, LH, LL}")
  
  save_gt(gt_counts, "Table_B2_confusion_counts")
  
} else {
  message("[B2] Skipped: `confusion` not found.")
}

# =============================================================================
# HMM transition matrix and share in High state
# Place: Section C tables; will auto-build Pmat/share_hi if missing.
# =============================================================================

# Build Pmat / share_hi (HMM-based or Markov fallback) if not already present.
ensure_Pmat_sharehi()

if (exists("Pmat")) {
  Pdf <- as.data.frame(Pmat)
  Pdf$from <- rownames(Pdf)
  
  # Long format for gt (from-state, to-state, probability).
  P_long <- tidyr::pivot_longer(Pdf, -from, names_to = "to", values_to = "p")
  
  # Nicely formatted 2x2 transition matrix (Low/High).
  # Rows: current state, columns: next state.
  gtP <- P_long %>%
    tidyr::pivot_wider(names_from = to, values_from = p) %>%
    dplyr::select(from, dplyr::everything()) %>%
    gt::gt(rowname_col = "from") %>%
    gt::fmt_number(columns = -from, decimals = 3) %>%
    gt::tab_header(title = "HMM: transition probabilities (rows: current, cols: next)")
  save_gt(gtP, "Table_C2_HMM_transitions")
  
  # Separate table for occupation measure (share of time in High state).
  if (exists("share_hi")) {
    sh <- tibble::tibble(quantity = "Share of rounds in High state", value = as.numeric(share_hi))
    gtSH <- sh %>%
      gt::gt() %>%
      gt::fmt_number(columns = "value", decimals = 3) %>%
      gt::tab_header(title = "HMM: occupation measure")
    save_gt(gtSH, "Table_C2b_HMM_share_high")
  }
} else {
  message("[C2] Skipped: could not build `Pmat`.")
}

# =============================================================================
# Behavioral humility distributions (β_peer, β_self)
# Place: after lmer heterogeneity section (has `mod_rs_learn_full`), or this
#        will fit a minimal model if needed.
# =============================================================================

# Obtain per-player β_peer (imitation of peers) and β_self (habit / own-history),
# using whichever source is available (pre-computed results or internal model).
# The object `beh` is a tibble with one row per player.
beh <- tryCatch(ensure_behav_hum(), error = function(e) NULL)

if (!is.null(beh) && all(c("beta_peer","beta_self") %in% names(beh))) {
  # Deciles / median of imitation and habit slopes (P10, median, P90).
  q_peer <- stats::quantile(beh$beta_peer, c(.1,.5,.9), na.rm = TRUE)
  q_self <- stats::quantile(beh$beta_self, c(.1,.5,.9), na.rm = TRUE)
  
  # Shares of anti-imitators (β_peer < 0) and "strong copiers" (> P90).
  share_anti   <- mean(beh$beta_peer < 0, na.rm = TRUE)
  share_strong <- mean(beh$beta_peer > q_peer["90%"], na.rm = TRUE)
  
  tab <- tibble::tibble(
    quantity = c("β_peer P10","β_peer Median","β_peer P90",
                 "β_self P10","β_self Median","β_self P90",
                 "Share anti-imitators (β_peer<0)",
                 "Share strong copiers (β_peer>P90)"),
    value = c(q_peer[1], q_peer[2], q_peer[3],
              q_self[1], q_self[2], q_self[3],
              share_anti, share_strong)
  )
  
  gt_beh <- tab %>%
    gt::gt() %>%
    gt::fmt_number(columns = "value", decimals = 3) %>%
    gt::tab_header(title = "Behavioral humility: distributional summaries")
  
  save_gt(gt_beh, "Table_C3_behav_hum_summaries")
} else {
  message("[C3] Skipped: could not compute β_peer/β_self summaries.")
}

# =============================================================================
# Structural parameters (d_i, φ_i, c^H at cap)
# Place: after you run `results <- run_extensions(pgg_data, ...)`.
# =============================================================================

# Structural section:
# - Expects results$struct_params, containing player-level:
#     * d_hat   : altruism parameter (d_i),
#     * phi_hat : friction or decay parameter (φ_i),
#     * cH_pred : predicted high-state contribution c^H_i, etc.
# - Produces a table summarizing the distribution of d_i and notes about caps
#   and small φ_i, as used in the structural part of the paper.
if (exists("results") && is.list(results) && "struct_params" %in% names(results)) {
  sp <- results$struct_params
  
  # Share of players at contribution cap (endowment = 12) in predicted high state,
  # and share of players with low φ (≤ 0.1).
  share_at_cap <- mean(sp$cH_pred >= 12 - 1e-9, na.rm = TRUE)
  phi_near0    <- mean(sp$phi_hat <= 0.1,        na.rm = TRUE)
  
  # Summary statistics for d_i (altruism parameter) over players.
  ds <- summary(sp$d_hat)
  d_tab <- tibble::tibble(stat = names(ds), d_hat = as.numeric(ds))
  
  gt_d <- d_tab %>%
    gt::gt() %>%
    gt::fmt_number(columns = "d_hat", decimals = 3) %>%
    gt::tab_header(title = "Structural altruism d_i — summary statistics") %>%
    gt::tab_source_note(
      source_note = sprintf("c^H at cap (12L): %0.1f%%;  φ_i ≤ 0.1: %0.1f%%",
                            100*share_at_cap, 100*phi_near0)
    )
  save_gt(gt_d, "Table_C4_structural_dhat_summary")
} else {
  # Proxy/placeholder if structural estimation has not been run.
  # Uses final-round contributions as a crude proxy for structural primitives,
  # so that some descriptive information is available even without structural
  # estimation routines.
  finals <- pgg_data %>%
    dplyr::filter(round_n == max(round_n, na.rm = TRUE)) %>%
    dplyr::summarise(
      N_players = dplyr::n_distinct(player),
      mean_final = mean(contributing, na.rm = TRUE),
      sd_final   = sd(contributing,   na.rm = TRUE),
      share_at_cap = mean(contributing >= 12 - 1e-9, na.rm = TRUE)
    )
  gt_proxy <- finals %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "quantity", values_to = "value") %>%
    gt::gt() %>%
    gt::fmt_number(columns = "value", decimals = 3) %>%
    gt::tab_header(title = "Structural parameters — proxies (run `run_extensions()` for d_i, φ_i)")
  save_gt(gt_proxy, "Table_C4_structural_proxies")
}

# =============================================================================
# Early-warning model — odds ratios + AUC line
# Place: anywhere after `pgg_data`; this recomputes the GLM/AUC if missing.
# =============================================================================

# The early-warning model predicts whether a player will be in the High state
# in the final round (round 10) using:
#   - summary statistics of early contributions (rounds 1–3),
#   - baseline covariates observed at round 1.
# If `table2_coefs` already exists (from a previous run), we reuse it.
if (!exists("table2_coefs")) {
  set.seed(123)
  early_rounds <- 1:3    # Early observation window.
  final_round  <- 10     # Horizon for "High" outcome.
  c_star <- get_c_star(pgg_data)
  
  # Early-round features per player (mean, SD, and slope of contributions).
  feat_early <- pgg_data %>%
    dplyr::filter(round_n %in% early_rounds) %>%
    dplyr::group_by(player) %>%
    dplyr::summarise(
      c_mean_e = mean(contributing, na.rm = TRUE),
      c_sd_e   = stats::sd(contributing, na.rm = TRUE),
      c_r1     = contributing[round_n == min(early_rounds)][1],
      c_r3     = contributing[round_n == max(early_rounds)][1],
      c_slope  = (c_r3 - c_r1)/(max(early_rounds) - min(early_rounds)),
      .groups = "drop"
    )
  
  # Baseline covariates from round 1 (gender, age, social ties, etc.).
  # These variables are used in the form present in the dataset.
  base1 <- pgg_data %>%
    dplyr::filter(round_n == 1) %>%
    dplyr::distinct(player, village_code, friends, finauto_q_2_temp, gender, age) #, IH)
  
  # Outcome: indicator of High contribution in final round (relative to c*).
  ydat <- pgg_data %>%
    dplyr::filter(round_n == final_round) %>%
    dplyr::transmute(player, final_high = contributing > c_star)
  
  # Merge features, baseline covariates, and outcome at the player level.
  dat_ew <- feat_early %>%
    dplyr::left_join(base1, by = "player") %>%
    dplyr::left_join(ydat, by = "player") %>%
    dplyr::filter(!is.na(final_high))
  
  # 70/30 train–test split at the player level.
  idx <- sample.int(nrow(dat_ew), floor(0.7*nrow(dat_ew)))
  train <- dat_ew[idx,]; test <- dat_ew[-idx,]
  
  # Logistic regression predicting final High from early contributions and covariates.
  # Coefficients will later be exponentiated to yield odds ratios.
  ew_glm <- stats::glm(final_high ~ c_mean_e + c_sd_e + c_slope + friends + finauto_q_2_temp + gender + age, #+ IH,
                       data = train, family = binomial)
  test$pred <- stats::predict(ew_glm, newdata = test, type = "response")
  
  # ROC and AUC on holdout set (out-of-sample performance).
  roc_obj <- pROC::roc(test$final_high, test$pred)
  auc_val <<- as.numeric(pROC::auc(roc_obj))
  opt     <<- pROC::coords(roc_obj, "best", best.method = "youden",
                           ret = c("threshold","sensitivity","specificity"))
  
  # Tidy odds ratios with confidence intervals for the main table (Table C5).
  table2_coefs <<- broom::tidy(ew_glm, conf.int = TRUE, exponentiate = TRUE) %>%
    dplyr::mutate(term = dplyr::recode(term,
                                       "(Intercept)"        = "Intercept",
                                       "c_mean_e"           = "Early mean",
                                       "c_sd_e"             = "Early SD",
                                       "c_slope"            = "Early slope",
                                       "friends"            = "Friends",
                                       "finauto_q_2_temp"   = "Financial autonomy",
                                       "gender"             = "Gender (1=male)",
                                       "age"                = "Age")) # "IH" = "Intellectual humility"
}

if (exists("table2_coefs")) {
  # Format odds ratio table: OR, 95% CI, p-values.
  tab <- table2_coefs %>%
    dplyr::mutate(estimate = as.numeric(estimate),
                  conf.low = as.numeric(conf.low),
                  conf.high = as.numeric(conf.high)) %>%
    dplyr::rename(OR = estimate, lo95 = conf.low, hi95 = conf.high, p = p.value) %>%
    dplyr::mutate(across(c(OR, lo95, hi95), round, 3),
                  p = signif(p, 3))
  
  gt_ew <- tab %>%
    gt::gt() %>%
    gt::tab_header(title = "Early-warning model (GLM on final High): odds ratios") %>%
    gt::cols_label(term = "Covariate", OR = "OR",
                   lo95 = "95% CI (low)", hi95 = "95% CI (high)", p = "p-value")
  save_gt(gt_ew, "Table_C5_early_warning_ORs")
  
  # Supplemental table with AUC and associated diagnostic metrics on holdout data.
  if (exists("auc_val") || exists("opt")) {
    auc_line <- tibble::tibble(
      metric = c("AUC (holdout)", "Best threshold (Youden J)", "Sensitivity@best", "Specificity@best"),
      value  = c(if (exists("auc_val")) auc_val else NA_real_,
                 if (exists("opt"))     as.numeric(opt["threshold"])   else NA_real_,
                 if (exists("opt"))     as.numeric(opt["sensitivity"]) else NA_real_,
                 if (exists("opt"))     as.numeric(opt["specificity"]) else NA_real_)
    )
    gt_auc <- auc_line %>%
      gt::gt() %>%
      gt::fmt_number(columns = "value", decimals = 3) %>%
      gt::tab_header(title = "Early-warning diagnostic metrics")
    save_gt(gt_auc, "Table_C5b_early_warning_AUC_thresholds")
  }
} else {
  message("[C5] Skipped: `table2_coefs` not found and could not be built.")
}

# =============================================================================
# OPTIONAL: Additional models — dynamic state logit (M1)
# Place: optional; will fit models if missing.
# =============================================================================

# Dynamic state logit for binary state s_t (High vs Low) with:
#   - Lagged state s_(t-1),
#   - Lagged group-level peer mean contributions,
#   - Baseline covariates and random intercepts at village, group, and player.
# This block estimates:
#   - m1_state    : baseline random-intercept model,
#   - m1_state_rs : optional model with player-level random slopes.
if (!exists("m1_state") && !exists("m1_state_rs")) {
  c_star <- get_c_star(pgg_data)
  
  # Base panel: one row per player–round with finite contribution.
  panel <- pgg_data %>%
    dplyr::filter(is.finite(contributing)) %>%
    dplyr::mutate(player = as.character(player),
                  village_code = as.character(village_code),
                  group = as.character(group))
  
  # Within-group contemporaneous peer summary: mean contribution and share above c*.
  # After computing at time t, construct lagged versions to use as predictors.
  panel <- panel %>%
    dplyr::group_by(group, round_n) %>%
    dplyr::mutate(
      G = dplyr::n(),
      peer_mean_c_t     = ifelse(G > 1, (sum(contributing) - contributing)/(G - 1), NA_real_),
      peer_share_high_t = ifelse(G > 1, (sum(contributing >= c_star) - (contributing >= c_star))/(G - 1), NA_real_)
    ) %>% dplyr::ungroup() %>%
    dplyr::arrange(group, round_n) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      peer_mean_c_tm1     = dplyr::lag(peer_mean_c_t),
      peer_share_high_tm1 = dplyr::lag(peer_share_high_t)
    ) %>% dplyr::ungroup()
  
  # Binary state indicator s_t (High vs Low) with lag, and round index r_t.
  panel <- panel %>%
    dplyr::arrange(player, round_n) %>%
    dplyr::group_by(player) %>%
    dplyr::mutate(s = as.integer(contributing >= c_star), s_lag = dplyr::lag(s)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(r_t = as.integer(round_n))
  
  # Initial state s1 (round 1) per player (baseline state).
  s1_df <- panel %>% dplyr::filter(round_n == 1) %>% dplyr::select(player, s1 = s)
  
  # Player-level average of lagged peer means (used as baseline peer covariate).
  peerbar_df <- panel %>%
    dplyr::group_by(player) %>%
    dplyr::summarise(peer_mean_c_tm1_bar = mean(peer_mean_c_tm1, na.rm = TRUE), .groups = "drop")
  
  # Baseline player/village covariates from round 1, standardized where appropriate.
  # Many covariates capture social network structure and demographics.
  base <- pgg_data %>%
    dplyr::filter(round_n == 1) %>%
    dplyr::distinct(player, village_code, .keep_all = TRUE) %>%
    dplyr::mutate(
      player       = as.character(player),
      village_code = as.character(village_code),
      male  = as.integer(tolower(as.character(gender)) %in% c("1","male","man","m")),
      age   = suppressWarnings(as.numeric(age)),
      friends      = suppressWarnings(as.numeric(friends)),
      adversaries  = suppressWarnings(as.numeric(adversaries)),
      network_density_fr  = suppressWarnings(as.numeric(network_density_fr)),
      network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
      network_size        = suppressWarnings(as.numeric(network_size)),
      b0100               = suppressWarnings(as.numeric(b0100)),
      FI                  = suppressWarnings(as.numeric(FI)),
      access_routes       = suppressWarnings(as.numeric(access_routes)),
      b0600 = as.factor(b0600),
      b0200 = as.factor(b0200),
      marital_status = as.factor(marital_status)
    ) %>%
    # Standardize numeric covariates (z-scores) for better interpretation
    # and to avoid scaling issues in logistic mixed models.
    dplyr::mutate(dplyr::across(
      c(age, friends, adversaries, network_density_fr, network_density_adv,
        network_size, b0100, FI, access_routes),
      ~ as.numeric(scale(.x)),
      .names = "{.col}_z"
    ))
  
  # Candidate standardized baseline covariates; we intersect to keep only
  # those actually present in the dataset.
  zvars <- intersect(
    c("age_z","friends_z","adversaries_z","network_density_fr_z","network_density_adv_z",
      "network_size_z","b0100_z","FI_z","access_routes_z"),
    names(base)
  )
  
  # Final dataset for dynamic state logit.
  # Contains:
  #   - current and lagged states (s, s_lag, s1),
  #   - lagged peer mean contributions,
  #   - baseline covariates and random-effect identifiers.
  m1_df <- panel %>%
    dplyr::left_join(base %>% dplyr::select(player, village_code, dplyr::all_of(zvars),
                                            male, b0600, b0200, marital_status),
                     by = c("player","village_code")) %>%
    dplyr::left_join(s1_df,     by = "player") %>%
    dplyr::left_join(peerbar_df,by = "player") %>%
    dplyr::mutate(peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12) %>%
    dplyr::filter(is.finite(s), is.finite(s_lag),
                  is.finite(peer_mean_c_tm1_scaled),
                  is.finite(s1), is.finite(peer_mean_c_tm1_bar)) %>%
    dplyr::mutate(
      village_code = factor(village_code),
      group        = factor(group),
      player       = factor(player)
    )
  
  # Baseline dynamic state logit with random intercepts at village, group, player.
  # Dependent variable: s (state High vs Low at time t).
  m1_state <<- lme4::glmer(
    s ~ s_lag + peer_mean_c_tm1_scaled + r_t + s1 + peer_mean_c_tm1_bar +
      male + age_z + friends_z + adversaries_z +
      network_density_fr_z + network_density_adv_z + network_size_z +
      b0100_z + FI_z + marital_status.x + b0200.x + b0600.x + access_routes_z +
      (1 | village_code) + (1 | group) + (1 | player),
    data = m1_df, family = binomial,
    control = lme4::glmerControl(optimizer = "nloptwrap",
                                 calc.derivs = FALSE,
                                 optCtrl = list(maxfun = 2e4))
  )
  
  # Optional random slopes version with richer player heterogeneity.
  # Adds player-level random slopes for r_t, s_lag, and peer_mean_c_tm1_scaled.
  m1_state_rs <<- tryCatch(
    lme4::glmer(
      s ~ s_lag + peer_mean_c_tm1_scaled + r_t + s1 + peer_mean_c_tm1_bar +
        male + age_z + friends_z + adversaries_z +
        network_density_fr_z + network_density_adv_z + network_size_z +
        b0100_z + FI_z + marital_status + b0200 + b0600 + access_routes_z +
        (1 | village_code) + (1 | group) + (1 + r_t + s_lag + peer_mean_c_tm1_scaled || player),
      data = m1_df, family = binomial,
      control = lme4::glmerControl(optimizer = "nloptwrap",
                                   calc.derivs = FALSE,
                                   optCtrl = list(maxfun = 2e4))
    ),
    error = function(e) NULL
  )
}

# Export dynamic state logit models in OR form if at least one was estimated.
if (exists("m1_state") || exists("m1_state_rs")) {
  mods <- list()
  if (exists("m1_state"))    mods[["(1) Dynamic state logit (ORs)"]] <- m1_state
  if (exists("m1_state_rs")) mods[["(2) + random slopes (ORs)"]]      <- m1_state_rs
  
  if (requireNamespace("modelsummary", quietly = TRUE)) {
    # Preferred route: modelsummary formatting of logistic OR tables.
    tblD <- modelsummary::msummary(
      models = mods, output = "gt", statistic = "({std.error})",
      stars = TRUE, exponentiate = TRUE
    )
    save_gt(tblD, "Table_D1_dynamic_state_logit_ORs")
  } else {
    # Fallback: tidy fixed effects from first model only (m1_state or m1_state_rs).
    md <- mods[[1]]
    td <- broom.mixed::tidy(md, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
      dplyr::transmute(term, OR = estimate, lo95 = conf.low, hi95 = conf.high, p = p.value)
    gtD <- td %>%
      gt::gt() %>%
      gt::fmt_number(columns = c("OR","lo95","hi95","p"), decimals = 3) %>%
      gt::tab_header(title = "Dynamic state logit (ORs)")
    save_gt(gtD, "Table_D1_dynamic_state_logit_ORs")
  }
}

# =============================================================================
# α–ROBUSTNESS (Sα.1–Sα.4): Structural re-estimates at α ∈ {0.3, 0.5, 0.7}
# Paste this AFTER your Helpers and AFTER `pgg_data` exists.
# Requires your structural function `run_extensions()` to be available.
# =============================================================================

# --- small saver: both HTML and LaTeX into "tables/" -------------------------
# Wrapper to save gt tables as both HTML and LaTeX (.tex) in a ./tables folder.
# - This is used in the α–robustness section below.
if (!exists("save_gt_both")) {
  save_gt_both <- function(gt_tbl, stem, outdir = file.path(getwd(), "tables")) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    html_fn <- file.path(outdir, paste0(stem, ".html"))
    tex_fn  <- file.path(outdir, paste0(stem,  ".tex"))
    gt::gtsave(gt_tbl, filename = html_fn)
    gt::gtsave(gt_tbl, filename = tex_fn)
    message("Saved: ", basename(html_fn), " and ", basename(tex_fn))
    invisible(list(html = html_fn, tex = tex_fn))
  }
}

# --- utility: detect column by possible names --------------------------------
# Picks the first existing column among candidate names, to allow for
# flexible naming in outputs from run_extensions().
.pick_col <- function(df, candidates) {
  x <- intersect(candidates, names(df))
  if (length(x) == 0) return(NULL)
  x[1]
}

# --- utility: find an ID column (player) for alignment -----------------------
# Used when matching structural parameters across α values (e.g. 0.3 vs 0.5).
.find_idcol <- function(df) {
  .pick_col(df, c("player","player_id","id","ID","subject","Subject","participant"))
}

# --- helper: summary quantiles for a numeric vector --------------------------
# Produces common quantiles + mean and count for a vector x.
summ_vec <- function(x) {
  x <- as.numeric(x)
  qs <- stats::quantile(x, probs = c(0, .10, .25, .50, .75, .90, 1), na.rm = TRUE)
  out <- tibble::tibble(
    stat = c("Min","P10","Q1","Median","Q3","P90","Max","Mean","N"),
    value = c(unname(qs[c(1,2,3,4,5,6,7)]), mean(x, na.rm = TRUE), sum(is.finite(x)))
  )
  out
}

# --- wrapper: call your structural routine with a given α --------------------
# This assumes run_extensions(data, alpha=..., ...) exists and returns a list
# with element struct_params containing (d_hat, phi_hat, cH_pred, ...).
# The wrapper is defensive about the argument name used for α and other params,
# so it can accommodate different naming conventions in `run_extensions()`.
run_structural_for_alpha <- function(alpha,
                                     b_default = 2, kappa_default = 1,
                                     data = pgg_data) {
  if (!exists("run_extensions")) stop("`run_extensions()` not found in your session.")
  fmls <- names(formals(run_extensions))
  
  # First (unnamed) argument is the data.
  args <- list(data)
  
  # Try to pass α under any sensible argument name used by run_extensions().
  alpha_names <- c("alpha", "alpha_param", "a_param", "alpha_w", "alpha0")
  an <- intersect(alpha_names, fmls)
  if (length(an)) args[[an[1]]] <- alpha
  
  # Pass defaults for other structural parameters if they exist in the signature.
  if ("b_param"       %in% fmls) args[["b_param"]]       <- b_default
  if ("kappa_param"   %in% fmls) args[["kappa_param"]]   <- kappa_default
  if ("seed"          %in% fmls) args[["seed"]]          <- 123L
  
  res <- try(do.call(run_extensions, args), silent = TRUE)
  if (inherits(res, "try-error")) {
    stop("Could not run `run_extensions()` for alpha = ", alpha,
         ". Try exposing the α argument name in `run_extensions()`.")
  }
  res
}

# --- tidy structural primitives for one α ------------------------------------
# Extracts for a given struct_params object:
#   - Distribution of d_i and φ_i (quantiles),
#   - Share at cap (c^H = endowment) and share with φ_i ≤ 0.1 (if cH available).
# Returns:
#   - A list with tab_primitives (tidy summary table).
summarise_struct_alpha <- function(sp, data_final_round, endowment = 12) {
  # Detect key columns by flexible naming.
  dcol   <- .pick_col(sp, c("d_hat","d","d_i","dhat"))
  phicol <- .pick_col(sp, c("phi_hat","phi","phi_i","phihat"))
  cHcol  <- .pick_col(sp, c("cH_pred","cH","c_H_pred","cHhat"))
  
  if (is.null(dcol) || is.null(phicol)) {
    stop("Could not find d_hat/phi_hat in `struct_params`.")
  }
  
  # Sα.1: distribution tables for d_i and φ_i.
  tab_d   <- summ_vec(sp[[dcol]]) |>
    dplyr::mutate(quantity = "d_i") |>
    dplyr::select(quantity, stat, value)
  tab_phi <- summ_vec(sp[[phicol]]) |>
    dplyr::mutate(quantity = "phi_i") |>
    dplyr::select(quantity, stat, value)
  
  # Share measures (if predicted c^H is available).
  share_rows <- NULL
  if (!is.null(cHcol)) {
    share_at_cap <- mean(sp[[cHcol]] >= (endowment - 1e-9), na.rm = TRUE)
    phi_near0    <- mean(sp[[phicol]] <= 0.1, na.rm = TRUE)
    share_rows   <- tibble::tibble(
      quantity = c("share_at_cap_cH", "share_phi_le_0.1"),
      stat     = "Share",
      value    = c(share_at_cap, phi_near0)
    )
  }
  
  list(
    tab_primitives = dplyr::bind_rows(tab_d, tab_phi, share_rows)
  )
}

# --- rank concordance vs α = 0.5 ---------------------------------------
# Spearman / Kendall rank correlations between parameter estimates across α.
# Compares d_i and φ_i under α = 0.5 vs α ∈ {0.3, 0.7}.
rank_concord_two <- function(sp_base, sp_alt) {
  # Detect ID columns for inner join; else fall back to row order.
  idb <- .find_idcol(sp_base); ida <- .find_idcol(sp_alt)
  if (!is.null(idb) && !is.null(ida)) {
    x <- dplyr::inner_join(
      sp_base |> dplyr::select(!!idb, d_base  = dplyr::any_of(c("d_hat","d","d_i","dhat")),
                               phi_base= dplyr::any_of(c("phi_hat","phi","phi_i","phihat"))),
      sp_alt  |> dplyr::select(!!ida, d_alt   = dplyr::any_of(c("d_hat","d","d_i","dhat")),
                               phi_alt = dplyr::any_of(c("phi_hat","phi","phi_i","phihat"))),
      by = setNames(ida, idb)
    )
  } else {
    # Fallback: assume same ordering across α-runs.
    n <- min(nrow(sp_base), nrow(sp_alt))
    x <- tibble::tibble(
      d_base   = sp_base[[ .pick_col(sp_base, c("d_hat","d","d_i","dhat")) ]][seq_len(n)],
      phi_base = sp_base[[ .pick_col(sp_base, c("phi_hat","phi","phi_i","phihat")) ]][seq_len(n)],
      d_alt    = sp_alt [[ .pick_col(sp_alt,  c("d_hat","d","d_i","dhat")) ]][seq_len(n)],
      phi_alt  = sp_alt [[ .pick_col(sp_alt,  c("phi_hat","phi","phi_i","phihat")) ]][seq_len(n)]
    )
  }
  
  out <- tibble::tibble(
    parameter = c("d_i","phi_i"),
    spearman  = c(stats::cor(x$d_base,  x$d_alt,  method = "spearman", use = "complete.obs"),
                  stats::cor(x$phi_base,x$phi_alt,method = "spearman", use = "complete.obs")),
    kendall   = c(stats::cor(x$d_base,  x$d_alt,  method = "kendall",  use = "complete.obs"),
                  stats::cor(x$phi_base,x$phi_alt,method = "kendall",  use = "complete.obs"))
  )
  out
}

# --- Sα.3: light PPC for a given α -------------------------------------------
# Computes simple posterior predictive diagnostics for structural predictions:
#   - median absolute error of final c vs predicted c^H at player level,
#   - difference in means (final mean - mean c^H),
#   - RMSE of round-wise mean contributions vs predicted round means (if available).
# These are robust, low-dimensional checks on the fit of the structural model.
ppc_light <- function(res_alpha, data, endowment = 12) {
  sp <- res_alpha$struct_params
  if (is.null(sp)) stop("`struct_params` missing in results.")
  
  cHcol <- .pick_col(sp, c("cH_pred","cH","c_H_pred","cHhat"))
  idcol <- .find_idcol(sp)
  
  # Final-round empirical contributions.
  df_final <- data |>
    dplyr::filter(round_n == max(round_n, na.rm = TRUE)) |>
    dplyr::select(player, contributing)
  
  # Player-level median absolute error if we can match IDs.
  med_abs_err <- NA_real_
  if (!is.null(cHcol) && !is.null(idcol) &&
      idcol %in% names(sp) && "player" %in% names(df_final)) {
    
    tmp <- dplyr::left_join(
      df_final, sp |> dplyr::select(player = !!idcol, cH = !!cHcol),
      by = "player"
    )
    med_abs_err <- stats::median(abs(tmp$contributing - tmp$cH), na.rm = TRUE)
  }
  
  # Δ mean between final-round mean and mean predicted c^H.
  delta_final_mean <- NA_real_
  if (!is.null(cHcol)) {
    mu_final <- mean(df_final$contributing, na.rm = TRUE)
    mu_cH    <- mean(sp[[cHcol]],          na.rm = TRUE)
    delta_final_mean <- mu_final - mu_cH
  }
  
  # RMSE of round-wise mean contributions vs predicted round means (if available).
  rmse_round_means <- NA_real_
  cand <- c("ppc_round_means","pred_round_means","sim_round_means")
  fld  <- intersect(cand, names(res_alpha))
  if (length(fld) && is.data.frame(res_alpha[[fld[1]]])) {
    pr <- res_alpha[[fld[1]]]
    rn <- .pick_col(pr, c("round","round_n"))
    mp <- .pick_col(pr, c("mean_pred","mean","mu_pred"))
    if (!is.null(rn) && !is.null(mp)) {
      emp <- data |>
        dplyr::group_by(round_n) |>
        dplyr::summarise(mean_emp = mean(contributing, na.rm = TRUE), .groups = "drop")
      comp <- dplyr::inner_join(emp, pr |> dplyr::rename(round_n = !!rn, mean_pred = !!mp),
                                by = "round_n")
      rmse_round_means <- sqrt(mean((comp$mean_emp - comp$mean_pred)^2, na.rm = TRUE))
    }
  }
  
  tibble::tibble(
    metric = c("Median |c_final - c_H| (per player)", "Δ mean (final - mean c_H)", "RMSE of round means"),
    value  = c(med_abs_err, delta_final_mean, rmse_round_means)
  )
}

# --- policy invariance grid (requires a callable in results) -----------
# Evaluates how predicted c^H responds to a subsidy m, using any available
# policy evaluation function stored in the results object.
#
# The policy function is expected to accept m (or m and struct_params) and
# return either a vector or data frame with predicted c^H under that policy.
policy_grid <- function(res_alpha, m_grid = c(0, 0.25, 0.5)) {
  sp <- res_alpha$struct_params
  if (is.null(sp)) return(NULL)
  
  # Search for a policy hook in results: function that returns cH under policy m.
  f_candidates <- c("policy_eval", "predict_cH", "cH_of_m", "policy", "simulate_policy")
  f_name <- intersect(f_candidates, names(res_alpha))
  policy_fun <- NULL
  if (length(f_name) && is.function(res_alpha[[f_name[1]]])) {
    policy_fun <- res_alpha[[f_name[1]]]
  }
  
  # Normalize policy_fn outputs to a numeric vector of predicted cH.
  pull_cH <- function(m) {
    if (is.null(policy_fun)) return(NULL)
    # Try different signatures commonly used in policy functions.
    out <- try(policy_fun(m), silent = TRUE)
    if (inherits(out, "try-error")) {
      out <- try(policy_fun(m = m), silent = TRUE)
    }
    if (inherits(out, "try-error")) {
      out <- try(policy_fun(m = m, struct_params = sp), silent = TRUE)
    }
    if (inherits(out, "try-error")) return(NULL)
    
    if (is.data.frame(out)) {
      cc <- .pick_col(out, c("cH_pred","cH","c_H_pred"))
      if (!is.null(cc)) out <- out[[cc]] else out <- NULL
    }
    if (!is.null(out)) as.numeric(out) else NULL
  }
  
  # Build a small grid of policy scenarios m ∈ m_grid.
  # For each m, compute share at cap and mean predicted contributions.
  rows <- lapply(m_grid, function(m) {
    cH <- pull_cH(m)
    if (is.null(cH)) return(NULL)
    tibble::tibble(
      m = m,
      share_high = mean(cH >= 12 - 1e-9, na.rm = TRUE),
      mean_cH    = mean(cH, na.rm = TRUE)
      # Welfare measures could be added here if returned by policy_fun.
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) return(NULL)
  dplyr::bind_rows(rows)
}

# -------------------- RUN STRUCTURAL FOR α = 0.3, 0.5, 0.7 -------------------
# This block only runs if a structural routine run_extensions() is defined.
if (exists("run_extensions")) {
  
  set.seed(123)
  alphas <- c(0.3, 0.5, 0.7)
  res_by_alpha <- list()
  
  # If baseline results exist, treat them as α = 0.5.
  if (exists("results")) {
    res_by_alpha[["0.5"]] <- results
  }
  
  # Run structural model for each α in {0.3, 0.5, 0.7} not already in cache.
  for (a in alphas) {
    key <- sprintf("%.1f", a)
    if (!is.null(res_by_alpha[[key]])) next
    res_by_alpha[[key]] <- run_structural_for_alpha(a)
  }
  
  # Convenience labels: alpha03 / alpha05 / alpha07 for table names.
  alpha_lab <- function(a) sprintf("alpha%02d", round(a * 100))
  
  # -------------------- Sα.1: Structural primitives for α = .3 and .7 ----------
  # Summarize d_i, φ_i, and share at cap for α = 0.3 and α = 0.7.
  for (a in c(0.3, 0.7)) {
    lab <- alpha_lab(a)
    res <- res_by_alpha[[sprintf("%.1f", a)]]
    sp  <- res$struct_params
    if (is.null(sp)) {
      message("Skipping Sα.1 for α = ", a, ": `struct_params` missing.")
      next
    }
    
    df_final <- pgg_data |>
      dplyr::filter(round_n == max(round_n, na.rm = TRUE)) |>
      dplyr::select(player, contributing)
    
    S1 <- summarise_struct_alpha(sp, df_final)
    
    gt_S1 <- S1$tab_primitives |>
      gt::gt() |>
      gt::fmt_number(columns = "value", decimals = 3) |>
      gt::tab_header(title = paste0("Structural primitives — ", lab))
    
    save_gt_both(gt_S1, paste0("Table_Sa1_structural_primitives_", lab))
  }
  
  # -------------------- Rank concordance vs α = 0.5 ----------------------
  # Compute Spearman and Kendall rank correlations between parameters under
  # α = 0.5 and α ∈ {0.3, 0.7}.
  sp_05 <- res_by_alpha[["0.5"]]$struct_params
  if (!is.null(sp_05)) {
    tabs <- list()
    for (a in c(0.3, 0.7)) {
      sp_alt <- res_by_alpha[[sprintf("%.1f", a)]]$struct_params
      if (is.null(sp_alt)) next
      rc <- rank_concord_two(sp_05, sp_alt) |>
        dplyr::mutate(compare = paste0(alpha_lab(a), " vs alpha05"))
      tabs[[alpha_lab(a)]] <- rc
    }
    if (length(tabs)) {
      rc_all <- dplyr::bind_rows(tabs) |>
        dplyr::select(compare, parameter, spearman, kendall)
      
      gt_rc <- rc_all |>
        gt::gt() |>
        gt::fmt_number(columns = c("spearman","kendall"), decimals = 3) |>
        gt::tab_header(title = "Rank concordance across α (Spearman ρ, Kendall τ)")
      
      save_gt_both(gt_rc, "Table_Sa2_rank_concordance_alpha_pairs")
    } else {
      message("Sα.2 skipped: could not compute any rank concordance.")
    }
  } else {
    message("Sα.2 skipped: baseline α = 0.5 structural parameters not found.")
  }
  
  # -------------------- Sα.3: PPC light (and full if available) ----------------
  # Compute simple PPC metrics for α = 0.3 and α = 0.7.
  for (a in c(0.3, 0.7)) {
    lab <- alpha_lab(a)
    res <- res_by_alpha[[sprintf("%.1f", a)]]
    if (is.null(res) || is.null(res$struct_params)) {
      message("Skipping Sα.3 for α = ", a, ": results/struct_params missing.")
      next
    }
    ppc <- ppc_light(res, pgg_data)
    
    gt_ppc <- ppc |>
      gt::gt() |>
      gt::fmt_number(columns = "value", decimals = 3) |>
      gt::tab_header(title = paste0("Posterior predictive checks — ", lab))
    
    save_gt_both(gt_ppc, paste0("Table_Sa3_ppc_metrics_", lab))
  }
  
  # -------------------- Sα.4: Policy invariance on m-grid ----------------------
  # For α = 0.3 and α = 0.7, evaluate how c^H responds to a subsidy m grid.
  for (a in c(0.3, 0.7)) {
    lab <- alpha_lab(a)
    res <- res_by_alpha[[sprintf("%.1f", a)]]
    if (is.null(res) || is.null(res$struct_params)) next
    
    pol <- policy_grid(res, m_grid = c(0, 0.25, 0.5))
    if (is.null(pol)) {
      message("Sα.4 skipped for ", lab, ": no callable policy method found in `results`.")
      next
    }
    
    gt_pol <- pol |>
      gt::gt() |>
      gt::fmt_number(columns = dplyr::everything(), decimals = 3) |>
      gt::tab_header(title = paste0("Policy invariance (share High, mean c^H) — ", lab))
    
    save_gt_both(gt_pol, paste0("Table_Sa4_policy_grid_", lab))
  }
  
} else {
  
  message("Skipping structural α–ROBUSTNESS block: `run_extensions()` is not defined.")
  
}

# =============================================================================
# End α–ROBUSTNESS block
# =============================================================================


# ============================================================================
# α–ROBUSTNESS without gt / run_extensions / HMM / humility
# Produces Sα.1–Sα.4 for α ∈ {0.3, 0.7} into ./tables/ (LaTeX + HTML).
# Requires: pgg_data with columns (player, round_n, contributing).
# This is a closed-form approximation that does not call run_extensions().
# ============================================================================
#
# This second α–robustness block is a stand-alone fallback that:
#   - Uses simple closed-form expressions for c^H and d_i (with φ_i ≡ 0),
#   - Does not require run_extensions(), HMMs, or gt,
#   - Writes LaTeX and HTML tables using xtable/htmlTable instead.
# It replicates the structure of Sα.1–Sα.4 in a minimal way for environments
# where the full structural machinery is not available.

safe_lib <- function(pkgs) {
  for (p in pkgs) if (!suppressWarnings(require(p, character.only = TRUE))) {
    install.packages(p, dependencies = TRUE); library(p, character.only = TRUE)
  }
}
safe_lib(c("dplyr","tidyr","tibble","xtable","htmlTable","stats"))

# ---------- model constants (matching the paper specification) ----------
# Public-good game parameters:
#   - b      : marginal per-capita return from the public account.
#   - N      : group size.
#   - kappa  : private return from the outside option.
#   - endow  : individual endowment.
#   - Delta  : κ − b/N, net return differential used in closed-form formulas.
b <- 2
N <- 5
kappa <- 1
endow <- 12
Delta <- kappa - b/N   # = 1 - 2/5 = 0.6

# ---------- I/O helpers (LaTeX + HTML, no gt dependency) ----------
# Minimal table writers to ./tables without relying on gt.

.save_tex <- function(df, stem, caption = NULL, outdir = file.path(getwd(),"tables")) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  fn <- file.path(outdir, paste0(stem, ".tex"))
  tb <- xtable::xtable(df, caption = caption, align = rep("r", ncol(df)+1))
  # keep content as-is, remove rownames
  sink(fn); print(tb, include.rownames=FALSE, sanitize.text.function=identity, comment=FALSE); sink()
  message("Saved LaTeX: ", basename(fn))
}

.save_html <- function(df, stem, caption = NULL, outdir = file.path(getwd(),"tables")) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  fn <- file.path(outdir, paste0(stem, ".html"))
  # Simple HTML wrapper with minimal styling.
  html <- paste0(
    "<!DOCTYPE html><html><head><meta charset='utf-8'><title>", stem, "</title>",
    "<style>table,th,td{border:1px solid #ddd;border-collapse:collapse;padding:6px;}th{background:#f5f5f5;}</style>",
    "</head><body>",
    if (!is.null(caption)) paste0("<h3>", caption, "</h3>") else "",
    htmlTable::htmlTable(df),
    "</body></html>"
  )
  writeLines(html, fn, useBytes = TRUE)
  message("Saved HTML: ", basename(fn))
}

# Save a data frame as both LaTeX and HTML (no gt).
save_df_both <- function(df, stem, caption=NULL, outdir=file.path(getwd(),"tables")) {
  .save_tex(df, stem, caption, outdir)
  .save_html(df, stem, caption, outdir)
}

# ---------- utilities ----------
# Simple column picker by candidate names.
.pick <- function(df, opts) { x <- intersect(opts, names(df)); if (length(x)) x[1] else NULL }

# Summary vector: quantiles + mean + N, rounded.
summ_vec <- function(x, digits=3) {
  x <- as.numeric(x)
  qs <- stats::quantile(x, probs = c(0,.10,.25,.50,.75,.90,1), na.rm=TRUE)
  out <- tibble::tibble(
    stat  = c("Min","P10","Q1","Median","Q3","P90","Max","Mean","N"),
    value = c(unname(qs), mean(x,na.rm=TRUE), sum(is.finite(x)))
  )
  out$value <- round(out$value, digits)
  out
}

# Label α as alpha03, alpha07, etc. (matches naming in other block).
alpha_lab <- function(a) sprintf("alpha%02d", round(a*100))

# ---------- 1) steady contributions per player ----------
# Approximate c_i^H using average contributions over late rounds (8–10),
# and keep final-round contributions for PPC diagnostics.
steady <- pgg_data %>%
  dplyr::filter(round_n >= 8) %>%
  dplyr::group_by(player) %>%
  dplyr::summarise(cH_obs = mean(contributing, na.rm = TRUE),
                   c_final = contributing[which.max(round_n)],
                   .groups = "drop") %>%
  dplyr::mutate(cH_obs = pmin(pmax(cH_obs, 0), endow))

# ---------- 2) recover d_i(α) from c_i^H with φ_i ≡ 0 ----------
# Structural relationship (closed-form):
#   c^H = [Δ / (α d_i)]^{1/(α-1)}  ⇒  d_i = Δ / (α (c^H)^{α-1}).
# We invert this formula treating observed cH_obs as c^H.
alpha_grid <- c(0.3, 0.7)  # (α=0.5 is baseline in paper; here focus on 0.3, 0.7)

est_d_for_alpha <- function(cH, alpha) {
  cH <- pmin(pmax(as.numeric(cH), 1e-6), endow)   # avoid 0^{α-1} for α<1
  d <- Delta / (alpha * (cH)^(alpha - 1))
  d[!is.finite(d)] <- NA_real_
  d
}

# Construct structural parameter tables for each α in alpha_grid.
# For each α, we store:
#   - cH_obs  : observed proxy for c^H,
#   - c_final : final-round contribution,
#   - d_hat   : altruism parameter from closed form,
#   - phi_hat : set to 0 in this simplified approximation.
params_by_alpha <- lapply(alpha_grid, function(a) {
  steady %>% dplyr::transmute(player, cH_obs, c_final,
                              d_hat = est_d_for_alpha(cH_obs, a),
                              phi_hat = 0)
})
names(params_by_alpha) <- sprintf("%.1f", alpha_grid)

# ---------- 3) Sα.1 Structural primitives tables ----------
# For α = 0.3 and 0.7, compute distribution of d_i and φ_i,
# along with share at cap and share φ_i ≤ 0.1 (trivially 1 here).
for (a in alpha_grid) {
  lab <- alpha_lab(a)
  sp  <- params_by_alpha[[sprintf("%.1f", a)]]
  tab_d   <- summ_vec(sp$d_hat)   %>% dplyr::mutate(quantity="d_i")   %>% dplyr::select(quantity, stat, value)
  tab_phi <- summ_vec(sp$phi_hat) %>% dplyr::mutate(quantity="phi_i") %>% dplyr::select(quantity, stat, value)
  
  share_at_cap <- mean(sp$cH_obs >= (endow - 1e-9), na.rm = TRUE)
  share_phi    <- 1.0  # phi=0 ⇒ Pr(phi≤0.1)=1 by construction.
  shares <- tibble::tibble(quantity=c("share_at_cap_cH", "share_phi_le_0.1"),
                           stat="Share", value=c(round(share_at_cap,3), round(share_phi,3)))
  
  out <- dplyr::bind_rows(tab_d, tab_phi, shares)
  save_df_both(out,
               stem    = paste0("Table_Sa1_structural_primitives_", lab),
               caption = paste0("Structural primitives — ", lab))
}

# ---------- 4) Sα.2 Rank concordance (vs α = 0.5) ----------
# If α=0.5 baseline parameters are available elsewhere, they can be loaded instead.
# Here α=0.5 parameters are constructed using the same closed-form formula.
params_alpha05 <- steady %>% dplyr::transmute(player, cH_obs, c_final,
                                              d_hat = Delta/(0.5*(pmin(pmax(cH_obs,1e-6),endow))^(0.5-1)),
                                              phi_hat = 0)

# Helper: rank concordance between d_i under two α values (closed-form).
rank_concord <- function(base, alt, base_lab, alt_lab) {
  x <- dplyr::inner_join(base %>% dplyr::select(player, d_base=d_hat),
                         alt  %>% dplyr::select(player, d_alt =d_hat), by="player")
  tibble::tibble(
    compare  = paste0(alt_lab, " vs ", base_lab),
    parameter= "d_i",
    spearman = round(stats::cor(x$d_base, x$d_alt, method="spearman", use="complete.obs"), 3),
    kendall  = round(stats::cor(x$d_base, x$d_alt, method="kendall",  use="complete.obs"), 3)
  )
}

rc <- dplyr::bind_rows(
  rank_concord(params_alpha05, params_by_alpha[["0.3"]], "alpha05", "alpha03"),
  rank_concord(params_alpha05, params_by_alpha[["0.7"]], "alpha05", "alpha07")
)
save_df_both(rc,
             stem    = "Table_Sa2_rank_concordance_alpha_pairs",
             caption = "Rank concordance across α (Spearman ρ, Kendall τ)")

# ---------- 5) Sα.3 Posterior predictive checks ----------
# Using cH_obs as steady prediction:
#   - Median |c_final - cH| per player,
#   - Δ mean (final - mean cH),
#   - RMSE of round means, where predicted path is constant at mean(cH_obs).
ppc_for_alpha <- function(sp) {
  tmp <- dplyr::inner_join(sp %>% dplyr::select(player, cH = cH_obs),
                           steady %>% dplyr::select(player, c_final), by="player")
  med_abs <- stats::median(abs(tmp$c_final - tmp$cH), na.rm = TRUE)
  delta_mu <- mean(tmp$c_final, na.rm = TRUE) - mean(tmp$cH, na.rm = TRUE)
  
  emp_means <- pgg_data %>%
    dplyr::group_by(round_n) %>%
    dplyr::summarise(mean_emp = mean(contributing, na.rm = TRUE), .groups = "drop")
  mu_pred <- mean(sp$cH_obs, na.rm = TRUE)
  rmse_round <- sqrt(mean((emp_means$mean_emp - mu_pred)^2, na.rm = TRUE))
  
  tibble::tibble(
    metric = c("Median |c_final - c_H| (per player)",
               "Δ mean (final - mean c_H)",
               "RMSE of round means"),
    value  = round(c(med_abs, delta_mu, rmse_round), 3)
  )
}

for (a in alpha_grid) {
  lab <- alpha_lab(a)
  sp  <- params_by_alpha[[sprintf("%.1f", a)]]
  df  <- ppc_for_alpha(sp)
  save_df_both(df,
               stem    = paste0("Table_Sa3_ppc_metrics_", lab),
               caption = paste0("Posterior predictive checks — ", lab))
}

# ---------- 6) Sα.4 Policy invariance (subsidy m ∈ {0, .25, .5}) ----------
# Under this closed-form approximation, a policy m reduces κ to κ-m (bounded below by 0),
# which adjusts Δ(m) and hence c^H(m) for each d_i:
#   c^H(m) = [ Δ(m) / (α d_i) ]^{1/(α-1)}, clipped to [0, endow].
cH_of_m <- function(d_i, alpha, m) {
  kappa_p <- max(kappa - m, 0)
  Delta_m <- kappa_p - b/N
  cH <- ifelse(Delta_m <= 0, 0, (Delta_m / (alpha * d_i))^(1/(alpha - 1)))
  pmin(pmax(cH, 0), endow)
}

policy_grid_alpha <- function(sp, alpha, m_grid = c(0,0.25,0.5)) {
  rows <- lapply(m_grid, function(m) {
    d <- sp$d_hat; d <- d[is.finite(d) & d>0]
    cH_m <- cH_of_m(d, alpha, m)
    tibble::tibble(
      m          = m,
      share_high = round(mean(cH_m >= (endow - 1e-9), na.rm = TRUE), 3),
      mean_cH    = round(mean(cH_m, na.rm = TRUE), 3)
    )
  })
  dplyr::bind_rows(rows)
}

for (a in alpha_grid) {
  lab <- alpha_lab(a)
  sp  <- params_by_alpha[[sprintf("%.1f", a)]]
  pol <- policy_grid_alpha(sp, alpha = a, m_grid = c(0, 0.25, 0.5))
  save_df_both(pol,
               stem    = paste0("Table_Sa4_policy_grid_", lab),
               caption = paste0("Policy invariance (share High, mean c^H) — ", lab))
}

message("\nα-robustness tables written to ./tables/ :\n",
        "  - Table_Sa1_structural_primitives_alpha03/07.(html|tex)\n",
        "  - Table_Sa2_rank_concordance_alpha_pairs.(html|tex)\n",
        "  - Table_Sa3_ppc_metrics_alpha03/07.(html|tex)\n",
        "  - Table_Sa4_policy_grid_alpha03/07.(html|tex)\n")
