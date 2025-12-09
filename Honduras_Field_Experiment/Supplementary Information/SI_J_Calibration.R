# ================================
# Model calibration & visualization
# ================================
# This script:
#   (i)  constructs an empirical 2×2 High/Low transition matrix between round 1 and round 10,
#   (ii) calibrates a binary Fermi–Moran imitation process (with optional IH heterogeneity)
#        so that its simulated transition matrix matches the empirical one,
#   (iii) explores the loss surface over (d, k),
#   (iv) compares empirical vs. simulated transition matrices, and
#   (v) contrasts empirical vs. simulated dynamics in the share of High contributors over rounds.
#
# Notation follows the main text and SI:
#   - "H" / "L" are High/Low states defined by a threshold on contributions in round 1.
#   - d     : selection tilt in favor of High contributors (d < 0 favors Low).
#   - k     : imitation intensity (Fermi temperature).
#   - N     : group size (N = 5, as in the experiment).
#   - rounds: macroscopic horizon (here, 10 rounds, matching the experiment).

# ---- Packages ----
need <- c("dplyr","tidyr","tibble","ggplot2","pbapply","patchwork")
new  <- setdiff(need, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble)
  library(ggplot2); library(pbapply); library(patchwork)
})

# Use a clean minimalist theme for all ggplot output
theme_set(theme_minimal(base_size = 12))

# -------------------------------------------------------------------
# Input data: long panel of PGG decisions at the individual-round level
# -------------------------------------------------------------------
# Expected minimal columns in data_set.csv:
#   - player       : individual identifier
#   - round_n      : round index (1–10)
#   - contributing : amount contributed in that round (0–endowment, in Lempiras)
# Optional:
#   - IH           : individual heterogeneity index used to weight imitation (if present)
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# ---- Sanity checks on required columns ----
stopifnot(all(c("player","round_n","contributing") %in% names(pgg_data)))

# ===========================================
# 1) Empirical threshold & 2x2 transition matrix (R1 -> R10)
# ===========================================
rounds <- 10

# Threshold for High vs Low:
# We use the empirical mean contribution in round 1, as in the binary
# High/Low classification in the SI (ĉ in the notation), i.e. H if c > ĉ.
thresh <- pgg_data %>%
  filter(round_n == 1) %>%
  summarise(m = mean(contributing, na.rm = TRUE)) %>%
  pull(m)

# Empirical 2×2 transition matrix T_emp:
#   - Rows: initial state in round 1 (H/L)
#   - Cols: final state in round 10 (H/L)
#   - Entries: probabilities Pr(state_10 | state_1)
T_emp <- pgg_data %>%
  filter(round_n %in% c(1, rounds)) %>%
  mutate(state = if_else(contributing > thresh, "H", "L")) %>%
  dplyr::select(player, round_n, state) %>%
  tidyr::pivot_wider(
    names_from   = round_n,
    names_prefix = "r",
    values_from  = state
  ) %>%
  count(r1, r10, name = "n") %>%
  complete(r1 = c("H","L"), r10 = c("H","L"), fill = list(n = 0)) %>%
  tidyr::pivot_wider(names_from = r10, values_from = n) %>%
  tibble::column_to_rownames("r1") %>%
  as.matrix()

# Row-normalize to obtain transition probabilities T_emp(row | row state at R1).
# If a row has zero mass (no observations), keep a uniform row.
rs <- rowSums(T_emp)
rs[rs == 0] <- 1
T_emp <- T_emp / rs

# ===========================================
# 2) Simulators (homogeneous & IH-weighted)
# ===========================================
# These functions simulate a binary Fermi–Moran imitation process.
# Each run:
#   - starts from an i.i.d. Bernoulli(0.5) configuration of N players (0=L, 1=H),
#   - applies 'rounds' birth–death updates with Fermi-type selection,
#   - records the initial and final states for all N players,
#   - aggregates across n_sims runs into a 2×2 transition matrix.
# Output:
#   - 2×2 row-stochastic matrix over {H,L} × {H,L}.

fermi_moran_matrix <- function(d, k, rounds = 10, N = 5, n_sims = 4000){
  # d : payoff advantage of High vs Low (fitness 1+d vs 1)
  # k : intensity of selection in the Fermi weights (k = 0 is random drift)
  # rounds : number of microscopic birth–death updates
  # N : group size
  # n_sims : number of independent realizations used to estimate the matrix
  
  trans <- matrix(0, 2, 2, dimnames = list(c("H","L"), c("H","L")))
  
  for(sim in seq_len(n_sims)){
    # Initial state: Bernoulli(0.5) for each of the N players; 1=H, 0=L
    s0 <- sample(0:1, N, TRUE); s <- s0
    
    # Moran birth–death process over 'rounds' microscopic updates
    for(t in seq_len(rounds)){
      # Fitness: High types get 1 + d, Low types get 1
      fit <- ifelse(s == 1, 1 + d, 1)
      
      # Fermi weights: probability of being chosen as role model
      p   <- exp(k * fit); p <- p / sum(p)
      if (any(!is.finite(p))) p <- rep(1/N, N)  # numerical guard
      
      # Birth–death step: one random individual copies a random role model
      s[sample(N, 1)] <- s[sample(N, 1, prob = p)]
    }
    
    # Update empirical counts for the 2×2 transition matrix
    for(idx in seq_len(N)){
      trans[ ifelse(s0[idx], "H","L"), ifelse(s[idx], "H","L") ] <-
        trans[ ifelse(s0[idx], "H","L"), ifelse(s[idx], "H","L") ] + 1
    }
  }
  
  # Convert counts to row-stochastic transition probabilities
  prop.table(trans, 1)
}

moran_IH_matrix <- function(d, k, IH_vec, rounds = 10, N = 5, n_sims = 4000){
  # Same Fermi–Moran process as above, but with heterogeneous imitation weights:
  #   - For each simulation, each player draws an IH "type" h from IH_vec.
  #   - The Fermi weight becomes exp(k * h * fit).
  #   - Larger h implies stronger effective selection / norm salience.
  stopifnot(length(IH_vec) > 0)
  
  trans <- matrix(0, 2, 2, dimnames = list(c("H","L"), c("H","L")))
  
  for(sim in seq_len(n_sims)){
    # Draw heterogeneous 'h' types and initial H/L states
    h  <- sample(IH_vec, N, TRUE)
    s0 <- sample(0:1, N, TRUE); s <- s0
    
    for(t in seq_len(rounds)){
      fit <- ifelse(s == 1, 1 + d, 1)
      
      # Heterogeneous Fermi weights: exp(k * h_i * fit_i)
      p   <- exp(k * h * fit); p <- p / sum(p)
      if (any(!is.finite(p))) p <- rep(1/N, N)  # numerical guard
      
      s[sample(N, 1)] <- s[sample(N, 1, prob = p)]
    }
    
    for(idx in seq_len(N)){
      trans[ ifelse(s0[idx], "H","L"), ifelse(s[idx], "H","L") ] <-
        trans[ ifelse(s0[idx], "H","L"), ifelse(s[idx], "H","L") ] + 1
    }
  }
  
  prop.table(trans, 1)
}

# ===========================================
# 3) Calibration helpers
# ===========================================
# Loss function: squared Frobenius distance between simulated and empirical 2×2 matrices
rss_loss <- function(M) sum((M - T_emp)^2)

calibrate_moran <- function(use_IH = FALSE,
                            rounds = 10, N = 5, n_sims = 3000,
                            d_grid = seq(-2, 3, .5),
                            k_grid = seq(0, 1.5, .25),
                            IH_vec = NULL){
  
  # Wrapper that:
  #   1. chooses the appropriate simulator (homogeneous vs IH-weighted),
  #   2. performs a coarse grid search over (d, k),
  #   3. refines the best grid point via L-BFGS-B optimization,
  #   4. returns the optimal (d, k), RSS, and the corresponding transition matrix.
  
  simfun <- if (use_IH) {
    function(d,k) moran_IH_matrix(d,k,IH_vec,rounds,N,n_sims)
  } else {
    function(d,k) fermi_moran_matrix(d,k,rounds,N,n_sims)
  }
  
  # Coarse grid over parameter space
  grid   <- expand.grid(d = d_grid, k = k_grid)
  rss    <- pbapply::pbsapply(
    seq_len(nrow(grid)),
    \(i) rss_loss(simfun(grid$d[i], grid$k[i]))
  )
  
  # Starting point for local optimization: grid point with smallest RSS
  start  <- grid[which.min(rss), ]
  
  # Continuous refinement of (d, k) with box constraints
  opt <- optim(
    par    = c(start$d, start$k),
    fn     = \(par) rss_loss(simfun(par[1], par[2])),
    method = "L-BFGS-B",
    lower  = c(min(d_grid), 0),          # enforce k >= 0
    upper  = c(max(d_grid), max(k_grid)),
    control = list(maxit = 80)
  )
  
  # Return calibrated parameters, loss, and fitted transition matrix
  list(
    d = opt$par[1],
    k = opt$par[2],
    rss = opt$value,
    T = simfun(opt$par[1], opt$par[2])
  )
}

compute_loss_surface <- function(d_grid, k_grid, rounds = 10, N = 5, n_sims = 500) {
  # Compute RSS(d,k) on a grid using the homogeneous Fermi–Moran model.
  # This is used only for visualization of the calibration loss surface.
  grid <- expand.grid(d = d_grid, k = k_grid)
  rss  <- pbapply::pbsapply(seq_len(nrow(grid)), function(i) {
    M <- fermi_moran_matrix(
      grid$d[i], grid$k[i],
      rounds = rounds, N = N, n_sims = n_sims
    )
    sum((M - T_emp)^2)
  })
  cbind(grid, rss = rss)
}

plot_heat <- function(M, title){
  # Plot a 2×2 transition matrix M with cell values and heat shading.
  # Var1: row (initial state), Var2: column (final state).
  as.data.frame(as.table(M)) %>%
    ggplot(aes(Var1, Var2, fill = Freq)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", Freq)), size = 5) +
    labs(x = "Initial", y = "Final", title = title) +
    theme_minimal() +
    scale_fill_gradient(low = "white", high = "steelblue")
}

# ===========================================
# 4) Optional IH vector (min–max to [0,1]); skip if missing
# ===========================================
# If the dataset contains a column "IH" (e.g., an individual-level heterogeneity or
# norm-salience index), construct a normalized vector in [0,1] to use in moran_IH_matrix.
IH_vec <- NULL
if ("IH" %in% names(pgg_data)) {
  tmp <- pgg_data %>%
    filter(!is.na(IH)) %>%
    dplyr::select(IH)
  if (nrow(tmp) > 0 && diff(range(tmp$IH)) > 0) {
    IH_vec <- ((tmp$IH - min(tmp$IH)) / (max(tmp$IH) - min(tmp$IH))) %>%
      as.numeric()
  }
}

# ===========================================
# 5) Fit models (plain and, if available, IH-weighted)
# ===========================================
set.seed(2025)

# Baseline homogeneous Fermi–Moran calibration
fit_plain <- calibrate_moran(FALSE, rounds = rounds, N = 5, n_sims = 1500)

# IH-weighted calibration (only if IH_vec is available and non-degenerate)
if (!is.null(IH_vec)) {
  fit_IH <- calibrate_moran(TRUE,  rounds = rounds, N = 5, n_sims = 1500, IH_vec = IH_vec)
} else {
  fit_IH <- NULL
}

# Choose the better fit by RSS
best_fit <- if (!is.null(fit_IH) && fit_IH$rss < fit_plain$rss) fit_IH else fit_plain
label    <- if (!is.null(fit_IH) && fit_IH$rss < fit_plain$rss) "(heterog. IH)" else "(homogeneous)"

# Print a concise calibration summary
cat("\nCalibration summary:\n")
print(
  rbind(
    c(
      model = "No-IH",
      d  = round(fit_plain$d,3),
      k  = round(fit_plain$k,3),
      RSS = round(fit_plain$rss,4)
    ),
    if (!is.null(fit_IH)) c(
      model = "With IH",
      d  = round(fit_IH$d,3),
      k  = round(fit_IH$k,3),
      RSS = round(fit_IH$rss,4)
    ) else NULL
  ),
  row.names = FALSE
)

# ===========================================
# 6) Loss surface (coarse grid, quick)
# ===========================================
# Map the RSS(d,k) surface around the empirically relevant region.
set.seed(2025)
loss_df <- compute_loss_surface(
  d_grid = seq(-1.5, 1.5, by = 0.15),
  k_grid = seq(0,    1.20, by = 0.05),
  rounds = rounds, N = 5, n_sims = 500
)
opt_row <- loss_df[which.min(loss_df$rss), ]

# Heat map + contours of the loss surface,
# with both the grid minimum and the final calibrated (best_fit) point marked.
p_loss <- ggplot(loss_df, aes(d, k, fill = rss)) +
  geom_tile() +
  stat_contour(aes(z = rss), bins = 10, size = 0.25) +
  geom_point(data = opt_row, aes(d, k), shape = 4, size = 3, stroke = 1.2) +
  geom_point(aes(x = best_fit$d, y = best_fit$k),
             shape = 21, size = 3, stroke = 1) +
  labs(
    title = "Calibration loss surface  (RSS vs. d, k)",
    x = expression(d), y = expression(k),
    subtitle = sprintf(
      "Grid min: d=%.2f, k=%.2f  |  Optimized: d=%.2f, k=%.2f",
      opt_row$d, opt_row$k, best_fit$d, best_fit$k
    )
  ) +
  theme_minimal()

# ===========================================
# 7) Transition matrices: empirical vs calibrated
# ===========================================
# Left panel: empirical R1→R10 transition matrix.
p_emp <- plot_heat(T_emp, "Empirical  (R1 \u2192 R10)")

# Right panel: transition matrix from the calibrated Fermi–Moran model.
p_sim <- plot_heat(
  best_fit$T,
  sprintf("Calibrated %s\n(d = %.2f, k = %.2f)",
          label, best_fit$d, best_fit$k)
)

p_transitions <- p_emp + p_sim +
  plot_annotation(title = "Transition matrices")

# ===========================================
# 8) Share High over rounds: data vs simulation
# ===========================================
# Empirical share of High contributors (c > thresh) in each round.
emp_share <- pgg_data %>%
  filter(round_n >= 1, round_n <= rounds) %>%
  mutate(state = if_else(contributing > thresh, "H", "L")) %>%
  group_by(round_n) %>%
  summarise(
    share_H = mean(state == "H", na.rm = TRUE),
    n       = dplyr::n(),
    .groups = "drop"
  )

# Empirical initial share of High in round 1
init_share_H <- emp_share$share_H[emp_share$round_n == 1]

simulate_moran_shares <- function(d, k,
                                  rounds = 10, N = 5, R = 400,
                                  init_share = 0.5, updates_per_round = 5*N,
                                  seed = NULL) {
  # Simulate the time path of the share of High contributors under a
  # pairwise-comparison (Fermi) imitation process:
  #   - Population: N players, state theta ∈ {0=L,1=H}
  #   - Fitness: w = 1 + d * theta
  #   - Imitation probability: logistic in fitness difference with intensity k.
  #
  # Args:
  #   d, k            : same interpretation as above.
  #   rounds          : number of macro rounds to record (T).
  #   N               : population size in each simulation.
  #   R               : number of independent simulation runs.
  #   init_share      : initial fraction of High in round 1.
  #   updates_per_round: number of pairwise comparison updates per macro round.
  #   seed            : optional RNG seed for reproducibility.
  #
  # Returns:
  #   tibble with round_n, sim_mean, and 10–90% quantile band for the share of High.
  
  if (!is.null(seed)) set.seed(seed)
  
  shares <- matrix(NA_real_, nrow = R, ncol = rounds)
  
  for (r in seq_len(R)) {
    # Initialize states: Bernoulli(init_share) for each of N players
    theta <- as.integer(runif(N) < init_share)  # 1=H, 0=L
    shares[r, 1] <- mean(theta)
    
    for (t in 2:rounds) {
      # Perform multiple pairwise comparison updates within each macro round
      for (u in seq_len(updates_per_round)) {
        i <- sample.int(N, 1)
        j <- sample.int(N, 1); while (j == i) j <- sample.int(N, 1)
        
        w_i <- 1 + d * theta[i]
        w_j <- 1 + d * theta[j]
        
        # Fermi update: i overwrites j with probability logistic(k*(w_i - w_j))
        p_ij <- 1 / (1 + exp(-k * (w_i - w_j)))
        if (runif(1) < p_ij) theta[j] <- theta[i]
      }
      shares[r, t] <- mean(theta)
    }
  }
  
  tibble(
    round_n = rep(seq_len(rounds), each = R),
    share   = as.numeric(shares)
  ) %>%
    group_by(round_n) %>%
    summarise(
      sim_mean = mean(share),
      lo       = quantile(share, 0.10),
      hi       = quantile(share, 0.90),
      .groups  = "drop"
    )
}

# Simulated share of High using calibrated (d, k), with same N and horizon as data.
sim_share <- simulate_moran_shares(
  best_fit$d, best_fit$k,
  rounds = rounds, N = 5, R = 600,
  init_share = init_share_H, seed = 2025
)

# Plot empirical vs simulated share of High over rounds.
p_share <- ggplot() +
  geom_ribbon(
    data = sim_share,
    aes(round_n, ymin = lo, ymax = hi),
    alpha = 0.25
  ) +
  geom_line(
    data = sim_share,
    aes(round_n, sim_mean),
    linewidth = 1
  ) +
  geom_point(
    data = emp_share,
    aes(round_n, share_H)
  ) +
  geom_line(
    data = emp_share,
    aes(round_n, share_H),
    linetype = 2
  ) +
  labs(
    title = "Share of High contributors over rounds",
    subtitle = "Data (points/dashed) vs. calibrated simulation (line & 10–90% band)",
    x = "Round", y = "Share High"
  ) +
  theme_minimal()

# ===========================================
# 9) Combine panels & export
# ===========================================
# Final panel:
#   - top-left : loss surface over (d, k),
#   - top-right: empirical vs calibrated 2×2 transition matrices,
#   - bottom  : empirical vs simulated time path of the High share.
calibration_panel <- (p_loss | p_transitions) / p_share +
  plot_annotation(
    title = "Model calibration: loss, transitions, and dynamics",
    theme = theme(plot.title = element_text(face = "bold"))
  )

# Save as PDF for inclusion in the SI
ggsave("calibration_panel.pdf", calibration_panel, width = 12, height = 9, dpi = 300)

cat("\nSaved figure: calibration_panel.pdf\n")
