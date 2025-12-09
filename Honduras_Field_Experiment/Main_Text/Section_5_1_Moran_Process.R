###############################################################################
#  MORAN CALIBRATION  – homogeneous vs. rich heterogeneous fitness
#  --------------------------------------------------------------------------
#  •  homogeneous (“plain”):           fitness_i  = 1 + d⋅1{H_i}
#  •  heterogeneous (15 covariates):   fitness_i  = 1 + d⋅1{H_i} + γ·HET_i
#       – HET_i  =   scaled composite of the 15 baseline variables
#       – γ      =   single scale parameter that turns covariate dispersion
#                    into fitness dispersion (avoids estimating 15 slopes)
#
#  The imitation step uses a Fermi rule  p_i ∝ exp(k·fitness_i).
#  d, k (and γ when heterogeneous=TRUE) are calibrated by
#      argmin ‖T_sim(d,k,γ) – T_emp‖²_F   (grid → L-BFGS-B).
#
#  High‑level pipeline:
#    1. Build an empirical 2×2 “High vs Low” transition matrix T_emp
#       from round 1 to round 10 of the PGG.
#    2. Construct a 1D heterogeneity index HET_i from 15 baseline covariates.
#    3. Simulate a binary Moran process (simulate_moran) with parameters
#       (d, k, γ) and compute its 2×2 transition matrix.
#    4. Calibrate (d, k, γ) so that the simulated matrix matches T_emp
#       as closely as possible in squared Frobenius norm.
#    5. Fit both homogeneous (γ = 0) and heterogeneous (γ ≠ 0) versions.
#    6. Plot empirical vs calibrated transition matrices for visual comparison.
###############################################################################

suppressPackageStartupMessages({
  # dplyr / tidyr / tibble: data wrangling
  # pbapply: pb*sapply with progress bar (for grid search)
  # ggplot2: plotting; patchwork: combining plots side by side
  library(dplyr); library(tidyr); library(tibble)
  library(pbapply); library(ggplot2); library(patchwork)
})

# Load the dataset of repeated public goods games.
# The CSV is expected to contain at least: player, round_n, contributing.
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)


# Baseline covariates (round 1 only); ensure types & z-scores
baseline <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::distinct(player, .keep_all = TRUE) %>%
  dplyr::transmute(
    player        = as.character(player),
    village_code  = as.character(village_code),
    gender        = factor(gender),
    marital_status= factor(marital_status),
    b0600         = factor(b0600),     # religion (relevel later)
    access_routes = factor(access_routes),
    age           = suppressWarnings(as.numeric(age)),
    friends       = suppressWarnings(as.numeric(friends)),
    adversaries   = suppressWarnings(as.numeric(adversaries)),
    network_density_fr  = suppressWarnings(as.numeric(network_density_fr)),
    network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
    network_size        = suppressWarnings(as.numeric(network_size)),
    b0100              = suppressWarnings(as.numeric(b0100)),  # education
    b0200              = as.numeric(b0200),  # indigenous
    FI                 = suppressWarnings(as.numeric(FI))
  ) %>%
  # Standardize continuous covariates
  dplyr::mutate(dplyr::across(
    c(age,friends,adversaries,network_density_fr,network_density_adv,network_size,b0100,FI),
    ~ as.numeric(scale(.x)),
    .names = "{.col}_z"
  )) %>%
  # Set religion ref = "2" (Catholic) if present; adjust as needed
  dplyr::mutate(b0600 = stats::relevel(b0600, ref = "2"))


## 1.  Empirical 2×2 transition matrix  ---------------------------------------
rounds  <- 10  # number of rounds in the experiment (R1 → R10)

# Empirical threshold used to define “High” (H) vs “Low” (L) contributors:
# mean contribution in round 1 across all players.
thresh  <- pgg_data |>
  filter(round_n == 1) |>
  summarise(m = mean(contributing, na.rm = TRUE)) |>
  pull()

# Construct empirical 2×2 matrix:
#   rows: initial state in round 1 (H or L)
#   cols: final state in round 10 (H or L)
#   entries: row-normalized probabilities of ending in each state.
T_emp <- pgg_data |>
  filter(round_n %in% c(1, rounds)) |>                                 # keep R1 and R10
  mutate(state = if_else(contributing > thresh, "H", "L")) |>          # binarize contribution
  dplyr::select(player, round_n, state) |>                                    # one row per (player, round)
  pivot_wider(names_from = round_n, names_prefix = "r",
              values_from = state) |>                                  # columns r1, r10
  count(r1, r10, name = "n") |>                                        # counts in each (r1, r10) cell
  complete(r1 = c("H","L"), r10 = c("H","L"),
           fill = list(n = 0)) |>                                      # ensure all 4 cells present
  pivot_wider(names_from = r10, values_from = n) |>                    # rows=r1, cols=r10
  column_to_rownames("r1") |>
  as.matrix() |>                                                       # convert to numeric matrix
  prop.table(1)                                                        # row-normalize to probabilities

## 2.  Build a *single* heterogeneity score  (z-scored and averaged) ----------
# List of baseline covariates that enter the composite heterogeneity index HET_i.
# Names must match columns in the `baseline` data frame.
covars <- c("age", "gender", "friends", "adversaries", "FI",
            "marital_status", "network_density_fr", "network_density_adv",
            "network_size", "b0100", "b0600", "b0200",
            "access_routes")

# Construct HET_i:
#   1. Coerce selected columns to numeric.
#   2. z-score each covariate (mean 0, sd 1).
#   3. Take the row-wise mean across the 15 z-scored variables.
#   4. Pull out as a numeric vector to be used in the simulator.
heter_vec <- baseline |>
  mutate(across(all_of(covars), as.numeric)) |>
  mutate(across(all_of(covars), ~ as.numeric(scale(.x)))) |>
  mutate(HET = rowMeans(across(all_of(covars)), na.rm = TRUE)) |>
  pull(HET)

## 3.  Simulator --------------------------------------------------------------
# Simulate a Moran birth–death process with Fermi imitation rule
# in a group of size N, for a given (d, k, gamma) and heterogeneity vector.
#
# Arguments:
#   d        : selection tilt for High vs Low (enters fitness_i = 1 + d⋅s_i + γ⋅h_i)
#   k        : intensity of selection / imitation (Fermi “temperature”)
#   gamma    : scale of heterogeneity (γ = 0 ⇒ homogeneous case)
#   heter_vec: vector of HET_i values; if NULL, no heterogeneity is used
#   rounds   : number of Moran updates per simulation (here aligned with game rounds)
#   N        : group size (N = 5 as in the experiment)
#   n_sims   : number of independent groups to simulate to estimate the transition matrix
simulate_moran <- function(d, k, gamma   = 0,      # γ = 0  ⇒ homogeneous
                           heter_vec    = NULL,
                           rounds       = 10,
                           N            = 5,
                           n_sims       = 4000) {
  
  # Accumulator for the 2×2 transition counts (initial H/L → final H/L).
  trans <- matrix(0, 2, 2, dimnames = list(c("H","L"), c("H","L")))
  
  # Repeat many independent simulations of a size‑N group.
  for(sim in seq_len(n_sims)){
    # draw a group (with heterogeneity if requested)
    # h_j is the HET value for individual j in this simulated group.
    h <- if(is.null(heter_vec)) rep(0, N) else sample(heter_vec, N, TRUE)
    
    # Random initial strategies: 0 = Low, 1 = High (binary contribution state).
    s0 <- sample(0:1, N, TRUE)       # 0 = L , 1 = H
    s  <- s0                         # current state vector, updated over rounds
    
    # Moran/Fermi dynamics for the given number of rounds.
    for(t in seq_len(rounds)){
      # Individual fitness: baseline 1 + d⋅s_j + γ⋅h_j
      fit <- 1 + d * s + gamma * h    # intrinsic + heterogeneity
      
      # Fermi imitation probabilities: p_j ∝ exp(k⋅fit_j)
      p   <- exp(k * fit);   p <- p / sum(p)
      
      # Numerical safeguard: if overflow/NaN, fall back to uniform copying.
      if(any(!is.finite(p))) p <- rep(1/N, N)
      
      # Moran step:
      #   1. randomly pick a “death” site (uniform over N indices)
      #   2. pick a “parent” to copy, chosen with probability ∝ p
      #   3. new state at the death site equals parent’s state.
      s[sample(N, 1)] <- s[sample(N, 1, prob = p)]
    }
    
    # After all rounds, record initial vs final H/L for each individual
    # and update the 2×2 transition count matrix.
    for(idx in seq_len(N)){
      trans[ ifelse(s0[idx], "H","L"),
             ifelse(s [idx], "H","L") ] <-
        trans[ ifelse(s0[idx], "H","L"),
               ifelse(s [idx], "H","L") ] + 1
    }
  }
  
  # Convert counts to a row-stochastic matrix (empirical transition probabilities).
  prop.table(trans, 1)
}

## 4.  Calibration wrapper ----------------------------------------------------
# Loss function: squared Frobenius norm between a candidate transition matrix M
# and the empirical transition matrix T_emp.
rss_loss <- function(M) sum((M - T_emp)^2)

# Wrapper that:
#   • chooses a coarse grid over (d, k, γ) or (d, k),
#   • evaluates rss_loss at each grid point via simulate_moran,
#   • uses the best grid point as starting value for L-BFGS-B optimization,
#   • returns the optimal parameters and the corresponding simulated matrix.
calibrate_moran <- function(use_heter = FALSE,
                            rounds = 10, N = 5, n_sims = 3000,
                            d_grid = seq(-2, 3, .5),
                            k_grid = seq(0, 1.5, .25),
                            g_grid = seq(0, 1.5, .25),       # γ-grid
                            heter_vec = NULL) {
  
  # choose parameter grids ----------------------------------------------------
  if(use_heter){
    # Full 3D grid for (d, k, γ) when using heterogeneity.
    grid <- expand.grid(d = d_grid, k = k_grid, g = g_grid)
    
    # Helper: simulation function with all non‑parameter arguments “frozen”.
    simfun <- function(d,k,g) simulate_moran(d,k,g,heter_vec,rounds,N,n_sims)
    
    # Evaluate the loss at each grid point (with progress bar).
    rss <- pbsapply(seq_len(nrow(grid)), \(i)
                    rss_loss(simfun(grid$d[i], grid$k[i], grid$g[i])))
    
    # Choose the grid point with minimal RSS as starting guess for optimization.
    start <- grid[which.min(rss), ]
    
    # Local continuous optimization (bounded) starting from the best grid point.
    opt   <- optim(par     = c(start$d, start$k, start$g),
                   fn      = \(par) rss_loss(simfun(par[1], par[2], par[3])),
                   method  = "L-BFGS-B",
                   lower   = c(min(d_grid), 0, 0),                     # k, γ ≥ 0
                   upper   = c(max(d_grid), max(k_grid), max(g_grid)), # upper bounds from grid
                   control = list(maxit = 80))
    
    # Return a list summarizing the calibrated heterogeneous model.
    list(model = "Heterogeneous",
         d = opt$par[1], k = opt$par[2], gamma = opt$par[3],
         rss = opt$value,
         T = simfun(opt$par[1], opt$par[2], opt$par[3]))
    
  } else {
    # 2D grid for homogeneous model (γ fixed at 0, heter_vec ignored).
    grid <- expand.grid(d = d_grid, k = k_grid)
    
    # Helper simulation function for homogeneous case (gamma=0, heter_vec=NULL).
    simfun <- function(d,k) simulate_moran(d,k,0,NULL,rounds,N,n_sims)
    
    # RSS over grid for (d, k).
    rss <- pbsapply(seq_len(nrow(grid)), \(i)
                    rss_loss(simfun(grid$d[i], grid$k[i])))
    
    # Best starting point.
    start <- grid[which.min(rss), ]
    
    # Bounded optimization over (d, k).
    opt   <- optim(par     = c(start$d, start$k),
                   fn      = \(par) rss_loss(simfun(par[1], par[2])),
                   method  = "L-BFGS-B",
                   lower   = c(min(d_grid), 0),                      # k ≥ 0
                   upper   = c(max(d_grid), max(k_grid)),
                   control = list(maxit = 80))
    
    # Return a list summarizing the calibrated homogeneous model.
    list(model = "Homogeneous",
         d = opt$par[1], k = opt$par[2], gamma = 0,
         rss = opt$value,
         T = simfun(opt$par[1], opt$par[2]))
  }
}

## 5.  Fit both specifications -----------------------------------------------
set.seed(2025)  # reproducible grid search + simulation randomness

# Fit homogeneous Moran calibration (γ = 0).
fit_plain <- calibrate_moran(FALSE, rounds, 5, 1500)

# Fit heterogeneous Moran calibration (γ free, using heter_vec).
fit_heter <- calibrate_moran(TRUE,  rounds, 5, 1500,
                             heter_vec = heter_vec)

# Collect and print a small summary table comparing the two fits:
#   model type, calibrated (d, k, γ), and RSS.
fits <- rbind(
  c(model = fit_plain$model, d = round(fit_plain$d,3),
    k = round(fit_plain$k,3), γ = 0,
    RSS = round(fit_plain$rss,4)),
  c(model = fit_heter$model, d = round(fit_heter$d,3),
    k = round(fit_heter$k,3), γ = round(fit_heter$gamma,3),
    RSS = round(fit_heter$rss,4))
)
print(fits, row.names = FALSE)

## 6.  Visual inspection -------------------------------------------------------
# Choose the specification with the smaller RSS as the “best fit”.
best_fit <- if(fit_heter$rss < fit_plain$rss) fit_heter else fit_plain

# Text used in the plot title to indicate whether heterogeneity is used.
label    <- if(best_fit$model == "Heterogeneous") "(heterogeneity)" else "(homogeneous)"

# Helper for plotting a 2×2 transition matrix as a heatmap with cell labels.
plot_heat <- function(M, title){
  as.data.frame(as.table(M)) |>
    ggplot(aes(Var1, Var2, fill = Freq)) +
    geom_tile() +                                                # colored tiles for probabilities
    geom_text(aes(label = sprintf("%.2f", Freq)), size = 5) +    # numeric labels in each cell
    labs(x = "Initial", y = "Final", title = title) +
    theme_minimal() +
    scale_fill_gradient(low = "white", high = "steelblue")
}

# Plot empirical round‑1→round‑10 transition matrix.
p_emp <- plot_heat(T_emp, "Empirical  (R1 → R10)")

# Plot calibrated transition matrix from the best Moran model.
p_sim <- plot_heat(best_fit$T,
                   sprintf("Calibrated Moran %s\n(d = %.2f, k = %.2f, γ = %.2f)",
                           label, best_fit$d, best_fit$k, best_fit$gamma))

# Display the two heatmaps side by side (using patchwork).
print(p_emp + p_sim)
