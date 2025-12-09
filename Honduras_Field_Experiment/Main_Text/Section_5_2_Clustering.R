###############################################################################
# FULL-TRAJECTORY HIERARCHICAL CLUSTERING (10 rounds) – “BIFURCATION” TEST
# + Dynamic High/Low state logit (Eq. dyn_state_logit)
###############################################################################

## ------------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------------

pkgs <- c(
  "dplyr", "tidyr", "tibble",
  "fastcluster", "cluster", "factoextra",
  "RColorBrewer", "ggplot2", "ggrepel",
  "conflicted", "lme4", "broom.mixed"
)

new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new) > 0L) {
  install.packages(new, dependencies = TRUE)
}

invisible(lapply(pkgs, library, character.only = TRUE))

## Resolve naming conflicts in favour of dplyr
conflict_prefer("filter",      "dplyr")
conflict_prefer("select",      "dplyr")
conflict_prefer("lag",         "dplyr")
conflict_prefer("rename",      "dplyr")
conflict_prefer("summarise",   "dplyr")
conflict_prefer("count",       "dplyr")
conflict_prefer("inner_join",  "dplyr")
conflict_prefer("left_join",   "dplyr")

options(contrasts = c("contr.treatment", "contr.poly"))

## ------------------------------------------------------------------
## 1. Load data
## ------------------------------------------------------------------

pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

###############################################################################
## A. FULL-TRAJECTORY HIERARCHICAL CLUSTERING (10 rounds)
###############################################################################

###############################################################################
## A1. Build wide player × round matrix
###############################################################################

# Collapse accidental duplicates and average contributions per player × round
aggregated_contributions <- pgg_data %>%
  group_by(player, round_n) %>%
  summarise(
    contrib = mean(contributing, na.rm = TRUE),
    .groups = "drop"
  )

# Wide matrix: rows = players, columns = rounds r1–r10
contrib_mat <- aggregated_contributions %>%
  mutate(player = as.character(player)) %>%        # keep player as character ID
  pivot_wider(
    names_from   = round_n,
    values_from  = contrib,
    names_prefix = "r"
  ) %>%
  drop_na() %>%                                    # retain only complete 10-round trajectories
  column_to_rownames("player") %>%
  as.matrix()

# Z-score each round (Ward distances are scale-sensitive)
contrib_sc <- scale(contrib_mat)

###############################################################################
## A2. Ward-D2 hierarchical clustering
###############################################################################

# Euclidean distances between 10-round trajectories
dist_mat <- dist(contrib_sc, method = "euclidean")

# Ward.D2 agglomerative clustering (fast implementation)
hc_tree <- fastcluster::hclust(dist_mat, method = "ward.D2")

###############################################################################
## A3. Silhouette widths & choice of k
###############################################################################

# Average silhouette width for a given number of clusters k
get_sil <- function(k) {
  sil <- cluster::silhouette(cutree(hc_tree, k), dist_mat)
  mean(sil[, 3])
}

# Evaluate k in {2, …, 6}
ks         <- 2:6
sil_widths <- setNames(vapply(ks, get_sil, numeric(1)), ks)

# Report average silhouette widths
print(round(sil_widths, 3))

# Select best k by maximum average silhouette
best_k <- ks[which.max(sil_widths)]    # typically best_k = 2

###############################################################################
## A4. Attach HC cluster tag to every player
###############################################################################

# Cluster labels C1, C2, … for each player
hc_tag <- factor(
  cutree(hc_tree, k = best_k),
  labels = paste0("C", seq_len(best_k))
)

# Data frame with player ID and assigned hierarchical cluster
contrib_df <- tibble(
  player     = rownames(contrib_sc),
  hc_cluster = hc_tag
)

###############################################################################
## A5. End-state labels (LL/LH/HL/HH) & confusion matrix
###############################################################################

# Threshold for High/Low based on mean round-1 contribution
thresh <- mean(
  pgg_data$contributing[pgg_data$round_n == 1],
  na.rm = TRUE
)

# Classify each player as H/L in rounds 1 and 10, then define transition type
states <- pgg_data %>%
  filter(round_n %in% c(1, 10)) %>%
  mutate(
    state = if_else(contributing > thresh, "H", "L")
  ) %>%
  dplyr::select(player, round_n, state) %>%
  pivot_wider(
    names_from   = round_n,
    names_prefix = "r",
    values_from  = state
  ) %>%
  mutate(
    trans  = paste0(r1, r10),               # LL, LH, HL, HH
    player = as.character(player)           # same type as in contrib_df
  )

# Cross-tabulate HC clusters with end-state transitions
confusion <- states %>%
  inner_join(contrib_df, by = "player") %>%
  count(hc_cluster, trans) %>%
  pivot_wider(
    names_from  = trans,
    values_from = n,
    values_fill = 0
  )

cat("\n--- Confusion matrix: End-state vs HC clusters (k =", best_k, ") ---\n")
print(confusion)

###############################################################################
## A6. Diagnostics plots (optional)
###############################################################################

# 6-A. Silhouette plot for the selected k
fviz_silhouette(
  cluster::silhouette(cutree(hc_tree, best_k), dist_mat)
)

# 6-B. Dendrogram for a random sample (≤ 2 000 trajectories)
set.seed(1)
sample_idx <- sample(
  seq_len(nrow(contrib_sc)),
  size = min(2000, nrow(contrib_sc))
)

# Palette for cluster colours (handles k < 3 without warnings)
palette_k <- RColorBrewer::brewer.pal(max(best_k, 3), "Set2")[seq_len(best_k)]

fviz_dend(
  hc_tree,
  k           = best_k,
  k_colors    = palette_k,
  rect        = TRUE,
  show_labels = FALSE,
  select_rows = sample_idx,
  main        = paste("Ward hierarchical clustering – k =", best_k)
)

###############################################################################
## A7. Mean trajectory per HC cluster
###############################################################################

# Mean z-scored contribution path by hierarchical cluster
traj_df <- contrib_df %>%
  left_join(as_tibble(contrib_sc, rownames = "player"), by = "player") %>%
  pivot_longer(
    cols      = starts_with("r"),
    names_to  = "round",
    values_to = "contrib"
  ) %>%
  group_by(hc_cluster, round) %>%
  summarise(
    mean_contrib = mean(contrib),
    .groups      = "drop"
  ) %>%
  mutate(
    round = as.integer(sub("r", "", round))
  )

# Simple mean trajectory plot by HC cluster
ggplot(
  traj_df,
  aes(x = round, y = mean_contrib,
      colour = hc_cluster, group = hc_cluster)
) +
  geom_line(size = 1.1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Average contribution path by HC cluster",
    x     = "Round",
    y     = "Mean contribution (z-scored)"
  ) +
  theme_minimal()

# Map cluster codes to descriptive labels and colour-blind-safe palette
lab_map <- c(C1 = "High contributors", C2 = "Low contributors")
pal     <- c(C1 = "#E07A5F", C2 = "#2A9D8F")  # coral / teal

# Points at last round for direct labels
last_pts <- traj_df %>%
  group_by(hc_cluster) %>%
  filter(round == max(round)) %>%
  ungroup() %>%
  mutate(label = lab_map[hc_cluster])

# Mean trajectory plot with direct labels at the end of each line
ggplot(
  traj_df,
  aes(x = round, y = mean_contrib,
      colour = hc_cluster, group = hc_cluster)
) +
  geom_line(size = 1.4, lineend = "round") +
  geom_point(size = 2.6) +
  scale_x_continuous(
    breaks = 1:10,
    expand = expansion(add = c(0.1, 0.6))
  ) +
  scale_color_manual(values = pal, guide = "none") +
  labs(
    title    = "Average contribution path by cluster",
    subtitle = "Ten rounds; z-scored contributions",
    x        = "Round",
    y        = "Mean contribution (z-scored)"
  ) +
  geom_text_repel(
    data          = last_pts,
    aes(label     = label),
    nudge_x       = 0.25,
    segment.color = NA,
    size          = 4.2,
    fontface      = "bold"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(colour = "#5A5A5A"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line          = element_line(colour = "grey60", linewidth = 0.3),
    axis.ticks         = element_line(colour = "grey70")
  )

###############################################################################
## A8. Append cluster label to the main data set
###############################################################################

pgg_data <- pgg_data %>%
  mutate(player = as.character(player)) %>%
  left_join(contrib_df, by = "player")

############################  QUICK COMMENTARY  ############################
# 1) Silhouette widths
#    Higher values indicate better average separation between clusters.
#    Example (from one run):
#      k = 2 : 0.362  (best solution)
#      k = 3 : 0.258
#      k = 4 : 0.218
#      k = 5 : 0.176
#      k = 6 : 0.152
#    The two-cluster solution provides the clearest structure in the data.
#
# 2) Confusion matrix (HC clusters vs end-state LL/LH/HL/HH)
#    Based on the printed table, one cluster is predominantly LL (persistently
#    low) while the other contains almost all HH plus most movers (HL/LH):
#    consistent with a bifurcation into High/mixed vs Low regimes.
###########################################################################


###############################################################################
## B. Dynamic High/Low state logit (Eq. dyn_state_logit)
###############################################################################
## logit Pr(s_igt = 1 | s_{i,t-1}, m_{g,t-1}, X_i, α_v, δ_g, a_i)
##   = β0 + ρ s_{i,t-1} + λ m_{g,t-1} + θ t + X_i'γ + α_{v(i)} + δ_{g(i)} + a_i
##
## m_{g,t-1} = lagged leave‑one‑out peer mean / 12
## s_igt = 1{c_igt >= c*}, where c* is the round‑1 mean contribution.
## Wooldridge-style initial conditions: include s1 and peer_mean_c_tm1_bar.
###############################################################################

## B1. High/Low state and base panel -----------------------------------------

# Round‑1 threshold c*
c_star <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::summarise(c_star = mean(contributing, na.rm = TRUE)) %>%
  dplyr::pull(c_star)

# Panel with only variables needed to build dynamics and peers
panel <- pgg_data %>%
  dplyr::filter(is.finite(contributing)) %>%
  dplyr::transmute(
    player       = as.character(player),
    village_code = as.character(village_code),
    group        = as.character(group),
    round_n      = round_n,
    contributing = contributing,
    s            = as.integer(contributing >= c_star)  # High/Low
  )

## B2. Leave‑one‑out peer means and lags -------------------------------------

panel <- panel %>%
  dplyr::group_by(group, round_n) %>%
  dplyr::mutate(
    G = dplyr::n(),
    peer_mean_c_t = dplyr::if_else(
      G > 1,
      (sum(contributing, na.rm = TRUE) - contributing) / (G - 1),
      NA_real_
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, round_n) %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(peer_mean_c_tm1 = dplyr::lag(peer_mean_c_t)) %>%
  dplyr::ungroup()

## B3. Own‑state lag and round counter ---------------------------------------

panel <- panel %>%
  dplyr::arrange(player, round_n) %>%
  dplyr::group_by(player) %>%
  dplyr::mutate(s_lag = dplyr::lag(s)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(r_t = as.integer(round_n))

## B4. Wooldridge initial‑conditions regressors ------------------------------

# s1: initial state (round 1)
s1_df <- panel %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::select(player, s1 = s)

# peer_mean_c_tm1_bar: player‑average lagged peer mean (level scale)
peerbar_df <- panel %>%
  dplyr::group_by(player) %>%
  dplyr::summarise(
    peer_mean_c_tm1_bar = mean(peer_mean_c_tm1, na.rm = TRUE),
    .groups = "drop"
  )

## B5. Baseline covariates X_i (round 1) -------------------------------------

base <- pgg_data %>%
  dplyr::filter(round_n == 1) %>%
  dplyr::distinct(player, village_code, .keep_all = TRUE) %>%
  dplyr::mutate(
    player       = as.character(player),
    village_code = as.character(village_code),
    
    # Simple male dummy; adjust to your coding if needed
    male  = as.integer(tolower(as.character(gender)) %in%
                         c("1", "male", "man", "m")),
    
    # Continuous covariates (kept in original units)
    age                 = suppressWarnings(as.numeric(age)),
    friends             = suppressWarnings(as.numeric(friends)),
    adversaries         = suppressWarnings(as.numeric(adversaries)),
    network_density_fr  = suppressWarnings(as.numeric(network_density_fr)),
    network_density_adv = suppressWarnings(as.numeric(network_density_adv)),
    network_size        = suppressWarnings(as.numeric(network_size)),
    b0100               = suppressWarnings(as.numeric(b0100)),  # education 0–13
    FI                  = suppressWarnings(as.numeric(FI)),
    
    # Discrete covariates
    marital_status = factor(marital_status),
    b0200          = factor(b0200),  # Indigenous dummy
    b0600          = factor(b0600),  # religion dummies
    access_routes  = factor(access_routes)
  )

## B6. Merge everything into analysis frame ----------------------------------

panel <- panel %>%
  dplyr::left_join(
    base %>%
      dplyr::select(
        player, village_code,
        male, age, friends, adversaries,
        network_density_fr, network_density_adv, network_size,
        b0100, FI,
        marital_status, b0200, b0600, access_routes
      ),
    by = c("player", "village_code")
  ) %>%
  dplyr::left_join(s1_df,     by = "player") %>%
  dplyr::left_join(peerbar_df, by = "player") %>%
  dplyr::mutate(
    # m_{g,t-1}: lagged LOO peer mean scaled by 12L (one endowment)
    peer_mean_c_tm1_scaled = peer_mean_c_tm1 / 12
  )

# Keep usable rows: need s, s_lag, m_{g,t-1}, s1, peer_mean_c_tm1_bar, and X_i
m1_df <- panel %>%
  dplyr::filter(
    is.finite(s),
    is.finite(s_lag),
    is.finite(peer_mean_c_tm1_scaled),
    is.finite(s1),
    is.finite(peer_mean_c_tm1_bar)
  ) %>%
  tidyr::drop_na(
    male, age, friends, adversaries,
    network_density_fr, network_density_adv, network_size,
    b0100, FI,
    marital_status, b0200, b0600, access_routes
  ) %>%
  dplyr::mutate(
    village_code = factor(village_code),
    group        = factor(group),
    player       = factor(player)
  )

## B7. Dynamic state logit GLMM (Eq. dyn_state_logit) ------------------------


m1_df$male <- as.factor(m1_df$male)
m1_df$access_routes <- as.factor(m1_df$access_routes)
m1_df$FI <- as.factor(m1_df$FI)
m1_df$b0600 <- as.factor(m1_df$b0600)
m1_df$b0200 <- as.factor(m1_df$b0200)
m1_df$marital_status <- as.factor(m1_df$marital_status)

m_dyn_state <- glmer(
  s ~
    # dynamics
    s_lag +
    peer_mean_c_tm1_scaled +  # m_{g,t-1} (per 12L)
    r_t +
    # initial-conditions adjustments (Wooldridge)
    s1 + peer_mean_c_tm1_bar +
    # baseline covariates X_i
    male + age + friends + adversaries +
    network_density_fr + network_density_adv + network_size +
    b0100 + FI +
    marital_status + b0200 + b0600 + as.numeric(access_routes) +
    # random intercepts
    (1 | village_code) +
    (1 | group) +
    (1 | player),
  data    = m1_df,
  family  = binomial(link = "logit"),
  control = glmerControl(
    optimizer   = "nloptwrap",
    calc.derivs = FALSE,
    optCtrl     = list(maxfun = 2e4)
  )
)

summary(m_dyn_state)

## B8. Odds ratios (including peer effect) ------------------------------------

cat("\n--- Dynamic state logit: fixed effects as odds ratios ---\n")
or_table <- broom.mixed::tidy(
  m_dyn_state,
  effects      = "fixed",
  conf.int     = TRUE,
  exponentiate = TRUE
) %>%
  dplyr::select(term, estimate, conf.low, conf.high, p.value)

print(or_table, row.names = FALSE)

# Peer effect: OR per 12L and per 1L of lagged peer mean
beta_peer  <- lme4::fixef(m_dyn_state)["peer_mean_c_tm1_scaled"]
or_per_12L <- exp(beta_peer)        # one full endowment (12L)
or_per_1L  <- or_per_12L^(1 / 12)   # per 1L change

cat(sprintf(
  "\nOR for lagged peer mean m_{g,t-1} (per 12L): %.3f\n", or_per_12L
))
cat(sprintf(
  "OR for lagged peer mean m_{g,t-1} (per 1L):  %.3f\n", or_per_1L
))

# Optional: singularity check
if (lme4::isSingular(m_dyn_state, tol = 1e-4)) {
  message("Note: singular fit detected (some random-effect variance ~ 0).")
}
