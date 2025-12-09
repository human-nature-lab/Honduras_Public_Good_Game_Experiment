# Initial data analysis for Section 3.3 (Honduras repeated public goods game data)

## ------------------------------------------------------------------
## 0. Packages
## ------------------------------------------------------------------
pkgs <- c(
  "dplyr", "tidyr", "moments", "ggplot2",
  "lme4", "lmerTest", "reshape2","patchwork"
)

new <- setdiff(pkgs, rownames(installed.packages()))
if (length(new)) install.packages(new, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))



# Load required packages
library(dplyr)       # data manipulation
library(tidyr)       # reshaping data
library(moments)     # skewness and kurtosis
library(ggplot2)     # plotting
library(reshape2)    # melting matrices for heatmaps
library(patchwork)   # combining ggplot objects

#------------------------------------------------------------
# 0. Load and prepare data
#------------------------------------------------------------

# Load the dataset
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# Ensure key variables are numeric / integer
pgg_data <- pgg_data %>%
  mutate(
    round_n      = as.integer(round_n),
    contributing = as.numeric(contributing)
  )

#------------------------------------------------------------
# 1. Round-level mean and variance (Figure 1a)
#------------------------------------------------------------

# Mean and variance by round
stats_by_round <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    mean_contrib = mean(contributing, na.rm = TRUE),
    var_contrib  = var(contributing,  na.rm = TRUE),
    .groups      = "drop"
  )

# Rescale variance so it can be plotted on the same axis as the mean
m_min <- min(stats_by_round$mean_contrib, na.rm = TRUE)
m_max <- max(stats_by_round$mean_contrib, na.rm = TRUE)
v_min <- min(stats_by_round$var_contrib,  na.rm = TRUE)
v_max <- max(stats_by_round$var_contrib,  na.rm = TRUE)

scale_to_mean <- function(v) {
  (v - v_min) / (v_max - v_min) * (m_max - m_min) + m_min
}

inv_scale_to_var <- function(y) {
  (y - m_min) / (m_max - m_min) * (v_max - v_min) + v_min
}

# Mean and variance on the same plot with a secondary y-axis (Figure 1a)
p <- ggplot(stats_by_round, aes(x = round_n)) +
  geom_line(aes(y = mean_contrib, color = "Mean")) +
  geom_point(aes(y = mean_contrib, color = "Mean")) +
  geom_line(aes(y = scale_to_mean(var_contrib), color = "Variance")) +
  geom_point(aes(y = scale_to_mean(var_contrib), color = "Variance")) +
  scale_x_continuous(breaks = stats_by_round$round_n) +
  scale_y_continuous(
    name = "Mean contribution",
    sec.axis = sec_axis(~ inv_scale_to_var(.), name = "Variance")
  ) +
  scale_color_manual(values = c("Mean" = "blue", "Variance" = "red")) +
  labs(
    title = "Mean and Variance of Contributions per Round",
    x     = "Round",
    color = ""
  ) +
  theme_minimal() +
  theme(legend.position = "top")

# Separate mean and variance plots (additional diagnostics)
mean_plot <- ggplot(stats_by_round, aes(x = round_n, y = mean_contrib)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = stats_by_round$round_n) +
  labs(
    title = "Mean Contribution per Round",
    x     = "Round",
    y     = "Mean contribution"
  ) +
  theme_minimal()

variance_plot <- ggplot(stats_by_round, aes(x = round_n, y = var_contrib)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = stats_by_round$round_n) +
  labs(
    title = "Variance of Contributions per Round",
    x     = "Round",
    y     = "Variance"
  ) +
  theme_minimal()

#------------------------------------------------------------
# 2. Detailed round-level descriptive statistics
#------------------------------------------------------------

desc_round <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    n          = n(),                                   # sample size per round
    mean       = mean(contributing,  na.rm = TRUE),     # average contribution
    sd         = sd(contributing,     na.rm = TRUE),    # standard deviation
    var        = var(contributing,    na.rm = TRUE),    # variance
    cv         = sd / mean,                             # coefficient of variation
    median     = median(contributing, na.rm = TRUE),    # median
    IQR        = IQR(contributing,    na.rm = TRUE),    # interquartile range
    min        = min(contributing,    na.rm = TRUE),    # minimum
    max        = max(contributing,    na.rm = TRUE),    # maximum
    skew       = skewness(contributing, na.rm = TRUE),  # skewness
    kurt       = kurtosis(contributing, na.rm = TRUE),  # kurtosis
    prop_zero  = mean(contributing == 0,  na.rm = TRUE),# share contributing 0
    prop_full  = mean(contributing == 12, na.rm = TRUE),# share contributing 12
    .groups    = "drop"
  )

print(desc_round, width = Inf)

#------------------------------------------------------------
# 3. Distributions by round (Figure 1b and boxplots)
#------------------------------------------------------------

# Boxplot of contributions by round
boxplot_round <- ggplot(pgg_data, aes(x = factor(round_n), y = contributing)) +
  geom_boxplot(outlier.alpha = 0.2) +
  labs(
    x     = "Round",
    y     = "Contribution",
    title = "Boxplot of Contributions by Round"
  ) +
  theme_minimal()

# Kernel density of contributions across rounds (Figure 1b)
distribution_plot <- ggplot(pgg_data, aes(x = contributing, fill = factor(round_n))) +
  geom_density(alpha = 0.3) +
  labs(
    x     = "Contribution",
    y     = "Density",
    fill  = "Round",
    title = "Density of Contributions Across Rounds"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

#------------------------------------------------------------
# 4. Round-to-round correlations
#------------------------------------------------------------

# Player-by-round matrix of contributions (columns are rounds)
wide_mat <- pgg_data %>%
  dplyr::select(player, round_n, contributing) %>%
  pivot_wider(
    names_from   = round_n,
    values_from  = contributing,
    names_prefix = "r"
  ) %>%
  dplyr::select(-player)

# Correlation matrix of contributions across rounds
round_cor <- cor(wide_mat, use = "pairwise.complete.obs")
print(round_cor)
# Correlation matrix of contributions across rounds
round_cor <- cor(wide_mat, use = "pairwise.complete.obs")
print(round_cor)

# Convert matrix to long format using reshape2::melt
round_cor_long <- reshape2::melt(round_cor)

# (Optional) Make sure the order of rounds is preserved on axes
round_cor_long$Var1 <- factor(round_cor_long$Var1, levels = colnames(round_cor))
round_cor_long$Var2 <- factor(round_cor_long$Var2, levels = rownames(round_cor))

# Heatmap of round-to-round correlations
correlation_rounds <- ggplot(round_cor_long, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(
    low      = "white",
    mid      = "lightblue",
    high     = "steelblue",
    midpoint = 0.5,
    limits   = c(0, 1)
  ) +
  labs(
    x     = "Round",
    y     = "Round",
    fill  = "Correlation",
    title = "Correlation Matrix: Contributions Across Rounds"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

print(correlation_rounds)


#------------------------------------------------------------
# 5. Counts of contributor types per round (Figures 1c & 1d)
#------------------------------------------------------------

# Low (< 6) vs High (≥ 6) contributors per round (Figure 1c)
summary_tbl <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    low_contributors  = sum(contributing < 6,  na.rm = TRUE),
    high_contributors = sum(contributing >= 6, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  pivot_longer(
    cols      = c(low_contributors, high_contributors),
    names_to  = "type",
    values_to = "count"
  )

plot1 <- ggplot(summary_tbl, aes(x = factor(round_n), y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.9)) +
  labs(
    title = "High vs Low Contributors by Round",
    x     = "Round",
    y     = "Count",
    fill  = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Free riders (0) vs full contributors (12) per round (Figure 1d)
summary_tbl <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    free_riders       = sum(contributing == 0,  na.rm = TRUE),
    full_contributors = sum(contributing == 12, na.rm = TRUE),
    .groups           = "drop"
  ) %>%
  pivot_longer(
    cols      = c(free_riders, full_contributors),
    names_to  = "type",
    values_to = "count"
  )

plot2 <- ggplot(summary_tbl, aes(x = factor(round_n), y = count, fill = type)) +
  geom_col(position = position_dodge(width = 0.9)) +
  labs(
    title = "Free Riders vs Full Contributors by Round",
    x     = "Round",
    y     = "Count",
    fill  = "Type"
  ) +
  theme_minimal() +
  theme(legend.position = "right")

# Four-panel figure for initial data analysis (mean/variance, density, counts)
combined_plot_rounds <- (p | distribution_plot) /
  (plot1 | plot2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Contributions by Round: Mean/Variance, Densities, and Headcounts"
  )

print(combined_plot_rounds)

# Vertical combination of contributor-type bar plots (types only)
combined_plot <- plot1 / plot2 +
  plot_annotation(
    title = "Contributor Types Across Rounds"
  )

print(combined_plot)

#------------------------------------------------------------
# 6. Individual-level distributions
#------------------------------------------------------------

round_summary <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    mean_contribution     = mean(contributing, na.rm = TRUE),
    variance_contribution = var(contributing,  na.rm = TRUE),
    .groups               = "drop"
  )

print(round_summary)

# Histograms of contributions by round
hist_individual_rounds <- ggplot(pgg_data, aes(x = contributing)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ round_n) +
  labs(
    title = "Contributions per Round",
    x     = "Contribution",
    y     = "Frequency"
  ) +
  theme_minimal()

#------------------------------------------------------------
# 7. Group-level distributions
#------------------------------------------------------------

# Group totals per round and their round-level mean and variance
group_totals_summary <- pgg_data %>%
  group_by(round_n, group) %>%
  summarise(
    total_group_contribution = sum(contributing, na.rm = TRUE),
    .groups                  = "drop"
  ) %>%
  group_by(round_n) %>%
  summarise(
    mean_contribution     = mean(total_group_contribution, na.rm = TRUE),
    variance_contribution = var(total_group_contribution,  na.rm = TRUE),
    .groups               = "drop"
  )

print(group_totals_summary)

# Histograms of total group contributions by round
hist_group_totals <- pgg_data %>%
  group_by(round_n, group) %>%
  summarise(
    total_group_contribution = sum(contributing, na.rm = TRUE),
    .groups                  = "drop"
  ) %>%
  ggplot(aes(x = total_group_contribution)) +
  geom_histogram(binwidth = 1, fill = "purple", color = "black", alpha = 0.7) +
  facet_wrap(~ round_n) +
  labs(
    title = "Total Group Contributions per Round",
    x     = "Total group contribution",
    y     = "Frequency"
  ) +
  theme_minimal() +
  theme(strip.text = element_blank())

# Time series of mean and variance of total group contributions
group_contribution <- ggplot(group_totals_summary, aes(x = round_n, y = mean_contribution)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Mean Total Group Contribution per Round",
    x     = "Round",
    y     = "Mean group contribution"
  ) +
  theme_minimal()

group_variance <- ggplot(group_totals_summary, aes(x = round_n, y = variance_contribution)) +
  geom_line() +
  geom_point() +
  labs(
    title = "Variance of Total Group Contribution per Round",
    x     = "Round",
    y     = "Variance of group contribution"
  ) +
  theme_minimal()

#------------------------------------------------------------
# 8. Changes in contributions across rounds
#------------------------------------------------------------

# Change in contributions from one round to the next for each player
contribution_changes <- pgg_data %>%
  arrange(player, round_n) %>%
  group_by(player) %>%
  mutate(
    contribution_change = contributing - lag(contributing)
  ) %>%
  ungroup() %>%
  filter(!is.na(contribution_change))  # remove first round per player

# Mean and variance of changes per round (rounds 2–10)
change_summary <- contribution_changes %>%
  group_by(round_n) %>%
  summarise(
    mean_change     = mean(contribution_change, na.rm = TRUE),
    variance_change = var(contribution_change,  na.rm = TRUE),
    .groups         = "drop"
  )

print(change_summary)

# Histograms of contribution changes by round
hist_changes <- ggplot(contribution_changes, aes(x = contribution_change)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black", alpha = 0.7) +
  facet_wrap(~ round_n, scales = "free_y", ncol = 3) +
  labs(
    x = "Contribution change",
    y = "Frequency"
  ) +
  theme_minimal() +
  theme(strip.text = element_blank())

#------------------------------------------------------------
# 9. Per-round counts of key contribution categories (tables)
#------------------------------------------------------------

# Free riders (0 contribution) per round
free_riders_summary <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    free_riders = sum(contributing == 0, na.rm = TRUE),
    .groups     = "drop"
  )

print(free_riders_summary)

# Full contributors (12 contribution) per round
Contributing_All_summary <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    Contributing_All = sum(contributing == 12, na.rm = TRUE),
    .groups          = "drop"
  )

print(Contributing_All_summary)

# Contributors with contribution ≤ 6 per round
below_6_summary <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    free_riders = sum(contributing <= 6, na.rm = TRUE),
    .groups     = "drop"
  )

print(below_6_summary)

# Contributors with contribution > 6 per round
above_6_summary <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    Contributing_All = sum(contributing > 6, na.rm = TRUE),
    .groups          = "drop"
  )

print(above_6_summary)

#------------------------------------------------------------
# 10. Optional: display additional plots one by one
#------------------------------------------------------------

# Example: print additional plots as needed
# print(boxplot_round)
# print(distribution_plot)
# print(correlation_rounds)
# print(hist_individual_rounds)
# print(hist_group_totals)
# print(group_contribution)
# print(group_variance)
# print(hist_changes)
