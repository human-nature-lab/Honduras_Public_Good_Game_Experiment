###########################################################################
## Public goods game analysis: contribution heterogeneity and trajectories
##
## This script works with the main round-by-round dataset `pgg_data`
## containing (at least) the variables:
##   - player        : individual identifier
##   - round_n       : round index (1–10)
##   - contributing  : tokens contributed in that round (0–12)
##   - x             : endowment per round (maximum possible contribution)
##
## The code:
##   1. Describes variance and free-riding over rounds (Figure 1–style plots)
##   2. Plots counts at each exact contribution level 0–12 by round
##   3. Classifies players into behavioral types over the 10 rounds
##   4. Tracks how many players are in each category or contribution band
##      in each round (e.g., 0–3, 4–7, 8–10, 11–12 units).
###########################################################################

# -------------------------------------------------------------------------
# 0. Libraries and data
# -------------------------------------------------------------------------

# Load necessary libraries for data manipulation and plotting
library(dplyr)
library(ggplot2)
library(tidyr)

# `pgg_data` is expected to be a wide panel: one row per player–round.
# Columns used below:
#   - contributing: contribution in that round
#   - round_n     : round number
#   - x           : endowment / maximum possible contribution in that round

# Read the public goods game dataset from CSV.
# The file `data_set.csv` should include at least the variables listed above.
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# For the main paper, the experimental design fixes the endowment per round
# for all players; here we assume `x` is constant across players and rounds.

# -------------------------------------------------------------------------
# 1. Variance of contributions and free-riding over rounds
# -------------------------------------------------------------------------

# Round-level moments: variance and mean of individual contributions.
variance_per_round <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    variance_contribution = var(contributing, na.rm = TRUE),
    mean_contribution     = mean(contributing, na.rm = TRUE)
  )

# We now compute, for each round:
#   - the number of “free riders” (contributing exactly 0)
#   - the number of “full contributors” (contributing the full endowment x)
#   - their proportions within the group of players observed that round.

# Get the maximum possible contribution (assuming a constant endowment `x`).
# In the Honduras data this should equal 12 for all observations; if it does
# not, the code stops so that varying endowments can be handled explicitly.
max_contribution <- unique(pgg_data$x)
if (length(max_contribution) > 1) {
  stop("Endowment 'x' varies among players. Adjust the code to handle variable endowments.")
}

# Count free riders and full contributors by round.
contributor_counts <- pgg_data %>%
  group_by(round_n) %>%
  summarise(
    total_players        = n(),
    num_free_riders      = sum(contributing == 0,            na.rm = TRUE),
    num_full_contributors = sum(contributing == max_contribution, na.rm = TRUE)
  ) %>%
  mutate(
    prop_free_riders      = num_free_riders      / total_players,
    prop_full_contributors = num_full_contributors / total_players
  )

# Merge the round-level variance/mean with free-rider/full-contributor counts.
variance_contributor_data <- variance_per_round %>%
  inner_join(contributor_counts, by = "round_n")

# Plot variance of contributions over rounds
# (comparable to the variance trajectories reported in the paper).
ggplot(variance_contributor_data, aes(x = round_n, y = variance_contribution)) +
  geom_line(color = "blue") +
  geom_point(color = "blue") +
  labs(
    title = "Variance of Contributions Over Rounds",
    x     = "Round Number",
    y     = "Variance of Contributions"
  ) +
  theme_minimal()

# Plot the number of free riders over rounds
# (parallels the evolution of pure free-riders over the 10 rounds).
ggplot(variance_contributor_data, aes(x = round_n, y = num_free_riders)) +
  geom_line(color = "red") +
  geom_point(color = "red") +
  labs(
    title = "Number of Free Riders Over Rounds",
    x     = "Round Number",
    y     = "Number of Free Riders"
  ) +
  theme_minimal()

# Correlation between round-level variance and free-riding.
cor_variance_free_riders <- cor(
  variance_contributor_data$variance_contribution,
  variance_contributor_data$num_free_riders
)

cat("Correlation between variance and number of free riders:", cor_variance_free_riders, "\n")

# Simple linear regression of variance on the number of free riders per round.
# (Can be extended to include additional explanatory variables if desired.)
model_variance <- lm(variance_contribution ~ num_free_riders, data = variance_contributor_data)
summary(model_variance)

# -------------------------------------------------------------------------
# 2. Distribution of exact contribution levels 0–12 by round
# -------------------------------------------------------------------------

# Load necessary libraries (repeated here so this block can run standalone).
library(dplyr)
library(ggplot2)
library(gridExtra)

# We now describe how many players choose each exact contribution amount
# in {0, 1, …, 12} in each round, and then combine these with a panel showing
# the variance of contributions by round.

# Contribution amounts from 0 to 12 (consistent with a 12-unit endowment).
contribution_levels <- 0:12

# List to collect one ggplot object per contribution amount, plus the variance.
plot_list <- list()

# Grayscale palette: lighter grays for low amounts, darker for higher.
color_palette <- gray.colors(13, start = 0.9, end = 0.2)

for (amount in contribution_levels) {
  # For each contribution amount, count the number of players who chose
  # exactly that amount in each round.
  data_amount <- pgg_data %>%
    group_by(round_n) %>%
    summarise(count = sum(contributing == amount, na.rm = TRUE)) %>%
    ungroup()
  
  # Line plot for a single exact amount over rounds.
  p <- ggplot(data_amount, aes(x = round_n, y = count)) +
    geom_line(color = color_palette[amount + 1], size = 1) +
    geom_point(color = color_palette[amount + 1], size = 2) +
    labs(
      title = paste("Amount Contributed:", amount),
      x     = "Round",
      y     = "Number of Players"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 8),
      axis.text  = element_text(size = 7)
    )
  
  # Store plot in the list (index shift because 0-based amount).
  plot_list[[amount + 1]] <- p
}

# Round-level variance of contributions (same object as above but recomputed
# here so this block can run independently).
variance_data <- pgg_data %>%
  group_by(round_n) %>%
  summarise(variance = var(contributing, na.rm = TRUE)) %>%
  ungroup()

p_variance <- ggplot(variance_data, aes(x = round_n, y = variance)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "blue", size = 2) +
  labs(
    title = "Variance of Contributions per Round",
    x     = "Round",
    y     = "Variance"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 8),
    axis.text  = element_text(size = 7)
  )

# Add the variance panel as the 14th plot.
plot_list[[14]] <- p_variance

# Fill the remaining cells in the 4x4 grid with empty plots so that
# grid.arrange has exactly 16 grobs to place.
empty_plot <- ggplot() + theme_void()
plot_list[[15]] <- empty_plot
plot_list[[16]] <- empty_plot

# Arrange the 16 panels (13 contribution amounts + variance + 2 blanks)
# in a 4x4 grid.
library(gridExtra)
grid.arrange(grobs = plot_list, ncol = 4, nrow = 4)

# -------------------------------------------------------------------------
# 3. Player-level behavioral types over the 10 rounds
# -------------------------------------------------------------------------

# Load necessary libraries again so this block can be executed in isolation.
library(dplyr)
library(ggplot2)

# Constants for this game:
#   - max_contribution: full endowment per round (12 units)
#   - num_rounds      : total number of rounds observed (10)
max_contribution <- 12  # Maximum possible contribution
num_rounds       <- 10  # Total number of rounds

# Ensure `player`, `round_n`, and `contributing` are in appropriate formats.
pgg_data <- pgg_data %>%
  mutate(
    player       = as.factor(player),
    round_n      = as.integer(round_n),
    contributing = as.numeric(contributing)
  )

# ----- 3.1 Players with consistently very low or very high contributions -----

# Free Riders: players who contribute 0 in at least 80% of the rounds
# they played (persistent full free-riders).
free_riders <- pgg_data %>%
  group_by(player) %>%
  summarise(
    num_zero_contributions = sum(contributing == 0),
    total_rounds           = n()
  ) %>%
  mutate(
    proportion_zero = num_zero_contributions / total_rounds
  ) %>%
  filter(proportion_zero >= 0.8)  # ≥ 80% of rounds with zero contribution

# Full Contributors: players who contribute the full endowment (12 units)
# in at least 80% of the rounds they played.
full_contributors <- pgg_data %>%
  group_by(player) %>%
  summarise(
    num_full_contributions = sum(contributing == max_contribution),
    total_rounds           = n()
  ) %>%
  mutate(
    proportion_full = num_full_contributions / total_rounds
  ) %>%
  filter(proportion_full >= 0.8)  # ≥ 80% of rounds with full contribution

# ----- 3.2 Players with consistently high or medium contributions -----

# Consistent Cooperators: players who contribute between 8 and 11 units
# (high but not always full) in at least 80% of the rounds they played.
# We exclude full contributors so these remain a distinct category.
consistent_cooperators <- pgg_data %>%
  group_by(player) %>%
  summarise(
    num_high_contributions = sum(contributing >= 8 & contributing <= 11),
    total_rounds           = n()
  ) %>%
  mutate(
    proportion_high = num_high_contributions / total_rounds
  ) %>%
  filter(proportion_high >= 0.8) %>%  # ≥ 80% of rounds with high contributions
  anti_join(full_contributors, by = "player")  # Exclude full contributors

# Intermediate Contributors: players who mostly choose a medium range
# (4–7 units) in at least 80% of the rounds.
intermediate_contributors <- pgg_data %>%
  group_by(player) %>%
  summarise(
    num_medium_contributions = sum(contributing >= 4 & contributing <= 7),
    total_rounds             = n()
  ) %>%
  mutate(
    proportion_medium = num_medium_contributions / total_rounds
  ) %>%
  filter(proportion_medium >= 0.8)  # ≥ 80% of rounds with medium contributions

# ----- 3.3 Defectors (late drop in contributions) -----

# Defectors (strategic free riders): players who start with at least moderately
# high contributions but substantially reduce them in the latter half of the game.
defectors <- pgg_data %>%
  arrange(player, round_n) %>%
  group_by(player) %>%
  mutate(
    avg_first_half  = mean(contributing[round_n <= num_rounds / 2], na.rm = TRUE),
    avg_second_half = mean(contributing[round_n >  num_rounds / 2], na.rm = TRUE)
  ) %>%
  summarise(
    avg_first_half  = first(avg_first_half),
    avg_second_half = first(avg_second_half)
  ) %>%
  mutate(
    reduction = avg_first_half - avg_second_half
  ) %>%
  filter(
    avg_first_half >= 4,  # Started at least medium/high
    reduction     >= 4    # Large drop in the second half
  )

# ----- 3.4 Conditional cooperators and drifters -----

# Conditional Cooperators: players with substantial variation in contributions
# (high within-player standard deviation), who are not already in any of the
# previous categories. This is a simple variance-based proxy for conditional
# cooperation in response to others or to past outcomes.
conditional_cooperators <- pgg_data %>%
  group_by(player) %>%
  summarise(
    sd_contribution   = sd(contributing, na.rm = TRUE),
    mean_contribution = mean(contributing, na.rm = TRUE),
    total_rounds      = n()
  ) %>%
  filter(
    sd_contribution >= 3,               # High within-player variability
    total_rounds    >= num_rounds * 0.8 # Played at least 80% of rounds
  ) %>%
  anti_join(free_riders,              by = "player") %>%
  anti_join(full_contributors,        by = "player") %>%
  anti_join(consistent_cooperators,   by = "player") %>%
  anti_join(intermediate_contributors, by = "player") %>%
  anti_join(defectors,                by = "player")

# Drifters: players whose contributions are very erratic (even higher
# standard deviation) and who do not fit any of the structured patterns above.
drifters <- pgg_data %>%
  group_by(player) %>%
  summarise(
    sd_contribution = sd(contributing, na.rm = TRUE),
    total_rounds    = n()
  ) %>%
  filter(
    sd_contribution >= 4,               # Very high within-player variability
    total_rounds    >= num_rounds * 0.8 # Played at least 80% of rounds
  ) %>%
  anti_join(free_riders,              by = "player") %>%
  anti_join(full_contributors,        by = "player") %>%
  anti_join(consistent_cooperators,   by = "player") %>%
  anti_join(intermediate_contributors, by = "player") %>%
  anti_join(defectors,                by = "player") %>%
  anti_join(conditional_cooperators,  by = "player")

# ----- 3.5 Combine all player types and summarise -----

# All distinct players in the dataset.
all_players <- pgg_data %>%
  dplyr::select(player) %>%
  distinct()

# Assign each player to a unique behavioral category based on the above sets.
player_categories <- all_players %>%
  mutate(
    category = case_when(
      player %in% free_riders$player             ~ "Free Rider",
      player %in% full_contributors$player       ~ "Full Contributor",
      player %in% consistent_cooperators$player  ~ "Consistent Cooperator",
      player %in% intermediate_contributors$player ~ "Intermediate Contributor",
      player %in% defectors$player               ~ "Defector",
      player %in% conditional_cooperators$player ~ "Conditional Cooperator",
      player %in% drifters$player                ~ "Drifter",
      TRUE                                       ~ "Unclassified"
    )
  )

# Tabulate how many players fall in each behavioral category.
table(player_categories$category)

# ----- 3.6 Trajectories by player type -----

# Compute mean contribution by round within each player category.
plot_data <- pgg_data %>%
  inner_join(player_categories, by = "player") %>%
  group_by(category, round_n) %>%
  summarise(
    avg_contribution = mean(contributing, na.rm = TRUE),
    count            = n()
  ) %>%
  ungroup()

# Plot average contributions over rounds by category.
ggplot(plot_data, aes(x = round_n, y = avg_contribution, color = category)) +
  geom_line(size = 1) +
  labs(
    title = "Average Contributions Over Rounds by Player Category",
    x     = "Round Number",
    y     = "Average Contribution",
    color = "Category"
  ) +
  theme_minimal()

# Optionally export the classified player types to disk for further use.
# write.csv(player_categories, "player_categories.csv", row.names = FALSE)

# -------------------------------------------------------------------------
# 4. Round-level category counts and contribution bands
# -------------------------------------------------------------------------

# Load necessary libraries again so this block can be executed independently.
library(dplyr)
library(ggplot2)
library(gridExtra)

# ----- 4.1 Per-round type labels (cross-sectional categories) -----

# Here we re-classify *each player–round* (not the whole trajectory) into
# categories based on that round's contribution. This is a cross-sectional
# snapshot of the distribution of strategies in each round.
pgg_data_per_round <- pgg_data %>%
  mutate(
    category = case_when(
      contributing == 0                      ~ "Free Rider",
      contributing == max(contributing)     ~ "Full Contributor",  # full endowment in that round
      contributing >= 8  & contributing <= 11 ~ "Consistent Cooperator",
      contributing >= 4  & contributing <= 7  ~ "Intermediate Contributor",
      contributing >= 1  & contributing <= 3  ~ "Low Contributor",
      TRUE                                   ~ "Other"
    )
  )

# Count distinct players in each category by round.
category_counts_per_round <- pgg_data_per_round %>%
  group_by(round_n, category) %>%
  summarise(
    num_players = n_distinct(player)
  ) %>%
  ungroup()

# Order the factor levels to appear logically in the legend and plots.
category_counts_per_round$category <- factor(
  category_counts_per_round$category,
  levels = c(
    "Free Rider",
    "Low Contributor",
    "Intermediate Contributor",
    "Consistent Cooperator",
    "Full Contributor"
  )
)

# Plot the number of players in each category per round,
# showing how the composition of the population changes over time.
p_categories <- ggplot(category_counts_per_round, aes(x = round_n, y = num_players, color = category)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Number of Players in Each Category per Round",
    x     = "Round Number",
    y     = "Number of Players",
    color = "Category"
  ) +
  theme_minimal()

print(p_categories)

# ----- 4.2 Contribution bands (0–3, 4–7, 8–10, 11–12) -----

# Assign each player–round to one of the four contribution bands
# used in the descriptive analysis (0–3, 4–7, 8–10, 11–12 units).
pgg_data_ranges <- pgg_data %>%
  mutate(
    contribution_range = case_when(
      contributing >= 0  & contributing <= 3  ~ "0-3",
      contributing >= 4  & contributing <= 7  ~ "4-7",
      contributing >= 8  & contributing <= 10 ~ "8-10",
      contributing >= 11 & contributing <= 12 ~ "11-12",
      TRUE                                   ~ "Other"
    )
  )

# For each band and round, count how many distinct players fall into that band.
range_counts_per_round <- pgg_data_ranges %>%
  group_by(round_n, contribution_range) %>%
  summarise(
    num_players = n_distinct(player)
  ) %>%
  ungroup()

# Order contribution ranges for plotting.
range_counts_per_round$contribution_range <- factor(
  range_counts_per_round$contribution_range,
  levels = c("0-3", "4-7", "8-10", "11-12")
)

# Create a separate panel for each band showing its headcount over time.
plots <- list()
for (range in levels(range_counts_per_round$contribution_range)) {
  data_subset <- range_counts_per_round %>%
    filter(contribution_range == range)
  
  p <- ggplot(data_subset, aes(x = round_n, y = num_players)) +
    geom_line(color = "blue", size = 1) +
    geom_point(color = "blue", size = 2) +
    labs(
      title = paste("Number of Players Contributing", range, "per Round"),
      x     = "Round Number",
      y     = "Number of Players"
    ) +
    theme_minimal()
  
  plots[[range]] <- p
}

# Arrange the four band-specific plots in a 2x2 grid, matching the
# presentation of contribution bands over rounds.
grid.arrange(
  plots[["0-3"]], plots[["4-7"]],
  plots[["8-10"]], plots[["11-12"]],
  ncol = 2,
  top  = "Number of Players Contributing in Each Range per Round"
)
