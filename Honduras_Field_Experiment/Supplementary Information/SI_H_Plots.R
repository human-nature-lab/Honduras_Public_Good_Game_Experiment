### APPENDIX H ###############################################################
### Illustrative village examples: gender, friends, religion, autonomy
###
### This script uses the main panel dataset from the Honduras public goods
### experiment (one row per player × round). It produces the descriptive
### plots and tables for Appendix H, focusing on:
###   - Gender differences (Villages 25 and 4),
###   - Friendship degree and contributions (Village 25),
###   - Religion (village-level shares and Village 4),
###   - Financial autonomy (Village 48, Group XG94A3).
###
### Assumptions:
###   * 'data_set.csv' is a cleaned panel with variables:
###       village_code : village identifier
###       group        : PGG group ID (within-village)
###       player       : unique player ID
###       round_n      : round index (1–10)
###       contributing : contribution in Lempiras (0–12)
###       gender       : 1 = man, 0 = woman
###       friends      : baseline friend degree (within-village network)
###       b0600        : religion (0 = none, 1 = Protestant, 2 = Catholic)
###       finauto_q_2_temp : 1 = financially autonomous, 0 = not autonomous
###   * tidyr::pivot_wider() is available (via library(tidyr) or tidyverse).
###############################################################################

## GENDER #####################################################################

library(dplyr)
library(ggplot2)
library(patchwork)

## Import main panel dataset (one row per player–round).
## This is the same structure described in Section 3 of the main text.
pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# ── helper: one horizontal strip of sub-plots for a given village ──────────
# For a given village_code:
#   * constructs one time-series panel per group in that village,
#   * plots each player's contribution path over the 10 rounds,
#   * colors each line/point by gender (1 = man, 0 = woman).
# Used below for visualizing gender patterns in Villages 25 and 4
# (see Appendix H gender examples).
plot_village <- function(v_code, df = pgg_data) {
  # Subset to focal village
  dat     <- df %>% filter(village_code == v_code)        # restrict to chosen village
  grp_ids <- sort(unique(dat$group))                     # all groups in that village
  
  # Build a list of ggplot objects, one panel per group in the village
  plot_list <- lapply(grp_ids, function(g) {
    ggplot(dat %>% filter(group == g),
           aes(x = round_n, y = contributing,
               group  = player,
               colour = factor(gender))) +
      geom_line(size = .7, show.legend = FALSE) +        # trajectory over rounds
      geom_point(size = 1, show.legend = FALSE) +        # per-round observations
      scale_x_continuous(breaks = 1:10) +                # rounds 1–10 on x-axis
      # Map gender codes to colors:
      #   gender = 1 (men)   → blue
      #   gender = 0 (women) → pink
      scale_colour_manual(values = c("1" = "blue", "0" = "pink")) +
      labs(
        title = paste("Village", v_code, "| Group", g),
        x     = "Round",
        y     = "Contribution (Lempiras)"
      ) +
      theme_minimal() +
      theme(
        plot.title   = element_text(size = 10, hjust = .5),
        axis.title.y = element_text(vjust = 1.2)
      )
  })
  
  # Arrange group-specific panels horizontally (one row per village)
  wrap_plots(plot_list, nrow = 1)
}

# ── Examples: visual inspection for two villages ───────────────────────────
# These plots underpin the gender illustrations for Villages 25 and 4
# used in Appendix H.
plot_village(25)
plot_village(4)

# Village-level summary: average contribution by gender, within village × group.
# Output: one row per (village_code, group), with columns Men and Women
# giving the mean contribution by gender. This is used for the gender tables
# (e.g., Table H.1 and Table H.3).
pgg_data %>%
  filter(village_code %in% c(25, 4)) %>%                 # restrict to two villages
  group_by(village_code, group, gender) %>%
  summarize(
    avg_contrib = mean(contributing, na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  mutate(gender = if_else(gender == 1, "Men", "Women")) %>%  # relabel gender codes
  pivot_wider(
    names_from  = gender,                                  # columns = Men / Women
    values_from = avg_contrib
  ) %>%
  arrange(village_code, group)                             # tidy ordering by village, then group

## FRIENDS ####################################################################
## Village 25: link between contribution level and friendship degree.
##
## Steps:
##   1. Compute each player's mean contribution across rounds (player-level average).
##   2. Classify players into two bins:
##        "Above6"  = avg_contrib > 6 Lempiras,
##        "AtMost6" = avg_contrib ≤ 6 Lempiras.
##   3. Within each group in Village 25, compute the mean number of friends
##      separately for Above6 vs AtMost6.
## The resulting table corresponds to the friendship summary reported in Appendix H
## (Village 25, Table H.2).

df_v25 <- pgg_data %>%
  filter(village_code == 25) %>%                    # keep only village 25
  group_by(player, village_code, group) %>%
  summarize(
    avg_contrib = mean(contributing, na.rm = TRUE), # player-level mean contribution
    friends     = first(friends),                   # baseline "friends" covariate (network degree)
    .groups     = "drop"
  ) %>%
  # Contribution bin relative to 6 Lempiras
  mutate(bin = if_else(avg_contrib > 6, "Above6", "AtMost6")) %>%   # contribution bin
  group_by(village_code, group, bin) %>%
  summarize(
    avg_friends = mean(friends, na.rm = TRUE),      # mean friends by bin and group
    .groups     = "drop"
  ) %>%
  pivot_wider(
    names_from  = bin,                              # columns: Above6 / AtMost6
    values_from = avg_friends
  ) %>%
  arrange(group)                                    # village_code is constant (25)

## RELIGION ###################################################################
## Religion variable: b0600
##   0 = no religion,
##   1 = Protestant,
##   2 = Catholic.
## This section:
##   1. Computes village-level shares by religion.
##   2. Identifies the villages with the highest Protestant and Catholic shares
##      (used to motivate the focus on Village 4, a strongly Catholic village).
##   3. Produces plots and group-level summaries by religion for Village 4
##      (Appendix H, religion examples).

library(dplyr)

# 1. Compute village-level religion shares.
village_religion <- pgg_data %>%
  filter(!is.na(b0600)) %>%                       # drop missing religion
  group_by(village_code) %>%
  summarise(
    prop_protestant = mean(b0600 == 1, na.rm = TRUE),   # share Protestant
    prop_catholic   = mean(b0600 == 2, na.rm = TRUE),   # share Catholic
    .groups         = "drop"
  )

# 2. Identify the village with the highest Protestant share.
max_protestant_village <- village_religion %>%
  filter(prop_protestant == max(prop_protestant)) %>%
  pull(village_code)

# 3. Identify the village with the highest Catholic share.
max_catholic_village <- village_religion %>%
  filter(prop_catholic == max(prop_catholic)) %>%
  pull(village_code)

# 4. Display results in the console (diagnostic / documentation).
cat("Village with highest proportion of Protestants: ", max_protestant_village, "\n")
cat("Village with highest proportion of Catholics:   ", max_catholic_village,   "\n")

library(dplyr)
library(ggplot2)
library(patchwork)

# ── helper: one horizontal strip of sub-plots for a given village ──────────
# As in plot_village(), but coloring by religion instead of gender:
#   b0600 = 0: no religion (gray),
#           1: Protestant (steelblue),
#           2: Catholic (dark green).
# Used below for Village 4, which has a high Catholic share.
plot_village_religion <- function(v_code, df = pgg_data) {
  
  dat     <- df %>% filter(village_code == v_code)  # restrict to chosen village
  grp_ids <- sort(unique(dat$group))               # all groups in that village
  
  plot_list <- lapply(grp_ids, function(g) {
    
    ggplot(dat %>% filter(group == g),
           aes(x = round_n,
               y = contributing,
               group  = player,
               colour = factor(b0600))) +          # use religion code as color
      geom_line(size = .7) +                       # individual time path
      geom_point(size = 1) +                       # per-round contributions
      scale_x_continuous(breaks = 1:10) +
      scale_colour_manual(
        name   = "Religion",
        values = c("0" = "gray50",    # no religion
                   "1" = "steelblue", # Protestant
                   "2" = "darkgreen"),
        labels = c("No religion", "Protestant", "Catholic")
      ) +
      labs(title = paste("Village", v_code, "| Group", g),
           x = "Round", y = "Contribution (Lempiras)") +
      theme_minimal() +
      theme(
        plot.title      = element_text(size = 10, hjust = .5),
        axis.title.y    = element_text(vjust = 1.2),
        legend.position = "bottom",
        legend.title    = element_text(size = 8),
        legend.text     = element_text(size = 7)
      )
  })
  
  # Arrange group-specific panels horizontally
  wrap_plots(plot_list, nrow = 1)
}

# Example usage: visual inspection in a strongly Catholic village (Village 4).
# This plot underlies the religion-based trajectory examples in Appendix H.
plot_village_religion(4)

# Group-level comparison within Village 4:
#   For each group, compute the average contribution for Catholics (b0600 == 2)
#   and for non-Catholics (b0600 != 2). This table corresponds to the
#   Catholic vs non-Catholic comparison reported in Appendix H (Table H.4).
pgg_data %>%
  filter(village_code == 4) %>%
  group_by(group) %>%
  summarize(
    avg_catholic     = mean(contributing[b0600 == 2], na.rm = TRUE),
    avg_non_catholic = mean(contributing[b0600 != 2], na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  arrange(group)

## FINANCIAL AUTONOMY ########################################################
## Village 48, Group XG94A3: contributions by financial autonomy.
##
## finauto_q_2_temp:
##   1 = financially autonomous,
##   0 = not autonomous.
## This section:
##   1. Plots individual trajectories for Group XG94A3, colored by autonomy.
##   2. Computes mean contributions separately for autonomous vs non-autonomous
##      players (Table H.5).

library(dplyr)
library(ggplot2)

# 1. Plot Village 48 · Group XG94A3
# Colors individual contribution paths by financial autonomy indicator
# (finauto_q_2_temp: 1 = financially autonomous, 0 = not autonomous).
pgg_data %>%
  filter(village_code == 48, group == "XG94A3") %>%
  ggplot(aes(
    x      = round_n,
    y      = contributing,
    group  = player,
    colour = factor(finauto_q_2_temp)           # 1 = autonomous, 0 = not autonomous
  )) +
  geom_line(size = 0.7, show.legend = FALSE) +  # trajectories; legend suppressed in figure
  geom_point(size = 1,  show.legend = FALSE) +
  scale_x_continuous(breaks = 1:10) +
  # Manual color mapping (gray vs black). The 'labels' argument is unused here
  # because show.legend = FALSE; colors only differentiate lines visually.
  scale_colour_manual(
    values = c("0" = "grey", "1" = "black"),
    labels = c("0" = "Women", "1" = "Men")
  ) +
  labs(
    title = "Village 48 · Group XG94A3",
    x     = "Round",
    y     = "Contribution (Lempiras)"
  ) +
  theme_minimal()

# 2. Table: average contribution by financial autonomy status
# Within Village 48, Group XG94A3:
#   compute mean contributions separately for autonomous vs non-autonomous.
# This reproduces the averages reported in Appendix H (Table H.5).
pgg_data %>%
  filter(
    village_code    == 48,
    group           == "XG94A3",
    !is.na(finauto_q_2_temp)
  ) %>%
  group_by(finauto_q_2_temp) %>%
  summarise(
    avg_contribution = mean(contributing, na.rm = TRUE),
    .groups          = "drop"
  ) %>%
  mutate(
    status = if_else(finauto_q_2_temp == 1,
                     "Autonomous", "Not Autonomous")
  ) %>%
  dplyr::select(status, avg_contribution)
