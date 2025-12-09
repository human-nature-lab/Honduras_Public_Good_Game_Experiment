###############################################################################
# Figure 3: Within-village social networks by gender (Villages 25 and 4)
#
# This script constructs gender-annotated friendship networks for two villages
# and plots them side by side (Village 25 on the left, Village 4 on the right).
# Node color encodes gender, node size encodes total contributions in the
# public good game, and edges represent reported social ties.
#
# Note: The network adjacency data (gr3_h_in) are not part of the public
# replication package and are available from the authors upon request.
# For this paper, only the panel data from the public good games (pgg_data)
# are released for submission.
###############################################################################

library(igraph)      # low-level network representation and layouts
library(tidygraph)   # tidy interface to igraph objects
library(ggraph)      # grammar of graphics for network plotting
library(dplyr)       # data manipulation (filter, mutate, summarize, etc.)
library(tibble)      # tibbles and deframe()
library(ggplot2)     # general plotting infrastructure
library(grid)        # grid graphics (margins, units)
library(patchwork)   # easy composition of multiple ggplots

pgg_data <- read.csv("data_set.csv", stringsAsFactors = FALSE)

# helper: linear rescale for sizes (robust to all-NA/constant vectors)
# - takes an input numeric vector (e.g. total contributions),
# - maps it linearly into a target size range "to",
# - handles the degenerate cases of all NA or zero variance by returning
#   a constant vector at the midpoint of "to".
.rescale_linear <- function(x, to = c(3, 10)) {
  if (all(is.na(x))) return(rep(mean(to), length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) == 0) return(rep(mean(to), length(x)))
  (x - rng[1]) / diff(rng) * diff(to) + to[1]
}

# helper: separate nodes of a target group if they are too close
# - coords: matrix of 2D layout coordinates (one row per node)
# - groups: vector of group labels (here, gender "M"/"F")
# - target: group whose nodes should be gently pushed apart ("M" = men)
# - min_sep: minimum allowed pairwise distance among target nodes
# - step: step size for each repulsive adjustment
# - n_steps: max number of iterations
# This is used only for Village 25 to improve readability of overlapping
# male nodes in the Fruchterman–Reingold layout.
.separate_group <- function(coords, groups, target = "M",
                            min_sep = 0.03, step = 0.02, n_steps = 50) {
  idx <- which(groups == target)
  if (length(idx) < 2) return(coords)
  for (iter in seq_len(n_steps)) {
    d <- as.matrix(dist(coords[idx, , drop = FALSE]))
    diag(d) <- Inf
    close_pairs <- which(d < min_sep, arr.ind = TRUE)
    if (!nrow(close_pairs)) break
    for (k in seq_len(nrow(close_pairs))) {
      i <- idx[close_pairs[k, 1]]
      j <- idx[close_pairs[k, 2]]
      v <- coords[i, ] - coords[j, ]
      if (!any(is.finite(v)) || all(v == 0)) v <- runif(2, -1, 1)
      v <- v / sqrt(sum(v^2))
      coords[i, ] <- coords[i, ] + step * v
      coords[j, ] <- coords[j, ] - step * v
    }
  }
  coords
}

# Main plotting function for one village network
# v_code: village_code (numeric index into network_list and pgg_data)
# network_list: list of igraph objects, indexed by village_code (e.g., gr3_h_in)
# df_players: player-level panel data (must contain village_code,
#             respondent_master_id, gender, total_contributions)
# node_alpha: transparency of nodes
# edge_alpha: transparency of edges
#
# The function:
#   1. extracts the igraph object for the target village,
#   2. removes isolates (degree 0),
#   3. merges in gender and total contributions,
#   4. sets node colors (blue = men, magenta = women),
#   5. sets node sizes proportional to total contributions,
#   6. sets edge colors depending on whether any endpoint is male/female,
#   7. computes a Fruchterman–Reingold layout (with a male-separation tweak
#      for Village 25),
#   8. renders the network with ggraph using a fixed manual layout.
plot_village_network <- function(
    v_code,
    network_list = gr3_h_in,
    df_players   = pgg_data,
    node_alpha   = 0.70,   # more transparent nodes
    edge_alpha   = 0.30    # slightly lighter edges
) {
  # --- 1) fetch igraph for this village ---
  g_raw <- network_list[[v_code]]
  if (is.null(g_raw)) stop("No network for village ", v_code)
  
  # --- 1b) drop isolates (degree == 0) ---
  # We focus on nodes with at least one reported tie to avoid plotting
  # isolated individuals as unconnected points.
  keep_vids <- which(igraph::degree(g_raw, mode = "all") > 0)
  if (length(keep_vids) == 0) stop("All nodes are isolates; nothing to plot for village ", v_code)
  g <- igraph::induced_subgraph(g_raw, vids = keep_vids)
  
  # --- 2) player metadata for this village (gender + total contributions) ---
  # respondent_master_id is used as the key to merge panel and network data.
  meta <- df_players %>%
    filter(village_code == v_code) %>%
    distinct(respondent_master_id, gender, total_contributions)
  
  # Map gender codes {0,1} (stored as strings) into "M" / "F" labels
  player_genders <- meta %>%
    mutate(gender = ifelse(gender == "1", "M", "F")) %>%
    select(respondent_master_id, gender) %>%
    deframe()
  
  # Total contributions across rounds (or pre-computed aggregate) per player
  player_totals <- meta %>%
    select(respondent_master_id, total_contributions) %>%
    deframe()
  
  # --- 3) node attributes ---
  # Assign gender to each vertex using its name (assumed to be respondent_master_id).
  V(g)$gender <- ifelse(V(g)$name %in% names(player_genders),
                        player_genders[V(g)$name], NA_character_)
  # Node color: blue for men, magenta for women, grey for missing/unknown.
  V(g)$col <- case_when(
    V(g)$gender == "M" ~ "blue",
    V(g)$gender == "F" ~ "magenta",
    TRUE               ~ "grey80"
  )
  # Node size proportional to total contributions (capped via rescaling).
  totals_vec     <- ifelse(V(g)$name %in% names(player_totals),
                           player_totals[V(g)$name], NA_real_)
  sizes_rescaled <- .rescale_linear(totals_vec, to = c(3, 10))
  # Only gender-identified nodes receive contribution-based sizing;
  # others get a small default radius.
  V(g)$size <- ifelse(!is.na(V(g)$gender), sizes_rescaled, 1.5)
  
  # --- 4) edge colours ---
  # Edge color heuristics:
  #   - blue if any endpoint is male,
  #   - magenta if no endpoint is male but at least one is female,
  #   - light grey otherwise (missing gender).
  E(g)$col <- sapply(E(g), function(e) {
    ends2 <- ends(g, e)
    g1 <- V(g)[ends2[1]]$gender
    g2 <- V(g)[ends2[2]]$gender
    if ((!is.na(g1) && g1 == "M") || (!is.na(g2) && g2 == "M")) {
      "blue"
    } else if ((!is.na(g1) && g1 == "F") || (!is.na(g2) && g2 == "F")) {
      "magenta"
    } else {
      "grey85"
    }
  })
  
  # --- 5) layout: start with FR; for Village 25 push male nodes apart slightly ---
  # Fruchterman–Reingold layout provides a force-directed embedding.
  # For Village 25, we apply the .separate_group() routine to reduce
  # overlap among male nodes, as described in the discussion of Figure 3.
  set.seed(42)
  fr <- igraph::layout_with_fr(g, niter = 1000, start.temp = 0.05)
  if (v_code == 25) {
    fr <- .separate_group(fr, V(g)$gender, target = "M",
                          min_sep = 0.035, step = 0.02, n_steps = 80)
  }
  
  # --- 6) plot with manual coordinates; zero expansions and margins ---
  # Convert igraph to a tidygraph object and use ggraph with a manual layout
  # (so that the layout is fixed and reproducible across runs).
  tbl <- as_tbl_graph(g)
  ggraph(tbl, layout = "manual", x = fr[, 1], y = fr[, 2]) +
    geom_edge_link(aes(colour = I(col)), alpha = edge_alpha) +
    geom_node_point(aes(colour = I(col), size = I(size)), alpha = node_alpha) +
    scale_colour_identity() +
    scale_size_identity() +
    scale_x_continuous(expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_equal(expand = FALSE, clip = "off") +
    theme_graph() +
    theme(
      plot.margin     = margin(10, 10, 10, 10),
      panel.spacing   = unit(10, "pt"),
      plot.background = element_rect(fill = NA, colour = NA)
    )
}

# Left: Village 25 with male-node de-overlap; Right: Village 4 unchanged
# This composition reproduces Figure 3 in the paper: a side-by-side comparison
# of within-village networks, highlighting gender composition and spatial
# segregation of men vs women in the friendship graph.
p25 <- plot_village_network(25)
p4  <- plot_village_network(4)

combined <- ((p25 | p4) +
               plot_layout(widths = c(1, 1)) +
               plot_annotation(theme = theme(plot.margin = margin(10, 10, 10, 10)))) &
  theme(panel.spacing = unit(10, "pt"),
        plot.margin   = margin(10, 10, 10, 10))

print(combined)
