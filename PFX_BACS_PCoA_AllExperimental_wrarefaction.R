# Script to perform PCoA ordination using ILR Aitchison distance 
# on 16S SSU rRNA amplicon sequencing data
# Analyzing all samples but plotting only selected groups
# Using black/orange color scheme with different shapes/shades

setwd("/home/pacifica/R/antarctica/S03/big_experiment/")

# Load required libraries
library(tidyverse)     # For data manipulation and visualization
library(vegan)         # For PERMANOVA
library(zCompositions) # for cmultRepl Bayesian-multiplicative replacement of count zeros
library(compositions)  # For ILR transformation and Aitchison distance
library(ggplot2)       # For plotting
library(scales)        # For color manipulation

# Define input file paths
count_table_path <- "ESV_table_PFX_2025.csv"
metadata_path <- "Metadata_PFX_2025.csv"

# Read in the data
count_data <- read.csv(count_table_path, row.names = 1)
metadata <- read.csv(metadata_path, row.names = 1)

# Add sample names as a column in metadata
metadata$Sample <- rownames(metadata)

# Create a new column in metadata that combines Order and Timing
metadata$Group <- paste(metadata$Order, metadata$Timing, sep = ".")

# Define the groups to be plotted (but analyze all samples)
specified_groups <- c("BB.0", "BO.1", "BO.2", "OO.0", "OB.1", "OB.2")

# Step 1: Reduce zeros by removing ESVs with fewer than 5 reads on average
# Use all samples for this calculation
esv_means <- rowMeans(count_data)
filtered_esvs <- count_data[esv_means >= 5, ]

# Step 2: Work directly with abundance data (no percent conversion)
abundance_data <- filtered_esvs

cat("Using abundance data for Aitchison distance calculation\n")
cat("Data dimensions:", dim(abundance_data), "\n")
cat("Data range - Min:", min(abundance_data), "Max:", max(abundance_data), "\n")

# Step 3: Calculate ILR Aitchison distance for all samples using abundance data
# First, handle zeros for the ILR transformation
abundance_replaced <- cmultRepl(abundance_data, label=0, method="CZM")
# Note: this gets rid of samples without notable abundance, including all the negative controls

# REMEMBER: IF samples are columns, margin = 2, if samples are ROWS, margin = 1

# Transpose to get samples as rows
abundance_t <- t(abundance_replaced)

# Convert to a compositions object for ILR transformation
comp_data <- acomp(abundance_t)

# Perform ILR transformation
ilr_data <- ilr(comp_data)

# Calculate Euclidean distance on ILR-transformed data (equivalent to Aitchison distance)
aitchison_dist <- dist(ilr_data, method = "euclidean")

# Step 4: Run PCoA on Aitchison distance for all samples
pcoa_result <- cmdscale(aitchison_dist, k = 2, eig = TRUE)

# Calculate percent variance explained by each axis
# eigenvalues can include negative values (from non-Euclidean spaces);
# use only positive eigenvalues for percent variance calculation
positive_eigs <- pcoa_result$eig[pcoa_result$eig > 0]
total_variance <- sum(positive_eigs)
pct_var <- round(100 * positive_eigs / total_variance, 1)

pct_axis1 <- pct_var[1]
pct_axis2 <- pct_var[2]

cat("PCoA Axis 1 variance explained:", pct_axis1, "%\n")
cat("PCoA Axis 2 variance explained:", pct_axis2, "%\n")

# Create axis labels including percent variance
axis1_label <- paste0("PCoA1 (", pct_axis1, "%)")
axis2_label <- paste0("PCoA2 (", pct_axis2, "%)")

# Extract PCoA scores for all samples
pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
pcoa_scores$SampleID <- rownames(pcoa_scores)

# Join with metadata
pcoa_all_data <- left_join(pcoa_scores, metadata %>% mutate(SampleID=Sample), by = "SampleID")

# Filter to only include the specified groups for plotting
pcoa_plot_data <- pcoa_all_data %>%
  filter(Group %in% specified_groups)

# Create a named vector for group labels (for plot legend)
group_labels <- c(
  "BB.0" = "Black mat before experiment",
  "BO.1" = "Black then orange mat after 1 season",
  "BO.2" = "Black then orange mat after 2 seasons",
  "OO.0" = "Orange mat before experiment",
  "OB.1" = "Orange then black mat after 1 season",
  "OB.2" = "Orange then black mat after 2 seasons"
)

# Convert Group to a factor with the specified order for the plot data
pcoa_plot_data$Group <- factor(pcoa_plot_data$Group, levels = specified_groups)

# Define custom color and shape scheme
custom_colors <- c(
  "BB.0" = "black",
  "BO.1" = alpha("black", 0.7),
  "BO.2" = alpha("black", 0.4),
  "OO.0" = "#FF8C00",
  "OB.1" = "#FFA500",
  "OB.2" = "#FFD700"
)

# Custom shapes for each group
custom_shapes <- c(
  "BB.0" = 16,
  "BO.1" = 17,
  "BO.2" = 15,
  "OO.0" = 16,
  "OB.1" = 17,
  "OB.2" = 15
)

# Function to create convex hull data frames for each group
find_hull <- function(df) {
  hull_points <- chull(df$PCoA1, df$PCoA2)
  hull_points <- c(hull_points, hull_points[1])
  return(df[hull_points, ])
}

# Create convex hull data for each group
hulls <- pcoa_plot_data %>%
  group_by(Group) %>%
  do(find_hull(.)) %>%
  ungroup()

# Calculate group centroids for arrow placement AND centroid markers
group_centroids <- pcoa_plot_data %>%
  group_by(Group) %>%
  summarise(
    centroid_x = mean(PCoA1),
    centroid_y = mean(PCoA2),
    .groups = 'drop'
  )

# Function to calculate arrow endpoints on circle circumference
calculate_arrow_endpoints <- function(x1, y1, x2, y2, radius) {
  angle <- atan2(y2 - y1, x2 - x1)
  start_x <- x1 + radius * cos(angle)
  start_y <- y1 + radius * sin(angle)
  end_x <- x2 - radius * cos(angle)
  end_y <- y2 - radius * sin(angle)
  return(list(start_x = start_x, start_y = start_y, end_x = end_x, end_y = end_y))
}

# Define circle radius (adjust based on your plot scale)
circle_radius <- 1.8

# Create arrow data frames with adjusted endpoints
# Black arrows (BB.0 -> BO.1 -> BO.2)
bb0_centroid <- group_centroids[group_centroids$Group == "BB.0", ]
bo1_centroid <- group_centroids[group_centroids$Group == "BO.1", ]
bo2_centroid <- group_centroids[group_centroids$Group == "BO.2", ]

arrow1 <- calculate_arrow_endpoints(bb0_centroid$centroid_x, bb0_centroid$centroid_y,
                                    bo1_centroid$centroid_x, bo1_centroid$centroid_y,
                                    circle_radius)
arrow2 <- calculate_arrow_endpoints(bo1_centroid$centroid_x, bo1_centroid$centroid_y,
                                    bo2_centroid$centroid_x, bo2_centroid$centroid_y,
                                    circle_radius)

black_arrows <- data.frame(
  x = c(arrow1$start_x, arrow2$start_x),
  y = c(arrow1$start_y, arrow2$start_y),
  xend = c(arrow1$end_x, arrow2$end_x),
  yend = c(arrow1$end_y, arrow2$end_y),
  arrow_type = c("BB.0_to_BO.1", "BO.1_to_BO.2")
)

# Orange arrows (OO.0 -> OB.1 -> OB.2)
oo0_centroid <- group_centroids[group_centroids$Group == "OO.0", ]
ob1_centroid <- group_centroids[group_centroids$Group == "OB.1", ]
ob2_centroid <- group_centroids[group_centroids$Group == "OB.2", ]

arrow3 <- calculate_arrow_endpoints(oo0_centroid$centroid_x, oo0_centroid$centroid_y,
                                    ob1_centroid$centroid_x, ob1_centroid$centroid_y,
                                    circle_radius)
arrow4 <- calculate_arrow_endpoints(ob1_centroid$centroid_x, ob1_centroid$centroid_y,
                                    ob2_centroid$centroid_x, ob2_centroid$centroid_y,
                                    circle_radius)

orange_arrows <- data.frame(
  x = c(arrow3$start_x, arrow4$start_x),
  y = c(arrow3$start_y, arrow4$start_y),
  xend = c(arrow3$end_x, arrow4$end_x),
  yend = c(arrow3$end_y, arrow4$end_y),
  arrow_type = c("OO.0_to_OB.1", "OB.1_to_OB.2")
)

# Step 5: Create the PCoA plot with custom colors and shapes
pcoa_plot <- ggplot() +
  # Add gray points for all samples (background)
  geom_point(data = pcoa_all_data, 
             aes(x = PCoA1, y = PCoA2), 
             color = "gray90", 
             size = 1.5, 
             alpha = 0.2) +
  # Add convex hulls for each group
  geom_polygon(data = hulls, 
               aes(x = PCoA1, y = PCoA2, 
                   fill = Group, 
                   color = Group, 
                   group = Group), 
               alpha = 0.1,
               linetype = 2) +
  # Add black arrows for BB transitions
  geom_segment(data = black_arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black",
               size = 1.2,
               alpha = 0.8) +
  # Add orange arrows for OO transitions
  geom_segment(data = orange_arrows,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#FF8C00",
               size = 1.2,
               alpha = 0.8) +
  # Add colored points for the selected groups with custom shapes
  geom_point(data = pcoa_plot_data, 
             aes(x = PCoA1, y = PCoA2, 
                 color = Group, 
                 shape = Group), 
             size = 2, 
             alpha = 0.9) +
  # Add open circles at group centroids
  geom_point(data = group_centroids,
             aes(x = centroid_x, y = centroid_y, 
                 color = Group),
             shape = 1,
             size = 4,
             stroke = 1.5) +
  # Use custom colors and shapes
  scale_color_manual(values = custom_colors, labels = group_labels) +
  scale_shape_manual(values = custom_shapes, labels = group_labels) +
  scale_fill_manual(values = custom_colors, labels = group_labels) +
  # Apply minimal theme
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_line(color = "gray98")
  ) +
  # Add axis labels with percent variance explained
  labs(
    title = "",
    subtitle = "",
    x = axis1_label,
    y = axis2_label,
    color = "Sample Group",
    shape = "Sample Group",
    fill = "Sample Group"
  ) +
  # Hide the fill from the legend since it's redundant with color
  guides(fill = "none")

# Print the PCoA plot
print(pcoa_plot)

# Save the PCoA plot
ggsave(
  "pcoa_ilr_aitchison_abundance_black_orange_convexhull_arrows.pdf", 
  pcoa_plot, 
  width = 8, 
  height = 5
)

# Print information about the PCoA result
cat("PCoA Axis 1:", pct_axis1, "% variance explained\n")
cat("PCoA Axis 2:", pct_axis2, "% variance explained\n")
cat("Total variance explained by first 2 axes:", pct_axis1 + pct_axis2, "%\n")
cat("Total samples in analysis:", nrow(pcoa_all_data), "\n")
cat("Selected groups samples in plot:", nrow(pcoa_plot_data), "\n")


# ============================================================================
# SECTION: Plot all data groups (comprehensive view)
# ============================================================================

cat("\n=== Creating comprehensive PCoA plot with all groups ===\n")

# Get all unique groups from the metadata
all_groups <- unique(pcoa_all_data$Group)
cat("All groups found in data:", paste(all_groups, collapse = ", "), "\n")

# Remove any NA groups
all_groups <- all_groups[!is.na(all_groups)]

# Filter to include all groups
pcoa_all_groups_data <- pcoa_all_data %>%
  filter(!is.na(Group))

n_groups <- length(all_groups)
cat("Total number of groups:", n_groups, "\n")

# Function to assign colors based on first letter of group name
assign_color_by_letter <- function(group_names) {
  colors <- character(length(group_names))
  
  black_grays <- c("black", "#2F2F2F", "#404040", "#555555", "#696969", "#808080", "#A9A9A9", "#C0C0C0")
  oranges_yellows <- c("#FF8C00", "#FFA500", "#FFB347", "#FFCC33", "#FFD700", "#FFFF00", "#FFFF66", "#FFFF99")
  reds <- c("#8B0000", "#B22222", "#DC143C", "#FF0000", "#FF4500", "#FF6347", "#FF7F7F", "#FFA07A")
  blues <- c("#000080", "#0000CD", "#0000FF", "#4169E1", "#6495ED", "#87CEEB", "#ADD8E6", "#B0E0E6")
  
  b_groups <- group_names[startsWith(group_names, "B")]
  o_groups <- group_names[startsWith(group_names, "O")]
  a_groups <- group_names[startsWith(group_names, "A")]
  n_groups_vec <- group_names[startsWith(group_names, "N")]
  other_groups <- group_names[!startsWith(group_names, "B") & 
                                !startsWith(group_names, "O") & 
                                !startsWith(group_names, "A") & 
                                !startsWith(group_names, "N")]
  
  for (i in seq_along(group_names)) {
    group <- group_names[i]
    first_letter <- substr(group, 1, 1)
    
    if (first_letter == "B") {
      b_index <- which(b_groups == group)
      colors[i] <- black_grays[((b_index - 1) %% length(black_grays)) + 1]
    } else if (first_letter == "O") {
      o_index <- which(o_groups == group)
      colors[i] <- oranges_yellows[((o_index - 1) %% length(oranges_yellows)) + 1]
    } else if (first_letter == "A") {
      a_index <- which(a_groups == group)
      colors[i] <- reds[((a_index - 1) %% length(reds)) + 1]
    } else if (first_letter == "N") {
      n_index <- which(n_groups_vec == group)
      colors[i] <- blues[((n_index - 1) %% length(blues)) + 1]
    } else {
      other_colors <- c("#800080", "#9932CC", "#DA70D6", "#228B22", "#32CD32", "#90EE90")
      other_index <- which(other_groups == group)
      colors[i] <- other_colors[((other_index - 1) %% length(other_colors)) + 1]
    }
  }
  
  return(colors)
}

all_colors <- assign_color_by_letter(all_groups)
names(all_colors) <- all_groups

cat("\nColor assignments by group:\n")
for (i in seq_along(all_groups)) {
  cat(sprintf("  %s: %s\n", all_groups[i], all_colors[i]))
}

available_shapes <- c(16, 17, 15, 18, 19, 20, 8, 11, 12, 13, 14, 21, 22, 23, 24, 25)
all_shapes <- rep(available_shapes, length.out = n_groups)
names(all_shapes) <- all_groups

# Define explicit legend order for the comprehensive plot
legend_order <- c("BB.0", "BB.1", "BB.2",
                  "BO.0", "BO.1", "BO.2",
                  "AL.0", "AL.1", "AL.2",
                  "OB.0", "OB.1", "OB.2",
                  "OO.0", "OO.1", "OO.2")

# Keep only levels that are actually present in the data, preserving the desired order
ordered_levels <- legend_order[legend_order %in% all_groups]
# Append any groups present in the data but not listed above (failsafe)
ordered_levels <- c(ordered_levels, setdiff(all_groups, ordered_levels))

all_group_labels <- setNames(ordered_levels, ordered_levels)

pcoa_all_groups_data$Group <- factor(pcoa_all_groups_data$Group, levels = ordered_levels)

# Create convex hulls for all groups
hulls_all <- pcoa_all_groups_data %>%
  group_by(Group) %>%
  do(find_hull(.)) %>%
  ungroup()

# Calculate centroids for all groups
group_centroids_all <- pcoa_all_groups_data %>%
  group_by(Group) %>%
  summarise(
    centroid_x = mean(PCoA1),
    centroid_y = mean(PCoA2),
    .groups = 'drop'
  )

# Create the comprehensive PCoA plot
pcoa_all_plot <- ggplot() +
  # Add convex hulls for each group
  geom_polygon(data = hulls_all, 
               aes(x = PCoA1, y = PCoA2, 
                   fill = Group, 
                   color = Group, 
                   group = Group), 
               alpha = 0.1,
               linetype = 2) +
  # Add colored points for all groups
  geom_point(data = pcoa_all_groups_data, 
             aes(x = PCoA1, y = PCoA2, 
                 color = Group, 
                 shape = Group), 
             size = 3, 
             alpha = 0.8) +
  # Add open circles at group centroids
  geom_point(data = group_centroids_all,
             aes(x = centroid_x, y = centroid_y, 
                 color = Group),
             shape = 1,
             size = 5,
             stroke = 1.5) +
  scale_color_manual(values = all_colors, labels = all_group_labels, breaks = ordered_levels) +
  scale_shape_manual(values = all_shapes, labels = all_group_labels, breaks = ordered_levels) +
  scale_fill_manual(values = all_colors, labels = all_group_labels, breaks = ordered_levels) +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.8, "cm"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_line(color = "gray98")
  ) +
  # Add axis labels with percent variance explained
  labs(
    title = "PCoA Ordination - All Sample Groups",
    subtitle = "Based on ILR Aitchison distance of 16S rRNA amplicon data",
    x = axis1_label,
    y = axis2_label,
    color = "Sample Group",
    shape = "Sample Group",
    fill = "Sample Group"
  ) +
  guides(fill = "none")

# Adjust legend if there are many groups
if (n_groups > 8) {
  pcoa_all_plot <- pcoa_all_plot +
    theme(legend.position = "right") +
    guides(
      color = guide_legend(ncol = ceiling(n_groups/15), override.aes = list(size = 2)),
      shape = guide_legend(ncol = ceiling(n_groups/15), override.aes = list(size = 2))
    )
}

# Print the comprehensive PCoA plot
print(pcoa_all_plot)

# Save the comprehensive PCoA plot
ggsave(
  "FigS2_PCoA_all_experimental_groups.pdf", 
  pcoa_all_plot, 
  width = 9,
  height = 8
)


# ==========================================
# PERMANOVA ANALYSIS SECTION
# ==========================================

cat("\n=== PERMANOVA ANALYSIS (ABUNDANCE DATA) ===\n")

# Prepare data for PERMANOVA analysis
permanova_metadata <- pcoa_all_data[match(labels(aitchison_dist), pcoa_all_data$Sample), ]

if(!all(labels(aitchison_dist) == permanova_metadata$Sample)) {
  stop("Sample order mismatch between distance matrix and metadata!")
}

# 1. Overall PERMANOVA test for all groups
cat("\n1. Overall PERMANOVA test (all groups):\n")
cat("Testing: Group (Order.Timing combination)\n")

set.seed(1523)
permanova_all <- adonis2(aitchison_dist ~ Group, 
                         data = permanova_metadata, 
                         permutations = 999,
                         method = "euclidean")

print(permanova_all)


# Pairwise PERMANOVA tests for all groups above abundance thresholds
cat("\n Pairwise PERMANOVA tests for all groups above abundance thresholds:\n")

cat("\nPerforming manual pairwise PERMANOVA tests...\n")
metadata_plot <- metadata[metadata$Sample %in% colnames(abundance_replaced),]
group_combinations <- combn(unique(metadata_plot$Group), 2, simplify = FALSE)
pairwise_results <- data.frame(
  Group1 = character(),
  Group2 = character(),
  R2 = numeric(),
  F_statistic = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(i in seq_along(group_combinations)) {
  group_pair <- group_combinations[[i]]
  
  pair_metadata <- metadata_plot %>%
    filter(Group %in% group_pair)
  
  if(nrow(pair_metadata) < 4) {
    cat("Skipping", paste(group_pair, collapse = " vs "), 
        "- insufficient samples\n")
    next
  }
  
  pair_samples <- pair_metadata$Sample
  pair_dist <- as.dist(as.matrix(aitchison_dist)[pair_samples, pair_samples])
  
  set.seed(1523)
  pair_result <- adonis2(pair_dist ~ Group, 
                         data = pair_metadata, 
                         permutations = 999,
                         method = "euclidean")
  
  pairwise_results <- rbind(pairwise_results, data.frame(
    Group1 = group_pair[1],
    Group2 = group_pair[2],
    R2 = pair_result$R2[1],
    F_statistic = pair_result$F[1],
    p_value = pair_result$`Pr(>F)`[1]
  ))
}

# Apply false discovery rate correction
pairwise_results$p_adjusted <- p.adjust(pairwise_results$p_value, method = "fdr")
pairwise_results <- pairwise_results %>% 
  mutate(R2 = round(R2,2),
         F_statistic = round(F_statistic,2),
         p_adjusted = round(p_adjusted,4)) 
write_csv(pairwise_results, "pairwise_permanova_fdr.csv")


# =============================================================================
# RAREFACTION SECTION
# Rarefy count_data to even read depth, then repeat the full PCoA pipeline
# (ILR/Aitchison distance, both plots, PERMANOVA) and compare to unrarefied.
# =============================================================================

# -----------------------------------------------------------------------------
# R1. CHOOSE RAREFACTION DEPTH & RAREFY
# -----------------------------------------------------------------------------

# count_data is ESVs x samples; transpose to samples x ESVs for vegan
counts_t_pcoa <- as.data.frame(t(count_data))

sample_depths_pcoa <- rowSums(counts_t_pcoa)
cat("\n========== READ DEPTH SUMMARY (unrarefied, PCoA script) ==========\n")
print(summary(sample_depths_pcoa))
cat("\nPer-sample depths (sorted):\n")
print(sort(sample_depths_pcoa))

# Rarefy to minimum depth (most conservative / no samples lost).
# Change rare_depth_pcoa manually if the minimum is a clear outlier you want to drop.
rare_depth_pcoa <- 10000
cat("\nRarefying to:", rare_depth_pcoa, "reads per sample\n")

samples_kept_pcoa    <- names(sample_depths_pcoa[sample_depths_pcoa >= rare_depth_pcoa])
samples_dropped_pcoa <- names(sample_depths_pcoa[sample_depths_pcoa <  rare_depth_pcoa])
if (length(samples_dropped_pcoa) > 0) {
  cat("Samples dropped (below rarefaction depth):",
      paste(samples_dropped_pcoa, collapse = ", "), "\n")
} else {
  cat("All samples retained after rarefaction.\n")
}

set.seed(42)
counts_rare_pcoa <- as.data.frame(
  vegan::rrarefy(counts_t_pcoa[samples_kept_pcoa, ], sample = rare_depth_pcoa)
)
cat("Rarefied matrix:", nrow(counts_rare_pcoa), "samples x",
    ncol(counts_rare_pcoa), "ESVs\n")

# Transpose back to ESVs x samples to match the format expected by the rest
# of this pipeline (count_data orientation)
count_data_rare <- as.data.frame(t(counts_rare_pcoa))

# -----------------------------------------------------------------------------
# R2. FILTER LOW-ABUNDANCE ESVs & ILR / AITCHISON DISTANCE (RAREFIED)
# -----------------------------------------------------------------------------

# Apply the same ESV mean-abundance filter (>= 5 reads) used on unrarefied data
esv_means_rare      <- rowMeans(count_data_rare)
filtered_esvs_rare  <- count_data_rare[esv_means_rare >= 5, ]

cat("\nRarefied data dimensions after ESV filter:",
    dim(filtered_esvs_rare), "\n")
cat("Rarefied data range — Min:", min(filtered_esvs_rare),
    "Max:", max(filtered_esvs_rare), "\n")

# Zero replacement, ILR transformation, Aitchison distance
abundance_replaced_rare <- cmultRepl(filtered_esvs_rare, label = 0, method = "CZM")
abundance_t_rare        <- t(abundance_replaced_rare)            # samples x ESVs
comp_data_rare          <- acomp(abundance_t_rare)
ilr_data_rare           <- ilr(comp_data_rare)
aitchison_dist_rare     <- dist(ilr_data_rare, method = "euclidean")

# -----------------------------------------------------------------------------
# R3. PCoA ON RAREFIED AITCHISON DISTANCE
# -----------------------------------------------------------------------------

pcoa_result_rare <- cmdscale(aitchison_dist_rare, k = 2, eig = TRUE)

positive_eigs_rare <- pcoa_result_rare$eig[pcoa_result_rare$eig > 0]
total_variance_rare <- sum(positive_eigs_rare)
pct_var_rare        <- round(100 * positive_eigs_rare / total_variance_rare, 1)

pct_axis1_rare <- pct_var_rare[1]
pct_axis2_rare <- pct_var_rare[2]

cat("\nRarefied PCoA Axis 1 variance explained:", pct_axis1_rare, "%\n")
cat("Rarefied PCoA Axis 2 variance explained:", pct_axis2_rare, "%\n")

axis1_label_rare <- paste0("PCoA1 (", pct_axis1_rare, "%)")
axis2_label_rare <- paste0("PCoA2 (", pct_axis2_rare, "%)")

# Extract scores and join with metadata
pcoa_scores_rare <- as.data.frame(pcoa_result_rare$points)
colnames(pcoa_scores_rare) <- c("PCoA1", "PCoA2")
pcoa_scores_rare$SampleID  <- rownames(pcoa_scores_rare)

# Only retain metadata rows for samples that survived rarefaction
metadata_rare <- metadata %>%
  filter(Sample %in% samples_kept_pcoa)

pcoa_all_data_rare <- left_join(
  pcoa_scores_rare,
  metadata_rare %>% mutate(SampleID = Sample),
  by = "SampleID"
)

# Filter to selected groups
pcoa_plot_data_rare <- pcoa_all_data_rare %>%
  filter(Group %in% specified_groups) %>%
  mutate(Group = factor(Group, levels = specified_groups))

# -----------------------------------------------------------------------------
# R4. SELECTED-GROUPS PCoA PLOT (RAREFIED)
# Same visual style as the unrarefied selected-groups plot
# -----------------------------------------------------------------------------

rare_subtitle_pcoa <- paste0("Rarefied to ", rare_depth_pcoa, " reads/sample")

# Convex hulls
hulls_rare <- pcoa_plot_data_rare %>%
  group_by(Group) %>%
  do(find_hull(.)) %>%
  ungroup()

# Centroids
group_centroids_rare <- pcoa_plot_data_rare %>%
  group_by(Group) %>%
  summarise(centroid_x = mean(PCoA1), centroid_y = mean(PCoA2), .groups = "drop")

# Arrows between centroids (BB.0 -> BO.1 -> BO.2 and OO.0 -> OB.1 -> OB.2)
get_centroid <- function(df, grp) df[df$Group == grp, ]

bb0r  <- get_centroid(group_centroids_rare, "BB.0")
bo1r  <- get_centroid(group_centroids_rare, "BO.1")
bo2r  <- get_centroid(group_centroids_rare, "BO.2")
oo0r  <- get_centroid(group_centroids_rare, "OO.0")
ob1r  <- get_centroid(group_centroids_rare, "OB.1")
ob2r  <- get_centroid(group_centroids_rare, "OB.2")

arrow1r <- calculate_arrow_endpoints(bb0r$centroid_x, bb0r$centroid_y,
                                     bo1r$centroid_x, bo1r$centroid_y, circle_radius)
arrow2r <- calculate_arrow_endpoints(bo1r$centroid_x, bo1r$centroid_y,
                                     bo2r$centroid_x, bo2r$centroid_y, circle_radius)
arrow3r <- calculate_arrow_endpoints(oo0r$centroid_x, oo0r$centroid_y,
                                     ob1r$centroid_x, ob1r$centroid_y, circle_radius)
arrow4r <- calculate_arrow_endpoints(ob1r$centroid_x, ob1r$centroid_y,
                                     ob2r$centroid_x, ob2r$centroid_y, circle_radius)

black_arrows_rare <- data.frame(
  x    = c(arrow1r$start_x, arrow2r$start_x),
  y    = c(arrow1r$start_y, arrow2r$start_y),
  xend = c(arrow1r$end_x,   arrow2r$end_x),
  yend = c(arrow1r$end_y,   arrow2r$end_y),
  arrow_type = c("BB.0_to_BO.1", "BO.1_to_BO.2")
)

orange_arrows_rare <- data.frame(
  x    = c(arrow3r$start_x, arrow4r$start_x),
  y    = c(arrow3r$start_y, arrow4r$start_y),
  xend = c(arrow3r$end_x,   arrow4r$end_x),
  yend = c(arrow3r$end_y,   arrow4r$end_y),
  arrow_type = c("OO.0_to_OB.1", "OB.1_to_OB.2")
)

pcoa_plot_rare <- ggplot() +
  geom_point(data = pcoa_all_data_rare,
             aes(x = PCoA1, y = PCoA2),
             color = "gray90", size = 1.5, alpha = 0.2) +
  geom_polygon(data = hulls_rare,
               aes(x = PCoA1, y = PCoA2, fill = Group, color = Group, group = Group),
               alpha = 0.1, linetype = 2) +
  geom_segment(data = black_arrows_rare,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "black", size = 1.2, alpha = 0.8) +
  geom_segment(data = orange_arrows_rare,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               color = "#FF8C00", size = 1.2, alpha = 0.8) +
  geom_point(data = pcoa_plot_data_rare,
             aes(x = PCoA1, y = PCoA2, color = Group, shape = Group),
             size = 2, alpha = 0.9) +
  geom_point(data = group_centroids_rare,
             aes(x = centroid_x, y = centroid_y, color = Group),
             shape = 1, size = 4, stroke = 1.5) +
  scale_color_manual(values = custom_colors, labels = group_labels) +
  scale_shape_manual(values = custom_shapes, labels = group_labels) +
  scale_fill_manual(values  = custom_colors, labels = group_labels) +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title   = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_line(color = "gray98")
  ) +
  labs(
    title    = paste0("PCoA (Rarefied) — Selected Groups"),
    subtitle = rare_subtitle_pcoa,
    x = axis1_label_rare, y = axis2_label_rare,
    color = "Sample Group", shape = "Sample Group", fill = "Sample Group"
  ) +
  guides(fill = "none")

print(pcoa_plot_rare)
ggsave(
  "pcoa_rarefied_selected_groups.pdf",
  pcoa_plot_rare, width = 8, height = 5
)

# -----------------------------------------------------------------------------
# R5. ALL-GROUPS PCoA PLOT (RAREFIED)
# -----------------------------------------------------------------------------

cat("\n=== Creating rarefied comprehensive PCoA plot (all groups) ===\n")

all_groups_rare <- unique(pcoa_all_data_rare$Group)
all_groups_rare <- all_groups_rare[!is.na(all_groups_rare)]
cat("All groups in rarefied data:", paste(all_groups_rare, collapse = ", "), "\n")

pcoa_all_groups_data_rare <- pcoa_all_data_rare %>% filter(!is.na(Group))

n_groups_rare <- length(all_groups_rare)

# Reuse the same color-assignment function from the unrarefied section
all_colors_rare <- assign_color_by_letter(all_groups_rare)
names(all_colors_rare) <- all_groups_rare

all_shapes_rare <- rep(available_shapes, length.out = n_groups_rare)
names(all_shapes_rare) <- all_groups_rare

# Preserve the same legend order as the unrarefied all-groups plot
ordered_levels_rare <- legend_order[legend_order %in% all_groups_rare]
ordered_levels_rare <- c(ordered_levels_rare,
                         setdiff(all_groups_rare, ordered_levels_rare))
all_group_labels_rare <- setNames(ordered_levels_rare, ordered_levels_rare)

pcoa_all_groups_data_rare$Group <- factor(pcoa_all_groups_data_rare$Group,
                                          levels = ordered_levels_rare)

hulls_all_rare <- pcoa_all_groups_data_rare %>%
  group_by(Group) %>%
  do(find_hull(.)) %>%
  ungroup()

group_centroids_all_rare <- pcoa_all_groups_data_rare %>%
  group_by(Group) %>%
  summarise(centroid_x = mean(PCoA1), centroid_y = mean(PCoA2), .groups = "drop")

pcoa_all_plot_rare <- ggplot() +
  geom_polygon(data = hulls_all_rare,
               aes(x = PCoA1, y = PCoA2, fill = Group, color = Group, group = Group),
               alpha = 0.1, linetype = 2) +
  geom_point(data = pcoa_all_groups_data_rare,
             aes(x = PCoA1, y = PCoA2, color = Group, shape = Group),
             size = 3, alpha = 0.8) +
  geom_point(data = group_centroids_all_rare,
             aes(x = centroid_x, y = centroid_y, color = Group),
             shape = 1, size = 5, stroke = 1.5) +
  scale_color_manual(values = all_colors_rare, labels = all_group_labels_rare,
                     breaks = ordered_levels_rare) +
  scale_shape_manual(values = all_shapes_rare, labels = all_group_labels_rare,
                     breaks = ordered_levels_rare) +
  scale_fill_manual(values  = all_colors_rare, labels = all_group_labels_rare,
                    breaks = ordered_levels_rare) +
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(size = 8),
    legend.key.size  = unit(0.8, "cm"),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.title   = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_line(color = "gray98")
  ) +
  labs(
    title    = "PCoA Ordination (Rarefied) — All Sample Groups",
    subtitle = paste0("ILR Aitchison distance | ", rare_subtitle_pcoa),
    x = axis1_label_rare, y = axis2_label_rare,
    color = "Sample Group", shape = "Sample Group", fill = "Sample Group"
  ) +
  guides(fill = "none")

if (n_groups_rare > 8) {
  pcoa_all_plot_rare <- pcoa_all_plot_rare +
    theme(legend.position = "right") +
    guides(
      color = guide_legend(ncol = ceiling(n_groups_rare / 15),
                           override.aes = list(size = 2)),
      shape = guide_legend(ncol = ceiling(n_groups_rare / 15),
                           override.aes = list(size = 2))
    )
}

print(pcoa_all_plot_rare)
ggsave(
  "pcoa_rarefied_all_experimental_groups.pdf",
  pcoa_all_plot_rare, width = 9, height = 8
)

# -----------------------------------------------------------------------------
# R6. PERMANOVA ON RAREFIED DATA
# -----------------------------------------------------------------------------

cat("\n=== PERMANOVA ANALYSIS (RAREFIED DATA) ===\n")

permanova_metadata_rare <- pcoa_all_data_rare[
  match(labels(aitchison_dist_rare), pcoa_all_data_rare$Sample), ]

if (!all(labels(aitchison_dist_rare) == permanova_metadata_rare$Sample)) {
  stop("Sample order mismatch between rarefied distance matrix and metadata!")
}

# Overall PERMANOVA
cat("\n1. Overall PERMANOVA test (rarefied, all groups):\n")
set.seed(1523)
permanova_all_rare <- adonis2(
  aitchison_dist_rare ~ Group,
  data         = permanova_metadata_rare,
  permutations = 999,
  method       = "euclidean"
)
print(permanova_all_rare)

# Pairwise PERMANOVA
cat("\nPerforming pairwise PERMANOVA tests (rarefied)...\n")
metadata_plot_rare <- metadata_rare[metadata_rare$Sample %in% colnames(abundance_replaced_rare), ]
group_combinations_rare <- combn(unique(metadata_plot_rare$Group), 2, simplify = FALSE)

pairwise_results_rare <- data.frame(
  Group1 = character(), Group2 = character(),
  R2 = numeric(), F_statistic = numeric(), p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_along(group_combinations_rare)) {
  group_pair <- group_combinations_rare[[i]]
  pair_metadata <- metadata_plot_rare %>% filter(Group %in% group_pair)

  if (nrow(pair_metadata) < 4) {
    cat("Skipping", paste(group_pair, collapse = " vs "),
        "- insufficient samples\n")
    next
  }

  pair_samples  <- pair_metadata$Sample
  pair_dist_r   <- as.dist(as.matrix(aitchison_dist_rare)[pair_samples, pair_samples])

  set.seed(1523)
  pair_result_r <- adonis2(pair_dist_r ~ Group,
                            data         = pair_metadata,
                            permutations = 999,
                            method       = "euclidean")

  pairwise_results_rare <- rbind(pairwise_results_rare, data.frame(
    Group1      = group_pair[1],
    Group2      = group_pair[2],
    R2          = pair_result_r$R2[1],
    F_statistic = pair_result_r$F[1],
    p_value     = pair_result_r$`Pr(>F)`[1]
  ))
}

pairwise_results_rare$p_adjusted <- p.adjust(pairwise_results_rare$p_value, method = "fdr")
pairwise_results_rare <- pairwise_results_rare %>%
  mutate(
    R2          = round(R2, 2),
    F_statistic = round(F_statistic, 2),
    p_adjusted  = round(p_adjusted, 4)
  )

write_csv(pairwise_results_rare, "pairwise_permanova_rarefied_fdr.csv")
cat("Wrote pairwise_permanova_rarefied_fdr.csv\n")

# -----------------------------------------------------------------------------
# R7. COMPARISON SUMMARY: UNRAREFIED vs RAREFIED
# -----------------------------------------------------------------------------

cat("\n")
cat("=============================================================\n")
cat("  COMPARISON SUMMARY: UNRAREFIED vs RAREFIED PCoA / PERMANOVA\n")
cat("=============================================================\n")

# --- Variance explained by PCoA axes ---
cat("\n--- PCoA variance explained ---\n")
pcoa_var_comparison <- data.frame(
  Dataset = c("Unrarefied", "Rarefied"),
  Axis1_pct = c(pct_axis1,      pct_axis1_rare),
  Axis2_pct = c(pct_axis2,      pct_axis2_rare),
  Total_2ax  = c(pct_axis1      + pct_axis2,
                 pct_axis1_rare + pct_axis2_rare)
)
print(pcoa_var_comparison, row.names = FALSE)

# --- Overall PERMANOVA comparison ---
cat("\n--- Overall PERMANOVA results ---\n")
overall_perm_comparison <- data.frame(
  Dataset     = c("Unrarefied", "Rarefied"),
  R2          = c(round(permanova_all$R2[1],      4),
                  round(permanova_all_rare$R2[1], 4)),
  F_statistic = c(round(permanova_all$F[1],       2),
                  round(permanova_all_rare$F[1],  2)),
  p_value     = c(permanova_all$`Pr(>F)`[1],
                  permanova_all_rare$`Pr(>F)`[1])
)
print(overall_perm_comparison, row.names = FALSE)

# --- Pairwise PERMANOVA: classify each pair by significance change ---
cat("\n--- Pairwise PERMANOVA significance comparison ---\n")

# Create a common key for merging
pw_unrar <- pairwise_results %>%
  mutate(.key = paste(pmin(Group1, Group2), pmax(Group1, Group2), sep = " | ")) %>%
  dplyr::select(.key, R2_unrar = R2, F_unrar = F_statistic,
         p_unrar = p_value, padj_unrar = p_adjusted)

pw_rar <- pairwise_results_rare %>%
  mutate(.key = paste(pmin(Group1, Group2), pmax(Group1, Group2), sep = " | ")) %>%
  dplyr::select(.key, R2_rar = R2, F_rar = F_statistic,
         p_rar = p_value, padj_rar = p_adjusted)

pw_comparison <- full_join(pw_unrar, pw_rar, by = ".key") %>%
  rename(Pair = .key) %>%
  mutate(
    sig_unrar = !is.na(padj_unrar) & padj_unrar < 0.05,
    sig_rar   = !is.na(padj_rar)   & padj_rar   < 0.05,
    Status    = dplyr::case_when(
      sig_unrar & sig_rar   ~ "Significant in both",
      sig_unrar & !sig_rar  ~ "Lost significance after rarefaction",
      !sig_unrar & sig_rar  ~ "Gained significance after rarefaction",
      TRUE                  ~ "Non-significant in both"
    )
  )

# Print a concise summary table
cat("\nStatus counts:\n")
print(table(pw_comparison$Status))

cat("\nPairs that LOST significance after rarefaction:\n")
lost  <- pw_comparison %>% filter(Status == "Lost significance after rarefaction")
if (nrow(lost)  > 0) print(lost  %>% select(Pair, padj_unrar, padj_rar, Status),
                            row.names = FALSE) else cat("  (none)\n")

cat("\nPairs that GAINED significance after rarefaction:\n")
gained <- pw_comparison %>% filter(Status == "Gained significance after rarefaction")
if (nrow(gained) > 0) print(gained %>% select(Pair, padj_unrar, padj_rar, Status),
                             row.names = FALSE) else cat("  (none)\n")

cat("\nPairs significant in BOTH datasets:\n")
both  <- pw_comparison %>% filter(Status == "Significant in both")
if (nrow(both)  > 0) print(both  %>% dplyr::select(Pair, padj_unrar, padj_rar, R2_unrar, R2_rar),
                            row.names = FALSE) else cat("  (none)\n")

# --- R2 shift for pairs significant in both ---
if (nrow(both) > 0) {
  cat("\nR2 change (rarefied - unrarefied) for pairs significant in both:\n")
  r2_shift <- both %>%
    mutate(R2_change = round(R2_rar - R2_unrar, 3)) %>%
    dplyr::select(Pair, R2_unrar, R2_rar, R2_change) %>%
    arrange(desc(abs(R2_change)))
  print(r2_shift, row.names = FALSE)
}

# --- Samples retained ---
cat("\n--- Sample retention ---\n")
cat("Unrarefied: ", nrow(pcoa_all_data),      "samples in PCoA\n")
cat("Rarefied:   ", nrow(pcoa_all_data_rare), "samples in PCoA\n")
if (length(samples_dropped_pcoa) > 0) {
  cat("Dropped samples:", paste(samples_dropped_pcoa, collapse = ", "), "\n")
}

# --- Write comparison CSVs ---
write_csv(pcoa_var_comparison,
          "comparison_pcoa_variance_unrarefied_vs_rarefied.csv")
cat("\nWrote comparison_pcoa_variance_unrarefied_vs_rarefied.csv\n")

write_csv(overall_perm_comparison,
          "comparison_overall_permanova_unrarefied_vs_rarefied.csv")
cat("Wrote comparison_overall_permanova_unrarefied_vs_rarefied.csv\n")

write_csv(pw_comparison,
          "comparison_pairwise_permanova_unrarefied_vs_rarefied.csv")
cat("Wrote comparison_pairwise_permanova_unrarefied_vs_rarefied.csv\n")

cat("\n=============================================================\n")
cat("  END OF RAREFACTION SECTION\n")
cat("=============================================================\n")
