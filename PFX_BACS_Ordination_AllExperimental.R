# Script to perform NMDS ordination using ILR Aitchison distance 
# on 16S SSU rRNA amplicon sequencing data
# Analyzing all samples but plotting only selected groups
# Using black/orange color scheme with different shapes/shades

# Script written with the assistance of Claude Sonnet 4 LLM by Anthropic.

setwd("/home/pacifica/R/antarctica/S03/big_experiment/")

# Load required libraries
library(tidyverse)     # For data manipulation and visualization
library(vegan)         # For NMDS and PERMANOVA
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

# Step 2: MODIFIED - Work directly with abundance data (no percent conversion)
# The filtered_esvs already contains the abundance data we want to use
abundance_data <- filtered_esvs

cat("Using abundance data for Aitchison distance calculation\n")
cat("Data dimensions:", dim(abundance_data), "\n")
cat("Data range - Min:", min(abundance_data), "Max:", max(abundance_data), "\n")

# Step 3: Calculate ILR Aitchison distance for all samples using abundance data
# First, we need to handle zeros for the ILR transformation

abundance_replaced <- cmultRepl(abundance_data,  label=0, method="CZM")
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

# Step 4: Run NMDS on Aitchison distance for all samples
set.seed(1523)  # Set seed for reproducibility
nmds_result <- metaMDS(aitchison_dist, k = 2, trymax = 100)

# Extract NMDS scores for all samples
nmds_scores <- as.data.frame(scores(nmds_result))
nmds_scores$SampleID <- rownames(nmds_scores)

# Join with metadata
nmds_all_data <- left_join(nmds_scores, metadata %>% mutate(SampleID=Sample), by = "SampleID")

# Filter to only include the specified groups for plotting
nmds_plot_data <- nmds_all_data %>%
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
nmds_plot_data$Group <- factor(nmds_plot_data$Group, levels = specified_groups)

# Define custom color and shape scheme
# Create a custom color vector (black shades for BB groups, orange shades for OO groups)
custom_colors <- c(
  "BB.0" = "black",
  "BO.1" = alpha("black", 0.7),  # Black with transparency
  "BO.2" = alpha("black", 0.4),  # Lighter black with more transparency
  "OO.0" = "#FF8C00",            # Dark orange 
  "OB.1" = "#FFA500",            # Medium orange
  "OB.2" = "#FFD700"             # Light orange/gold
)

# Custom shapes for each group
custom_shapes <- c(
  "BB.0" = 16,  # Filled circle
  "BO.1" = 17,  # Filled triangle
  "BO.2" = 15,  # Filled square
  "OO.0" = 16,  # Filled circle
  "OB.1" = 17,  # Filled triangle
  "OB.2" = 15   # Filled square
)

# Function to create convex hull data frames for each group
find_hull <- function(df) {
  hull_points <- chull(df$NMDS1, df$NMDS2)
  # Add the first point at the end to close the polygon
  hull_points <- c(hull_points, hull_points[1])
  return(df[hull_points, ])
}

# Create a list to store convex hull data for each group
hulls <- nmds_plot_data %>%
  group_by(Group) %>%
  do(find_hull(.)) %>%
  ungroup()

# Calculate group centroids for arrow placement AND centroid markers
group_centroids <- nmds_plot_data %>%
  group_by(Group) %>%
  summarise(
    centroid_x = mean(NMDS1),
    centroid_y = mean(NMDS2),
    .groups = 'drop'
  )

# Function to calculate arrow endpoints on circle circumference
calculate_arrow_endpoints <- function(x1, y1, x2, y2, radius) {
  # Calculate the angle from point 1 to point 2
  angle <- atan2(y2 - y1, x2 - x1)
  
  # Calculate start point (on circumference of first circle, pointing toward second)
  start_x <- x1 + radius * cos(angle)
  start_y <- y1 + radius * sin(angle)
  
  # Calculate end point (on circumference of second circle, coming from first)
  end_x <- x2 - radius * cos(angle)
  end_y <- y2 - radius * sin(angle)
  
  return(list(start_x = start_x, start_y = start_y, end_x = end_x, end_y = end_y))
}

# Define circle radius (should match the visual size of centroid circles)
# Size 4 with stroke 1.5 in ggplot roughly corresponds to this radius in data units
circle_radius <- 1.8  # Adjust this value based on your plot scale

# Create arrow data frames with adjusted endpoints
# Black arrows (BB.0 -> BO.1 -> BO.2)
bb0_centroid <- group_centroids[group_centroids$Group == "BB.0", ]
bo1_centroid <- group_centroids[group_centroids$Group == "BO.1", ]
bo2_centroid <- group_centroids[group_centroids$Group == "BO.2", ]

# Calculate first black arrow (BB.0 -> BO.1)
arrow1 <- calculate_arrow_endpoints(bb0_centroid$centroid_x, bb0_centroid$centroid_y,
                                    bo1_centroid$centroid_x, bo1_centroid$centroid_y,
                                    circle_radius)

# Calculate second black arrow (BO.1 -> BO.2)
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

# Calculate first orange arrow (OO.0 -> OB.1)
arrow3 <- calculate_arrow_endpoints(oo0_centroid$centroid_x, oo0_centroid$centroid_y,
                                    ob1_centroid$centroid_x, ob1_centroid$centroid_y,
                                    circle_radius)

# Calculate second orange arrow (OB.1 -> OB.2)
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

# Step 5: Create the NMDS plot with custom colors and shapes
nmds_plot <- ggplot() +
  # Add gray points for all samples (background)
  geom_point(data = nmds_all_data, 
             aes(x = NMDS1, y = NMDS2), 
             color = "gray90", 
             size = 1.5, 
             alpha = 0.2) +
  # Add convex hulls for each group
  geom_polygon(data = hulls, 
               aes(x = NMDS1, y = NMDS2, 
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
  geom_point(data = nmds_plot_data, 
             aes(x = NMDS1, y = NMDS2, 
                 color = Group, 
                 shape = Group), 
             size = 2, 
             alpha = 0.9) +
  # Add open circles at group centroids
  geom_point(data = group_centroids,
             aes(x = centroid_x, y = centroid_y, 
                 color = Group),
             shape = 1,  # Open circle
             size = 4,   # Larger than regular points
             stroke = 1.5) +  # Thicker border
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
  # Add labels
  labs(
    title = "",
    subtitle = "",
    x = "NMDS1",
    y = "NMDS2",
    color = "Sample Group",
    shape = "Sample Group",
    fill = "Sample Group"
  ) +
  # # Limit x-axis to the range of plotted points with some padding
  # xlim(min(nmds_plot_data$NMDS1) - 0.1 * diff(range(nmds_plot_data$NMDS1)),
  #      max(nmds_plot_data$NMDS1) + 0.1 * diff(range(nmds_plot_data$NMDS1))) +
  # Add stress value to the plot (using plot data range for positioning)
  annotate(
    "text",
    x = min(nmds_plot_data$NMDS1),
    y = max(nmds_plot_data$NMDS2),
    label = paste("Stress =", round(nmds_result$stress, 3)),
    hjust = -4,
    vjust = 40
  ) +
  # Hide the fill from the legend since it's redundant with color
  guides(fill = "none")

# Print the NMDS plot
print(nmds_plot)

# Save the NMDS plot
ggsave(
  "nmds_ilr_aitchison_abundance_black_orange_convexhull_arrows.pdf", 
  nmds_plot, 
  width = 8, 
  height = 5
)

# Print information about the NMDS result
cat("NMDS Stress:", nmds_result$stress, "\n")
cat("NMDS Convergence:", nmds_result$converged, "\n")
cat("Number of iterations:", nmds_result$iterations, "\n")
cat("Total samples in analysis:", nrow(nmds_all_data), "\n")
cat("Selected groups samples in plot:", nrow(nmds_plot_data), "\n")


# ============================================================================
# SECTION: Plot all data groups (comprehensive view)
# ============================================================================

cat("\n=== Creating comprehensive NMDS plot with all groups ===\n")

# Get all unique groups from the metadata
all_groups <- unique(nmds_all_data$Group)
cat("All groups found in data:", paste(all_groups, collapse = ", "), "\n")

# Remove any NA groups
all_groups <- all_groups[!is.na(all_groups)]

# Filter to include all groups (removing any samples without group assignment)
nmds_all_groups_data <- nmds_all_data %>%
  filter(!is.na(Group))

# Create a comprehensive color palette for all groups based on first letter
n_groups <- length(all_groups)
cat("Total number of groups:", n_groups, "\n")

# Function to assign colors based on first letter of group name
assign_color_by_letter <- function(group_names) {
  colors <- character(length(group_names))
  
  # Define color palettes for each letter
  black_grays <- c("black", "#2F2F2F", "#404040", "#555555", "#696969", "#808080", "#A9A9A9", "#C0C0C0")
  oranges_yellows <- c("#FF8C00", "#FFA500", "#FFB347", "#FFCC33", "#FFD700", "#FFFF00", "#FFFF66", "#FFFF99")
  reds <- c("#8B0000", "#B22222", "#DC143C", "#FF0000", "#FF4500", "#FF6347", "#FF7F7F", "#FFA07A")
  blues <- c("#000080", "#0000CD", "#0000FF", "#4169E1", "#6495ED", "#87CEEB", "#ADD8E6", "#B0E0E6")
  
  # Count groups by first letter to manage palette distribution
  b_groups <- group_names[startsWith(group_names, "B")]
  o_groups <- group_names[startsWith(group_names, "O")]
  a_groups <- group_names[startsWith(group_names, "A")]
  n_groups <- group_names[startsWith(group_names, "N")]
  other_groups <- group_names[!startsWith(group_names, "B") & 
                                !startsWith(group_names, "O") & 
                                !startsWith(group_names, "A") & 
                                !startsWith(group_names, "N")]
  
  # Assign colors
  for (i in seq_along(group_names)) {
    group <- group_names[i]
    first_letter <- substr(group, 1, 1)
    
    if (first_letter == "B") {
      # Black/gray shades
      b_index <- which(b_groups == group)
      colors[i] <- black_grays[((b_index - 1) %% length(black_grays)) + 1]
    } else if (first_letter == "O") {
      # Orange/yellow shades
      o_index <- which(o_groups == group)
      colors[i] <- oranges_yellows[((o_index - 1) %% length(oranges_yellows)) + 1]
    } else if (first_letter == "A") {
      # Red shades
      a_index <- which(a_groups == group)
      colors[i] <- reds[((a_index - 1) %% length(reds)) + 1]
    } else if (first_letter == "N") {
      # Blue shades
      n_index <- which(n_groups == group)
      colors[i] <- blues[((n_index - 1) %% length(blues)) + 1]
    } else {
      # Other groups - use purple/green shades
      other_colors <- c("#800080", "#9932CC", "#DA70D6", "#228B22", "#32CD32", "#90EE90")
      other_index <- which(other_groups == group)
      colors[i] <- other_colors[((other_index - 1) %% length(other_colors)) + 1]
    }
  }
  
  return(colors)
}

# Apply the color assignment function
all_colors <- assign_color_by_letter(all_groups)

# Create named vector for colors
names(all_colors) <- all_groups

# Print color assignments for verification
cat("\nColor assignments by group:\n")
for (i in seq_along(all_groups)) {
  cat(sprintf("  %s: %s\n", all_groups[i], all_colors[i]))
}

# Create shape palette - cycle through available shapes
available_shapes <- c(16, 17, 15, 18, 19, 20, 8, 11, 12, 13, 14, 21, 22, 23, 24, 25)
all_shapes <- rep(available_shapes, length.out = n_groups)
names(all_shapes) <- all_groups

# Create labels for all groups (use the group name as label)
all_group_labels <- setNames(all_groups, all_groups)

# Convert Group to a factor for consistent plotting
nmds_all_groups_data$Group <- factor(nmds_all_groups_data$Group, levels = all_groups)

# Create convex hulls for all groups
hulls_all <- nmds_all_groups_data %>%
  group_by(Group) %>%
  do(find_hull(.)) %>%
  ungroup()

# Calculate centroids for all groups
group_centroids_all <- nmds_all_groups_data %>%
  group_by(Group) %>%
  summarise(
    centroid_x = mean(NMDS1),
    centroid_y = mean(NMDS2),
    .groups = 'drop'
  )

# Create the comprehensive NMDS plot
nmds_all_plot <- ggplot() +
  # Add convex hulls for each group
  geom_polygon(data = hulls_all, 
               aes(x = NMDS1, y = NMDS2, 
                   fill = Group, 
                   color = Group, 
                   group = Group), 
               alpha = 0.1,
               linetype = 2) +
  # Add colored points for all groups
  geom_point(data = nmds_all_groups_data, 
             aes(x = NMDS1, y = NMDS2, 
                 color = Group, 
                 shape = Group), 
             size = 3, 
             alpha = 0.8) +
  # Add open circles at group centroids for all groups
  geom_point(data = group_centroids_all,
             aes(x = centroid_x, y = centroid_y, 
                 color = Group),
             shape = 1,  # Open circle
             size = 5,   # Larger than regular points
             stroke = 1.5) +  # Thicker border
  # Use the comprehensive color and shape palettes
  scale_color_manual(values = all_colors, labels = all_group_labels) +
  scale_shape_manual(values = all_shapes, labels = all_group_labels) +
  scale_fill_manual(values = all_colors, labels = all_group_labels) +
  # Apply minimal theme
  theme_minimal() +
  theme(
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 8),  # Smaller text for many groups
    legend.key.size = unit(0.8, "cm"),     # Smaller legend keys
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_line(color = "gray95"),
    panel.grid.minor = element_line(color = "gray98")
  ) +
  # Add labels
  labs(
    title = "NMDS Ordination - All Sample Groups",
    subtitle = "Based on ILR Aitchison distance of 16S rRNA amplicon data",
    x = "NMDS1",
    y = "NMDS2",
    color = "Sample Group",
    shape = "Sample Group",
    fill = "Sample Group"
  ) +
  # Add stress value to the plot
  annotate(
    "text",
    x = min(nmds_all_groups_data$NMDS1),
    y = max(nmds_all_groups_data$NMDS2),
    label = paste("Stress =", round(nmds_result$stress, 3)),
    hjust = 0,
    vjust = 1,
    fontface = "bold"
  ) +
  # Hide the fill from the legend since it's redundant with color
  guides(fill = "none")

# Adjust legend if there are many groups
if (n_groups > 8) {
  nmds_all_plot <- nmds_all_plot +
    theme(legend.position = "right") +
    guides(
      color = guide_legend(ncol = ceiling(n_groups/15), override.aes = list(size = 2)),
      shape = guide_legend(ncol = ceiling(n_groups/15), override.aes = list(size = 2))
    )
}

# Print the comprehensive NMDS plot
print(nmds_all_plot)

# Save the comprehensive NMDS plot
ggsave(
  "FigS2_Ordination_all_experimental_groups.pdf", 
  nmds_all_plot, 
  width = 9,  # Wider to accommodate legend
  height = 8
)



# ==========================================
# PERMANOVA ANALYSIS SECTION
# ==========================================

cat("\n=== PERMANOVA ANALYSIS (ABUNDANCE DATA) ===\n")

# Prepare data for PERMANOVA analysis
# Use the same Aitchison distance matrix calculated above
# Create metadata dataframe that matches the distance matrix sample order
permanova_metadata <- nmds_all_data[match(labels(aitchison_dist), nmds_all_data$Sample), ]

# Verify that sample names match between distance matrix and metadata
if(!all(labels(aitchison_dist) == permanova_metadata$Sample)) {
  stop("Sample order mismatch between distance matrix and metadata!")
}

# 1. Overall PERMANOVA test for all groups
cat("\n1. Overall PERMANOVA test (all groups):\n")
cat("Testing: Group (Order.Timing combination)\n")

set.seed(1523)  # Set seed for reproducibility
permanova_all <- adonis2(aitchison_dist ~ Group, 
                         data = permanova_metadata, 
                         permutations = 999,
                         method = "euclidean")

print(permanova_all)


# Pairwise PERMANOVA tests for all groups above abundance thresholds
cat("\n Pairwise PERMANOVA tests for all groups above abundance thresholds:\n")

# Manual pairwise tests using vegan::adonis2
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
  
  # Filter data for this pair
  pair_metadata <- metadata_plot %>%
    filter(Group %in% group_pair)
  
  if(nrow(pair_metadata) < 4) {
    cat("Skipping", paste(group_pair, collapse = " vs "), 
        "- insufficient samples\n")
    next
  }
  
  # Filter distance matrix for this pair
  pair_samples <- pair_metadata$Sample
  pair_dist <- as.dist(as.matrix(aitchison_dist)[pair_samples, pair_samples])
  
  # Run PERMANOVA for this pair
  set.seed(1523)
  pair_result <- adonis2(pair_dist ~ Group, 
                         data = pair_metadata, 
                         permutations = 999,
                         method = "euclidean")
  
  # Store results
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
write_csv(pairwise_results,"pairwise_permanova_fdr.csv")
