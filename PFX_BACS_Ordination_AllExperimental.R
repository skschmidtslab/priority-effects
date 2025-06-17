# Script to perform NMDS ordination using ILR Aitchison distance 
# on 16S SSU rRNA amplicon sequencing data
# Analyzing all samples but plotting only selected groups
# Using black/orange color scheme with different shapes/shades

setwd("/home/pacifica/R/antarctica/S03/big_experiment/")

# Load required libraries
library(tidyverse)     # For data manipulation and visualization
library(vegan)         # For NMDS and PERMANOVA
library(compositions)  # For ILR transformation and Aitchison distance
library(ggplot2)       # For plotting
library(scales)        # For color manipulation

# Define input file paths
count_table_path <- "ESV_table_PFX_Science_2023.csv"
metadata_path <- "Metadata_mapping_table_PFX_Science_2023.csv"

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
# Using a simple multiplicative replacement strategy for zeros
# zero_replacement <- function(x, delta = 0.01) {
#   # Find zeros
#   zeros <- x == 0
#   # Count zeros
#   n_zeros <- sum(zeros)
#   if (n_zeros == 0) return(x)
#   
#   # Calculate the multiplicative replacement value
#   delta_star <- delta * min(x[x > 0]) / n_zeros
#   
#   # Replace zeros
#   x[zeros] <- delta_star
#   
#   # Return the modified vector
#   return(x)
# }

# Apply zero replacement to each sample (abundance data)
#abundance_replaced <- apply(abundance_data, 2, zero_replacement)
abundance_replaced <- cmultRepl(abundance_data,  label=0, method="CZM")
# REMEMBER: IF samples are columns, margin = 2, if samples are ROWS, margin = 1
#field.bacs.nonat.ilr <- ilr(field.bacs.nonat.czm)

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
nmds_scores$Sample <- rownames(nmds_scores)

# Join with metadata
nmds_all_data <- left_join(nmds_scores, metadata, by = "Sample")

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

# Calculate group centroids for arrow placement
group_centroids <- nmds_plot_data %>%
  group_by(Group) %>%
  summarise(
    centroid_x = mean(NMDS1),
    centroid_y = mean(NMDS2),
    .groups = 'drop'
  )

# Create arrow data frames
# Black arrows (BB.0 -> BO.1 -> BO.2)
black_arrows <- data.frame(
  x = c(group_centroids$centroid_x[group_centroids$Group == "BB.0"],
        group_centroids$centroid_x[group_centroids$Group == "BO.1"]),
  y = c(group_centroids$centroid_y[group_centroids$Group == "BB.0"],
        group_centroids$centroid_y[group_centroids$Group == "BO.1"]),
  xend = c(group_centroids$centroid_x[group_centroids$Group == "BO.1"],
           group_centroids$centroid_x[group_centroids$Group == "BO.2"]),
  yend = c(group_centroids$centroid_y[group_centroids$Group == "BO.1"],
           group_centroids$centroid_y[group_centroids$Group == "BO.2"]),
  arrow_type = c("BB.0_to_BO.1", "BO.1_to_BO.2")
)

# Orange arrows (OO.0 -> OB.1 -> OB.2)
orange_arrows <- data.frame(
  x = c(group_centroids$centroid_x[group_centroids$Group == "OO.0"],
        group_centroids$centroid_x[group_centroids$Group == "OB.1"]),
  y = c(group_centroids$centroid_y[group_centroids$Group == "OO.0"],
        group_centroids$centroid_y[group_centroids$Group == "OB.1"]),
  xend = c(group_centroids$centroid_x[group_centroids$Group == "OB.1"],
           group_centroids$centroid_x[group_centroids$Group == "OB.2"]),
  yend = c(group_centroids$centroid_y[group_centroids$Group == "OB.1"],
           group_centroids$centroid_y[group_centroids$Group == "OB.2"]),
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
             size = 4, 
             alpha = 0.9) +
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
  # Limit x-axis to the range of plotted points with some padding
  xlim(min(nmds_plot_data$NMDS1) - 0.1 * diff(range(nmds_plot_data$NMDS1)),
       max(nmds_plot_data$NMDS1) + 0.1 * diff(range(nmds_plot_data$NMDS1))) +
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
  "nmds_ilr_aitchison_abundance_black_orange_convexhull.pdf", 
  nmds_plot, 
  width = 8, 
  height = 6
)

# Print information about the NMDS result
cat("NMDS Stress:", nmds_result$stress, "\n")
cat("NMDS Convergence:", nmds_result$converged, "\n")
cat("Number of iterations:", nmds_result$iterations, "\n")
cat("Total samples in analysis:", nrow(nmds_all_data), "\n")
cat("Selected groups samples in plot:", nrow(nmds_plot_data), "\n")

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

# 2. PERMANOVA test for individual factors
cat("\n2. PERMANOVA tests for individual factors:\n")

# Test Order effect
cat("\nTesting: Order (Black/Orange initial mat type)\n")
set.seed(1523)
permanova_order <- adonis2(aitchison_dist ~ Order, 
                           data = permanova_metadata, 
                           permutations = 999,
                           method = "euclidean")
print(permanova_order)

# Test Timing effect
cat("\nTesting: Timing (0, 1, 2 seasons)\n")
set.seed(1523)
permanova_timing <- adonis2(aitchison_dist ~ Timing, 
                            data = permanova_metadata, 
                            permutations = 999,
                            method = "euclidean")
print(permanova_timing)

# Test interaction effect
cat("\nTesting: Order * Timing interaction\n")
set.seed(1523)
permanova_interaction <- adonis2(aitchison_dist ~ Order * Timing, 
                                 data = permanova_metadata, 
                                 permutations = 999,
                                 method = "euclidean")
print(permanova_interaction)

# 3. PERMANOVA test for only the plotted groups (specified_groups)
cat("\n3. PERMANOVA test for plotted groups only:\n")

# Filter metadata to only include specified groups
plot_metadata <- permanova_metadata %>%
  filter(Group %in% specified_groups)

# Filter distance matrix to only include specified groups
plot_samples <- plot_metadata$Sample
aitchison_dist_plot <- as.dist(as.matrix(aitchison_dist)[plot_samples, plot_samples])

cat("Testing plotted groups:", paste(specified_groups, collapse = ", "), "\n")
set.seed(1523)
permanova_plot <- adonis2(aitchison_dist_plot ~ Group, 
                          data = plot_metadata, 
                          permutations = 999,
                          method = "euclidean")
print(permanova_plot)

# 4. Pairwise PERMANOVA tests for plotted groups
cat("\n4. Pairwise PERMANOVA tests for plotted groups:\n")

  # Manual pairwise tests using vegan::adonis2
  cat("\nPerforming manual pairwise PERMANOVA tests...\n")
  group_combinations <- combn(unique(metadata$Group), 2, simplify = FALSE)
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
    pair_metadata <- metadata %>%
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
  
