# Script to plot rank abundance of the top 10 most abundant ESVs
# from 16S SSU rRNA amplicon sequencing data
setwd("/home/pacifica/R/antarctica/S03/big_experiment/")

# Load required libraries
library(tidyverse)  # For data manipulation and visualization
library(gridExtra)  # For arranging multiple plots
library(grid)       # For grid graphics functions

# Define input file paths (update these to match your file locations)
count_table_path <- "ESV_table_PFX_2025.csv"
taxonomy_table_path <- "ESV_taxonomy_PFX_2025.csv"
metadata_path <- "Metadata_PFX_2025.csv"


# Read in the data
count_data <- read.csv(count_table_path, row.names = 1)
taxonomy_data <- read.csv(taxonomy_table_path, row.names = 1)
metadata <- read.csv(metadata_path, row.names = 1)

# Add sample names as a column in metadata
metadata$Sample <- rownames(metadata)

# 1. Filter for only the specified groups (in alphabetical order)
specified_groups <- c("BB.0", "BO.1", "BO.2", "OO.0", "OB.1", "OB.2")

# Create a new column in metadata that combines Order and Timing
metadata$Group <- paste(metadata$Order, metadata$Timing, sep = ".")

# Filter metadata to include only the specified groups
filtered_metadata <- metadata %>%
  filter(Group %in% specified_groups)

# Get sample names for the filtered groups
filtered_samples <- filtered_metadata$Sample

# Filter count data to include only those samples
filtered_counts <- count_data[, filtered_samples]

# 2. Calculate percent abundance for each ESV in each sample
calculate_percent_abundance <- function(count_matrix) {
  sample_totals <- colSums(count_matrix)
  percent_matrix <- count_matrix / rep(sample_totals, each = nrow(count_matrix)) * 100
  return(percent_matrix)
}

percent_abundance <- calculate_percent_abundance(filtered_counts)

# 3. Function to format genus names by removing text after first "_"
format_genus <- function(genus_name) {
  if(is.na(genus_name) || genus_name == "") {
    return("Unknown")
  }
  # Split by underscore and take the first part
  parts <- strsplit(genus_name, "_")[[1]]
  return(parts[1])
}

# 4. Prepare the dataset for all groups
all_groups_data <- data.frame()

for (group in specified_groups) {
  # Get samples for this group
  group_samples <- filtered_metadata %>%
    filter(Group == group) %>%
    pull(Sample)
  
  # Subset count data for this group
  group_counts <- filtered_counts[, group_samples]
  
  # Calculate percent abundance for this group
  group_percent <- calculate_percent_abundance(group_counts)
  
  # Calculate mean percent abundance for each ESV in this group
  mean_percent_by_group <- rowMeans(group_percent)
  
  # Get top 10 ESVs for this group
  top_esv_indices <- order(mean_percent_by_group, decreasing = TRUE)[1:10]
  top_esv_names <- rownames(filtered_counts)[top_esv_indices]
  
  # Create a mapping of ESV to rank based on mean abundance
  esv_rank_map <- data.frame(
    ESV = top_esv_names,
    MeanAbundance = mean_percent_by_group[top_esv_indices],
    Rank = 1:10
  )
  
  # Get sample-level data
  for (sample in group_samples) {
    sample_percent <- group_percent[top_esv_names, sample]
    
    # Create sample dataframe with proper ranks based on the mapping
    sample_df <- data.frame(
      ESV = top_esv_names,
      Sample = sample,
      Group = group,
      PercentAbundance = sample_percent
    )
    
    # Join with the rank mapping to ensure consistent ranks
    sample_df <- left_join(sample_df, esv_rank_map, by = "ESV")
    
    # Add taxonomic information
    sample_df$Phylum <- sapply(
      sample_df$ESV,
      function(esv) {
        if(esv %in% rownames(taxonomy_data)) {
          return(taxonomy_data[esv, "Phylum"])
        } else {
          return("Unknown")
        }
      }
    )
    
    sample_df$Genus <- sapply(
      sample_df$ESV,
      function(esv) {
        if(esv %in% rownames(taxonomy_data)) {
          genus <- taxonomy_data[esv, "Genus"]
          return(format_genus(genus))  # Format genus name
        } else {
          return("Unknown")
        }
      }
    )
    
    # Create new x-axis label format: "ESV_number: Genus"
    sample_df$ESV_Genus <- paste(sample_df$ESV, sample_df$Genus, sep = ": ")
    
    # Add to the combined dataset
    all_groups_data <- rbind(all_groups_data, sample_df)
  }
}

# 5. Set the correct order for the Group factor and create descriptive labels
# Define the desired order of groups
group_order <- c("BB.0", "BO.1", "BO.2", "OO.0", "OB.1", "OB.2")

# Create a named vector for group labels
group_labels <- c(
  "BB.0" = "A. Black mat before experiment",
  "BO.1" = "B. Black then orange mat after 1 season",
  "BO.2" = "C. Black then orange mat after 2 seasons",
  "OO.0" = "D. Orange mat before experiment",
  "OB.1" = "E. Orange then black mat after 1 season",
  "OB.2" = "F. Orange then black mat after 2 seasons"
)

# Convert Group to a factor with the specified order
all_groups_data$Group <- factor(all_groups_data$Group, levels = group_order)

# Make Phylum a factor with all possible levels
all_groups_data$Phylum <- factor(all_groups_data$Phylum)
# 6. Create individual plots for each group with proper x-axis ordering

# First, determine the overall y-axis range across all groups
y_range <- range(all_groups_data$PercentAbundance, na.rm = TRUE)
y_limits <- c(0, max(y_range) * 1.05)  # Add 5% padding at the top

# Get all unique phyla across all groups for consistent color mapping
all_phyla <- unique(all_groups_data$Phylum)
all_phyla <- sort(all_phyla[!is.na(all_phyla)])  # Sort for consistency

# Create consistent color palette for all phyla
paul_tol_colors <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#DDCC77", "#CC6677")
if(length(all_phyla) <= 6) {
  phylum_colors <- paul_tol_colors[1:length(all_phyla)]
} else {
  # If more than 6 phyla, use viridis to get more colors
  phylum_colors <- viridis::viridis(length(all_phyla), option = "D")
}
names(phylum_colors) <- all_phyla
create_group_plot <- function(group_name, data, group_labels, y_limits, phylum_colors, show_legend = FALSE) {
  # Filter data for this specific group
  group_data <- data %>%
    filter(Group == group_name) %>%
    arrange(Rank)
  
  # Create factor with correct order based on rank
  group_data$ESV_Genus_ordered <- factor(group_data$ESV_Genus, 
                                         levels = unique(group_data$ESV_Genus))
  
  # Create the plot
  p <- ggplot(group_data, aes(x = ESV_Genus_ordered, y = PercentAbundance, 
                              fill = Phylum, color = Phylum)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
      plot.title = element_text(hjust = 0.5, size = 7, face = "bold"),
      legend.position = if(show_legend) "bottom" else "none",
      legend.title = element_text(size = 6, face = "bold"),
      panel.border = element_rect(color = "lightgrey", fill = NA, size = 0.5)
    ) +
    labs(
      x = "",
      y = "% Abundance",
      title = group_labels[group_name],
      fill = "Phylum",
      color = "Phylum"
    ) +
    scale_y_continuous(limits = y_limits, labels = scales::percent_format(scale = 1)) +
    scale_fill_manual(values = phylum_colors, limits = names(phylum_colors)) +
    scale_color_manual(values = phylum_colors, limits = names(phylum_colors)) +
    guides(
      fill = guide_legend(override.aes = list(alpha = 1)),
      color = guide_legend(override.aes = list(alpha = 1))
    )
  
  return(p)
}

# Create individual plots for each group
plot_list <- list()
for(i in seq_along(group_order)) {
  group_name <- group_order[i]
  # Don't show legend on any individual plot
  plot_list[[i]] <- create_group_plot(group_name, all_groups_data, group_labels, 
                                      y_limits, phylum_colors, show_legend = FALSE)
}

# Create a separate legend plot
# Create a dummy plot that contains all phyla for consistent legend
# Create dummy data with all phyla
dummy_data <- data.frame(
  ESV_Genus_ordered = rep("Dummy", length(all_phyla)),
  PercentAbundance = rep(0, length(all_phyla)),
  Phylum = factor(all_phyla, levels = all_phyla)
)

# Create dummy plot just for extracting complete legend
legend_plot <- ggplot(dummy_data, aes(x = ESV_Genus_ordered, y = PercentAbundance, 
                                      fill = Phylum, color = Phylum)) +
  geom_point() +  # Just need some geom to create the legend
  scale_fill_manual(values = phylum_colors, limits = names(phylum_colors)) +
  scale_color_manual(values = phylum_colors, limits = names(phylum_colors)) +
  theme(legend.position = "bottom") +
  labs(fill = "Phylum", color = "Phylum") +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, shape = 22, size = 4)),
    color = guide_legend(override.aes = list(alpha = 1, shape = 22, size = 4))
  )

legend_only <- ggpubr::get_legend(legend_plot)

# Arrange plots in a 3x2 grid with legend at the bottom
combined_plot <- grid.arrange(
  plot_list[[1]], plot_list[[2]], plot_list[[3]],
  plot_list[[4]], plot_list[[5]], plot_list[[6]],
  legend_only,
  ncol = 3,
  nrow = 3,
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6),
                        c(7, 7, 7)),  # Legend spans all 3 columns
  heights = c(4, 4, 1)  # Make legend row shorter
)

# Print the combined plot
print(combined_plot)

# Save the combined plot
ggsave(
  "three_column_rank_abundance.pdf", 
  combined_plot, 
  width = 9, 
  height = 6
)

