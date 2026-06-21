# Script to plot relative abundance of pre-defined ESVs of interest
# across treatment groups/timepoints from 16S SSU rRNA amplicon sequencing data
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
# (computed once across all filtered samples -- a sample's percent abundance
# only depends on that sample's own column total, so it does not need to be
# recalculated separately within each group)
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

# 4. Define the pre-determined set of taxa to plot, in the desired display order
target_esvs <- c("ESV_1", "ESV_18", "ESV_5", "ESV_6", "ESV_4", "ESV_2")

# Check that all target ESVs are present in the count table; warn if not
missing_esvs <- setdiff(target_esvs, rownames(filtered_counts))
if (length(missing_esvs) > 0) {
  warning(paste("The following target ESVs were not found in the count table and will be skipped:",
                paste(missing_esvs, collapse = ", ")))
  target_esvs <- intersect(target_esvs, rownames(filtered_counts))
}

# Create a mapping of ESV to a fixed rank/order (used to control x-axis ordering)
esv_rank_map <- data.frame(
  ESV = target_esvs,
  Rank = seq_along(target_esvs)
)

# 5. Prepare the dataset for all groups, using the same fixed taxon list in every group
all_groups_data <- data.frame()

for (group in specified_groups) {
  # Get samples for this group
  group_samples <- filtered_metadata %>%
    filter(Group == group) %>%
    pull(Sample)
  
  # Get sample-level data for the target ESVs
  for (sample in group_samples) {
    sample_percent <- percent_abundance[target_esvs, sample]
    
    # Create sample dataframe with the fixed taxon order
    sample_df <- data.frame(
      ESV = target_esvs,
      Sample = sample,
      Group = group,
      PercentAbundance = sample_percent
    )
    
    # Join with the rank mapping to ensure consistent ordering
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

# 6. Set the correct order for the Group factor and create descriptive labels
# Define the desired order of groups (these encode timepoint within each
# treatment lineage: BB.0 -> BO.1 -> BO.2 is the black-origin lineage,
# OO.0 -> OB.1 -> OB.2 is the orange-origin lineage)
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

# 7. Create individual plots for each group with proper x-axis ordering

# First, determine the overall y-axis range across all groups, for consistent
# scaling so the panels are directly comparable
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
  
  # Create factor with correct order based on the fixed taxon order
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
      panel.border = element_rect(color = "lightgrey", fill = NA, linewidth = 0.5)
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
  "selected_taxa_abundance_by_group.pdf",
  combined_plot,
  width = 9,
  height = 6
)

# ============================================================
# FIGURE 2: Abundance of each taxon over time, colored by mat lineage
# Each panel is one taxon; x-axis is the timepoint (0, 1, 2); points are
# colored/filled by which mat color the plot originated from:
#   - Black mat origin (BB.0, BO.1, BO.2): black outline, gray fill
#   - Orange mat origin (OO.0, OB.1, OB.2): orange outline, orange fill
# ============================================================

# Derive Timing (0, 1, 2) and mat lineage directly from the Group labels
# already present in all_groups_data (e.g. "BO.1" -> Timing = 1, lineage
# determined by the first letter of the group code)
trajectory_data <- all_groups_data %>%
  mutate(
    Timing = as.numeric(sub("^[A-Za-z]+\\.", "", as.character(Group))),
    Lineage = ifelse(substr(as.character(Group), 1, 1) == "B",
                     "Black mat first", "Orange mat first")
  )

trajectory_data$Lineage <- factor(trajectory_data$Lineage,
                                  levels = c("Black mat first", "Orange mat first"))

# Keep taxon panels in the same fixed order used in Figure 1
trajectory_data <- trajectory_data %>% arrange(Rank)
trajectory_data$ESV_Genus <- factor(trajectory_data$ESV_Genus,
                                    levels = unique(trajectory_data$ESV_Genus))

# Define point colors (outline) and fills per lineage
lineage_colors <- c("Black mat first" = "black", "Orange mat first" = "darkorange3")
lineage_fills  <- c("Black mat first" = "black", "Orange mat first" = "orange")

trajectory_plot <- ggplot(trajectory_data, aes(x = Timing, y = PercentAbundance,
                                               color = Lineage, fill = Lineage)) +
  geom_point(shape = 21, size = 2.5, alpha = 0.85, stroke = 0.6,
             position = position_jitter(width = 0.08, height = 0)) +
  stat_summary(fun = mean, geom = "line", aes(group = Lineage), linewidth = 0.6) +
  facet_wrap(~ ESV_Genus, ncol = 3, scales = "free_y") +
  scale_x_continuous(breaks = c(0, 1, 2)) +
  scale_color_manual(values = lineage_colors) +
  scale_fill_manual(values = lineage_fills) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 8, face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    x = "Timepoint (season)",
    y = "% Abundance",
    title = ""
  )

# Print the trajectory plot
print(trajectory_plot)

# Save the trajectory plot
ggsave(
  "selected_taxa_trajectories_by_lineage.pdf",
  trajectory_plot,
  width = 9,
  height = 6
)
