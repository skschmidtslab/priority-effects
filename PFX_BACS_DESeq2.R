# A script to analyze differential abundance of ESVs in treatments of 
# the Antarctic field priority effects experiment

# Script coded with the help of Anthropic's Claude Sonnet 4 LLM

# Load necessary libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ggrepel) # For volcano plot?
library(EnhancedVolcano)

# Define input files
setwd("/home/pacifica/R/antarctica/S03/big_experiment/")
esv_table_file <- "ESV_table_PFX_2025.csv"
taxonomy_file <- "ESV_taxonomy_PFX_2025.csv"
metadata_file <- "Metadata_PFX_2025.csv"

# Read input files
esv_table <- read.csv(esv_table_file, row.names = 1)
taxonomy_table <- read.csv(taxonomy_file, row.names = 1)
metadata <- read.csv(metadata_file, row.names = 1)

# Convert to matrices for phyloseq
otu_mat <- as.matrix(esv_table)
tax_mat <- as.matrix(taxonomy_table)

# Create phyloseq objects
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
META <- sample_data(metadata)

# Ensure sample names match across objects
sample_names(META) <- rownames(metadata)

# Create phyloseq object
physeq <- phyloseq(OTU, TAX, META)

# Create a combined variable for Order and Timing
sample_data(physeq)$OrderTiming <- paste0(sample_data(physeq)$Order, ".", sample_data(physeq)$Timing)

# Create subsets for the specified group comparisons
# 1. OO.0 vs BB.0
cat("\n--- ANALYSIS 1: OO.0 vs BB.0 ---\n")
ps_time0 <- subset_samples(physeq, Timing == 0 & Order %in% c("OO", "BB"))
sample_data(ps_time0)$OrderGroup <- factor(sample_data(ps_time0)$Order, levels = c("BB", "OO"))
# Check sample sizes
cat("Sample sizes for Time 0 comparison (OO.0 vs BB.0):\n")
print(table(sample_data(ps_time0)$OrderGroup))
# Filter taxa with low prevalence
ps_time0_filtered <- filter_taxa(ps_time0, function(x) sum(x > 0) >= (0.1 * length(x)), TRUE)
cat("Number of ESVs after filtering:", ntaxa(ps_time0_filtered), "\n")
# Convert phyloseq to DESeq2 object
dds.0 <- phyloseq_to_deseq2(ps_time0_filtered, ~ OrderGroup)
# Use poscounts for size factor estimation (crucial for microbiome data)
dds.0 <- DESeq(dds.0, sfType = "poscounts")
# Extract results
res.0 <- results(dds.0, contrast = c("OrderGroup","OO","BB"))

# VOLCANO PLOT - Overview of differential abundance
# Convert results to dataframe
res0_df <- as.data.frame(res.0)
res0_df$taxon <- rownames(res0_df)
res0_df$significant <- res0_df$padj < 0.05 & !is.na(res0_df$padj)

# Add abundance information to your results dataframe
# Calculate mean abundance across all samples for each taxon
abundance_data <- as.data.frame(otu_table(ps_time0_filtered))
mean_abundance <- rowMeans(abundance_data)

# Add abundance info to results
res0_df$mean_abundance <- mean_abundance[rownames(res0_df)]

# Extract taxonomy information
tax_info <- as.data.frame(tax_table(ps_time0_filtered))

# Create more informative labels (e.g., using Genus or Family)
res0_df$tax_label <- paste0(substr(rownames(res0_df), 1, 8), ": ",
                            sub("_.*", "", tax_info[rownames(res0_df), "Genus"]))

# Identify top abundant significant taxa (e.g., top 10)
top_abundant_sig <- res0_df %>%
  filter(significant == TRUE) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = 10) %>% 
  bind_rows(res0_df %>% filter(abs(log2FoldChange) > 10)) %>% 
  unique()

# create a new column that combines significance and direction
res0_df$color_group <- ifelse(!res0_df$significant, "non_significant",
                              ifelse(res0_df$log2FoldChange > 0, "sig_positive", "sig_negative"))

# Create plot with labels
p0 <- ggplot(res0_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color_group), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("non_significant" = "gray60", 
                                "sig_positive" = "orange", 
                                "sig_negative" = "black"),
                     labels = c("non_significant" = "NS",
                                "sig_positive" = "More in orange or OB",
                                "sig_negative" = "More in black or BO")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  # Add labels for top abundant significant taxa
  geom_text_repel(data = top_abundant_sig,
                  aes(label = tax_label),
                  size = 3,
                  max.overlaps = 15) +
  labs(x = "Log2 Fold Change", 
       y = "-log10(Adjusted P-value)",
       title = "A. Starting mats",
       subtitle = paste0("Significant taxa: ", sum(res0_df$significant, na.rm = TRUE)),
       color = "Significance") +  # Add legend title
  theme_minimal() +
  theme(legend.position = "bottom")
p0


# Now for second comparison: OB vs BO after first season  
# 2. OB.1 vs BO.1
cat("\n\n--- ANALYSIS 2: OB.1 vs BO.1 ---\n")
ps_time1 <- subset_samples(physeq, Timing == 1 & Order %in% c("OB", "BO"))
sample_data(ps_time1)$OrderGroup <- factor(sample_data(ps_time1)$Order, levels = c("BO", "OB"))

# Check sample sizes
cat("Sample sizes for Time 1 comparison (OB.1 vs BO.1):\n")
print(table(sample_data(ps_time1)$OrderGroup))

# Filter taxa with low prevalence
ps_time1_filtered <- filter_taxa(ps_time1, function(x) sum(x > 0) >= (0.1 * length(x)), TRUE)
cat("Number of ESVs after filtering:", ntaxa(ps_time1_filtered), "\n")

# Convert phyloseq to DESeq2 object
dds.1 <- phyloseq_to_deseq2(ps_time1_filtered, ~ OrderGroup)
# Use poscounts for size factor estimation (crucial for microbiome data)
dds.1 <- DESeq(dds.1, sfType = "poscounts")
# Extract results
res.1 <- results(dds.1, contrast = c("OrderGroup","OB","BO"))

# VOLCANO PLOT - Overview of differential abundance
# Convert results to dataframe
res1_df <- as.data.frame(res.1)
res1_df$taxon <- rownames(res1_df)
res1_df$significant <- res1_df$padj < 0.05 & !is.na(res1_df$padj)

# Add abundance information to your results dataframe
# Calculate mean abundance across all samples for each taxon
abundance_data_1 <- as.data.frame(otu_table(ps_time1_filtered))
mean_abundance_1 <- rowMeans(abundance_data_1)

# Add abundance info to results
res1_df$mean_abundance_1 <- mean_abundance_1[rownames(res1_df)]

# Extract taxonomy information
tax_info_1 <- as.data.frame(tax_table(ps_time1_filtered))

# Create more informative labels (e.g., using Genus or Family)
res1_df$tax_label <- paste0(substr(rownames(res1_df), 1, 8), ": ",
                            sub("_.*", "", tax_info[rownames(res1_df), "Genus"]))


# Identify top abundant significant taxa (e.g., top 10)
top_abundant_sig_1 <- res1_df %>%
  filter(significant == TRUE) %>%
  arrange(desc(mean_abundance_1)) %>%
  slice_head(n = 10) %>% 
  bind_rows(res1_df %>% filter(abs(log2FoldChange) > 10)) %>% 
  unique()

# create a new column that combines significance and direction
res1_df$color_group <- ifelse(!res1_df$significant, "non_significant",
                              ifelse(res1_df$log2FoldChange > 0, "sig_positive", "sig_negative"))

# Create plot with labels
p1 <- ggplot(res1_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color_group), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("non_significant" = "gray60", 
                                "sig_positive" = "orange", 
                                "sig_negative" = "black"),
                     labels = c("non_significant" = "NS",
                                "sig_positive" = "More in orange or OB",
                                "sig_negative" = "More in black or BO")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  # Add labels for top abundant significant taxa
  geom_text_repel(data = top_abundant_sig_1,
                  aes(label = tax_label),
                  size = 3,
                  max.overlaps = 15) +
  labs(x = "Log2 Fold Change", 
       y = "-log10(Adjusted P-value)",
       title = "B. After one season",
       subtitle = paste0("Significant taxa: ", sum(res1_df$significant, na.rm = TRUE)),
       color = "Significance") +  # Add legend title
  theme_minimal() +
  theme(legend.position = "bottom")
p1



# And for third comparison: OB vs BO after second season of growth
# 3. OB.2 vs BO.2
cat("\n\n--- ANALYSIS 3: OB.2 vs BO.2 ---\n")
ps_time2 <- subset_samples(physeq, Timing == 2 & Order %in% c("OB", "BO"))
sample_data(ps_time2)$OrderGroup <- factor(sample_data(ps_time2)$Order, levels = c("BO", "OB"))

# Check sample sizes
cat("Sample sizes for Time 2 comparison (OB.2 vs BO.2):\n")
print(table(sample_data(ps_time2)$OrderGroup))

# Filter taxa with low prevalence
ps_time2_filtered <- filter_taxa(ps_time2, function(x) sum(x > 0) >= (0.1 * length(x)), TRUE)
cat("Number of ESVs after filtering:", ntaxa(ps_time2_filtered), "\n")

# Convert phyloseq to DESeq2 object
dds.2 <- phyloseq_to_deseq2(ps_time2_filtered, ~ OrderGroup)
# Use poscounts for size factor estimation (crucial for microbiome data)
dds.2 <- DESeq(dds.2, sfType = "poscounts")
# Extract results
res.2 <- results(dds.2, contrast = c("OrderGroup","OB","BO"))

# VOLCANO PLOT - Overview of differential abundance
# Convert results to dataframe
res2_df <- as.data.frame(res.2)
res2_df$taxon <- rownames(res2_df)
res2_df$significant <- res2_df$padj < 0.05 & !is.na(res2_df$padj)

# Add abundance information to your results dataframe
# Calculate mean abundance across all samples for each taxon
abundance_data_2 <- as.data.frame(otu_table(ps_time2_filtered))
mean_abundance_2 <- rowMeans(abundance_data_2)

# Add abundance info to results
res2_df$mean_abundance_2 <- mean_abundance_2[rownames(res2_df)]

# Extract taxonomy information
tax_info_2 <- as.data.frame(tax_table(ps_time2_filtered))

# Create more informative labels (e.g., using Genus or Family)
res2_df$tax_label <- paste0(substr(rownames(res2_df), 1, 8), ": ",
                            sub("_.*", "", tax_info[rownames(res2_df), "Genus"]))

# Identify top abundant significant taxa (e.g., top 10)
top_abundant_sig_2 <- res2_df %>%
  filter(significant == TRUE) %>%
  arrange(desc(mean_abundance_2)) %>%
  slice_head(n = 10) %>% 
  bind_rows(res2_df %>% filter(abs(log2FoldChange) > 8)) %>% 
  unique()

# create a new column that combines significance and direction
res2_df$color_group <- ifelse(!res2_df$significant, "non_significant",
                              ifelse(res2_df$log2FoldChange > 0, "sig_positive", "sig_negative"))

# Create plot with labels
p2 <- ggplot(res2_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = color_group), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("non_significant" = "gray60", 
                                "sig_positive" = "orange", 
                                "sig_negative" = "black"),
                     labels = c("non_significant" = "NS",
                                "sig_positive" = "More in orange or OB",
                                "sig_negative" = "More in black or BO")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  # Add labels for top abundant significant taxa
  geom_text_repel(data = top_abundant_sig_2,
                  aes(label = tax_label),
                  size = 3,
                  max.overlaps = 15) +
  labs(x = "Log2 Fold Change", 
       y = "-log10(Adjusted P-value)",
       title = "C. After two seasons",
       subtitle = paste0("Significant taxa: ", sum(res2_df$significant, na.rm = TRUE)),
       color = "Significance") +  # Add legend title
  theme_minimal() +
  theme(legend.position = "bottom")
p2


### Combine the plots for a figure
# Load patchwork library for combining plots
library(patchwork)

# Remove individual legends from each plot
p0_no_legend <- p0 + theme(legend.position = "none")
p1_no_legend <- p1 + theme(legend.position = "none") 
p2_no_legend <- p2 + theme(legend.position = "none")

# Combine plots with patchwork and add a single legend at the bottom
combined_plot <- (p0_no_legend | p1_no_legend | p2_no_legend) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
combined_plot

# Optional: Save the plot
ggsave("combined_volcano_plots.png", combined_plot, 
       width = 15, height = 6, dpi = 300, bg = "white")

write_csv(res0_df,"Table_S2_Dif_abund_starting_mats.csv")
write_csv(res1_df,"Table_S3_Dif_abund_first_season.csv")
write_csv(res2_df,"Table_S4_Dif_abund_second_season.csv")
