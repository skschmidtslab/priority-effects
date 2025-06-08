
# Load necessary libraries
library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(ggplot2)
library(microbiome) # For additional visualization
library(patchwork) # For combining plots
library(vegan) # For PERMANOVA
library(ggrepel) # For volcano plot?

# Define input files
setwd("/home/pacifica/R/antarctica/S03/big_experiment/")
esv_table_file <- "ESV_table_PFX_Science_2023.csv"
taxonomy_file <- "ESV_taxonomy_PFX_Science_2023.csv"
metadata_file <- "Metadata_mapping_table_PFX_Science_2023.csv"

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
  
  
# 
# # Function to create a bar plot of the top differentially abundant taxa


# Complete script for ANCOMBC analysis with improved error handling

# 1. First, let's check the ancombc2 output structure with a simple example
# This will help understand what fields are available
check_ancombc_structure <- function(ps_obj, formula_str, group_var) {
  cat("\nRunning test ANCOMBC to check output structure\n")
  
  test_result <- ancombc2(
    data = ps_obj,
    fix_formula = formula_str,
    p_adj_method = "BH",
    prv_cut = 0.10,
    lib_cut = 1000,
    group = group_var,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05
  )
  
  cat("\nANCOMBC2 output structure:\n")
  print(names(test_result))
  
  # Check what's in res
  if ("res" %in% names(test_result)) {
    cat("\nFields in 'res':\n")
    print(names(test_result$res))
  }
  
  # Check what's in res_pair if it exists
  if ("res_pair" %in% names(test_result)) {
    cat("\nFields in 'res_pair':\n")
    print(names(test_result$res_pair))
    
    # Check if beta exists and its dimensions
    if (!is.null(test_result$res_pair$beta)) {
      cat("\nDimensions of 'res_pair$beta':", dim(test_result$res_pair$beta), "\n")
      cat("Column names of 'res_pair$beta':", colnames(test_result$res_pair$beta), "\n")
    } else {
      cat("\n'res_pair$beta' is NULL\n")
    }
  }
  
  return(test_result)
}

# 2. Improved function for ANCOMBC analysis
# Modified run_ancombc_analysis function to correctly extract results from ANCOMBC output
# Corrected run_ancombc_analysis function
run_ancombc_analysis <- function(ps_obj, formula_str, group_var, comparison_name) {
  cat("\nRunning ANCOMBC for", comparison_name, "\n")
  
  # Run ancombc2 with simplified parameters
  ancom_result <- ancombc2(
    data = ps_obj,
    fix_formula = formula_str,
    p_adj_method = "BH",
    prv_cut = 0.10,    # ESVs present in at least 10% of samples
    lib_cut = 1000,    # Minimum library size
    group = group_var,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    verbose = TRUE     # Add verbose output for debugging
  )
  
  # Print structure of result
  cat("\nStructure of ancom_result:\n")
  print(names(ancom_result))
  
  # Extract the results from ancom_result$res
  if ("res" %in% names(ancom_result) && !is.null(ancom_result$res)) {
    cat("\nUsing 'res' for two-group comparison\n")
    
    # Get column names from res
    res_cols <- names(ancom_result$res)
    
    # Find the column names for the group of interest
    group_levels <- levels(sample_data(ps_obj)[[group_var]])
    target_group <- group_levels[2]  # The second group (compared to first)
    
    # Find matching column names
    lfc_col <- grep(paste0("lfc_", group_var, target_group), res_cols, value = TRUE)
    se_col <- grep(paste0("se_", group_var, target_group), res_cols, value = TRUE)
    W_col <- grep(paste0("W_", group_var, target_group), res_cols, value = TRUE)
    p_col <- grep(paste0("p_", group_var, target_group), res_cols, value = TRUE)
    q_col <- grep(paste0("q_", group_var, target_group), res_cols, value = TRUE)
    diff_col <- grep(paste0("diff_", group_var, target_group), res_cols, value = TRUE)
    
    # If we can't find the exact pattern, try more general matching
    if (length(lfc_col) == 0) lfc_col <- grep("lfc_", res_cols, value = TRUE)[2]  # Skip intercept
    if (length(se_col) == 0) se_col <- grep("se_", res_cols, value = TRUE)[2]
    if (length(W_col) == 0) W_col <- grep("W_", res_cols, value = TRUE)[2]
    if (length(p_col) == 0) p_col <- grep("p_", res_cols, value = TRUE)[2]
    if (length(q_col) == 0) q_col <- grep("q_", res_cols, value = TRUE)[2]
    if (length(diff_col) == 0) diff_col <- grep("diff_", res_cols, value = TRUE)[2]
    
    cat("\nFound columns for extraction:", lfc_col, se_col, W_col, p_col, q_col, diff_col, "\n")
    
    # KEY FIX: Use the taxon names from the ancom_result$res directly
    # This ensures we only process taxa that were actually analyzed
    results_df <- data.frame(
      ESV = ancom_result$res$taxon
    )
    
    # Add columns if they exist
    if (length(lfc_col) > 0) results_df$lfc <- ancom_result$res[[lfc_col]]
    else results_df$lfc <- NA
    
    if (length(se_col) > 0) results_df$se <- ancom_result$res[[se_col]]
    else results_df$se <- NA
    
    if (length(W_col) > 0) results_df$W <- ancom_result$res[[W_col]]
    else results_df$W <- NA
    
    if (length(p_col) > 0) results_df$p_val <- ancom_result$res[[p_col]]
    else results_df$p_val <- NA
    
    if (length(q_col) > 0) results_df$q_val <- ancom_result$res[[q_col]]
    else results_df$q_val <- NA
    
    if (length(diff_col) > 0) results_df$diff_abn <- ancom_result$res[[diff_col]]
    else results_df$diff_abn <- NA
    
    # Add taxonomy - making sure we only look for ESVs that are in the results
    tax_subset <- tax_table(ps_obj)[results_df$ESV, ]
    results_df <- cbind(results_df, as.data.frame(tax_subset))
    
  } else {
    # Fallback: create empty results framework
    cat("\nNo standard result structure found - creating empty dataframe\n")
    
    results_df <- data.frame(
      ESV = character(0),
      lfc = numeric(0),
      se = numeric(0),
      W = numeric(0),
      p_val = numeric(0),
      q_val = numeric(0),
      diff_abn = logical(0)
    )
  }
  
  # Sort by p-value or q-value if available, handling NA values
  if (!all(is.na(results_df$q_val))) {
    results_df <- results_df[order(results_df$q_val, na.last = TRUE), ]
  } else if (!all(is.na(results_df$p_val))) {
    results_df <- results_df[order(results_df$p_val, na.last = TRUE), ]
  }
  
  # Return both the raw ANCOMBC result and processed results dataframe
  return(list(ancombc_obj = ancom_result, results = results_df))
}

# Modified volcano plot function with custom colors and specific taxa labeling
generate_volcano_plot <- function(ancom_results, title, sig_threshold = 0.05) {
  # Check if required columns exist
  if (!"lfc" %in% names(ancom_results$results) || 
      !"q_val" %in% names(ancom_results$results) ||
      all(is.na(ancom_results$results$lfc)) || 
      all(is.na(ancom_results$results$q_val))) {
    
    cat("Cannot generate volcano plot - missing required data columns.\n")
    return(NULL)
  }
  
  # Create a simplified taxonomy label
  ancom_results$results$tax_label <- paste0(
    ancom_results$results$Phylum, ";", 
    ancom_results$results$Class, ";",
    ancom_results$results$Genus
  )
  
  # Create custom significance and direction categories
  ancom_results$results$sig_direction <- ifelse(
    ancom_results$results$q_val >= sig_threshold, 
    "NS",
    ifelse(ancom_results$results$lfc < 0, "More in black", "More in orange")
  )
  
  # Define colors for the three categories
  color_palette <- c(
    "NS" = "grey",
    "More in black" = "black", 
    "More in orange" = "orange"
  )
  
  # Identify taxa of interest for labeling - ONLY SIGNIFICANT TAXA
  # 1. Top 10 most significant ESVs (only if they are significant)
  top_esvs <- data.frame()
  if (sum(!is.na(ancom_results$results$q_val) & ancom_results$results$q_val < sig_threshold) > 0) {
    significant_results <- ancom_results$results[ancom_results$results$q_val < sig_threshold & !is.na(ancom_results$results$q_val), ]
    top_esvs <- significant_results[order(significant_results$q_val), ][1:min(10, nrow(significant_results)), ]
  }
  
  # 2. Cyanobacteria (only if significant)
  cyanobacteria <- ancom_results$results[
    !is.na(ancom_results$results$Phylum) & 
      ancom_results$results$Phylum == "Cyanobacteria" &
      !is.na(ancom_results$results$q_val) &
      ancom_results$results$q_val < sig_threshold, 
  ]
  
  # 3. Specific genera of interest (only if significant)
  genera_of_interest <- c("Flavobacterium", "Polaromonas", "Cryobacterium")
  specific_genera <- ancom_results$results[
    !is.na(ancom_results$results$Genus) & 
      ancom_results$results$Genus %in% genera_of_interest &
      !is.na(ancom_results$results$q_val) &
      ancom_results$results$q_val < sig_threshold, 
  ]
  
  # Combine all taxa to be labeled (remove duplicates by ESV)
  taxa_to_label <- unique(rbind(
    top_esvs[, names(ancom_results$results)],
    cyanobacteria,
    specific_genera
  ))
  
  # Create labels - prioritize Genus, fall back to higher taxonomy if Genus is NA
  taxa_to_label$label <- ifelse(
    !is.na(taxa_to_label$Genus) & taxa_to_label$Genus != "",
    taxa_to_label$Genus,
    ifelse(
      !is.na(taxa_to_label$Family) & taxa_to_label$Family != "",
      paste0(taxa_to_label$Family, " (Family)"),
      paste0(taxa_to_label$Class, " (Class)")
    )
  )
  
  # Create the volcano plot
  p <- ggplot(ancom_results$results, aes(x = lfc, y = -log10(q_val), color = sig_direction)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    scale_color_manual(values = color_palette) +
    theme_bw() +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-log10(adjusted p-value)",
      color = "Significance"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10)
    )
  
  # Add labels for taxa of interest if ggrepel is available
  if (nrow(taxa_to_label) > 0 && !all(is.na(taxa_to_label$label))) {
    # Check if ggrepel is available, if not use regular geom_text
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = taxa_to_label,
        aes(label = label),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.3,
        show.legend = FALSE
      )
    } else {
      # Fallback to regular geom_text if ggrepel is not available
      warning("ggrepel package not available, using geom_text instead")
      p <- p + geom_text(
        data = taxa_to_label,
        aes(label = label),
        size = 3,
        vjust = -0.5,
        show.legend = FALSE
      )
    }
  }
  
  # Add annotation showing number of significant taxa in each direction
  n_sig_pos <- sum(ancom_results$results$sig_direction == "More in orange", na.rm = TRUE)
  n_sig_neg <- sum(ancom_results$results$sig_direction == "More in black", na.rm = TRUE)
  
  annotation_text <- paste0(
    "Significant taxa:\n",
    "Orange first: ", n_sig_pos, "\n",
    "Black first: ", n_sig_neg
  )
  
  # Add annotation to the plot
  p <- p + annotate(
    "text",
    x = Inf, y = -Inf,
    label = annotation_text,
    hjust = 1.1, vjust = -0.1,
    size = 3,
    color = "black",
    fontface = "italic"
  )
  
  return(p)
}

plot_top_taxa <- function(ancom_results, physeq_obj, group_var, n_taxa = 10, sig_threshold = 0.05) {
  # Get significant taxa
  sig_taxa <- ancom_results$results[ancom_results$results$q_val < sig_threshold, ]
  
  if(nrow(sig_taxa) == 0) {
    cat("No significant taxa found for this comparison.\n")
    return(NULL)
  }
  
  # Sort by absolute log fold change
  sig_taxa$abs_lfc <- abs(sig_taxa$lfc)
  sig_taxa <- sig_taxa[order(sig_taxa$abs_lfc, decreasing = TRUE), ]
  
  # Select top n taxa
  if(nrow(sig_taxa) > n_taxa) {
    top_taxa <- sig_taxa$ESV[1:n_taxa]
  } else {
    top_taxa <- sig_taxa$ESV
  }
  
  # Subset phyloseq object to these taxa
  ps_subset <- prune_taxa(top_taxa, physeq_obj)
  
  # Transform to relative abundance
  ps_rel <- transform_sample_counts(ps_subset, function(x) x / sum(x))
  
  # Extract abundance data
  abund_data <- psmelt(ps_rel)
  
  # Add taxonomy information
  abund_data$Taxon <- paste0(abund_data$Genus, " (", abund_data$Family, ")")
  
  # Create a more readable label
  abund_data[[group_var]] <- as.factor(abund_data[[group_var]])
  
  # Plot
  p <- ggplot(abund_data, aes(x = Taxon, y = Abundance, fill = get(group_var))) +
    geom_boxplot(outlier.size = 0.5) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    ) +
    labs(
      x = "Taxon",
      y = "Relative Abundance",
      fill = group_var
    ) +
    scale_fill_brewer(palette = "Set1")
  
  return(p)
}

# Function to create taxonomic composition barplots for different groups
plot_taxonomic_composition <- function(ps_obj, level = "Phylum", group_var = "OrderTiming") {
  # Agglomerate taxa at the specified level
  ps_glom <- tax_glom(ps_obj, taxrank = level)
  
  # Convert to relative abundance
  ps_rel <- transform_sample_counts(ps_glom, function(x) x / sum(x))
  
  # Get the data
  plot_data <- psmelt(ps_rel)
  
  # Organize taxa by abundance (group low abundance taxa as "Other")
  taxa_abund <- aggregate(Abundance ~ get(level), plot_data, mean)
  names(taxa_abund)[1] <- level
  taxa_abund <- taxa_abund[order(taxa_abund$Abundance, decreasing = TRUE), ]
  
  # Keep top 10 taxa, group the rest as "Other"
  top_taxa <- taxa_abund[[level]][1:min(10, nrow(taxa_abund))]
  plot_data[[level]] <- ifelse(plot_data[[level]] %in% top_taxa,
                               as.character(plot_data[[level]]),
                               "Other")
  
  # Summarize by group
  plot_data$group <- plot_data[[group_var]]
  
  # Plot
  p <- ggplot(plot_data, aes(x = group, y = Abundance, fill = get(level))) +
    geom_bar(stat = "identity", position = "fill") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    labs(
      x = "Group",
      y = "Relative Abundance",
      title = paste("Taxonomic Composition at", level, "Level")
    ) +
    scale_fill_brewer(palette = "Paired", name = level)
  
  return(p)
}

# Function to write significant results to file
write_sig_results <- function(ancom_results, filename, sig_threshold = 0.05) {
  sig_results <- ancom_results$results[ancom_results$results$q_val < sig_threshold, ]
  write.csv(sig_results, filename, row.names = FALSE)
  return(sig_results)
}





# 3. Using the functions

# First run a diagnostic check to understand the structure
ps_time0_check <- check_ancombc_structure(
  ps_obj = ps_time0_filtered,
  formula_str = "OrderGroup",
  group_var = "OrderGroup"
)

# Then run the actual analysis with the improved function
ancom_time0 <- run_ancombc_analysis(
  ps_obj = ps_time0_filtered,
  formula_str = "OrderGroup",
  group_var = "OrderGroup", 
  comparison_name = "OO.0 vs BB.0"
)

# Continue with visualization and result export
# Generate volcano plot if results are available
volcano_plot_time0 <- generate_volcano_plot(
  ancom_results = ancom_time0,
  title = "A. Starting mats",
  sig_threshold = 0.05
)
print(volcano_plot_time0)

# Generate taxa abundance plot for significant taxa
taxa_plot_time0 <- plot_top_taxa(
  ancom_results = ancom_time0,
  physeq_obj = ps_time0_filtered,
  group_var = "OrderGroup",
  n_taxa = 10,
  sig_threshold = 0.05
)
if (!is.null(taxa_plot_time0)) {
  print(taxa_plot_time0)
}

# Write significant results to file
sig_results_time0 <- write_sig_results(
  ancom_results = ancom_time0,
  filename = "significant_taxa_OO0_vs_BB0.csv",
  sig_threshold = 0.05
)



##### Time point 2 (ob vs bo after 1 season )
# First run a diagnostic check to understand the structure
ps_time1_check <- check_ancombc_structure(
  ps_obj = ps_time1_filtered,
  formula_str = "OrderGroup",
  group_var = "OrderGroup"
)

# Then run the actual analysis with the improved function
ancom_time1 <- run_ancombc_analysis(
  ps_obj = ps_time1_filtered,
  formula_str = "OrderGroup",
  group_var = "OrderGroup", 
  comparison_name = "OO.0 vs BB.0"
)

# Continue with visualization and result export
# Generate volcano plot if results are available
volcano_plot_time1 <- generate_volcano_plot(
  ancom_results = ancom_time1,
  title = "B. After one season",
  sig_threshold = 0.05
)
print(volcano_plot_time1)

# Generate taxa abundance plot for significant taxa
taxa_plot_time1 <- plot_top_taxa(
  ancom_results = ancom_time1,
  physeq_obj = ps_time1_filtered,
  group_var = "OrderGroup",
  n_taxa = 10,
  sig_threshold = 0.05
)
if (!is.null(taxa_plot_time1)) {
  print(taxa_plot_time1)
}

# Write significant results to file
sig_results_time1 <- write_sig_results(
  ancom_results = ancom_time1,
  filename = "significant_taxa_OO0_vs_BB0.csv",
  sig_threshold = 0.05
)


##### TIME POINT 3 (OB VS BO AFTER 2 SEASONS)
# First run a diagnostic check to understand the structure
ps_time2_check <- check_ancombc_structure(
  ps_obj = ps_time2_filtered,
  formula_str = "OrderGroup",
  group_var = "OrderGroup"
)

# Then run the actual analysis with the improved function
ancom_time2 <- run_ancombc_analysis(
  ps_obj = ps_time2_filtered,
  formula_str = "OrderGroup",
  group_var = "OrderGroup", 
  comparison_name = "After two seasons"
)

# Continue with visualization and result export
# Generate volcano plot if results are available
volcano_plot_time2 <- generate_volcano_plot(
  ancom_results = ancom_time2,
  title = "C. After two seasons",
  sig_threshold = 0.05
)
print(volcano_plot_time2)

# Generate taxa abundance plot for significant taxa
taxa_plot_time2 <- plot_top_taxa(
  ancom_results = ancom_time2,
  physeq_obj = ps_time2_filtered,
  group_var = "OrderGroup",
  n_taxa = 10,
  sig_threshold = 0.05
)
if (!is.null(taxa_plot_time2)) {
  print(taxa_plot_time2)
}

# Write significant results to file
sig_results_time2 <- write_sig_results(
  ancom_results = ancom_time2,
  filename = "significant_taxa_OO0_vs_BB0.csv",
  sig_threshold = 0.05
)


# Combine volcano plots
combined_volcano <- volcano_plot_time0 / volcano_plot_time1 / volcano_plot_time2 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save the combined plot
ggsave("ancombc_volcano_plots.pdf", combined_volcano, width = 6, height = 12)
ggsave("ancombc_volcano_plots.png", combined_volcano, width = 6, height = 12)
cat("Combined volcano plots saved to: ancombc_volcano_plots.pdf\n")















# Generate top taxa plots if there are significant results
if(nrow(sig_time0) > 0) {
  top_taxa_plot0 <- plot_top_taxa(ancom_time0, ps_time0_filtered, "OrderGroup")
  ggsave("top_taxa_OO0_vs_BB0.pdf", top_taxa_plot0, width = 10, height = 7)
  cat("Top taxa plot for OO.0 vs BB.0 saved to: top_taxa_OO0_vs_BB0.pdf\n")
}

if(nrow(sig_time1) > 0) {
  top_taxa_plot1 <- plot_top_taxa(ancom_time1, ps_time1_filtered, "OrderGroup")
  ggsave("top_taxa_OB1_vs_BO1.pdf", top_taxa_plot1, width = 10, height = 7)
  cat("Top taxa plot for OB.1 vs BO.1 saved to: top_taxa_OB1_vs_BO1.pdf\n")
}

if(nrow(sig_time2) > 0) {
  top_taxa_plot2 <- plot_top_taxa(ancom_time2, ps_time2_filtered, "OrderGroup")
  ggsave("top_taxa_OB2_vs_BO2.pdf", top_taxa_plot2, width = 10, height = 7)
  cat("Top taxa plot for OB.2 vs BO.2 saved to: top_taxa_OB2_vs_BO2.pdf\n")
}

# Create taxonomic composition plots
phylum_plot <- plot_taxonomic_composition(physeq, "Phylum", "OrderTiming")
ggsave("taxonomic_composition_phylum.pdf", phylum_plot, width = 10, height = 6)
cat("Phylum composition plot saved to: taxonomic_composition_phylum.pdf\n")

family_plot <- plot_taxonomic_composition(physeq, "Family", "OrderTiming") 
ggsave("taxonomic_composition_family.pdf", family_plot, width = 10, height = 6)
cat("Family composition plot saved to: taxonomic_composition_family.pdf\n")

# Calculate and plot alpha diversity metrics
plot_alpha_diversity <- function(ps_obj, group_var = "OrderTiming") {
  # Calculate alpha diversity metrics
  alpha_div <- data.frame(
    sample_data(ps_obj),
    Observed = phyloseq::estimate_richness(ps_obj, measures = "Observed")$Observed,
    Shannon = phyloseq::estimate_richness(ps_obj, measures = "Shannon")$Shannon,
    Simpson = phyloseq::estimate_richness(ps_obj, measures = "Simpson")$Simpson
  )
  
  # Melt for plotting
  alpha_div_long <- reshape2::melt(alpha_div, 
                                   id.vars = colnames(sample_data(ps_obj)), 
                                   variable.name = "Metric", 
                                   value.name = "Value")
  
  # Plot
  p <- ggplot(alpha_div_long, aes(x = get(group_var), y = Value, fill = get(group_var))) +
    geom_boxplot() +
    facet_wrap(~ Metric, scales = "free_y") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(
      x = group_var,
      y = "Alpha Diversity",
      title = "Alpha Diversity by Group"
    ) +
    scale_fill_brewer(palette = "Set1")
  
  return(p)
}

# Generate alpha diversity plot
alpha_div_plot <- plot_alpha_diversity(physeq)
ggsave("alpha_diversity.pdf", alpha_div_plot, width = 10, height = 6)
cat("Alpha diversity plot saved to: alpha_diversity.pdf\n")

# Perform beta diversity analysis
# Calculate Bray-Curtis distance
ps_rel_all <- transform_sample_counts(physeq, function(x) x / sum(x))
ord <- ordinate(ps_rel_all, method = "PCoA", distance = "bray")

# Plot PCoA
beta_plot <- plot_ordination(ps_rel_all, ord, color = "OrderTiming") +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    title = "PCoA based on Bray-Curtis Distance",
    color = "Group"
  ) +
  scale_color_brewer(palette = "Set1")

ggsave("beta_diversity_pcoa.pdf", beta_plot, width = 8, height = 6)
cat("Beta diversity PCoA plot saved to: beta_diversity_pcoa.pdf\n")

# Perform PERMANOVA test
# Extract OTU table and convert to relative abundance
otu_table_rel <- t(as(otu_table(ps_rel_all), "matrix"))
# Calculate distance matrix
dist_matrix <- vegdist(otu_table_rel, method = "bray")
# Get metadata
metadata_df <- as(sample_data(ps_rel_all), "data.frame")
# Run PERMANOVA
cat("\n\n--- PERMANOVA RESULTS ---\n")
permanova_result <- adonis2(dist_matrix ~ OrderTiming, data = metadata_df)
print(permanova_result)

# Save the permanova results
capture.output(permanova_result, file = "permanova_results.txt")
cat("PERMANOVA results saved to: permanova_results.txt\n")

# Create a summary of the analysis
cat("\n\n--- ANALYSIS SUMMARY ---\n")
cat("Total number of ESVs analyzed:", ntaxa(physeq), "\n")
cat("Total number of samples:", nsamples(physeq), "\n")
cat("Group comparisons performed: OO.0 vs BB.0, OB.1 vs BO.1, OB.2 vs BO.2\n\n")

cat("Differential abundance results:\n")
cat("- OO.0 vs BB.0:", nrow(sig_time0), "significantly different ESVs\n")
cat("- OB.1 vs BO.1:", nrow(sig_time1), "significantly different ESVs\n")
cat("- OB.2 vs BO.2:", nrow(sig_time2), "significantly different ESVs\n\n")

cat("Output files created:\n")
cat("1. Significant taxa tables: significant_taxa_*.csv\n")
cat("2. Volcano plots: ancombc_volcano_plots.pdf\n")
cat("3. Top differentially abundant taxa plots: top_taxa_*.pdf\n")
cat("4. Taxonomic composition plots: taxonomic_composition_*.pdf\n")
cat("5. Alpha diversity plot: alpha_diversity.pdf\n")
cat("6. Beta diversity plot: beta_diversity_pcoa.pdf\n")
cat("7. PERMANOVA results: permanova_results.txt\n")