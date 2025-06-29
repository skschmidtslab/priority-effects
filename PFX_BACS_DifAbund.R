
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
  #ps_time0_filtered <- filter_taxa(ps_time0, function(x) sum(x > 0) >= (0.1 * length(x)), TRUE)
  ps_time0_filtered <- filter_taxa(ps_time0, function(x) sum(x > 0) >= 2, TRUE)
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

# 2. ANCOMBC analysis
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



# Modified volcano plot function with custom colors and ESV-based labeling
generate_volcano_plot <- function(ancom_results, title, sig_threshold = 0.05) {
  # Check if required columns exist
  if (!"lfc" %in% names(ancom_results$results) || 
      !"q_val" %in% names(ancom_results$results) ||
      all(is.na(ancom_results$results$lfc)) || 
      all(is.na(ancom_results$results$q_val))) {
    
    cat("Cannot generate volcano plot - missing required data columns.\n")
    return(NULL)
  }
  
  # Debug information
  cat("\n=== VOLCANO PLOT DEBUG INFO ===\n")
  cat("Number of rows in results:", nrow(ancom_results$results), "\n")
  cat("LFC range:", range(ancom_results$results$lfc, na.rm = TRUE), "\n")
  cat("Q-value range:", range(ancom_results$results$q_val, na.rm = TRUE), "\n")
  cat("Number of non-NA LFC values:", sum(!is.na(ancom_results$results$lfc)), "\n")
  cat("Number of non-NA q-values:", sum(!is.na(ancom_results$results$q_val)), "\n")
  cat("Number of infinite LFC values:", sum(is.infinite(ancom_results$results$lfc)), "\n")
  cat("Number of infinite q-values:", sum(is.infinite(ancom_results$results$q_val)), "\n")
  
  # Remove rows with NA or infinite values in key columns
  valid_rows <- !is.na(ancom_results$results$lfc) & 
    !is.na(ancom_results$results$q_val) &
    !is.infinite(ancom_results$results$lfc) &
    !is.infinite(ancom_results$results$q_val) &
    ancom_results$results$q_val > 0  # q-values must be positive for log transformation
  
  if (sum(valid_rows) == 0) {
    cat("ERROR: No valid data points for plotting after removing NA/infinite values\n")
    return(NULL)
  }
  
  cat("Number of valid rows for plotting:", sum(valid_rows), "\n")
  
  # Subset to valid data
  plot_data <- ancom_results$results[valid_rows, ]
  cat("Plot data LFC range:", range(plot_data$lfc), "\n")
  cat("Plot data q-value range:", range(plot_data$q_val), "\n")
  cat("================================\n")
  
  # Create a simplified taxonomy label using the valid plot data
  plot_data$tax_label <- paste0(
    plot_data$Phylum, ";", 
    plot_data$Class, ";",
    plot_data$Genus
  )
  
  # Create custom significance and direction categories
  plot_data$sig_direction <- ifelse(
    plot_data$q_val >= sig_threshold, 
    "NS",
    ifelse(plot_data$lfc < 0, "More in black", "More in orange")
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
  if (sum(!is.na(plot_data$q_val) & plot_data$q_val < sig_threshold) > 0) {
    significant_results <- plot_data[plot_data$q_val < sig_threshold & !is.na(plot_data$q_val), ]
    top_esvs <- significant_results[order(significant_results$q_val), ][1:min(10, nrow(significant_results)), ]
  }
  
  # 2. Cyanobacteria (only if significant)
  cyanobacteria <- plot_data[
    !is.na(plot_data$Phylum) & 
      plot_data$Phylum == "Cyanobacteria" &
      !is.na(plot_data$q_val) &
      plot_data$q_val < sig_threshold, 
  ]
  
  # 3. Specific genera of interest (only if significant)
  genera_of_interest <- c("Flavobacterium", "Polaromonas", "Cryobacterium")
  specific_genera <- plot_data[
    !is.na(plot_data$Genus) & 
      plot_data$Genus %in% genera_of_interest &
      !is.na(plot_data$q_val) &
      plot_data$q_val < sig_threshold, 
  ]
  
  # Combine all taxa to be labeled (remove duplicates by ESV)
  taxa_to_label <- unique(rbind(
    top_esvs[, names(plot_data)],
    cyanobacteria,
    specific_genera
  ))
  
  # Create labels in "ESV: Taxonomy" format
  # Function to clean taxonomy names by removing "_" and everything after it
  clean_taxonomy <- function(tax_name) {
    # Handle vector input by applying to each element
    sapply(tax_name, function(x) {
      if (is.na(x) || x == "") {
        return(x)
      }
      # Split by "_" and take only the first part
      return(strsplit(x, "_")[[1]][1])
    }, USE.NAMES = FALSE)
  }
  
  # Prioritize Genus, fall back to higher taxonomy if Genus is NA
  # Clean the taxonomy names to remove suffixes
  taxa_to_label$label <- paste0(
    taxa_to_label$ESV, ": ",
    ifelse(
      !is.na(taxa_to_label$Genus) & taxa_to_label$Genus != "",
      clean_taxonomy(taxa_to_label$Genus),
      ifelse(
        !is.na(taxa_to_label$Family) & taxa_to_label$Family != "",
        clean_taxonomy(taxa_to_label$Family),
        ifelse(
          !is.na(taxa_to_label$Order) & taxa_to_label$Order != "",
          clean_taxonomy(taxa_to_label$Order),
          clean_taxonomy(taxa_to_label$Class)
        )
      )
    )
  )
  
  # Create the volcano plot using the cleaned plot_data
  p <- ggplot(plot_data, aes(x = lfc, y = -log10(q_val), color = sig_direction)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = -log10(sig_threshold), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.3) +
    scale_color_manual(values = color_palette) +
    theme_bw() +
    labs(
      title = title,
      x = "Log2 Fold Difference",
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
  n_sig_pos <- sum(plot_data$sig_direction == "More in orange", na.rm = TRUE)
  n_sig_neg <- sum(plot_data$sig_direction == "More in black", na.rm = TRUE)
  
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

# # Generate taxa abundance plot for significant taxa
# taxa_plot_time0 <- plot_top_taxa(
#   ancom_results = ancom_time0,
#   physeq_obj = ps_time0_filtered,
#   group_var = "OrderGroup",
#   n_taxa = 10,
#   sig_threshold = 0.05
# )
# if (!is.null(taxa_plot_time0)) {
#   print(taxa_plot_time0)
# }

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

# # Generate taxa abundance plot for significant taxa
# taxa_plot_time1 <- plot_top_taxa(
#   ancom_results = ancom_time1,
#   physeq_obj = ps_time1_filtered,
#   group_var = "OrderGroup",
#   n_taxa = 10,
#   sig_threshold = 0.05
# )
# if (!is.null(taxa_plot_time1)) {
#   print(taxa_plot_time1)
# }

# Write significant results to file
sig_results_time1 <- write_sig_results(
  ancom_results = ancom_time1,
  filename = "significant_taxa_OB1_vs_BO1.csv",
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



# Write significant results to file
sig_results_time2 <- write_sig_results(
  ancom_results = ancom_time2,
  filename = "significant_taxa_OB2_vs_BO2.csv",
  sig_threshold = 0.05
)


# Combine volcano plots
combined_volcano <- volcano_plot_time0 / volcano_plot_time1 / volcano_plot_time2 +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# Save the combined plot
ggsave("ancombc_volcano_plots.pdf", combined_volcano, width = 6, height = 12)
ggsave("ancombc_volcano_plots.png", combined_volcano, width = 6, height = 12)





######### -------- Figure S3: Relative abundance of ESV 18 (Nostoc) ------------

# Transform to relative abundances
physeq_rel <- transform_sample_counts(physeq, function(x) x / sum(x))

# Extract ESV_18 abundance data
esv18_data <- data.frame(
  Sample = sample_names(physeq_rel),
  ESV_18_RelAbund = as.numeric(otu_table(physeq_rel)["ESV_18", ]),
  OrderTiming = sample_data(physeq_rel)$OrderTiming,
  Order = sample_data(physeq_rel)$Order,
  Timing = sample_data(physeq_rel)$Timing
)

# Check if ESV_18 exists in the dataset
if(!"ESV_18" %in% taxa_names(physeq)) {
  warning("ESV_18 not found in the dataset. Available ESVs starting with 'ESV_18':")
  print(taxa_names(physeq)[grep("^ESV_18", taxa_names(physeq))])
  stop("Please check the ESV name")
}

# Create the boxplot
esv18_boxplot <- ggplot(esv18_data, aes(x = OrderTiming, y = ESV_18_RelAbund, fill = Order)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1.5) +
  labs(
    title = "Relative Abundance of ESV_18 across Order-Timing Groups",
    x = "Order.Timing",
    y = "Relative Abundance",
    fill = "Order"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  guides(fill = guide_legend(title = "Order"))

# Display the plot
print(esv18_boxplot)

# Optional: Save the plot
ggsave("Fig_S3_ESV_18_RelativeAbundance_Boxplot.png", plot = esv18_boxplot, 
       width = 12, height = 8, dpi = 300, bg = "white")



