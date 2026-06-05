# =============================================================================
# Alpha Diversity Analysis: Richness & Evenness
# Factors: Order, Timing, Order x Timing interaction
# =============================================================================
setwd("/home/pacifica/R/antarctica/S03/big_experiment/")

library(tidyverse)
library(vegan)       # for diversity metrics
library(emmeans)     # for pairwise contrasts
library(ggpubr)      # for stat_compare_means / publication-ready plots
library(patchwork)   # for combining plots

# Define input file paths
count_table_path <- "ESV_table_PFX_2025.csv"
metadata_path <- "Metadata_PFX_2025.csv"

# Read in the data
count_data <- read.csv(count_table_path, row.names = 1)
metadata <- read.csv(metadata_path, row.names = 1)


# -----------------------------------------------------------------------------
# 1. PREPARE DATA
# -----------------------------------------------------------------------------

# --- Transpose count_data so samples are rows, ESVs are columns ---
# Assumes: count_data rows = ESVs, columns = sample IDs
counts_t <- as.data.frame(t(count_data))  # samples x ESVs

# --- Compute alpha diversity metrics per sample ---
alpha_df <- data.frame(
  Sample         = rownames(counts_t),
  Richness       = specnumber(counts_t),              # observed ESV richness
  Shannon        = diversity(counts_t, index = "shannon"),
  Simpson        = diversity(counts_t, index = "simpson"),
  Pielou_J       = diversity(counts_t, index = "shannon") / log(specnumber(counts_t))
  # Pielou's J = Shannon / ln(S) — evenness metric; NaN when Richness = 0
)

# --- Join with metadata ---
alpha_meta <- alpha_df %>%
  left_join(
    metadata %>% select(Sample, Order, Timing),
    by = "Sample"
  ) %>%
  mutate(
    Order  = factor(Order, levels = c("BB", "BO", "AL", "OB", "OO", "NO")),
    # Timing is numeric in the CSV (0 / 1 / 2); recode explicitly before factoring
    # so that 0 is never silently dropped to NA by factor(levels = c("1","2"))
    Timing = case_when(
      Timing == 0 ~ "Pre-experiment",
      Timing == 1 ~ "One Season",
      Timing == 2 ~ "Two Seasons",
      TRUE        ~ NA_character_
    ),
    Timing = factor(Timing, levels = c("Pre-experiment", "One Season", "Two Seasons")),
    # Remove samples with zero richness (no valid evenness)
    Pielou_J = ifelse(is.nan(Pielou_J) | is.infinite(Pielou_J), NA, Pielou_J)
  )

# Quick check
glimpse(alpha_meta)
cat("\nSample counts per group:\n")
print(count(alpha_meta, Order, Timing))

# -----------------------------------------------------------------------------
# 2. STATISTICAL MODELS
# -----------------------------------------------------------------------------

# --- Linear models for each metric ---
# Richness
lm_rich  <- lm(Richness  ~ Order * Timing, data = alpha_meta)
# Shannon
lm_shan  <- lm(Shannon   ~ Order * Timing, data = alpha_meta)
# Pielou's J (evenness)
lm_even  <- lm(Pielou_J  ~ Order * Timing, data = alpha_meta)

cat("\n========== RICHNESS: ANOVA ==========\n")
print(anova(lm_rich))

cat("\n========== SHANNON: ANOVA ==========\n")
print(anova(lm_shan))

cat("\n========== EVENNESS (Pielou J): ANOVA ==========\n")
print(anova(lm_even))

# --- Estimated marginal means & pairwise contrasts ---

# Order main effect
cat("\n--- Richness: Order contrasts ---\n")
emm_rich_order <- emmeans(lm_rich, ~ Order)
print(pairs(emm_rich_order, adjust = "BH"))

cat("\n--- Shannon: Order contrasts ---\n")
emm_shan_order <- emmeans(lm_shan, ~ Order)
print(pairs(emm_shan_order, adjust = "BH"))

cat("\n--- Evenness: Order contrasts ---\n")
emm_even_order <- emmeans(lm_even, ~ Order)
print(pairs(emm_even_order, adjust = "BH"))

# Timing main effect
cat("\n--- Richness: Timing contrasts ---\n")
emm_rich_timing <- emmeans(lm_rich, ~ Timing)
print(pairs(emm_rich_timing, adjust = "BH"))

# Interaction contrasts (Order within each Timing level)
cat("\n--- Richness: Order | Timing interaction contrasts ---\n")
emm_rich_int <- emmeans(lm_rich, ~ Order | Timing)
print(pairs(emm_rich_int, adjust = "BH"))

cat("\n--- Shannon: Order | Timing interaction contrasts ---\n")
emm_shan_int <- emmeans(lm_shan, ~ Order | Timing)
print(pairs(emm_shan_int, adjust = "BH"))

cat("\n--- Evenness: Order | Timing interaction contrasts ---\n")
emm_even_int <- emmeans(lm_even, ~ Order | Timing)
print(pairs(emm_even_int, adjust = "BH"))

# -----------------------------------------------------------------------------
# 3. VISUALIZATIONS
# -----------------------------------------------------------------------------

# Shared theme
theme_alpha <- theme_bw(base_size = 12) +
  theme(
    strip.background  = element_rect(fill = "grey92"),
    legend.position   = "bottom",
    axis.text.x       = element_text(angle = 30, hjust = 1)
  )

# Colour palette for Order factor
# Named vector ordered BB, BO, AL, OB, OO, NO to match factor levels —
# ggplot2 draws the legend in the order the vector is defined, so this
# must match the factor levels set above for the legend to be correct.
order_colors <- c(BB = "#1A1A1A", BO = "#808080", AL = "#8B4513",
                  OB = "#DAA520", OO = "#E06C00", NO = "#0055A4")

# ----- 3a. Boxplots: Richness by Order, faceted by Timing -----
p_rich_box <- ggplot(alpha_meta, aes(x = Order, y = Richness, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey30") +
  facet_wrap(~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Observed Richness", x = "Order (Treatment)", y = "No. of ESVs") +
  theme_alpha

# ----- 3b. Boxplots: Shannon by Order, faceted by Timing -----
p_shan_box <- ggplot(alpha_meta, aes(x = Order, y = Shannon, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey30") +
  facet_wrap(~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Shannon Diversity", x = "Order (Treatment)", y = "Shannon H'") +
  theme_alpha

# ----- 3c. Boxplots: Pielou's J by Order, faceted by Timing -----
p_even_box <- ggplot(alpha_meta, aes(x = Order, y = Pielou_J, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, na.rm = TRUE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey30", na.rm = TRUE) +
  facet_wrap(~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Evenness (Pielou's J)", x = "Order (Treatment)", y = "Pielou's J") +
  theme_alpha

# ----- 3d. Interaction plot: mean ± SE for Richness -----
alpha_summary <- alpha_meta %>%
  group_by(Order, Timing) %>%
  summarise(
    mean_rich = mean(Richness, na.rm = TRUE),
    se_rich   = sd(Richness,   na.rm = TRUE) / sqrt(n()),
    mean_shan = mean(Shannon,  na.rm = TRUE),
    se_shan   = sd(Shannon,    na.rm = TRUE) / sqrt(n()),
    mean_even = mean(Pielou_J, na.rm = TRUE),
    se_even   = sd(Pielou_J,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_rich_int <- ggplot(alpha_summary,
                     aes(x = Timing, y = mean_rich, colour = Order, group = Order)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_rich - se_rich, ymax = mean_rich + se_rich),
                width = 0.1, linewidth = 0.7) +
  scale_colour_manual(values = order_colors) +
  labs(title = "Richness: Order × Timing Interaction",
       x = "Timing", y = "Mean Richness ± SE") +
  theme_alpha

p_shan_int <- ggplot(alpha_summary,
                     aes(x = Timing, y = mean_shan, colour = Order, group = Order)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_shan - se_shan, ymax = mean_shan + se_shan),
                width = 0.1, linewidth = 0.7) +
  scale_colour_manual(values = order_colors) +
  labs(title = "Shannon: Order × Timing Interaction",
       x = "Timing", y = "Mean Shannon H' ± SE") +
  theme_alpha

p_even_int <- ggplot(alpha_summary,
                     aes(x = Timing, y = mean_even, colour = Order, group = Order)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_even - se_even, ymax = mean_even + se_even),
                width = 0.1, linewidth = 0.7) +
  scale_colour_manual(values = order_colors) +
  labs(title = "Evenness: Order × Timing Interaction",
       x = "Timing", y = "Mean Pielou's J ± SE") +
  theme_alpha

# ----- 3e. Combined figure -----
combined_plot <- (p_rich_box | p_shan_box | p_even_box) /
  (p_rich_int | p_shan_int | p_even_int) +
  plot_annotation(
    title    = "Alpha Diversity by Treatment Order and Timing",
    subtitle = "Top row: distributions per group  |  Bottom row: Order × Timing interaction",
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

print(combined_plot)

# Save if desired
 ggsave("alpha_diversity_combined.pdf", combined_plot, width = 16, height = 10)

# -----------------------------------------------------------------------------
# 4. OPTIONAL: Kruskal-Wallis + Dunn post-hoc (non-parametric alternative)
# -----------------------------------------------------------------------------
# Use if residuals are clearly non-normal (check with shapiro.test on residuals)

cat("\n========== NORMALITY CHECK (Shapiro-Wilk on model residuals) ==========\n")
cat("Richness residuals:  p =", shapiro.test(residuals(lm_rich))$p.value, "\n")
cat("Shannon residuals:   p =", shapiro.test(residuals(lm_shan))$p.value, "\n")
cat("Evenness residuals:  p =", shapiro.test(residuals(lm_even))$p.value, "\n")
cat("(p < 0.05 suggests non-normality; consider non-parametric tests below)\n")

# Kruskal-Wallis for Order effect on Richness
cat("\n--- Kruskal-Wallis: Richness ~ Order ---\n")
print(kruskal.test(Richness ~ Order, data = alpha_meta))

cat("\n--- Kruskal-Wallis: Richness ~ Timing ---\n")
print(kruskal.test(Richness ~ Timing, data = alpha_meta))

# Dunn post-hoc (requires dunn.test or rstatix)
if (requireNamespace("rstatix", quietly = TRUE)) {
  cat("\n--- Dunn post-hoc: Richness ~ Order (BH-adjusted) ---\n")
  dunn_rich <- rstatix::dunn_test(alpha_meta, Richness ~ Order, p.adjust.method = "BH")
  print(dunn_rich)
  cat("\n--- Dunn post-hoc: Shannon ~ Order (BH-adjusted) ---\n")
  dunn_shan <- rstatix::dunn_test(alpha_meta, Shannon ~ Order, p.adjust.method = "BH")
  print(dunn_shan)
  cat("\n--- Dunn post-hoc: Evenness (Pielou J) ~ Order (BH-adjusted) ---\n")
  dunn_even <- rstatix::dunn_test(alpha_meta, Pielou_J ~ Order, p.adjust.method = "BH")
  print(dunn_even)
} else {
  message("Install 'rstatix' for Dunn post-hoc: install.packages('rstatix')")
  dunn_rich <- dunn_shan <- dunn_even <- NULL
}

# -----------------------------------------------------------------------------
# 5. EXPORT STATISTICAL RESULTS TO CSV
# -----------------------------------------------------------------------------

# Helper: convert anova() output to a tidy data frame
anova_to_df <- function(aov_obj, metric) {
  df <- as.data.frame(aov_obj)
  df$Term   <- rownames(df)
  df$Metric <- metric
  df$Test   <- "ANOVA"
  rownames(df) <- NULL
  df[, c("Metric", "Test", "Term",
         "Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")]
}

# Helper: convert pairs() / emmeans contrast output to a tidy data frame
emmeans_to_df <- function(pairs_obj, metric, contrast_label) {
  df <- as.data.frame(summary(pairs_obj))
  df$Metric   <- metric
  df$Test     <- "emmeans_pairs"
  df$Contrast_group <- contrast_label
  df
}

# Helper: convert kruskal.test output to a one-row data frame
kruskal_to_df <- function(kt, metric, predictor) {
  data.frame(
    Metric      = metric,
    Test        = "Kruskal-Wallis",
    Predictor   = predictor,
    statistic   = kt$statistic,
    df          = kt$parameter,
    p.value     = kt$p.value,
    stringsAsFactors = FALSE
  )
}

# 5a. ANOVA tables
anova_results <- bind_rows(
  anova_to_df(anova(lm_rich), "Richness"),
  anova_to_df(anova(lm_shan), "Shannon"),
  anova_to_df(anova(lm_even), "Pielou_J")
)
write.csv(anova_results, "stats_anova.csv", row.names = FALSE)
cat("\nWrote stats_anova.csv\n")

# 5b. emmeans pairwise contrasts
emmeans_results <- bind_rows(
  emmeans_to_df(pairs(emm_rich_order,  adjust = "BH"), "Richness", "Order"),
  emmeans_to_df(pairs(emm_shan_order,  adjust = "BH"), "Shannon",  "Order"),
  emmeans_to_df(pairs(emm_even_order,  adjust = "BH"), "Pielou_J", "Order"),
  emmeans_to_df(pairs(emm_rich_timing, adjust = "BH"), "Richness", "Timing"),
  emmeans_to_df(pairs(emm_rich_int,    adjust = "BH"), "Richness", "Order|Timing"),
  emmeans_to_df(pairs(emm_shan_int,    adjust = "BH"), "Shannon",  "Order|Timing"),
  emmeans_to_df(pairs(emm_even_int,    adjust = "BH"), "Pielou_J", "Order|Timing")
)
write.csv(emmeans_results, "stats_emmeans_contrasts.csv", row.names = FALSE)
cat("Wrote stats_emmeans_contrasts.csv\n")

# 5c. Kruskal-Wallis results
kw_results <- bind_rows(
  kruskal_to_df(kruskal.test(Richness ~ Order,  data = alpha_meta), "Richness", "Order"),
  kruskal_to_df(kruskal.test(Richness ~ Timing, data = alpha_meta), "Richness", "Timing"),
  kruskal_to_df(kruskal.test(Shannon  ~ Order,  data = alpha_meta), "Shannon",  "Order"),
  kruskal_to_df(kruskal.test(Pielou_J ~ Order,  data = alpha_meta), "Pielou_J", "Order")
)
write.csv(kw_results, "stats_kruskal_wallis.csv", row.names = FALSE)
cat("Wrote stats_kruskal_wallis.csv\n")

# 5d. Dunn post-hoc results (if rstatix was available)
if (!is.null(dunn_rich)) {
  dunn_results <- bind_rows(
    dunn_rich %>% mutate(Metric = "Richness"),
    dunn_shan %>% mutate(Metric = "Shannon"),
    dunn_even %>% mutate(Metric = "Pielou_J")
  )
  write.csv(dunn_results, "stats_dunn_posthoc.csv", row.names = FALSE)
  cat("Wrote stats_dunn_posthoc.csv\n")
} else {
  cat("Dunn results not available (rstatix not installed); skipping stats_dunn_posthoc.csv\n")
}

# 5e. Shapiro-Wilk normality check
shapiro_results <- data.frame(
  Metric  = c("Richness", "Shannon", "Pielou_J"),
  Test    = "Shapiro-Wilk",
  W       = c(shapiro.test(residuals(lm_rich))$statistic,
              shapiro.test(residuals(lm_shan))$statistic,
              shapiro.test(residuals(lm_even))$statistic),
  p.value = c(shapiro.test(residuals(lm_rich))$p.value,
              shapiro.test(residuals(lm_shan))$p.value,
              shapiro.test(residuals(lm_even))$p.value)
)
write.csv(shapiro_results, "stats_shapiro_normality.csv", row.names = FALSE)
cat("Wrote stats_shapiro_normality.csv\n")






# =============================================================================
# 6. RAREFACTION ANALYSIS
# Rarefy count_data to even depth, repeat all analyses, compare to unrarefied
# =============================================================================

# -----------------------------------------------------------------------------
# 6a. CHOOSE RAREFACTION DEPTH & RAREFY
# -----------------------------------------------------------------------------

# Inspect the per-sample read depth distribution to choose a threshold.
# Common practice: rarefy to the minimum depth, or a depth that retains
# most samples while excluding clear outlier-low samples.
sample_depths <- rowSums(counts_t)
cat("\n========== READ DEPTH SUMMARY (unrarefied) ==========\n")
print(summary(sample_depths))
cat("\nPer-sample depths (sorted):\n")
print(sort(sample_depths))

# Set rarefaction depth. Default = minimum sample depth (most conservative).
# Change rare_depth manually if the minimum is an extreme outlier.
rare_depth <- 10000
cat("\nRarefying to:", rare_depth, "reads per sample\n")

# Samples below rare_depth are automatically dropped by rrarefy / we filter them.
samples_kept   <- names(sample_depths[sample_depths >= rare_depth])
samples_dropped <- names(sample_depths[sample_depths < rare_depth])
if (length(samples_dropped) > 0) {
  cat("Samples dropped (below rarefaction depth):", paste(samples_dropped, collapse = ", "), "\n")
} else {
  cat("All samples retained after rarefaction.\n")
}

# Rarefy using vegan::rrarefy (random subsampling without replacement).
# Set a seed for reproducibility.
set.seed(42)
counts_rare <- as.data.frame(
  rrarefy(counts_t[samples_kept, ], sample = rare_depth)
)

cat("Rarefied matrix dimensions:", nrow(counts_rare), "samples x",
    ncol(counts_rare), "ESVs\n")

# -----------------------------------------------------------------------------
# 6b. COMPUTE ALPHA DIVERSITY ON RAREFIED COUNTS
# -----------------------------------------------------------------------------

alpha_rare_df <- data.frame(
  Sample   = rownames(counts_rare),
  Richness = specnumber(counts_rare),
  Shannon  = diversity(counts_rare, index = "shannon"),
  Simpson  = diversity(counts_rare, index = "simpson"),
  Pielou_J = diversity(counts_rare, index = "shannon") / log(specnumber(counts_rare))
)

# Join with metadata
alpha_rare_meta <- alpha_rare_df %>%
  left_join(
    metadata %>% select(Sample, Order, Timing),
    by = "Sample"
  ) %>%
  mutate(
    Order  = factor(Order, levels = c("BB", "BO", "AL", "OB", "OO", "NO")),
    Timing = case_when(
      Timing == 0 ~ "Pre-experiment",
      Timing == 1 ~ "One Season",
      Timing == 2 ~ "Two Seasons",
      TRUE        ~ NA_character_
    ),
    Timing   = factor(Timing, levels = c("Pre-experiment", "One Season", "Two Seasons")),
    Pielou_J = ifelse(is.nan(Pielou_J) | is.infinite(Pielou_J), NA, Pielou_J)
  )

cat("\nSample counts per group (rarefied):\n")
print(count(alpha_rare_meta, Order, Timing))

# -----------------------------------------------------------------------------
# 6c. STATISTICAL MODELS ON RAREFIED DATA
# -----------------------------------------------------------------------------

lm_rare_rich <- lm(Richness ~ Order * Timing, data = alpha_rare_meta)
lm_rare_shan <- lm(Shannon  ~ Order * Timing, data = alpha_rare_meta)
lm_rare_even <- lm(Pielou_J ~ Order * Timing, data = alpha_rare_meta)

cat("\n========== RAREFIED — RICHNESS: ANOVA ==========\n")
print(anova(lm_rare_rich))

cat("\n========== RAREFIED — SHANNON: ANOVA ==========\n")
print(anova(lm_rare_shan))

cat("\n========== RAREFIED — EVENNESS (Pielou J): ANOVA ==========\n")
print(anova(lm_rare_even))

# emmeans contrasts — Order main effect
cat("\n--- Rarefied Richness: Order contrasts ---\n")
emm_rare_rich_order <- emmeans(lm_rare_rich, ~ Order)
print(pairs(emm_rare_rich_order, adjust = "BH"))

cat("\n--- Rarefied Shannon: Order contrasts ---\n")
emm_rare_shan_order <- emmeans(lm_rare_shan, ~ Order)
print(pairs(emm_rare_shan_order, adjust = "BH"))

cat("\n--- Rarefied Evenness: Order contrasts ---\n")
emm_rare_even_order <- emmeans(lm_rare_even, ~ Order)
print(pairs(emm_rare_even_order, adjust = "BH"))

# Timing main effect
cat("\n--- Rarefied Richness: Timing contrasts ---\n")
emm_rare_rich_timing <- emmeans(lm_rare_rich, ~ Timing)
print(pairs(emm_rare_rich_timing, adjust = "BH"))

# Interaction contrasts
cat("\n--- Rarefied Richness: Order | Timing interaction contrasts ---\n")
emm_rare_rich_int <- emmeans(lm_rare_rich, ~ Order | Timing)
print(pairs(emm_rare_rich_int, adjust = "BH"))

cat("\n--- Rarefied Shannon: Order | Timing interaction contrasts ---\n")
emm_rare_shan_int <- emmeans(lm_rare_shan, ~ Order | Timing)
print(pairs(emm_rare_shan_int, adjust = "BH"))

cat("\n--- Rarefied Evenness: Order | Timing interaction contrasts ---\n")
emm_rare_even_int <- emmeans(lm_rare_even, ~ Order | Timing)
print(pairs(emm_rare_even_int, adjust = "BH"))

# Normality check on rarefied residuals
cat("\n========== RAREFIED — NORMALITY CHECK (Shapiro-Wilk) ==========\n")
cat("Richness residuals:  p =", shapiro.test(residuals(lm_rare_rich))$p.value, "\n")
cat("Shannon residuals:   p =", shapiro.test(residuals(lm_rare_shan))$p.value, "\n")
cat("Evenness residuals:  p =", shapiro.test(residuals(lm_rare_even))$p.value, "\n")

# Kruskal-Wallis
cat("\n--- Rarefied Kruskal-Wallis: Richness ~ Order ---\n")
print(kruskal.test(Richness ~ Order, data = alpha_rare_meta))
cat("\n--- Rarefied Kruskal-Wallis: Richness ~ Timing ---\n")
print(kruskal.test(Richness ~ Timing, data = alpha_rare_meta))

# Dunn post-hoc
if (requireNamespace("rstatix", quietly = TRUE)) {
  cat("\n--- Rarefied Dunn: Richness ~ Order ---\n")
  dunn_rare_rich <- rstatix::dunn_test(alpha_rare_meta, Richness ~ Order, p.adjust.method = "BH")
  print(dunn_rare_rich)
  cat("\n--- Rarefied Dunn: Shannon ~ Order ---\n")
  dunn_rare_shan <- rstatix::dunn_test(alpha_rare_meta, Shannon ~ Order, p.adjust.method = "BH")
  print(dunn_rare_shan)
  cat("\n--- Rarefied Dunn: Evenness ~ Order ---\n")
  dunn_rare_even <- rstatix::dunn_test(alpha_rare_meta, Pielou_J ~ Order, p.adjust.method = "BH")
  print(dunn_rare_even)
} else {
  dunn_rare_rich <- dunn_rare_shan <- dunn_rare_even <- NULL
}

# -----------------------------------------------------------------------------
# 6d. VISUALIZATIONS ON RAREFIED DATA
# -----------------------------------------------------------------------------

rare_subtitle <- paste0("Rarefied to ", rare_depth, " reads/sample")

# Boxplots
p_rare_rich_box <- ggplot(alpha_rare_meta, aes(x = Order, y = Richness, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey30") +
  facet_wrap(~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Observed Richness (Rarefied)", subtitle = rare_subtitle,
       x = "Order (Treatment)", y = "No. of ESVs") +
  theme_alpha

p_rare_shan_box <- ggplot(alpha_rare_meta, aes(x = Order, y = Shannon, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey30") +
  facet_wrap(~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Shannon Diversity (Rarefied)", subtitle = rare_subtitle,
       x = "Order (Treatment)", y = "Shannon H'") +
  theme_alpha

p_rare_even_box <- ggplot(alpha_rare_meta, aes(x = Order, y = Pielou_J, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, na.rm = TRUE) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, colour = "grey30", na.rm = TRUE) +
  facet_wrap(~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Evenness — Pielou's J (Rarefied)", subtitle = rare_subtitle,
       x = "Order (Treatment)", y = "Pielou's J") +
  theme_alpha

# Interaction plots
alpha_rare_summary <- alpha_rare_meta %>%
  group_by(Order, Timing) %>%
  summarise(
    mean_rich = mean(Richness, na.rm = TRUE),
    se_rich   = sd(Richness,   na.rm = TRUE) / sqrt(n()),
    mean_shan = mean(Shannon,  na.rm = TRUE),
    se_shan   = sd(Shannon,    na.rm = TRUE) / sqrt(n()),
    mean_even = mean(Pielou_J, na.rm = TRUE),
    se_even   = sd(Pielou_J,   na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

p_rare_rich_int <- ggplot(alpha_rare_summary,
                           aes(x = Timing, y = mean_rich, colour = Order, group = Order)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_rich - se_rich, ymax = mean_rich + se_rich),
                width = 0.1, linewidth = 0.7) +
  scale_colour_manual(values = order_colors) +
  labs(title = "Richness: Order × Timing (Rarefied)", subtitle = rare_subtitle,
       x = "Timing", y = "Mean Richness ± SE") +
  theme_alpha

p_rare_shan_int <- ggplot(alpha_rare_summary,
                           aes(x = Timing, y = mean_shan, colour = Order, group = Order)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_shan - se_shan, ymax = mean_shan + se_shan),
                width = 0.1, linewidth = 0.7) +
  scale_colour_manual(values = order_colors) +
  labs(title = "Shannon: Order × Timing (Rarefied)", subtitle = rare_subtitle,
       x = "Timing", y = "Mean Shannon H' ± SE") +
  theme_alpha

p_rare_even_int <- ggplot(alpha_rare_summary,
                           aes(x = Timing, y = mean_even, colour = Order, group = Order)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_even - se_even, ymax = mean_even + se_even),
                width = 0.1, linewidth = 0.7) +
  scale_colour_manual(values = order_colors) +
  labs(title = "Evenness: Order × Timing (Rarefied)", subtitle = rare_subtitle,
       x = "Timing", y = "Mean Pielou's J ± SE") +
  theme_alpha

# Combined rarefied figure
combined_rare_plot <- (p_rare_rich_box | p_rare_shan_box | p_rare_even_box) /
  (p_rare_rich_int | p_rare_shan_int | p_rare_even_int) +
  plot_annotation(
    title    = "Alpha Diversity by Treatment Order and Timing (Rarefied)",
    subtitle = paste0("Top row: distributions per group  |  Bottom row: Order × Timing interaction\n",
                      rare_subtitle),
    theme    = theme(plot.title    = element_text(size = 14, face = "bold"),
                     plot.subtitle = element_text(size = 11, colour = "grey40"))
  )

print(combined_rare_plot)
ggsave("alpha_diversity_rarefied_combined.pdf", combined_rare_plot, width = 16, height = 10)

# -----------------------------------------------------------------------------
# 6e. SIDE-BY-SIDE COMPARISON PLOTS (Unrarefied vs Rarefied)
# -----------------------------------------------------------------------------
# Stack the unrarefied and rarefied boxplots for each metric so differences
# in distribution / group patterns are immediately visible.

# Add a Dataset label to each data frame for faceting
alpha_meta_lab      <- alpha_meta      %>% mutate(Dataset = "Unrarefied")
alpha_rare_meta_lab <- alpha_rare_meta %>% mutate(Dataset = paste0("Rarefied (", rare_depth, " reads)"))
alpha_combined_lab  <- bind_rows(alpha_meta_lab, alpha_rare_meta_lab) %>%
  mutate(Dataset = factor(Dataset, levels = c("Unrarefied",
                                              paste0("Rarefied (", rare_depth, " reads)"))))

p_cmp_rich <- ggplot(alpha_combined_lab, aes(x = Order, y = Richness, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4, colour = "grey30") +
  facet_grid(Dataset ~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Richness: Unrarefied vs Rarefied",
       x = "Order (Treatment)", y = "No. of ESVs") +
  theme_alpha

p_cmp_shan <- ggplot(alpha_combined_lab, aes(x = Order, y = Shannon, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4, colour = "grey30") +
  facet_grid(Dataset ~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Shannon: Unrarefied vs Rarefied",
       x = "Order (Treatment)", y = "Shannon H'") +
  theme_alpha

p_cmp_even <- ggplot(alpha_combined_lab, aes(x = Order, y = Pielou_J, fill = Order)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8, na.rm = TRUE) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.4, colour = "grey30", na.rm = TRUE) +
  facet_grid(Dataset ~ Timing) +
  scale_fill_manual(values = order_colors) +
  labs(title = "Evenness (Pielou's J): Unrarefied vs Rarefied",
       x = "Order (Treatment)", y = "Pielou's J") +
  theme_alpha

combined_cmp_plot <- p_cmp_rich / p_cmp_shan / p_cmp_even +
  plot_annotation(
    title    = "Alpha Diversity: Unrarefied vs Rarefied Comparison",
    theme    = theme(plot.title = element_text(size = 14, face = "bold"))
  )

print(combined_cmp_plot)
ggsave("alpha_diversity_comparison_unrarefied_vs_rarefied.pdf",
       combined_cmp_plot, width = 16, height = 18)

# -----------------------------------------------------------------------------
# 6f. EXPORT RAREFIED STATISTICAL RESULTS TO CSV
# -----------------------------------------------------------------------------

# ANOVA
anova_rare_results <- bind_rows(
  anova_to_df(anova(lm_rare_rich), "Richness"),
  anova_to_df(anova(lm_rare_shan), "Shannon"),
  anova_to_df(anova(lm_rare_even), "Pielou_J")
)
write.csv(anova_rare_results, "stats_rarefied_anova.csv", row.names = FALSE)
cat("\nWrote stats_rarefied_anova.csv\n")

# emmeans contrasts
emmeans_rare_results <- bind_rows(
  emmeans_to_df(pairs(emm_rare_rich_order,  adjust = "BH"), "Richness", "Order"),
  emmeans_to_df(pairs(emm_rare_shan_order,  adjust = "BH"), "Shannon",  "Order"),
  emmeans_to_df(pairs(emm_rare_even_order,  adjust = "BH"), "Pielou_J", "Order"),
  emmeans_to_df(pairs(emm_rare_rich_timing, adjust = "BH"), "Richness", "Timing"),
  emmeans_to_df(pairs(emm_rare_rich_int,    adjust = "BH"), "Richness", "Order|Timing"),
  emmeans_to_df(pairs(emm_rare_shan_int,    adjust = "BH"), "Shannon",  "Order|Timing"),
  emmeans_to_df(pairs(emm_rare_even_int,    adjust = "BH"), "Pielou_J", "Order|Timing")
)
write.csv(emmeans_rare_results, "stats_rarefied_emmeans_contrasts.csv", row.names = FALSE)
cat("Wrote stats_rarefied_emmeans_contrasts.csv\n")

# Kruskal-Wallis
kw_rare_results <- bind_rows(
  kruskal_to_df(kruskal.test(Richness ~ Order,  data = alpha_rare_meta), "Richness", "Order"),
  kruskal_to_df(kruskal.test(Richness ~ Timing, data = alpha_rare_meta), "Richness", "Timing"),
  kruskal_to_df(kruskal.test(Shannon  ~ Order,  data = alpha_rare_meta), "Shannon",  "Order"),
  kruskal_to_df(kruskal.test(Pielou_J ~ Order,  data = alpha_rare_meta), "Pielou_J", "Order")
)
write.csv(kw_rare_results, "stats_rarefied_kruskal_wallis.csv", row.names = FALSE)
cat("Wrote stats_rarefied_kruskal_wallis.csv\n")

# Dunn post-hoc
if (!is.null(dunn_rare_rich)) {
  dunn_rare_results <- bind_rows(
    dunn_rare_rich %>% mutate(Metric = "Richness"),
    dunn_rare_shan %>% mutate(Metric = "Shannon"),
    dunn_rare_even %>% mutate(Metric = "Pielou_J")
  )
  write.csv(dunn_rare_results, "stats_rarefied_dunn_posthoc.csv", row.names = FALSE)
  cat("Wrote stats_rarefied_dunn_posthoc.csv\n")
} else {
  cat("Rarefied Dunn results not available (rstatix not installed).\n")
}

# Shapiro-Wilk
shapiro_rare_results <- data.frame(
  Metric  = c("Richness", "Shannon", "Pielou_J"),
  Test    = "Shapiro-Wilk",
  W       = c(shapiro.test(residuals(lm_rare_rich))$statistic,
              shapiro.test(residuals(lm_rare_shan))$statistic,
              shapiro.test(residuals(lm_rare_even))$statistic),
  p.value = c(shapiro.test(residuals(lm_rare_rich))$p.value,
              shapiro.test(residuals(lm_rare_shan))$p.value,
              shapiro.test(residuals(lm_rare_even))$p.value)
)
write.csv(shapiro_rare_results, "stats_rarefied_shapiro_normality.csv", row.names = FALSE)
cat("Wrote stats_rarefied_shapiro_normality.csv\n")

# -----------------------------------------------------------------------------
# 6g. PROGRAMMATIC COMPARISON SUMMARY
# Compare ANOVA p-values and significant emmeans contrasts between the two
# datasets and print a human-readable summary + write a comparison CSV.
# -----------------------------------------------------------------------------

cat("\n")
cat("=============================================================\n")
cat("  COMPARISON SUMMARY: UNRAREFIED vs RAREFIED\n")
cat("=============================================================\n")

# --- Helper: extract ANOVA p-values into a named vector ---
anova_pvals <- function(lm_obj) {
  av <- as.data.frame(anova(lm_obj))
  setNames(av[["Pr(>F)"]], rownames(av))
}

# Collect p-values for all terms x metrics x dataset
terms <- c("Order", "Timing", "Order:Timing")

compare_anova <- function(metric_label, lm_unrar, lm_rar) {
  p_u <- anova_pvals(lm_unrar)
  p_r <- anova_pvals(lm_rar)
  purrr::map_dfr(terms, function(trm) {
    pu <- if (trm %in% names(p_u)) p_u[[trm]] else NA_real_
    pr <- if (trm %in% names(p_r)) p_r[[trm]] else NA_real_
    sig_u  <- !is.na(pu) & pu < 0.05
    sig_r  <- !is.na(pr) & pr < 0.05
    status <- dplyr::case_when(
      sig_u  & sig_r  ~ "Significant in both",
      sig_u  & !sig_r ~ "Lost significance after rarefaction",
      !sig_u & sig_r  ~ "Gained significance after rarefaction",
      TRUE            ~ "Non-significant in both"
    )
    data.frame(Metric = metric_label, Term = trm,
               p_unrarefied = pu, p_rarefied = pr,
               Status = status, stringsAsFactors = FALSE)
  })
}

anova_comparison <- bind_rows(
  compare_anova("Richness", lm_rich,      lm_rare_rich),
  compare_anova("Shannon",  lm_shan,      lm_rare_shan),
  compare_anova("Pielou_J", lm_even,      lm_rare_even)
)

cat("\n--- ANOVA term significance: Unrarefied vs Rarefied ---\n")
print(anova_comparison, row.names = FALSE)

# --- Compare significant emmeans contrasts ---
sig_contrasts <- function(pairs_obj, dataset_label, metric, contrast_group) {
  df <- as.data.frame(summary(pairs_obj))
  p_col <- intersect(c("p.value", "p.value.BH", "p.adj"), names(df))[1]
  df %>%
    filter(.data[[p_col]] < 0.05) %>%
    mutate(Dataset = dataset_label, Metric = metric, Contrast_group = contrast_group) %>%
    select(Dataset, Metric, Contrast_group, contrast, estimate,
           any_of(c("SE", "df", "t.ratio")), !!p_col)
}

sig_emm_unrar <- bind_rows(
  sig_contrasts(pairs(emm_rich_order,  adjust = "BH"), "Unrarefied", "Richness", "Order"),
  sig_contrasts(pairs(emm_shan_order,  adjust = "BH"), "Unrarefied", "Shannon",  "Order"),
  sig_contrasts(pairs(emm_even_order,  adjust = "BH"), "Unrarefied", "Pielou_J", "Order"),
  sig_contrasts(pairs(emm_rich_timing, adjust = "BH"), "Unrarefied", "Richness", "Timing"),
  sig_contrasts(pairs(emm_rich_int,    adjust = "BH"), "Unrarefied", "Richness", "Order|Timing"),
  sig_contrasts(pairs(emm_shan_int,    adjust = "BH"), "Unrarefied", "Shannon",  "Order|Timing"),
  sig_contrasts(pairs(emm_even_int,    adjust = "BH"), "Unrarefied", "Pielou_J", "Order|Timing")
)

sig_emm_rar <- bind_rows(
  sig_contrasts(pairs(emm_rare_rich_order,  adjust = "BH"), "Rarefied", "Richness", "Order"),
  sig_contrasts(pairs(emm_rare_shan_order,  adjust = "BH"), "Rarefied", "Shannon",  "Order"),
  sig_contrasts(pairs(emm_rare_even_order,  adjust = "BH"), "Rarefied", "Pielou_J", "Order"),
  sig_contrasts(pairs(emm_rare_rich_timing, adjust = "BH"), "Rarefied", "Richness", "Timing"),
  sig_contrasts(pairs(emm_rare_rich_int,    adjust = "BH"), "Rarefied", "Richness", "Order|Timing"),
  sig_contrasts(pairs(emm_rare_shan_int,    adjust = "BH"), "Rarefied", "Shannon",  "Order|Timing"),
  sig_contrasts(pairs(emm_rare_even_int,    adjust = "BH"), "Rarefied", "Pielou_J", "Order|Timing")
)

# Identify contrasts that appear in only one dataset (changed significance)
key_u <- sig_emm_unrar %>% mutate(.key = paste(Metric, Contrast_group, contrast))
key_r <- sig_emm_rar   %>% mutate(.key = paste(Metric, Contrast_group, contrast))

lost_sig    <- key_u %>% filter(!(.key %in% key_r$.key)) %>% select(-.key)
gained_sig  <- key_r %>% filter(!(.key %in% key_u$.key)) %>% select(-.key)
kept_sig    <- key_u %>% filter(  .key %in% key_r$.key)  %>% select(-.key)

cat("\n--- emmeans contrasts significant in BOTH datasets ---\n")
if (nrow(kept_sig) > 0) print(kept_sig, row.names = FALSE) else cat("  (none)\n")

cat("\n--- emmeans contrasts that LOST significance after rarefaction ---\n")
if (nrow(lost_sig) > 0) print(lost_sig, row.names = FALSE) else cat("  (none)\n")

cat("\n--- emmeans contrasts that GAINED significance after rarefaction ---\n")
if (nrow(gained_sig) > 0) print(gained_sig, row.names = FALSE) else cat("  (none)\n")

# --- Metric-level mean comparison ---
metric_means <- bind_rows(
  alpha_meta      %>% summarise(across(c(Richness, Shannon, Pielou_J), mean, na.rm = TRUE)) %>%
    mutate(Dataset = "Unrarefied"),
  alpha_rare_meta %>% summarise(across(c(Richness, Shannon, Pielou_J), mean, na.rm = TRUE)) %>%
    mutate(Dataset = "Rarefied")
) %>% select(Dataset, everything())

cat("\n--- Overall mean alpha diversity: Unrarefied vs Rarefied ---\n")
print(metric_means, row.names = FALSE)

# --- Write comparison CSVs ---
write.csv(anova_comparison,
          "stats_comparison_anova_unrarefied_vs_rarefied.csv", row.names = FALSE)
cat("\nWrote stats_comparison_anova_unrarefied_vs_rarefied.csv\n")

write.csv(bind_rows(
    lost_sig   %>% mutate(Change = "Lost significance"),
    gained_sig %>% mutate(Change = "Gained significance"),
    kept_sig   %>% mutate(Change = "Remained significant")
  ),
  "stats_comparison_emmeans_unrarefied_vs_rarefied.csv", row.names = FALSE)
cat("Wrote stats_comparison_emmeans_unrarefied_vs_rarefied.csv\n")

write.csv(metric_means, "stats_comparison_mean_diversity_unrarefied_vs_rarefied.csv",
          row.names = FALSE)
cat("Wrote stats_comparison_mean_diversity_unrarefied_vs_rarefied.csv\n")

cat("\n=============================================================\n")
cat("  END OF RAREFACTION ANALYSIS\n")
cat("=============================================================\n")


# Create Table S1: 
sample_depths_df <- tibble(
  Sample = names(sample_depths_pcoa),
  SampleDepth = as.numeric(sample_depths_pcoa)
)

metadata_merged <- metadata %>%
  left_join(sample_depths_df, by = "Sample") 

metadata_fortable <- metadata_merged %>% 
  dplyr::select(MatOrder=Order,
         Season=Timing,
         SampleID=Sample,
         ReadCount=SampleDepth,
         Treatment_description)
write_csv(metadata_fortable,"Table_S1_ReadDepths.csv")
