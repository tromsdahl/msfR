# MSF R Data Analysis for Metabolomics v2.0 (updated 07/18/2025)
# ---------------------------------------------------------------
# Project: Project Name (PI Name -- Institution)
# Date: 07/18/2025
# Experiment: Data Analysis of samples -- Selected group example

# Obj: Determine the LM response in plasma

# Required columns in data file:
# sample_no
# sample_name
# sample_type
# rep_no
# component_name
# component_group <- needed for lipidomics
# is
# area
# rt <- needed for data quality check
# prot_conc

# load libraries
library(tidyverse)
library(data.table)
library(gplots)
library(ggrepel)
library(ComplexUpset)
library(rstatix)

# prefix and threshold values and other parameters
prefix <- "20250718_Example_Data_LM_"

l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adj_pval_max <- 0.05

is_normalize <- "LXA4_d5.IS"

group_comparisons <- "group_comparisons.txt"
fc_comparisons <- "FC_comparisons.txt"

# load data
met_data <- setDT(read.delim("met_data.txt"))

# select specific sample types
# Tip: change sample_type %in% c(...)  to !(sample_type %in% c(...)) if you want to exclude samples
met_data <- met_data[sample_type %in% c("Group 27", "Group 30", "Group 33", "Group 36", "Pooled QC", "Ext PQC"), ]

# normalize data to internal standard and sample starting amount (e.g. protein concentration, tissue weight, cell no., etc)
met_data_is <- met_data[component_name == is_normalize, ]
norm_met_data <- met_data[is == "n", ]

norm_met_data[met_data_is, ":="(norm_area = area / (i.area * prot_conc), is_used = i.component_name, is_area = i.area), on = .(sample_no)]

# write final table of normalized data
write.table(norm_met_data, paste0(prefix,"norm_data.txt"), quote = F, row.names = F, sep = "\t")


# Create Metaboanalyst-compatible table
norm_met_data_metabo <- dcast(norm_met_data[is == "n" & (is.finite(norm_area))], paste0(sample_type,"_",rep_no) + sample_type ~ component_name, value.var = "norm_area")
norm_met_data_metabo <- setnames(norm_met_data_metabo, old = c("sample_type", "sample_type_1"), new = c("sample_name","sample_type"))

# remove blanks and IS samples
norm_met_data_metabo <- norm_met_data_metabo[!(grep("Blank",sample_type)), ]

# remove columns where all are NA
norm_met_data_metabo <- norm_met_data_metabo[, .SD, .SDcols = function(x){!(all(is.na(x)))}]

write.table(norm_met_data_metabo, paste0(prefix,"norm_data_metabo.txt"), quote = F, row.names = F, sep = "\t")


#============================================================================
# STATISTICAL WORKFLOW
#============================================================================
met_data_stats <- melt(norm_met_data_metabo, id.vars = 1:2, variable.name = "component_name", value.name = "norm_area")
met_data_stats[norm_area == 0, norm_area := NA] # Skyline sometimes sets values to zero, but these should be NA

# remove lipids/metabolites with missing values
# Function to remove lipids/metabolites where 2/3 of replicates of one group are missing within a paired comparison
comparisons <- read.delim(group_comparisons)

filtering_table <- vector("list", length(comparisons[,1]) + 1) # +1 for the pooled QC
for (i in seq_along(comparisons[,1])) {
  
  # select data for groups within comparison
  group_data <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3])]
  
  # count how many lipids/metabolites per sample group are missing, and how many samples there are per group
  group_1_ct <- length(unique(group_data[sample_type == comparisons[i,2]]$sample_name)) # number of replicates in group 1
  group_2_ct <- length(unique(group_data[sample_type == comparisons[i,3]]$sample_name)) # number of replicates in group 2
  missing_ct <- group_data[is.na(norm_area), .N, by = .(sample_type, component_name)] # count number of missing values per metabolite/lipid
  met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_1_ct) | N >= ((2/3) * group_2_ct)]$component_name)) # identify which lipids/metabolites are missing in at least 2/3 of either sample group
  
  # remove the selected lipids/metabolites for the particular groups from the met_data_stats table and add to filtering table list
  filtering_table[[i]] <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3]) & !(component_name %in% met_remove), ]
  
}

# filter out 2/3 of missing lipids/metabolites from pooled QC
group_data <- met_data_stats[!(sample_type %in% c(comparisons$group_1, comparisons$group_2))]

group_ct <- length(unique(group_data$sample_name)) # number of replicates in (presumably) pooled QC
missing_ct <- group_data[is.na(norm_area), .N, by = .(sample_type, component_name)]
met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_ct)]$component_name))

filtering_table[[length(comparisons[,1]) + 1]] <- met_data_stats[!(sample_type %in% c(comparisons$group_1, comparisons$group_2)) & !(component_name %in% met_remove)]

# Remake met_data_stats table to continue analysis
met_data_stats <- rbindlist(filtering_table)

met_data_stats <- unique(met_data_stats)

met_data_subquant <- dcast(met_data_stats, component_name ~ sample_name, value.var = "norm_area")


#============================================================================
# CORRELATION HEATMAP
#============================================================================
met_data_corr <- sapply(met_data_subquant[, 2:ncol(met_data_subquant)], log10)
row.names(met_data_corr) <- met_data_subquant$component_name

cormat <- round(cor(met_data_corr, use = "complete.obs"), 2)

colorpalette <- colorRampPalette(c("blue","yellow","red"))

png(paste0(prefix, "sample_correlation_heatmap.png"), res = 300, units = "in", width = 16, height = 16)
heatmap.2(x = cormat, col = colorpalette(100), trace = "none", density.info = "none", margins = c(10,10), lhei = c(1,8), lwid = c(2,8))
dev.off()


#============================================================================
# METABOLITE/LIPID COUNT
#============================================================================
# count number of metabolites/lipids per sample
met_count <- met_data_stats[!(is.na(norm_area)), .N, by = .(sample_name, sample_type)]

# calculate mean and SD of metabolite counts
met_count_avgSD <- met_count[, .(avg_ct = round(mean(N), 1), sd_ct = round(sd(N), 1)), by = .(sample_type)]

# plot as bar graph with error bars
ggplot(met_count_avgSD, aes(x = sample_type, y = avg_ct)) +
  geom_col() +
  geom_errorbar(aes(ymin = avg_ct - sd_ct, ymax = avg_ct + sd_ct), width = 0.5) +
  geom_text(aes(label = avg_ct, y = avg_ct + sd_ct), vjust = -1, size = 3) +
  labs(x = "Sample Type", y = "Avg. Num. Metabolites Quantified") +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, max(met_count_avgSD$avg_ct) + 20), expand = c(0,0))
ggsave(paste0(prefix, "met_count.png"), device = "png", dpi = 300, units = "in", width = 10, height = 6)


#============================================================================
# INTENSITIES OF METABOLITES (NORMALIZED) PER SAMPLE
#============================================================================
subquant_medians <- met_data_stats[, .(med_log_10_int = round(log10(median(norm_area, na.rm = T)), 2)), by = .(sample_name, sample_type)]

# Plot normalized metabolite intensities
ggplot(met_data_stats, aes(x = sample_name, y = log10(norm_area))) +
  geom_boxplot(aes(fill = sample_type)) +
  geom_text(data = subquant_medians, aes(x = sample_name, y = med_log_10_int, label = med_log_10_int), vjust = -0.5, size = 2) +
  labs(x = NULL, y = "Log10NormalizedIntensity", fill = "Sample Group") +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    legend.position = "bottom"
  )
ggsave(paste0(prefix, "sample_intensities.png"), device = "png", dpi = 300, units = "in", width = 16, height = 6)


#============================================================================
# IMPUTATION AND LOG10 TRANSFORMATION
#============================================================================
# imputation by 1/5 of minimum value of metabolite across all samples
imputate <- met_data_stats[!(grep("QC", sample_type)), .(imput = min(norm_area, na.rm = T) / 5), by = .(component_name)]
met_data_stats[imputate, imput := i.imput, on = .(component_name)]
met_data_stats[is.na(norm_area), norm_area := imput]

# log10 transform
met_data_stats[, log10norm_area := log10(norm_area)]


#============================================================================
# PCA ANALYSIS
#============================================================================
# create transposed matrix for PCA
met_data_pca_mat <- na.omit(dcast(met_data_stats, component_name ~ sample_name, value.var = "log10norm_area"))
pca_mat_names <- met_data_pca_mat$component_name
met_data_pca_mat <- as.matrix(met_data_pca_mat[, 2:ncol(met_data_pca_mat)])
row.names(met_data_pca_mat) <- pca_mat_names
met_data_pca_mat <- t(met_data_pca_mat)

# calculate principle component analysis
pca <- prcomp(na.omit(met_data_pca_mat), center = T, scale = T)

# create table with values from PCA analysis for plotting
pcaplot <- data.table(sample_name = row.names(pca$x[,1:2]), PC1 = pca$x[,1], PC2 = pca$x[,2])
pcaplot[, sample_type := str_extract(sample_name, "[^_]+")]

# Plot PC1 vs PC2 and save to png
ggplot(pcaplot, aes(x = PC1, y = PC2, fill = sample_type)) +
  geom_point(color = "black", pch = 21, size = 4) +
  geom_text_repel(aes(label = sample_name), max.overlaps = 50) +
  labs(x = paste0("PC 1 (",round(summary(pca)[["importance"]][2,1] * 100,1),"%)"), y = paste0("PC 2 (",round(summary(pca)[["importance"]][2,2] * 100,1),"%)"), fill = "Sample Type") +
  theme(aspect.ratio = 1)
ggsave(paste0(prefix,"PCA_plot.png"), device = "png", dpi = 300, units = "in", width = 8, height = 8)

# Plot PC1 vs PC2 and save to png (without labels)
ggplot(pcaplot, aes(x = PC1, y = PC2, fill = sample_type)) +
  geom_point(color = "black", pch = 21, size = 4) +
  #geom_text_repel(aes(label = sample_name), max.overlaps = 50) +
  labs(x = paste0("PC 1 (",round(summary(pca)[["importance"]][2,1] * 100,1),"%)"), y = paste0("PC 2 (",round(summary(pca)[["importance"]][2,2] * 100,1),"%)"), fill = "Sample Type") +
  theme(aspect.ratio = 1)
ggsave(paste0(prefix,"PCA_plot_ul.png"), device = "png", dpi = 300, units = "in", width = 8, height = 8)


#============================================================================
# LOADINGS BIPLOT FOR PCA ANALYSIS
#============================================================================
# Create table with loadings data
loading_table <- data.table(component_name = row.names(pca$rotation[,1:2]),
                            PC1 = pca$rotation[,1],
                            PC2 = pca$rotation[,2])

# multiplication factor to scale loadings to PCA for biplot
mult <- min(
  (max(pcaplot$PC2) - min(pcaplot$PC2) / (max(loading_table$PC2) - min(loading_table$PC2))),
  (max(pcaplot$PC1) - min(pcaplot$PC1) / (max(loading_table$PC1) - min(loading_table$PC1)))
)

# apply scaling factor
loading_table[, ":="(PC1 = 0.8 * mult * PC1,
                     PC2 = 0.8 * mult * PC2)]

# create biplot
ggplot(pcaplot, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = sample_type), color = "black", pch = 21, size = 4) +
  labs(x = paste0("PC 1 (",round(summary(pca)[["importance"]][2,1] * 100,1),"%)"), y = paste0("PC 2 (",round(summary(pca)[["importance"]][2,2] * 100,1),"%)"), fill = "Sample Type") +
  geom_segment(data = loading_table, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(angle = 10,, length = unit(0.1, "inches"), type = "closed")) +
  geom_text(data = loading_table, aes(label = component_name), size = 3, color = "orange") +
  theme(aspect.ratio = 1)
ggsave(paste0(prefix,"PCA_biplot.png"), device = "png", dpi = 300, units = "in", width = 8, height = 8)


#============================================================================
# GROUP COMPARISON STATISTICS
#============================================================================
# Function to calculate p, adj p, and log2FC for each comparison within the comparisons table
# vectors to store values
comparison_table <- vector("list", length(comparisons[,1]))
names(comparison_table) <- comparisons[,1]

# Function to calculate p values, adj p values, and log2FC
for (i in seq_along(comparisons[,1])) {
  # to monitor progress
  print(paste0("Now analyzing ", comparisons[i,1], "..."))
  
  # autoscaling samples within comparison against each other
  autoscale <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3]), .(mean_val = mean(log10norm_area, na.rm = T), sd_val = sd(log10norm_area, na.rm = T)), by = .(component_name)]
  met_data_stats[autoscale, ":="(mean_val = i.mean_val, sd_val = i.sd_val), on = .(component_name)]
  met_data_stats[, trans_area := (log10norm_area - mean_val) / sd_val]
  
  # filter out components that might be missing from one sample group
  n_replicates <- length(unique(met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3])]$sample_name))
  n_components <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3]), .N, by = .(component_name)]
  components <- as.character(n_components[N == n_replicates]$component_name)
  
  # calculate p and p adj values
  p_table <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3]) & component_name %in% components] %>% 
    group_by(component_name) %>% 
    t_test(trans_area ~ sample_type, var.equal = T) %>% 
    adjust_pvalue(method = "BH") %>% 
    add_significance()
  
  # calculate log2FC
  group_1 <- met_data_stats[sample_type == comparisons[i,2] & component_name %in% components, .(avg_norm_area = mean(norm_area)), by = .(component_name)]
  group_2 <- met_data_stats[sample_type == comparisons[i,3] & component_name %in% components, .(avg_norm_area = mean(norm_area)), by = .(component_name)]
  
  log2FC <- group_1[group_2, .(component_name, fold_change = avg_norm_area / i.avg_norm_area, log2FC = log2(avg_norm_area / i.avg_norm_area)), on = .(component_name)]
  
  comparison_table[[i]] <- list(p_table,log2FC)
  
  # to monitor progress
  print(paste0(comparisons[i,1], " analysis complete!"))
}


unpack <- vector("list",length(comparison_table))
# unpack above list into readable table
for (i in seq_along(comparison_table)) {
  component_name <- comparison_table[[i]][[1]][[1]]
  comparison <- rep(names(comparison_table)[[i]],length(component_name))
  log2FC <- comparison_table[[i]][[2]][["log2FC"]]
  fold_change <- comparison_table[[i]][[2]][["fold_change"]]
  t_value <- comparison_table[[i]][[1]][["statistic"]]
  df <- comparison_table[[i]][[1]][["df"]]
  p_value <- comparison_table[[i]][[1]][["p"]]
  p_adj <- comparison_table[[i]][[1]][["p.adj"]]
  
  results <- data.table(component_name, comparison, log2FC, fold_change, t_value, df, p_value, p_adj)
  
  unpack[[i]] <- results
  
}

results <- rbindlist(unpack)

# write results to file
write.table(results, paste0(prefix,"results.txt"), quote = F, row.names = F, sep = "\t")


#============================================================================
# VOLCANO PLOTS (p < 0.05)
#============================================================================
# using unadjusted p values
results$volcano <- "NO CHANGE"

results[log2FC > l2fc_max & p_value < pval_max, volcano := "UP"]
results[log2FC < l2fc_min & p_value < pval_max, volcano := "DOWN"]

cols <- c("UP" = "red", "NO CHANGE" = "gray", "DOWN" = "blue")

# plot volcano plots
ggplot(results, aes(x = log2FC, y = -log10(p_value))) +
  geom_hline(yintercept = -log10(pval_max), color = "gray", linetype = 2) +
  geom_vline(xintercept = c(l2fc_min, l2fc_max), color = "gray", linetype = 2) +
  geom_point(aes(color = volcano), size = 0.25) +
  #geom_text_repel(data = results[volcano %in% c("UP","DOWN")], aes(label = component_name), size = 2, max.overlaps = 50) +
  scale_color_manual(values = cols) +
  facet_wrap(comparison ~ ., ncol = 4) +
  labs(x = "Log2 Fold Change", y = "-Log10(p)", color = NULL)
ggsave(paste0(prefix, "volcano_unadj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# UPSET PLOT (p < 0.05)
#============================================================================
# using unadjusted p values
# change volcano values to logical
results$volcano <- results$volcano %in% c("UP","DOWN")

# dcast to wide format
results_wide <- dcast(results, component_name ~ comparison, value.var = "volcano")

# remove metabolite names
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all FALSE
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))
ggsave(paste0(prefix, "upset_plot_unadj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# 2nd VOLCANO PLOT (adj p < 0.05)
#============================================================================
# using adj p value
results$volcano <- NULL
results$volcano <- "NO CHANGE"

results[log2FC > l2fc_max & p_adj < adj_pval_max, volcano := "UP"]
results[log2FC < l2fc_min & p_adj < adj_pval_max, volcano := "DOWN"]

cols <- c("UP" = "red", "NO CHANGE" = "gray", "DOWN" = "blue")

# plot volcano plots
ggplot(results, aes(x = log2FC, y = -log10(p_value))) +
  geom_hline(yintercept = -log10(pval_max), color = "gray", linetype = 2) +
  geom_vline(xintercept = c(l2fc_min, l2fc_max), color = "gray", linetype = 2) +
  geom_point(aes(color = volcano), size = 0.25) +
  #geom_text_repel(data = results[volcano %in% c("UP","DOWN")], aes(label = component_name), size = 2, max.overlaps = 50) +
  scale_color_manual(values = cols) +
  facet_wrap(comparison ~ ., ncol = 4) +
  labs(x = "Log2 Fold Change", y = "-Log10(adj p)", color = NULL)
ggsave(paste0(prefix, "volcano_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# UPSET PLOT (adj p < 0.05)
#============================================================================
# using adj p values
# change volcano values to logical
results$volcano <- results$volcano %in% c("UP","DOWN")

# dcast to wide format
results_wide <- dcast(results, component_name ~ comparison, value.var = "volcano")

# remove metabolite names
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all FALSE
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))
ggsave(paste0(prefix, "upset_plot_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# SIGNIFICANTLY DIFFERENT METABOLITES (p < 0.05)
#============================================================================
# separate metabolites with p < 0.05
results_sig <- results[p_value < pval_max, ]

# count number of significant metabolites per group comparison
results_sig_ct <- results_sig[, .N, by = .(comparison)]

results_sig_ct[, position := N + 1]

# plot number of significant metabolites per group
ggplot(results_sig_ct, aes(x = comparison, y = N)) +
  geom_col(color = "black", fill = "gray") +
  geom_text(aes(y = position, label = N), vjust = 1) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  labs(x = "Comparison", y = "Sig. Metabolites (p < 0.05)") +
  scale_y_continuous(limits = c(0, max(results_sig_ct$N + 10)), expand = c(0,0))
ggsave(paste0(prefix, "signif_met.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# UPSET PLOT OF SIGNIFICANT METABOLITES (p < 0.05)
#============================================================================
# which metabolites are significant (p < 0.05)
results$significant <- results$p_value < pval_max

# cast to wide
results_wide <- dcast(results, component_name ~ comparison, value.var = "significant")

# remove metabolites name
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all false
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))
ggsave(paste0(prefix, "upset_signif_met.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# SIGNIFICANTLY DIFFERENT METABOLITES (adj p < 0.05)
#============================================================================
# separate metabolites with adj p < 0.05
results_sig <- results[p_adj < adj_pval_max, ]

# count number of significant metabolites per group comparison
results_sig_ct <- results_sig[, .N, by = .(comparison)]

results_sig_ct[, position := N + 1]

# plot number of significant metabolites per group
ggplot(results_sig_ct, aes(x = comparison, y = N)) +
  geom_col(color = "black", fill = "gray") +
  geom_text(aes(y = position, label = N), vjust = 1) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  labs(x = "Comparison", y = "Sig. Metabolites (adj p < 0.05)") +
  scale_y_continuous(limits = c(0,max(results_sig_ct$N + 10)), expand = c(0,0))
ggsave(paste0(prefix, "signif_met_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# UPSET PLOT OF SIGNIFICANT METABOLITES (adj p < 0.05)
#============================================================================
# which metabolites are significant (adj p < 0.05)
results$significant <- results$p_adj < adj_pval_max

# cast to wide
results_wide <- dcast(results, component_name ~ comparison, value.var = "significant")

# remove metabolites name
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all false
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))
ggsave(paste0(prefix, "upset_signif_met_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# DIFFERENTIALLY ABUNDANT METABOLITES (p < 0.05)
#============================================================================
# using unadjusted p values
# upregulated or downregulated?
results[log2FC > l2fc_max & p_value < pval_max, direction := "Increased"][log2FC < l2fc_min & p_value < pval_max, direction := "Decreased"]

# count number of up/down regulated metabolites using table
upregulatedcounts <- results[direction == "Increased", .N, by = .(comparison)]
downregulatedcounts <- results[direction == "Decreased", .N, by = .(comparison)]

# label for change direction
upregulatedcounts$direction <- "Increased"
downregulatedcounts$direction <- "Decreased"
significantcounts <- rbind(upregulatedcounts, downregulatedcounts)

# flip sign of decreased metabolites
significantcounts[direction == "Decreased", N := N * -1]

# add position for plotting label
significantcounts$position <- 0
significantcounts[N < 0, position := N - 1][N > 0, position := N + 1]

# set colors for bar plots
cols <- c("Increased" = "red", "Decreased" = "blue")

# plot number of increased vs. decreased metabolites
ggplot(significantcounts, aes(x = comparison, y = N)) +
  geom_bar(stat = "identity", aes(fill = factor(direction, levels = c("Increased", "Decreased")))) +
  scale_fill_manual(values = cols, name = "Direction") +
  geom_text(aes(y = position, label = abs(N))) +
  geom_hline(yintercept = 0, color = "black") +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(x = "Comparison", y = "Num. of Sig. Changes")
ggsave(paste0(prefix, "diff_abund_unadj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# DIFFERENTIALLY ABUNDANT METABOLITES (adj p < 0.05)
#============================================================================
# using adj p values
results[, direction := NULL][log2FC > l2fc_max & p_adj < adj_pval_max, direction := "Increased"][log2FC < l2fc_min & p_adj < adj_pval_max, direction := "Decreased"]

# count number of up/down regulated metabolites using table
upregulatedcounts <- results[direction == "Increased", .N, by = .(comparison)]
downregulatedcounts <- results[direction == "Decreased", .N, by = .(comparison)]

# label for change direction
upregulatedcounts$direction <- "Increased"
downregulatedcounts$direction <- "Decreased"
significantcounts <- rbind(upregulatedcounts, downregulatedcounts)

# flip sign of decreased metabolites
significantcounts[direction == "Decreased", N := N * -1]

# add position for plotting label
significantcounts$position <- 0
significantcounts[N < 0, position := N - 1][N > 0, position := N + 1]

# set colors for bar plots
cols <- c("Increased" = "red", "Decreased" = "blue")

# plot number of increased vs. decreased metabolites
ggplot(significantcounts, aes(x = comparison, y = N)) +
  geom_bar(stat = "identity", aes(fill = factor(direction, levels = c("Increased", "Decreased")))) +
  scale_fill_manual(values = cols, name = "Direction") +
  geom_text(aes(y = position, label = abs(N))) +
  geom_hline(yintercept = 0, color = "black") +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(x = "Comparison", y = "Num. of Sig. Changes")
ggsave(paste0(prefix, "diff_abund_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)

#============================================================================
# HIERARCHICAL CLUSTERING HEATMAP
#============================================================================
# plotting heatmap for each comparison
# NOTE: This assumes equal number of replicates per sample group being compared
hm_cols <- colorRampPalette(c("blue","yellow","red"))(50)

for (i in seq_along(comparisons[,1])) {
  
  # autoscaling samples within comparison against each other
  autoscale <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3]), .(mean_val = mean(log10norm_area, na.rm = T), sd_val = sd(log10norm_area, na.rm = T)), by = .(component_name)]
  met_data_stats[autoscale, ":="(mean_val = i.mean_val, sd_val = i.sd_val), on = .(component_name)]
  met_data_stats[, trans_area := (log10norm_area - mean_val) / sd_val]
  
  met_select <- results[comparison == comparisons[i,1]]
  met_select <- met_select[order(met_select$p_value), ]
  heatmap_comp <- head(met_select, 25)
  heatmap_comp <- dcast(met_data_stats[component_name %in% heatmap_comp$component_name & sample_type %in% c(comparisons[i,2], comparisons[i,3])], component_name ~ sample_name, value.var = "trans_area")
  hm_row <- as.character(heatmap_comp$component_name)
  heatmap_comp <- as.matrix(heatmap_comp[, 2:ncol(heatmap_comp)])
  row.names(heatmap_comp) <- hm_row
  
  png(paste0(prefix,comparisons[i,1],"_heatmap.png"), width = 7, height = 7, units = "in", res = 300)
  heatmap.2(heatmap_comp, Colv = NA, dendrogram = "row", xlab = "Samples", density.info = "none", trace = "none", col = hm_cols,
            hclustfun = function(d){hclust(d, method = "ward.D2")}, key.title = "Transformed Area", key.xlab = NA,
            cexCol = 0.7, cexRow = 0.7, margins = c(8,10), lhei = c(1,7))#, ColSideColors = c(rep("orange",ncol(heatmap_comp)/2),rep("purple",ncol(heatmap_comp)/2)))
  dev.off()
  
}


#============================================================================
# PCA PER GROUP COMPARISON
#============================================================================
pca_plots <- vector("list", length(comparisons[,1]))

for (i in seq_along(comparisons[,1])) {
  # autoscaling samples within comparison against each other
  autoscale <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3]), .(mean_val = mean(log10norm_area, na.rm = T), sd_val = sd(log10norm_area, na.rm = T)), by = .(component_name)]
  met_data_stats[autoscale, ":="(mean_val = i.mean_val, sd_val = i.sd_val), on = .(component_name)]
  met_data_stats[, trans_area := (log10norm_area - mean_val) / sd_val]
  
  pca_table <- met_data_stats[sample_type %in% c(comparisons[i,2], comparisons[i,3])]
  pca_table <- na.omit(dcast(pca_table, component_name ~ sample_name, value.var = "trans_area"))
  pca_table_names <- pca_table$component_name
  pca_matrix <- as.matrix(pca_table[,2:ncol(pca_table)])
  row.names(pca_matrix) <- pca_table_names
  pca_matrix <- t(pca_matrix)
  
  pca <- prcomp(na.omit(pca_matrix), center = T, scale = T)
  
  pcaplot <- data.table(sample_name = row.names(pca$x[,1:2]), PC1 = pca$x[,1], PC2 = pca$x[,2])
  pcaplot[, sample_type := str_extract(sample_name, "[^_]+")]
  
  pca_plots[[i]] <- ggplot(pcaplot, aes(x = PC1, y = PC2, fill = sample_type)) +
    geom_point(color = "black", pch = 21, size = 4) +
    geom_text_repel(aes(label = sample_name), max.overlaps = 50) +
    labs(x = paste0("PC 1 (",round(summary(pca)[["importance"]][2,1] * 100,1),"%)"), y = paste0("PC 2 (",round(summary(pca)[["importance"]][2,2] * 100,1),"%)"), fill = "Sample Type") +
    theme(aspect.ratio = 1)
  
  print(pca_plots[i])
  
  ggsave(paste0(prefix,comparisons[i,1],"_pca.png"), device = "png", dpi = 300, units = "in", width = 8, height = 8)
}


#============================================================================
# QC PLOT OF CVs
#============================================================================
cv_data <- met_data_stats[, .(cv = (sd(norm_area, na.rm = T) / mean(norm_area, na.rm = T) * 100)), by = .(sample_type, component_name)]
cv_med <- cv_data[, .(med_cv = median(cv)), by = .(sample_type)]

ggplot(cv_data, aes(x = cv, fill = sample_type)) +
  geom_histogram(binwidth = 5, color = "black") +
  geom_vline(data = cv_med, aes(xintercept = med_cv), linetype = 2, color = "black") +
  geom_text(data = cv_med, aes(label = paste0("Median = ",format(round(med_cv, 2), nsmall = 2),"%"), x = med_cv + 60, y = length(unique(cv_data$component_name)) / 5)) +
  labs(x = "RSD", y = "Num. of Metabolites", fill = "Sample Group") +
  facet_wrap(vars(sample_type), ncol = 4)
ggsave(paste0(prefix, "cv_plot.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# QC HEATMAP OF MISSINGNESS PER SAMPLE
#============================================================================
met_data_missing <- sapply(met_data_subquant[, 2:ncol(met_data_subquant)], is.na)
row.names(met_data_missing) <- met_data_subquant$component_name
met_data_missing[met_data_missing == TRUE] <- 1

met_data_missing <- met_data_missing[which(rowSums(met_data_missing, na.rm = T) > 0), ]

png(paste0(prefix, "missingness_heatmap.png"), res = 300, units = "in", width = 16, height = 16)
heatmap.2(x = met_data_missing, col = c("black","white"), trace = "none", density.info = "none", margins = c(10,10), key = F)
legend(x = 0, y = 1, legend = c("Missing value", "Valid value"), fill = c("white", "black"), border = "black", bty = "n", xpd = T, cex = 2)
dev.off()


#============================================================================
# FILTERED METABOLITES
#============================================================================
# This is the same function that would filter out lipids/metabolites missing in 2/3 of a sample group, but in this case selects for those lipids that are missing
# Function to remove lipids/metabolites where 2/3 of replicates of one group are missing within a paired comparison
filtered_met <- melt(norm_met_data_metabo, id.vars = 1:2, variable.name = "component_name", value.name = "norm_area")
filtered_met <- filtered_met[norm_area == 0, norm_area := NA] # Skyline sometimes sets values to zero, but these should be NA

filtering_table <- vector("list", length(comparisons[,1]) + 1) # +1 for the pooled QC
for (i in seq_along(comparisons[,1])) {
  
  # select data for groups within comparison
  group_data <- filtered_met[sample_type %in% c(comparisons[i,2], comparisons[i,3])]
  
  # count how many lipids/metabolites per sample group are missing, and how many samples there are per group
  group_1_ct <- length(unique(group_data[sample_type == comparisons[i,2]]$sample_name)) # number of replicates in group 1
  group_2_ct <- length(unique(group_data[sample_type == comparisons[i,3]]$sample_name)) # number of replicates in group 2
  missing_ct <- group_data[is.na(norm_area), .N, by = .(sample_type, component_name)] # count number of missing values per metabolite/lipid
  met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_1_ct) | N >= ((2/3) * group_2_ct)]$component_name)) # identify which lipids/metabolites are missing in at least 2/3 of either sample group
  
  # select the lipids/metabolites that are missing for the particular groups from the filtered_met table and add to filtering table list
  filtering_table[[i]] <- filtered_met[sample_type %in% c(comparisons[i,2], comparisons[i,3]) & (component_name %in% met_remove), ]
  
}

# select out 2/3 of missing lipids/metabolites from pooled QC
group_data <- filtered_met[!(sample_type %in% c(comparisons$group_1, comparisons$group_2))]

group_ct <- length(unique(group_data$sample_name)) # number of replicates in (presumably) pooled QC
missing_ct <- group_data[is.na(norm_area), .N, by = .(sample_type, component_name)]
met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_ct)]$component_name))

filtering_table[[length(comparisons[,1]) + 1]] <- filtered_met[!(sample_type %in% c(comparisons$group_1, comparisons$group_2)) & (component_name %in% met_remove)]

# Remake filtered_met table to continue analysis
filtered_met <- rbindlist(filtering_table)

filtered_met <- unique(filtered_met)

# write filtered lipid/metabolite results to file
write.table(filtered_met, paste0(prefix,"filtered.txt"), quote = F, row.names = F, sep = "\t")


#=========================================================================
# 2D FC ANALYSIS PLOT
#=========================================================================
# load FC comparisons 
fc_groups <- read.delim(fc_comparisons)

fc_plots <- vector("list", length(fc_groups[,1]))

for (i in seq_along(fc_groups[,1])) {
  # select for lipids/metabolites that are found in each comparison
  fc_met_x <- as.character(unique(met_data_stats[sample_type %in% c(comparisons[comparisons[,1] == fc_groups[i,2], c(2,3)]), ]$component_name))
  fc_met_y <- as.character(unique(met_data_stats[sample_type %in% c(comparisons[comparisons[,1] == fc_groups[i,3], c(2,3)]), ]$component_name))
  
  fc_met <- intersect(fc_met_x, fc_met_y)
  
  fc_met <- met_data_stats[sample_type %in% c(comparisons[comparisons[,1] == fc_groups[i,2],c(2,3)], comparisons[comparisons[,1] == fc_groups[i,3],c(2,3)]) &
                                  component_name %in% fc_met]

  # remove metabolites not found in all four sample groups
  no_of_grps <- length(unique(fc_met$sample_name))
  fc_met_count <- fc_met[, .N, by = .(component_name)]
  fc_met_count_remove <- as.vector(fc_met_count[N < no_of_grps]$component_name)
  
  fc_met <- fc_met[!(component_name %in% fc_lipid_count_remove)]
  
  # one-way ANOVA between the two group comparison groups
  p_anova <- fc_met %>% 
    group_by(component_name) %>% 
    anova_test(log10norm_area ~ sample_type, type = 1) %>% 
    adjust_pvalue(method = "BH") %>% 
    add_significance()
  
  # create table of FC values
  fc_plot_table <- data.table(p_anova)[, .(component_name, p)]
  x_axis <- results[comparison == fc_groups[i,2] & component_name %in% unique(fc_met$component_name), ]
  y_axis <- results[comparison == fc_groups[i,3] & component_name %in% unique(fc_met$component_name), ]
  fc_plot_table[x_axis, x_axis := i.log2FC, on = .(component_name)]
  fc_plot_table[y_axis, y_axis := i.log2FC, on = .(component_name)]
  
  # plot 2D FC
  fc_plots[[i]] <- ggplot(na.omit(fc_plot_table), aes(x = x_axis, y = y_axis, color = -log10(p))) +
    geom_point(size = 3, pch = 21, color = "black") +
    geom_point(size = 2) +
    geom_vline(xintercept = c(-1,1), color = "red", linewidth = 0.5, linetype = 2) +
    geom_hline(yintercept = c(-1,1), color = "red", linewidth = 0.5, linetype = 2) +
    geom_text_repel(aes(label = component_name), color = "black") +
    scale_color_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 1.3) +
    labs(x = fc_groups[i,2], y = fc_groups[i,3], color = "-log10(p)")
  
  print(fc_plots[i])
  
  ggsave(paste0(prefix, "fc_plot_", fc_groups[i,1],".png"), device = "png", dpi = 300, units = "in", width = 9, height = 6)
  
  write.table(fc_plot_table, paste0(prefix,"fc_table_", fc_groups[i,1],".txt"), quote = F, row.names = F, sep = "\t")

}
