# MSF R Data Analysis for Proteomics v2.0 (updated 07/22/2025)
# -------------------------------------------------------------
# Project: Project Name (PI Name -- Institution)
# Date: 7/22/2025
# Experiment: Full Analysis of Proteomics

# Obj: Determine the response of proteins in plasma
# 1. Annotate full proteomics data (pg_matrix) with the experimental annotation file
# 2. Perform statistical analysis for each of the indicated paired groups
# 3. (If applicable) perform GO-SEA analysis


# Indicate the full proteomics pg matrix
prot_pg_matrix <- "report.pg_matrix.tsv"

# Indicate the full proteomics exoeriment annotation file
prot_annot <- "experiment_annotation.txt"

# Indicate the group comparisons file for statistical analysis
comparisons <- "group_comparisons.txt"

# Indicate the prefix to use when saving tables and figures, and indicate the log2 fold change (lg2fc) and p value (pval) thresholds to use in statistical analysis
prefix <- "20250722_Example_Data_prot_"
l2fc_max <- 1
l2fc_min <- -1
pval_max <- 0.05
adj_pval_max <- 0.05

# Indicate the file path for the GO term definitions file
GO_def <- "20210218_GO_definitions.txt"


#============================================================================
# LOADING DATA FROM FILES & ADDING SAMPLE INFORMATION TO DATA
#============================================================================
library(tidyverse)
library(data.table)

# load full proteomics annotations
full_annot <- setDT(read.delim(prot_annot))

# load full proteomics data
full_data <- setDT(read.delim(prot_pg_matrix))

# re-format data tables from wide to long format
full_data <- melt(full_data, 1:5, variable.name = "sample", value.name = "intensity")

# replace : and \ under the sample column of the annotation files with . so that it matches the data tables
full_annot[, sample := gsub(":|\\\\|-",".", sample)][, sample := gsub("..", "X..",sample, fixed = T)] # this second part is needed if the ip address is included

# Add condition and replicate information to the data tables
full_data[full_annot, ":="(condition = i.condition, replicate = i.replicate), on = .(sample)]

# exclude certain samples or groups (if desired or needed) select specific sample types
# Tip: change condition %in% c(...)  to !(condition %in% c(...)) if you want to exclude samples; or comment out this line to use all samples
full_data <- full_data[(condition %in% c("Group 27", "Group 30", "Group 33", "Group 36", "pool-01", "pool-02", "pool-03", "pool-04", "pool-05", "pool-06", "pool-07", "pool-08", "pool-09", "pool-10", "pool-11")), ]


#============================================================================
# STATISTICAL ANALYSIS OF FULL PROTEOMICS
#============================================================================
# load libraries
library(gplots)
library(ggrepel)
library(ComplexUpset)
library(rstatix)

# create Metaboanalyst-compatible table
full_data[, ":="(norm_int = log2(intensity), sample_name = paste0(condition,"_",replicate))]
full_data_metabo <- dcast(full_data, sample_name + condition ~ Protein.Group, value.var = "norm_int")

# remove columns where all are NA
full_data_metabo <- full_data_metabo[, .SD, .SDcols = function(x){!(all(is.na(x)))}]

#write.table(full_data_metabo, paste0(prefix,"norm_data_metabo.txt"), quote = F, row.names = F, sep = "\t") # Uncomment this line if you'd like to output a Metaboanalyst-compatible table

#============================================================================
# STATISTICAL WORKFLOW
#============================================================================
full_data_stats <- melt(full_data_metabo, id.vars = 1:2, variable.name = "component_name", value.name = "norm_area")
full_data_stats[norm_area == 0, norm_area := NA] # 0's should be treated as NA


# REMOVE PROTEINS WHERE 2/3 MISSING
# function to remove proteins where 2/3 of replicates of one group within a paired comparison are missing
comparisons <- read.delim(comparisons)

filtering_table <- vector("list", length(comparisons[,1]) + 1) # +1 for the pooled QC
for (i in seq_along(comparisons[,1])) {
  
  # select data for groups within comparison
  group_data <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3])]
  
  # count how many proteins per sample group are missing, and how many samples there are per group
  group_1_ct <- length(unique(group_data[condition == comparisons[i,2]]$sample_name)) # number of replicates in group 1
  group_2_ct <- length(unique(group_data[condition == comparisons[i,3]]$sample_name)) # number of replicates in group 2
  missing_ct <- group_data[is.na(norm_area), .N, by = .(condition, component_name)] # count number of missing values per protein
  met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_1_ct) | N >= ((2/3) * group_2_ct)]$component_name)) # identify which proteins are missing in at least 2/3 of either sample group
  
  # remove the selected proteins for the particular groups from the full_data_stats table and add to filtering table list
  filtering_table[[i]] <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3]) & !(component_name %in% met_remove), ]
  
}

# filter out 2/3 of missing proteins from pooled QC
group_data <- full_data_stats[!(condition %in% c(comparisons$group_1, comparisons$group_2))]

group_ct <- length(unique(group_data$sample_name)) # number of replicates in (presumably) pooled QC
missing_ct <- group_data[is.na(norm_area), .N, by = .(condition, component_name)]
met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_ct)]$component_name))

filtering_table[[length(comparisons[,1]) + 1]] <- full_data_stats[!(condition %in% c(comparisons$group_1, comparisons$group_2)) & !(component_name %in% met_remove)]

# Remake full_data_stats table to continue analysis
full_data_stats <- rbindlist(filtering_table)

full_data_stats <- unique(full_data_stats)

full_data_subquant <- dcast(full_data_stats, component_name ~ sample_name, value.var = "norm_area")


#============================================================================
# CORRELATION HEATMAP
#============================================================================
full_data_corr <- as.matrix(full_data_subquant[, 2:ncol(full_data_subquant)])
row.names(full_data_corr) <- full_data_subquant$component_name

cormat <- round(cor(full_data_corr, use = "complete.obs"), 4)

colorpalette <- colorRampPalette(c("blue","yellow","red"))

png(paste0(prefix, "sample_correlation_heatmap.png"), res = 300, units = "in", width = 16, height = 16)
heatmap.2(x = cormat, col = colorpalette(100), trace = "none", density.info = "none", margins = c(10,10), lhei = c(1,8), lwid = c(2,8))
dev.off()


#============================================================================
# PROTEIN COUNT
#============================================================================
# count number of proteins per sample
met_count <- full_data_stats[!(is.na(norm_area)), .N, by = .(sample_name, condition)]

# calculate mean and SD of protein counts
met_count_avgSD <- met_count[, .(avg_ct = round(mean(N), 1), sd_ct = round(sd(N), 1)), by = .(condition)]

# plot as bar graph with error bars
ggplot(met_count_avgSD, aes(x = condition, y = avg_ct)) +
  geom_col() +
  geom_errorbar(aes(ymin = avg_ct - sd_ct, ymax = avg_ct + sd_ct), width = 0.5) +
  geom_text(aes(label = avg_ct, y = avg_ct + sd_ct), vjust = 0.5, hjust = 1.2, size = 3, angle = -90) +
  labs(x = "Sample Type", y = "Avg. Num. Proteins") +
  theme(
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, max(met_count_avgSD$avg_ct) + 200), expand = c(0,0))
ggsave(paste0(prefix, "protein_count.png"), device = "png", dpi = 300, units = "in", width = 10, height = 6)


#============================================================================
# INTENSITIES OF PROTEINS (NORMALIZED) PER SAMPLE
#============================================================================
subquant_medians <- full_data_stats[, .(med_log_2_int = round(median(norm_area, na.rm = T), 2)), by = .(sample_name, condition)]

# Plot normalized protein intensities
ggplot(full_data_stats, aes(x = sample_name, y = norm_area)) +
  geom_boxplot(aes(fill = condition)) +
  geom_text(data = subquant_medians, aes(x = sample_name, y = med_log_2_int, label = med_log_2_int), vjust = -0.5, size = 2) +
  labs(x = NULL, y = "Log2NormalizedIntensity", fill = "Sample Group") +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    legend.position = "bottom"
  )
ggsave(paste0(prefix, "sample_intensities.png"), device = "png", dpi = 300, units = "in", width = 16, height = 6)


# IMPUTATION
# imputation by 1/5 of minimum value of proteins across all samples
imputate <- full_data_stats[!(grep("Pool", condition, ignore.case = T)), .(imput = min(norm_area, na.rm = T) / 5), by = .(component_name)]
full_data_stats[imputate, imput := i.imput, on = .(component_name)]
full_data_stats[is.na(norm_area), norm_area := imput]


#============================================================================
# PCA ANALYSIS
#============================================================================
# create transposed matrix for PCA
met_data_pca_mat <- na.omit(dcast(full_data_stats, component_name ~ sample_name, value.var = "norm_area"))
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
# LOADINGS BIPLOT FOR PCA ANALYSIS <- Too many proteins to be that useful; will comment out for now
#============================================================================
# Create table with loadings data
#loading_table <- data.table(component_name = row.names(pca$rotation[,1:2]),
#                            PC1 = pca$rotation[,1],
 #                           PC2 = pca$rotation[,2])

# multiplication factor to scale loadings to PCA for biplot
#mult <- min(
#  (max(pcaplot$PC2) - min(pcaplot$PC2) / (max(loading_table$PC2) - min(loading_table$PC2))),
#  (max(pcaplot$PC1) - min(pcaplot$PC1) / (max(loading_table$PC1) - min(loading_table$PC1)))
#)

# apply scaling factor
#loading_table[, ":="(PC1 = 1.5 * mult * PC1,
#                     PC2 = 1.5 * mult * PC2)]

# create biplot
#ggplot(pcaplot, aes(x = PC1, y = PC2)) +
#  geom_point(aes(fill = sample_type), color = "black", pch = 21, size = 4) +
#  labs(x = paste0("PC 1 (",round(summary(pca)[["importance"]][2,1] * 100,1),"%)"), y = paste0("PC 2 (",round(summary(pca)[["importance"]][2,2] * 100,1),"%)"), fill = "Sample Type") +
#  geom_segment(data = loading_table, aes(x = 0, y = 0, xend = PC1, yend = PC2), arrow = arrow(angle = 10,, length = unit(0.1, "inches"), type = "closed")) +
#  geom_text(data = loading_table, aes(label = component_name), size = 3, color = "orange") +
#  theme(aspect.ratio = 1)
#ggsave(paste0(prefix,"PCA_biplot.png"), device = "png", dpi = 300, units = "in", width = 8, height = 8)


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
  autoscale <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3]), .(mean_val = mean(norm_area, na.rm = T), sd_val = sd(norm_area, na.rm = T)), by = .(component_name)]
  full_data_stats[autoscale, ":="(mean_val = i.mean_val, sd_val = i.sd_val), on = .(component_name)]
  full_data_stats[, trans_area := (norm_area - mean_val) / sd_val]
  
  # filter out components that might be missing from one sample group
  n_replicates <- length(unique(full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3])]$sample_name))
  n_components <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3]), .N, by = .(component_name)]
  components <- as.character(n_components[N == n_replicates]$component_name)
  
  # calculate p and p adj values
  p_table <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3]) & component_name %in% components] %>% 
    group_by(component_name) %>% 
    t_test(norm_area ~ condition, var.equal = T) %>% 
    adjust_pvalue(method = "BH") %>% 
    add_significance()
  
  # calculate log2FC
  group_1 <- full_data_stats[condition == comparisons[i,2] & component_name %in% components, .(avg_norm_area = mean(2^norm_area)), by = .(component_name)]
  group_2 <- full_data_stats[condition == comparisons[i,3] & component_name %in% components, .(avg_norm_area = mean(2^norm_area)), by = .(component_name)]
  
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

# remove protein names
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all FALSE
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))#, min_size = 5) # min_size may need to change; uncomment and delete preceding parentheses to add min_size
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

# remove protein names
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all FALSE
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))#, min_size = 5) # min_size may need to change; uncomment and delete preceding parentheses to add min_size
ggsave(paste0(prefix, "upset_plot_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# SIGNIFICANTLY DIFFERENT PROTEINS (p < 0.05)
#============================================================================
# separate proteins with p < 0.05
results_sig <- results[p_value < pval_max, ]

# count number of significant proteins per group comparison
results_sig_ct <- results_sig[, .N, by = .(comparison)]

results_sig_ct[, position := N + 5]

# plot number of significant proteins per group
ggplot(results_sig_ct, aes(x = comparison, y = N)) +
  geom_col(color = "black", fill = "gray") +
  geom_text(aes(y = position, label = N), vjust = 0) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  labs(x = "Comparison", y = "Sig. Proteins (p < 0.05)") +
  scale_y_continuous(limits = c(0, max(results_sig_ct$N + 100)), expand = c(0,0))
ggsave(paste0(prefix, "signif.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# UPSET PLOT OF SIGNIFICANT PROTEINS (p < 0.05)
#============================================================================
# which proteins are significant (p < 0.05)
results$significant <- results$p_value < pval_max

# cast to wide
results_wide <- dcast(results, component_name ~ comparison, value.var = "significant")

# remove proteins name
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all false
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))#, min_size = 5) # min_size may need to change; uncomment and delete preceding parentheses to add min_size
ggsave(paste0(prefix, "upset_signif.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# SIGNIFICANTLY DIFFERENT PROTEINS (adj p < 0.05)
#============================================================================
# separate proteins with adj p < 0.05
results_sig <- results[p_adj < adj_pval_max, ]

# count number of significant proteins per group comparison
results_sig_ct <- results_sig[, .N, by = .(comparison)]

results_sig_ct[, position := N + 5]

# plot number of significant proteins per group
ggplot(results_sig_ct, aes(x = comparison, y = N)) +
  geom_col(color = "black", fill = "gray") +
  geom_text(aes(y = position, label = N), vjust = 1) +
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
    panel.border = element_rect(fill = NA, color = "black")
  ) +
  labs(x = "Comparison", y = "Sig. Proteins (adj p < 0.05)") +
  scale_y_continuous(limits = c(0,max(results_sig_ct$N + 20)), expand = c(0,0))
ggsave(paste0(prefix, "signif_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# UPSET PLOT OF SIGNIFICANT PROTEINS (adj p < 0.05)
#============================================================================
# which proteins are significant (adj p < 0.05)
results$significant <- results$p_adj < adj_pval_max

# cast to wide
results_wide <- dcast(results, component_name ~ comparison, value.var = "significant")

# remove proteins name
results_wide <- results_wide[, 2:ncol(results_wide)]

# remove rows that are all false
results_wide <- results_wide[which(rowSums(results_wide, na.rm = T) > 0), ]

# plot upset plot
upset(results_wide, intersect = unique(results$comparison))#, min_size = 5) # min_size may need to change; uncomment and delete preceding parentheses to add min_size
ggsave(paste0(prefix, "upset_signif_adj.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# DIFFERENTIALLY ABUNDANT PROTEINS (p < 0.05)
#============================================================================
# using unadjusted p values
# upregulated or downregulated?
results[log2FC > l2fc_max & p_value < pval_max, direction := "Increased"][log2FC < l2fc_min & p_value < pval_max, direction := "Decreased"]

# count number of up/down regulated proteins using table
upregulatedcounts <- results[direction == "Increased", .N, by = .(comparison)]
downregulatedcounts <- results[direction == "Decreased", .N, by = .(comparison)]

# label for change direction
upregulatedcounts$direction <- "Increased"
downregulatedcounts$direction <- "Decreased"
significantcounts <- rbind(upregulatedcounts, downregulatedcounts)

# flip sign of decreased proteins
significantcounts[direction == "Decreased", N := N * -1]

# add position for plotting label
significantcounts$position <- 0
significantcounts[N < 0, position := N - 3][N > 0, position := N + 3]

# set colors for bar plots
cols <- c("Increased" = "red", "Decreased" = "blue")

# plot number of increased vs. decreased proteins
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
# DIFFERENTIALLY ABUNDANT PROTEINS (adj p < 0.05)
#============================================================================
# using adj p values
results[, direction := NULL][log2FC > l2fc_max & p_adj < adj_pval_max, direction := "Increased"][log2FC < l2fc_min & p_adj < adj_pval_max, direction := "Decreased"]

# count number of up/down regulated proteins using table
upregulatedcounts <- results[direction == "Increased", .N, by = .(comparison)]
downregulatedcounts <- results[direction == "Decreased", .N, by = .(comparison)]

# label for change direction
upregulatedcounts$direction <- "Increased"
downregulatedcounts$direction <- "Decreased"
significantcounts <- rbind(upregulatedcounts, downregulatedcounts)

# flip sign of decreased proteins
significantcounts[direction == "Decreased", N := N * -1]

# add position for plotting label
significantcounts$position <- 0
significantcounts[N < 0, position := N - 3][N > 0, position := N + 3]

# set colors for bar plots
cols <- c("Increased" = "red", "Decreased" = "blue")

# plot number of increased vs. decreased proteins
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
  autoscale <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3]), .(mean_val = mean(norm_area, na.rm = T), sd_val = sd(norm_area, na.rm = T)), by = .(component_name)]
  full_data_stats[autoscale, ":="(mean_val = i.mean_val, sd_val = i.sd_val), on = .(component_name)]
  full_data_stats[, trans_area := (norm_area - mean_val) / sd_val]
  
  met_select <- results[comparison == comparisons[i,1]]
  met_select <- met_select[order(met_select$p_value), ]
  heatmap_comp <- head(met_select, 25)
  heatmap_comp <- dcast(full_data_stats[component_name %in% heatmap_comp$component_name & condition %in% c(comparisons[i,2], comparisons[i,3])], component_name ~ sample_name, value.var = "trans_area")
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
  autoscale <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3]), .(mean_val = mean(norm_area, na.rm = T), sd_val = sd(norm_area, na.rm = T)), by = .(component_name)]
  full_data_stats[autoscale, ":="(mean_val = i.mean_val, sd_val = i.sd_val), on = .(component_name)]
  full_data_stats[, trans_area := (norm_area - mean_val) / sd_val]
  
  pca_table <- full_data_stats[condition %in% c(comparisons[i,2], comparisons[i,3])]
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
cv_data <- full_data_stats[, .(cv = (sd(2^norm_area, na.rm = T) / mean(2^norm_area, na.rm = T) * 100)), by = .(condition, component_name)]
cv_med <- cv_data[, .(med_cv = median(cv, na.rm = T)), by = .(condition)]

ggplot(cv_data, aes(x = cv, fill = condition)) +
  geom_histogram(binwidth = 5, color = "black") +
  geom_vline(data = cv_med, aes(xintercept = med_cv), linetype = 2, color = "black") +
  geom_text(data = cv_med, aes(label = paste0("Median = ",format(round(med_cv, 2), nsmall = 2),"%"), x = med_cv + 60, y = length(unique(cv_data$component_name)) / 5)) +
  labs(x = "RSD", y = "Num. of Proteins", fill = "Sample Group") +
  facet_wrap(vars(condition), ncol = 4)
ggsave(paste0(prefix, "cv_plot.png"), device = "png", dpi = 300, units = "in", width = 12, height = 6)


#============================================================================
# QC HEATMAP OF MISSINGNESS PER SAMPLE
#============================================================================
full_data_missing <- sapply(full_data_subquant[, 2:ncol(full_data_subquant)], is.na)
row.names(full_data_missing) <- full_data_subquant$component_name
full_data_missing[full_data_missing == TRUE] <- 1

full_data_missing <- full_data_missing[which(rowSums(full_data_missing, na.rm = T) > 0), ]

png(paste0(prefix, "missingness_heatmap.png"), res = 300, units = "in", width = 16, height = 16)
heatmap.2(x = full_data_missing, col = c("black","white"), trace = "none", density.info = "none", margins = c(10,10), key = F)
legend(x = 0, y = 1, legend = c("Missing value", "Valid value"), fill = c("white", "black"), border = "black", bty = "n", xpd = T, cex = 2)
dev.off()


#============================================================================
# FILTERED PROTEINS
#============================================================================
# This is the same function that would filter out proteins missing in 2/3 of a sample group, but in this case selects for those proteins that are missing
# Function to remove proteins where 2/3 of replicates of one group are missing within a paired comparison
filtered_proteins <- melt(full_data_metabo, id.vars = 1:2, variable.name = "component_name", value.name = "norm_area")
filtered_proteins <- filtered_proteins[norm_area == 0, norm_area := NA] # Skyline sometimes sets values to zero, but these should be NA

filtering_table <- vector("list", length(comparisons[,1]) + 1) # +1 for the pooled QC
for (i in seq_along(comparisons[,1])) {
  
  # select data for groups within comparison
  group_data <- filtered_proteins[condition %in% c(comparisons[i,2], comparisons[i,3])]
  
  # count how many proteins per sample group are missing, and how many samples there are per group
  group_1_ct <- length(unique(group_data[condition == comparisons[i,2]]$sample_name)) # number of replicates in group 1
  group_2_ct <- length(unique(group_data[condition == comparisons[i,3]]$sample_name)) # number of replicates in group 2
  missing_ct <- group_data[is.na(norm_area), .N, by = .(condition, component_name)] # count number of missing values per protein
  met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_1_ct) | N >= ((2/3) * group_2_ct)]$component_name)) # identify which proteins are missing in at least 2/3 of either sample group
  
  # select the proteins that are missing for the particular groups from the filtered_proteins table and add to filtering table list
  filtering_table[[i]] <- filtered_proteins[condition %in% c(comparisons[i,2], comparisons[i,3]) & (component_name %in% met_remove), ]
  
}

# select out 2/3 of missing proteins from pooled QC
group_data <- filtered_proteins[!(condition %in% c(comparisons$group_1, comparisons$group_2))]

group_ct <- length(unique(group_data$sample_name)) # number of replicates in (presumably) pooled QC
missing_ct <- group_data[is.na(norm_area), .N, by = .(condition, component_name)]
met_remove <- unique(as.character(missing_ct[N >= ((2/3) * group_ct)]$component_name))

filtering_table[[length(comparisons[,1]) + 1]] <- filtered_proteins[!(condition %in% c(comparisons$group_1, comparisons$group_2)) & (component_name %in% met_remove)]

# Remake filtered_proteins table to continue analysis
filtered_proteins <- rbindlist(filtering_table)

filtered_proteins <- unique(filtered_proteins)

# write filtered protein results to file
write.table(filtered_proteins, paste0(prefix,"filtered.txt"), quote = F, row.names = F, sep = "\t")


#============================================================================
# ANNOTATING THE RESULTS FILE
#============================================================================
# select only the results
results_unannot <- results[,1:8]
ids_to_map <- paste(as.character(unique(results_unannot$component_name)), collapse = ",")

# library to handle URLs
library(httr)

# function to determine whether the submitted job (i.e. annotations) are ready to download
isJobReady <- function(jobId) {
  pollingInterval = 5
  nTries = 20
  for (i in 1:nTries) {
    url <- paste("https://rest.uniprot.org/idmapping/status/", jobId, sep = "")
    r <- GET(url = url, accept_json())
    status <- content(r, as = "parsed")
    if (!is.null(status[["results"]]) || !is.null(status[["failedIds"]])) {
      return(TRUE)
    }
    if (!is.null(status[["messages"]])) {
      print(status[["messages"]])
      return (FALSE)
    }
    Sys.sleep(pollingInterval)
  }
  return(FALSE)
}

# function to retrieve results from query
getResultsURL <- function(redirectURL) {
  if (grepl("/idmapping/results/", redirectURL, fixed = TRUE)) {
    url <- gsub("/idmapping/results/", "/idmapping/stream/", redirectURL)
  } else {
    url <- gsub("/results/", "/results/stream/", redirectURL)
  }
}

# Uniprot IDs to match to database for annotating the results
files = list(
  ids = ids_to_map,
  from = "UniProtKB_AC-ID",
  to = "UniProtKB"
)

r <- POST(url = "https://rest.uniprot.org/idmapping/run", body = files, encode = "multipart", accept_json())

submission <- content(r, as = "parsed")

if (isJobReady(submission[["jobId"]])) {
  url <- paste("https://rest.uniprot.org/idmapping/details/", submission[["jobId"]], sep = "")
  r <- GET(url = url, accept_json())
  details <- content(r, as = "parsed")
  url <- getResultsURL(details[["redirectURL"]])
  # Using TSV format see: https://www.uniprot.org/help/api_queries#what-formats-are-available
  url <- paste(url, "?format=tsv&fields=accession,protein_name,gene_names,length,go_id", sep = "")
  r <- GET(url = url, accept_json())
  resultsTable = setDT(read.table(text = content(r), sep = "\t", header=TRUE))
  #print(resultsTable)
}

# add protein name, gene name, and length to the unannotated results table
results_unannot[resultsTable, annotation := paste0("Protein: ", i.Protein.names, "; Gene name: ", i.Gene.Names), on = .(component_name = From)]

# write annotated results
write.table(results_unannot, paste0(prefix, "results_annot.txt"), quote = F, row.names = F, sep = "\t")


#=============================================================
# CREATING GO ID TABLE FOR GO-SEA ANALYSIS
#=============================================================
# retrieve GO IDs
go_annotations <- resultsTable[, .(GeneID = From, GO_ID = Gene.Ontology.IDs)]

# calculate length new table will need to be
go_column <- go_annotations$GO_ID
go_column <- lapply(go_column, function(x){as.vector(unlist(str_split(x, "; ")))})

go_column_2 <- unlist(go_column)

# determine number that each line needs to be repeated
n_rep <- vector("numeric", length(go_column))

for (i in seq_along(go_column)) {
  
  n_rep[[i]] <- length(go_column[[i]])
  
}

# create replicates for each item in GO terms table
go_table_final <- vector("list", length(go_column))

for (i in seq_along(go_annotations$GeneID)) {
  
  # replicate that line for each number of GO IDs associated with it
  rep_lines <- do.call("rbind", replicate(n_rep[i], go_annotations[i], simplify = F))
  
  # add in these new entries into the final go_table
  go_table_final[[i]] <- rep_lines
  
}

go_table_final <- rbindlist(go_table_final)

go_table_final <- cbind(go_table_final, go_ids = go_column_2)

go_table_final$GO_ID <- NULL

go_table_final <- go_table_final[go_ids != "", ]


#============================================================================
# GO-SEA ANALYSIS
#============================================================================
# load libraries
library(fgsea)

# Read in GO definitions (names of GO terms)
go_def <- setDT(read.delim(file = GO_def))

# Get list of all unique GO terms
go_terms <- unique(go_table_final$go_ids)

# Create a list to store GO annotations in GSEA format
go_gsea <- vector("list", length(go_terms))

for (i in seq_along(go_terms)) {
  
  go_members <- unique(go_table_final[go_ids == go_terms[i]]$GeneID)
  
  go_gsea[[i]] <- go_members
  
}

# Name each element of the go_gsea list with the GO term ID
names(go_gsea) <- go_terms

# Get a list of all comparisons in the results
comparisons <- unique(results$comparison)

# create list to contain the gsea results
all_gsea_results <- vector("list", length(comparisons))

# loop through each comparison to run GSEA analysis on it
for (i in seq_along(comparisons)) {
  
  # print current comparison to console
  print(paste0("Now performing GSEA on ", comparisons[i], " comparison..."))
  
  # current comparison's results
  current_comparison <- results[comparison == comparisons[i], ]
  
  # create named vector containing log2FC values in order
  current_ranking <- setorder(current_comparison, log2FC)$log2FC
  names(current_ranking) <- setorder(current_comparison, log2FC)$component_name
  
  # run gsea for this comparison
  gsea_result <- fgsea(pathways = go_gsea,
                       stats = current_ranking,
                       minSize = 15,
                       maxSize = 500)
  
  # add column with group comparison
  gsea_result <- as.data.table(gsea_result)[, comparison := comparisons[i]]
  
  # convert leading edge list to a text string of semicolon separated accessions
  gsea_result$leadingEdge <- sapply(gsea_result$leadingEdge, paste, collapse = ";")
  
  # add gsea results to larger list of results
  all_gsea_results[[i]] <- gsea_result
  
  # print completion of the gsea analysis for current comparison
  print(paste0("GSEA on ", comparisons[i], " complete!"))
  
}

all_gsea_results <- rbindlist(all_gsea_results)

# add GO term definitions
all_gsea_results[go_def, definition := i.definition, on = .(pathway = termid)]

# write GSEA results table
write.table(all_gsea_results, paste0(prefix, "GO_SEA.txt"), sep = "\t", quote = F, row.names = F)


#===================================================================
# TOP TEN GSEA RESULTS FROM GSEA ANALYSIS
#===================================================================
# create list to contain top ten GSEA results from each comparison
top_10_gsea <- vector("list", length(comparisons))

for (i in seq_along(comparisons)) {
  # Top 10 GSEA results for each comparison (if p value is < set pval threshold)
  top_10_gsea[[i]] <- head(setorder(all_gsea_results[comparison == comparisons[i] & pval < pval_max], pval)$pathway, 10)
  
}

# combine top 10 results per comparison and remove duplicates
top_10_gsea <- unique(unlist(top_10_gsea))

# add -log10(p) column for plotting
top_10_gsea <- all_gsea_results[pathway %in% top_10_gsea, ][, log_p := -log10(pval)]

# flip sign of -log10(p) value if NES is negative
top_10_gsea[NES < 0, log_p := log_p * -1][is.na(log_p), log_p := 0]

# cast to wide format
gsea_wide <- dcast(top_10_gsea, pathway + definition ~ comparison, value.var = "log_p")

# create matrix from wide format table
logp_mat <- as.matrix(gsea_wide[, 3:ncol(gsea_wide)])

# name rows with definition
row.names(logp_mat) <- gsea_wide$definition

# calculate distances (using euclidean method)
d <- dist(logp_mat)

# perform hierarchical clustering (to make consistent with the heatmap clustering, the same Ward method will be used)
clustered <- hclust(d, method = "ward.D2")

# plot results for top 10 GSEA results
ggplot(top_10_gsea, aes(x = comparison, y = definition)) +
  geom_point(aes(size = abs(NES), color = log_p)) +
  #geom_text(aes(label = round(abs(NES), 2)), hjust = -0.5) +
  #geom_text(aes(label = paste0("p = ", format(padj, digits = 3, scientific = T))), hjust = 1.2) +
  scale_color_gradient2(low = "blue", mid = "yellow", high = "red", name = "-LogP") +
  scale_y_discrete(limits = rev(gsea_wide$definition[clustered$order])) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5),
        axis.text = element_text(color = "black")) +
  scale_x_discrete(limits = comparisons) +
  scale_size(name = "NES") +
  labs(x = NULL, y = NULL, title = "GSEA (Top 10 per Comparison)")
ggsave(paste0(prefix,"GO_SEA_plot.png"), units = "in", device = "png", dpi = 300, width = 8, height = 8)
