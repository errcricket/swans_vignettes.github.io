# 02_cluster_stability_functions.R
#
# functions for subsampling-based cluster stability assessment via scclusteval
#
# changes from original:
# - `resolution` now REQUIRED arg throughout (no silent default). original
#   code defaulted run_stability_analysis()'s resolution to 0.2, but calling
#   wrapper (assess_stability(), see script 03) hardcoded resolution = 0.5
#   for every call regardless of resolution input clusters were actually
#   generated at. silently reclustered subsampled replicates at different
#   resolution than original clustering being evaluated, biasing stability
#   scores. making resolution required, threaded through explicitly from
#   top-level call, closes off this failure mode
# - removed unused `qs` + `viridis` library calls (qs2 + manual Zissou
#   colors are what's actually used)
# - output filename now reflects actual sampling_fractions used, rather
#   than hardcoded '.8.9.' suffix that goes stale if fractions change
# - removed plot_stability_2(), exact duplicate of plot_stability()
# - stability threshold (originally hardcoded 0.75 in 4 places) now single
#   parameter w/ default, passed through consistently
#
# Requires: Seurat, scclusteval, tidyverse, qs2, ggplot2

library(Seurat)
library(scclusteval)
library(tidyverse)
library(qs2)


#' run scclusteval subsampling-based stability analysis on a Seurat object
#'
#' @param seurat_obj Seurat object w/ normalized data + cluster identities
#'   set via Idents(). must contain raw counts (required by scclusteval to
#'   reprocess subsampled replicates)
#' @param output_dir dir to save results
#' @param resolution REQUIRED. clustering resolution used for ORIGINAL
#'   cluster identities in `seurat_obj`. subsampled replicates reclustered
#'   at this same resolution, so must match how Idents(seurat_obj) was
#'   generated -- passing different value than original clustering
#'   silently biases stability scores (see note above)
#' @param sampling_fractions fractions of cells to subsample (default c(0.8, 0.9))
#' @param n_reps num repetitions per fraction (default 10)
#' @param pc.use num PCs used for original clustering (must match)
#' @param seed random seed for reproducibility
#' @param verbose print progress messages
#' @return list w/ jaccard_results + output_dir
run_stability_analysis <- function(seurat_obj, output_dir, resolution, sampling_fractions = c(0.8, 0.9), n_reps = 10, pc.use = 14, seed = 42, verbose = TRUE)
{
  if (missing(resolution) || is.null(resolution))
  {
    stop('`resolution` required + must match resolution used to ',
         'generate original cluster identities in `seurat_obj`. passing ',
         'a mismatched resolution here silently biases stability scores by ',
         'reclustering subsampled replicates at wrong granularity.')
  }

  if (verbose) cat('starting scclusteval analysis (resolution =', resolution, ')...\n')

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (!('counts' %in% Layers(seurat_obj))) stop('raw counts required in seurat object for scclusteval')

  if (verbose) cat(sprintf('input: %d cells, %d clusters\n', ncol(seurat_obj), length(unique(Idents(seurat_obj)))))

  set.seed(seed)

  all_jaccard_results <- list()

  for (frac in sampling_fractions)
  {
    if (verbose) cat(sprintf('processing fraction: %.1f\n', frac))

    fraction_results <- list()

    for (i in seq_len(n_reps))
    {
      if (verbose) cat(sprintf('  rep %d/%d\n', i, n_reps))

      subset_obj <- RandomSubsetData(object = seurat_obj, rate = frac, random.subset.seed = seed + i)

      # preprocess + recluster at SAME resolution as original clustering
      subset_obj <- PreprocessSubsetDataV2(object = subset_obj, num.pc = 30, pc.use = pc.use, resolution = resolution)

      jaccard <- PairWiseJaccardSets(ident1 = Idents(seurat_obj)[colnames(subset_obj)], ident2 = Idents(subset_obj))
      fraction_results[[i]] <- jaccard
    }

    all_jaccard_results[[as.character(frac)]] <- fraction_results
  }

  if (verbose) cat('scclusteval completed.\n')

  fraction_tag <- paste(gsub('\\.', '', as.character(sampling_fractions)), collapse = '_')
  qs_save(all_jaccard_results, file.path(output_dir, paste0('jaccard_results.', fraction_tag, '.qs2')))

  return(list(jaccard_results = all_jaccard_results, output_dir = output_dir, resolution = resolution))
}


#' parse jaccard results + calculate per-cluster stability
#'
#' @param stability_output full output from run_stability_analysis()
#' @param stability_threshold Jaccard threshold above which a replicate
#'   considered "stable" (default 0.75)
#' @return data frame w/ stability metrics per cluster
calculate_cluster_stability <- function(stability_output, stability_threshold = 0.75)
{
  jaccard_results <- stability_output$jaccard_results

  all_rows <- list()

  for (frac in names(jaccard_results))
  {
    for (i in seq_along(jaccard_results[[frac]]))
    {
      jmat <- jaccard_results[[frac]][[i]]
      highest <- apply(jmat, 1, max)
      df <- data.frame(celltype = names(highest), jaccard = highest, fraction = as.numeric(frac), rep = i, row.names = NULL)
      all_rows[[length(all_rows) + 1]] <- df
    }
  }

  jaccard_combined <- do.call(rbind, all_rows)

  stability_metrics <- jaccard_combined %>%
    group_by(celltype) %>%
    summarise(
      mean_jaccard = mean(jaccard, na.rm = TRUE),
      sd_jaccard = sd(jaccard, na.rm = TRUE),
      min_jaccard = min(jaccard, na.rm = TRUE),
      max_jaccard = max(jaccard, na.rm = TRUE),
      n_stable = sum(jaccard > stability_threshold, na.rm = TRUE),
      pct_stable = 100 * mean(jaccard > stability_threshold, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean_jaccard))

  qs_save(stability_metrics, file.path(stability_output$output_dir, 'cluster_stability.qs2'))

  return(stability_metrics)
}


#' aggregate cluster-level stability to cell-type level
#'
#' @param cluster_stability_df output of calculate_cluster_stability()
#' @param cluster_to_celltype_mapping data frame: cluster, celltype
#' @param stability_threshold Jaccard threshold for STABLE/UNSTABLE status
#'   (default 0.75; should match threshold used in calculate_cluster_stability())
#' @return data frame w/ cell-type-level stability status
calculate_celltype_stability <- function(cluster_stability_df, cluster_to_celltype_mapping, stability_threshold = 0.75)
{
  celltype_stability <- cluster_stability_df %>%
    select(celltype, mean_jaccard, sd_jaccard, min_jaccard, max_jaccard, n_stable, pct_stable) %>%
    mutate(stability_status = ifelse(mean_jaccard > stability_threshold, 'STABLE', 'UNSTABLE')) %>%
    arrange(desc(mean_jaccard))

  return(celltype_stability)
}


#' create publication-ready stability plots (heatmap + summary barplot)
#'
#' @param cluster_stability_df output of calculate_cluster_stability() after
#'   assess_stability()'s relabeling step -- must already carry both
#'   `cluster` (number) + `celltype` (label) cols
#' @param output_dir dir to save plots
#' @param stability_threshold Jaccard threshold shown as reference line
#'   (default 0.75)
#' @param save_plots whether to write PDFs to output_dir (default TRUE)
#' @return list w/ heatmap + barplot ggplot objects
plot_stability <- function(cluster_stability_df, output_dir, stability_threshold = 0.75, save_plots = TRUE)
{
  library(ggplot2)

  zissou_colors <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')

  # no join needed here -- cluster_stability_df already has both cluster +
  # celltype natively (see assess_stability()). joining by celltype label
  # here previously caused a many-to-many blow-up whenever multiple
  # clusters shared the same label, since celltype isn't a unique key
  plot_data <- cluster_stability_df %>% filter(!is.na(celltype))

  p_heatmap <- ggplot(plot_data, aes(x = celltype, y = as.factor(cluster), fill = mean_jaccard)) +
    geom_tile(color = 'white', linewidth = 0.5) +
    scale_fill_gradientn(colors = zissou_colors, name = 'mean jaccard\nindex', limits = c(0, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_text(size = 12, face = 'bold'),
          plot.title = element_text(size = 14, face = 'bold')) +
    labs(title = 'cluster stability under subsampling (scclusteval)', x = 'cell type', y = 'cluster') +
    geom_text(aes(label = round(mean_jaccard, 2)), size = 2.5, color = 'black')

  celltype_stats <- plot_data %>%
    group_by(celltype) %>%
    summarise(mean_jaccard = mean(mean_jaccard), sd_jaccard = sd(mean_jaccard),
              stability_status = ifelse(mean_jaccard > stability_threshold, 'STABLE', 'UNSTABLE'),
              .groups = 'drop') %>%
    arrange(desc(mean_jaccard))

  p_barplot <- ggplot(celltype_stats, aes(x = factor(celltype, levels = sort(unique(celltype), decreasing = TRUE)), y = mean_jaccard, fill = stability_status)) +
    geom_col(color = 'black', linewidth = 0.5) +
    geom_errorbar(aes(ymin = mean_jaccard - sd_jaccard, ymax = mean_jaccard + sd_jaccard), width = 0.2, linewidth = 0.5) +
    geom_hline(yintercept = stability_threshold, linetype = 'dashed', color = 'red', linewidth = 1) +
    scale_fill_manual(values = c('STABLE' = zissou_colors[1], 'UNSTABLE' = zissou_colors[5])) +
    coord_flip() +
    theme_minimal() +
    theme(axis.title = element_text(size = 12, face = 'bold'), plot.title = element_text(size = 14, face = 'bold'), legend.position = 'bottom') +
    labs(title = 'cell type stability summary', subtitle = paste0('mean jaccard index (threshold = ', stability_threshold, ')'),
         x = 'cell type', y = 'mean jaccard index', fill = 'status') +
    ylim(0, 1)

  if (save_plots)
  {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    ggsave(file.path(output_dir, 'stability_heatmap.pdf'), p_heatmap, width = 8, height = 10)
    ggsave(file.path(output_dir, 'celltype_stability_barplot.pdf'), p_barplot, width = 8, height = 6)
    cat(sprintf('plots saved to %s\n', output_dir))
  }

  return(list(heatmap = p_heatmap, barplot = p_barplot))
}


#' save stability tables as CSVs for manuscript use
export_stability_results <- function(cluster_stability_df, celltype_stability_df, output_dir)
{
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(cluster_stability_df, file.path(output_dir, 'cluster_stability.csv'), row.names = FALSE)
  write.csv(celltype_stability_df, file.path(output_dir, 'celltype_stability.csv'), row.names = FALSE)
  cat(sprintf('results exported to %s\n', output_dir))
}


#' generate a summary sentence for manuscript use
summarize_stability <- function(celltype_stability_df, stability_threshold = 0.75)
{
  mean_jac <- mean(celltype_stability_df$mean_jaccard)
  sd_jac <- sd(celltype_stability_df$mean_jaccard)
  n_stable <- sum(celltype_stability_df$stability_status == 'STABLE')
  n_total <- nrow(celltype_stability_df)

  sprintf('%d of %d cell types exceed the stability threshold (jaccard > %.2f) under subsampling (mean jaccard = %.3f +/- %.3f).',
          n_stable, n_total, stability_threshold, mean_jac, sd_jac)
}
