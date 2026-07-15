# 01_gene_robustness.R
#
# gene signature robustness across analysis schemas, via pairwise Jaccard
# similarity of top marker genes per meta-cluster
#
# generalized version, consumes automated cross-schema cluster matching
# produced by 00_match_clusters_across_schemas.R, in place of hand-built
# per-schema annotation vectors. meta-cluster IDs anonymous
# (meta_cluster_1, meta_cluster_2, ...) unless renamed after matching
#
# Requires: tidyverse, ggplot2, patchwork

library(tidyverse)
library(ggplot2)
library(patchwork)

set.seed(42)

zissou_colors <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')


#' load per-schema z-score tables + tag each row w/ its schema
#'
#' expects one file per schema at:
#'   {z_score_dir}/{project}_z_scores.{schema}_snn_res.{resolution}.txt
#' matching SWANS's standard z-score table output naming
#'
#' @param z_score_dir dir containing z-score tables
#' @param project SWANS project name (used in file naming pattern)
#' @param schema_cols character vector of schema identifiers, e.g.
#'   c('sct.harmony', 'sct.cca', 'standard.rpca') -- WITHOUT the
#'   '_snn_res.X' suffix, since that's supplied separately
#' @param resolution clustering resolution, must match schema_cols +
#'   resolution used for cross-schema matching in 00_match_clusters_across_schemas.R
#' @return combined data frame of all schemas' z-score tables, w/ a `schema`
#'   col added (format matches meta.data col names, e.g. 'sct.harmony_snn_res.0.2')
load_zscore_tables <- function(z_score_dir, project, schema_cols, resolution)
{
  all_tables <- list()

  for (schema in schema_cols)
  {
    file_path <- file.path(z_score_dir, paste0(project, '_z_scores.', schema, '_snn_res.', resolution, '.txt'))

    if (!file.exists(file_path))
    {
      warning('z-score file not found, skipping: ', file_path)
      next
    }

    tbl <- read.table(file_path, header = TRUE, sep = '\t')
    tbl$schema <- paste0(schema, '_snn_res.', resolution)
    all_tables[[schema]] <- tbl
  }

  if (length(all_tables) == 0)
  {
    stop('no z-score tables found in ', z_score_dir, ' for requested schemas.')
  }

  return(do.call(rbind, all_tables))
}


#' join z-score data to meta-cluster assignments + select top genes per
#' meta-cluster, aggregating across clusters that map to same meta-cluster
#'
#' @param zscore_data output of load_zscore_tables()
#' @param matched_clusters output of match_clusters_across_schemas() --
#'   cols: schema, cluster, meta_cluster, is_singleton
#' @param top_n num top genes (by mean z-score) to retain per meta-cluster
#'   per schema (default 100)
#' @param exclude_singletons if TRUE (default), meta-clusters w/ no
#'   cross-schema match dropped before robustness comparison, since Jaccard
#'   comparison needs at least 2 schemas per meta-cluster to be meaningful
#' @return data frame: schema | meta_cluster | gene | mean_z (top_n rows per
#'   schema x meta_cluster group)
aggregate_genes_by_metacluster <- function(zscore_data, matched_clusters, top_n = 100, exclude_singletons = TRUE)
{
  matched_clusters$cluster <- as.character(matched_clusters$cluster)
  zscore_data$cluster <- as.character(zscore_data$cluster)

  if (exclude_singletons)
  {
    matched_clusters <- matched_clusters %>% filter(!is_singleton)
  }

  joined <- zscore_data %>%
    left_join(matched_clusters, by = c('schema', 'cluster'))

  n_unmatched <- sum(is.na(joined$meta_cluster))
  if (n_unmatched > 0)
  {
    warning(n_unmatched, ' z-score rows had no corresponding entry in ',
            'matched_clusters + were dropped. check that `schema` + ',
            '`cluster` values align between the two inputs.')
    joined <- joined %>% filter(!is.na(meta_cluster))
  }

  # aggregate genes within each (schema, meta_cluster): mean z-score across
  # any raw clusters mapping to same meta-cluster within that schema, then
  # take top_n unique genes
  aggregated <- joined %>%
    group_by(schema, meta_cluster, gene) %>%
    summarise(mean_z = mean(z.score, na.rm = TRUE), .groups = 'drop') %>%
    group_by(schema, meta_cluster) %>%
    arrange(desc(mean_z), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()

  return(aggregated)
}


#' compute pairwise Jaccard similarity of top genes across schemas, per
#' meta-cluster, + generate heatmap plots
#'
#' @param aggregated_genes output of aggregate_genes_by_metacluster()
#' @param out_dir dir to write per-meta-cluster Jaccard tables/plots +
#'   combined panel plot
#' @return list w/ combined Jaccard data frame + combined panel plot
compute_gene_jaccard <- function(aggregated_genes, out_dir)
{
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  all_jac_data <- list()
  all_plots <- list()

  meta_clusters <- unique(aggregated_genes$meta_cluster)

  for (mc in meta_clusters)
  {
    mc_data <- aggregated_genes %>% filter(meta_cluster == mc)
    schemas <- unique(mc_data$schema)

    if (length(schemas) < 2)
    {
      warning('meta-cluster ', mc, ' has fewer than 2 schemas represented, skipping Jaccard comparison.')
      next
    }

    gene_lists <- list()
    for (s in schemas)
    {
      gene_lists[[s]] <- mc_data %>% filter(schema == s) %>% pull(gene)
    }

    n <- length(gene_lists)
    jac_mat <- matrix(NA, nrow = n, ncol = n)
    rownames(jac_mat) <- names(gene_lists)
    colnames(jac_mat) <- names(gene_lists)

    for (i in 1:n)
    {
      for (j in 1:n)
      {
        intersection <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
        union_n <- length(union(gene_lists[[i]], gene_lists[[j]]))
        jac_mat[i, j] <- if (union_n > 0) intersection / union_n else NA
      }
    }

    jac_df <- as.data.frame(jac_mat) %>%
      rownames_to_column('row') %>%
      pivot_longer(-row, names_to = 'col', values_to = 'jaccard') %>%
      mutate(meta_cluster = mc)

    write.table(jac_df, file.path(out_dir, paste0(mc, '_gene_jaccard.txt')),
                sep = '\t', quote = FALSE, row.names = FALSE)
    all_jac_data[[mc]] <- jac_df

    p <- ggplot(jac_df, aes(x = col, y = row, fill = jaccard)) +
      geom_tile(color = 'white', linewidth = 0.5) +
      scale_fill_gradientn(colors = zissou_colors, limits = c(0, 1), name = 'Jaccard Index') +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_blank(),
            plot.title = element_text(size = 12, face = 'bold')) +
      geom_text(aes(label = round(jaccard, 2)), size = 3) +
      labs(title = paste0(mc, ' -- Gene Signature Jaccard Index'))

    ggsave(file.path(out_dir, paste0(mc, '_gene_jaccard.pdf')), p, width = 10, height = 10)
    all_plots[[mc]] <- p
  }

  combined_jac_data <- do.call(rbind, all_jac_data)
  write.table(combined_jac_data, file.path(out_dir, 'all_metacluster_gene_jaccard.txt'),
              sep = '\t', quote = FALSE, row.names = FALSE)

  panel_plot <- wrap_plots(all_plots, ncol = 3, guides = 'collect') & theme(legend.position = 'bottom')
  ggsave(file.path(out_dir, 'all_metaclusters_gene_jaccard_panel.pdf'), panel_plot,
         width = 16, height = ceiling(length(all_plots) / 3) * 5)

  return(list(jaccard_data = combined_jac_data, panel_plot = panel_plot))
}


#' collate gene-signature Jaccard results across every resolution scanned
#' by match+robustness loop into one combined table + summary plot -- helps
#' identify which resolution gives most reproducible gene signatures
#' overall, before committing to one schema for annotation
#'
#' @param gene_robustness_dir top-level dir w/ res_<X> subfolders, each w/
#'   an all_metacluster_gene_jaccard.txt from compute_gene_jaccard()
#' @return list w/ combined_txt (path), combined_pdf (path, or NULL if
#'   fewer than 2 resolutions found), combined_data (data frame)
combine_gene_robustness_across_resolutions <- function(gene_robustness_dir)
{
  res_dirs <- list.dirs(gene_robustness_dir, recursive = FALSE, full.names = TRUE)
  res_dirs <- res_dirs[grepl('^res_[0-9.]+$', basename(res_dirs))]

  if (length(res_dirs) == 0)
  {
    stop('no res_<X> subfolders found in ', gene_robustness_dir,
         '. run "match"+"robustness" steps first.')
  }

  all_data <- list()
  for (rd in res_dirs)
  {
    jac_file <- file.path(rd, 'all_metacluster_gene_jaccard.txt')
    if (!file.exists(jac_file))
    {
      warning('no all_metacluster_gene_jaccard.txt found in ', rd, ' -- skipping.')
      next
    }
    res_val <- sub('^res_', '', basename(rd))
    jac_data <- read.table(jac_file, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
    jac_data$resolution <- res_val
    all_data[[rd]] <- jac_data
  }

  if (length(all_data) == 0)
  {
    stop('no valid all_metacluster_gene_jaccard.txt files found under ', gene_robustness_dir)
  }

  combined_data <- do.call(rbind, all_data)
  combined_txt <- file.path(gene_robustness_dir, 'all_resolutions_gene_jaccard.txt')
  write.table(combined_data, combined_txt, sep = '\t', quote = FALSE, row.names = FALSE)

  combined_pdf <- NULL
  n_res <- length(unique(combined_data$resolution))
  if (n_res >= 2)
  {
    # drop self-comparisons (row == col, jaccard == 1 trivially). meta_cluster
    # IDs are NOT comparable across resolutions (matching reruns fresh per
    # resolution), so can't color/facet by meta_cluster here -- schema
    # identity (row/col) IS stable across resolutions, so classify each pair
    # as within- vs cross-normalization instead, matching the concordance
    # pattern already established in the manuscript's single-resolution
    # analysis (within-norm pairs consistently more concordant than
    # cross-norm pairs)
    summary_data <- combined_data %>%
      filter(row != col) %>%
      mutate(
        norm_row = sub('\\..*$', '', row),
        norm_col = sub('\\..*$', '', col),
        pair_type = ifelse(norm_row == norm_col, 'within-normalization', 'cross-normalization')
      )

    p <- ggplot(summary_data, aes(x = resolution, y = jaccard, fill = pair_type)) +
      geom_boxplot(outlier.alpha = 0.4, position = position_dodge(width = 0.75)) +
      scale_fill_manual(values = c('within-normalization' = zissou_colors[1],
                                    'cross-normalization' = zissou_colors[4]),
                         name = NULL) +
      theme_minimal() +
      labs(title = 'gene signature Jaccard concordance by resolution',
           subtitle = 'within-normalization: sct-sct or standard-standard pairs. cross-normalization: sct vs standard pairs.',
           x = 'resolution', y = 'Jaccard index (cross-schema pairs)') +
      theme(plot.title = element_text(size = 12, face = 'bold'),
            plot.subtitle = element_text(size = 9, color = 'grey30'),
            legend.position = 'bottom')

    combined_pdf <- file.path(gene_robustness_dir, 'all_resolutions_gene_jaccard_summary.pdf')
    ggsave(combined_pdf, p, width = 9, height = 6)
  }

  list(combined_txt = combined_txt, combined_pdf = combined_pdf, combined_data = combined_data)
}


# ============================================================================
# example usage
# ============================================================================
# source('00_match_clusters_across_schemas.R')
#
# seurat_obj <- qs2::qs_read('path/to/analyzed_seurat_object.qs2')
# schema_cols_full <- detect_schema_columns(seurat_obj, resolution = 0.2)
# matched_clusters <- match_clusters_across_schemas(seurat_obj, schema_cols_full)
#
# # schema identifiers WITHOUT _snn_res.X suffix, for z-score file naming
# schema_short <- c('sct.harmony', 'sct.cca', 'sct.rpca',
#                    'standard.harmony', 'standard.cca', 'standard.rpca')
#
# zscore_data <- load_zscore_tables(
#   z_score_dir = 'data/endpoints/PROJECT/analysis/report/tables',
#   project = 'PROJECT',
#   schema_cols = schema_short,
#   resolution = 0.2
# )
#
# aggregated <- aggregate_genes_by_metacluster(zscore_data, matched_clusters, top_n = 100)
#
# results <- compute_gene_jaccard(aggregated, out_dir = 'gene_robustness_output/')