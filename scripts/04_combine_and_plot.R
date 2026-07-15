# 04_combine_and_plot.R
#
# combines per-schema celltype_stability.csv outputs (produced by
# 03_cluster_stability_usage.R) into a single table + generates the
# faceted comparison plot across schemas. also combines per-resolution
# gene-robustness outputs (produced by compute_gene_jaccard() in
# 01_gene_robustness.R) into one table + one multi-page PDF across
# resolutions.
#
# changes from original:
# - schemas to combine discovered from `output_dir` subfolders rather
#   than hardcoded in an `idents` vector, so works regardless of how
#   many schemas were actually run
# - no dataset-specific paths or filenames
#
# Requires: tidyverse, ggplot2. pdftools ONLY required if actually using
# combine_gene_robustness_across_resolutions() (loaded lazily there, not
# up front) -- stability combine doesn't need it at all.

library(tidyverse)
library(ggplot2)

zissou_colors <- c('#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00')


#' Discover and combine celltype_stability.csv files across all schema
#' subfolders under output_dir.
#'
#' Expects the directory structure produced by assess_stability():
#'   {output_dir}/{schema}/scclusteval/stability/tables/celltype_stability.csv
#'
#' @param output_dir Base directory containing one subfolder per schema
#' @return Combined data frame with a `schema` column, duplicate cell types
#'   within a schema aggregated by mean
combine_stability_results <- function(output_dir)
{
  schema_dirs <- list.dirs(output_dir, recursive = FALSE, full.names = FALSE)

  all_stability <- list()

  for (schema in schema_dirs)
  {
    stability_file <- file.path(output_dir, schema, 'scclusteval', 'stability', 'tables', 'celltype_stability.csv')

    if (!file.exists(stability_file))
    {
      warning('No celltype_stability.csv found for schema "', schema, '" at expected path: ', stability_file, ' -- skipping.')
      next
    }

    tbl <- read.csv(stability_file, header = TRUE)
    tbl$schema <- schema
    all_stability[[schema]] <- tbl
  }

  if (length(all_stability) == 0)
  {
    stop('No celltype_stability.csv files found under any subfolder of ', output_dir)
  }

  combined <- do.call(rbind, all_stability)
  rownames(combined) <- NULL

  # aggregate duplicate cell types within the same schema (e.g. if a schema's
  # annotation collapsed multiple raw clusters into one celltype label)
  combined <- combined %>%
    group_by(schema, celltype) %>%
    summarise(
      mean_jaccard = mean(mean_jaccard, na.rm = TRUE),
      sd_jaccard = mean(sd_jaccard, na.rm = TRUE),
      min_jaccard = min(min_jaccard, na.rm = TRUE),
      max_jaccard = max(max_jaccard, na.rm = TRUE),
      n_stable = sum(n_stable, na.rm = TRUE),
      pct_stable = mean(pct_stable, na.rm = TRUE),
      .groups = 'drop'
    )

  return(combined)
}


#' Generate the faceted cell-type stability comparison plot across schemas.
#'
#' @param combined_stability Output of combine_stability_results()
#' @param out_path Path to save the plot PDF
#' @param stability_threshold Solid reference line (default 0.75, "stable")
#' @param moderate_threshold Dashed reference line (default 0.60, "moderately stable")
#' @param celltype_order Optional character vector specifying celltype
#'   display order (top to bottom on the flipped axis). If NULL, alphabetical.
#' @return The ggplot object
plot_combined_stability <- function(combined_stability, out_path, stability_threshold = 0.75, moderate_threshold = 0.60, celltype_order = NULL)
{
  if (is.null(celltype_order))
  {
    celltype_order <- sort(unique(combined_stability$celltype), decreasing = TRUE)
  }

  combined_stability$celltype <- factor(combined_stability$celltype, levels = celltype_order)

  n_schemas <- length(unique(combined_stability$schema))

  p <- ggplot(combined_stability, aes(x = celltype, y = mean_jaccard, fill = mean_jaccard)) +
    geom_col(color = 'black', linewidth = 0.3) +
    geom_errorbar(aes(ymin = mean_jaccard - sd_jaccard, ymax = mean_jaccard + sd_jaccard), width = 0.2, linewidth = 0.5) +
    geom_hline(yintercept = stability_threshold, linetype = 'dashed', color = 'black', linewidth = 0.8) +
    geom_hline(yintercept = moderate_threshold, linetype = 'dashed', color = 'grey50', linewidth = 0.5) +
    facet_wrap(~schema, ncol = 2) +
    coord_flip() +
    ylim(0, 1) +
    scale_fill_gradientn(colors = zissou_colors, limits = c(0, 1), name = 'mean Jaccard\nindex') +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 11, face = 'bold'),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 10, face = 'bold'),
      legend.position = 'right',
      panel.grid.major.y = element_blank()
    ) +
    labs(
      title = 'Cluster stability across schemas: mean Jaccard index under subsampling',
      x = 'Cell Type',
      y = 'Mean Jaccard Index'
    )

  ggsave(out_path, p, width = 12, height = max(6, ceiling(n_schemas / 2) * 3.5))

  return(p)
}


#' Discover + combine per-resolution gene-robustness outputs (produced by
#' compute_gene_jaccard() in 01_gene_robustness.R, one call per resolution
#' in the match/robustness loop in 05_run_validation.R) into one table +
#' one multi-page PDF spanning every resolution.
#'
#' Expects the directory structure produced by that loop:
#'   {gene_robustness_dir}/res_{resolution}/all_metacluster_gene_jaccard.txt
#'   {gene_robustness_dir}/res_{resolution}/all_metaclusters_gene_jaccard_panel.pdf
#'
#' No separate resolution column is added to the combined table -- schema
#' identifiers already carry resolution (e.g. sct.harmony_snn_res.0.2), so
#' it's recoverable from row/col directly. meta_cluster IDs, however, are
#' only unique WITHIN a resolution (community detection reruns per
#' resolution) -- do not group/filter by meta_cluster alone across
#' resolutions w/o also keying on row/col.
#'
#' @param gene_robustness_dir base dir containing one res_{resolution}
#'   subfolder per resolution (i.e. opt$`output-dir`/gene_robustness)
#' @return combined data frame (row, col, jaccard, meta_cluster), one
#'   res_{resolution} folder's rows stacked per resolution
combine_gene_robustness_across_resolutions <- function(gene_robustness_dir)
{
  res_dirs <- list.dirs(gene_robustness_dir, recursive = FALSE, full.names = TRUE)
  res_dirs <- res_dirs[grepl('^res_', basename(res_dirs))]

  if (length(res_dirs) == 0)
  {
    stop('no res_* subfolders found under ', gene_robustness_dir)
  }

  all_jac <- list()
  all_pdfs <- character(0)

  for (res_dir in res_dirs)
  {
    txt_file <- file.path(res_dir, 'all_metacluster_gene_jaccard.txt')
    pdf_file <- file.path(res_dir, 'all_metaclusters_gene_jaccard_panel.pdf')

    if (!file.exists(txt_file))
    {
      warning('no all_metacluster_gene_jaccard.txt found under ', res_dir, ' -- skipping.')
      next
    }

    all_jac[[res_dir]] <- read.table(txt_file, header = TRUE, sep = '\t')

    if (file.exists(pdf_file)) all_pdfs <- c(all_pdfs, pdf_file)
  }

  if (length(all_jac) == 0)
  {
    stop('no all_metacluster_gene_jaccard.txt files found under any res_* subfolder of ', gene_robustness_dir)
  }

  combined <- do.call(rbind, all_jac)
  rownames(combined) <- NULL

  combined_txt <- file.path(gene_robustness_dir, 'all_resolutions_gene_jaccard.txt')
  write.table(combined, combined_txt, sep = '\t', quote = FALSE, row.names = FALSE)

  if (length(all_pdfs) > 0)
  {
    if (!requireNamespace('pdftools', quietly = TRUE))
    {
      stop('package "pdftools" required to merge per-resolution PDF panels ',
           '(pdf_combine()) -- install it (install.packages("pdftools")) ',
           'or run install_dependencies.R. not needed for stability combine, ',
           'only for combine_gene_robustness_across_resolutions().')
    }
    combined_pdf <- file.path(gene_robustness_dir, 'all_resolutions_gene_jaccard_panel.pdf')
    pdftools::pdf_combine(all_pdfs, output = combined_pdf)
  } else
  {
    combined_pdf <- NULL
    warning('no per-resolution panel PDFs found -- combined_txt written, but no combined PDF.')
  }

  return(list(combined_data = combined, combined_txt = combined_txt, combined_pdf = combined_pdf))
}


# ============================================================================
# Example usage
# ============================================================================
# output_dir <- 'CLUSTER_STABILITY/PROJECT/'
#
# combined_stability <- combine_stability_results(output_dir)
# write.csv(combined_stability, file.path(output_dir, 'combined_celltype_stability.csv'), row.names = FALSE)
#
# p <- plot_combined_stability(
#   combined_stability,
#   out_path = file.path(output_dir, 'combined_stability_facets.pdf'),
#   celltype_order = rev(c('B_cells', 'Endothelial', 'Fibroblasts', 'Follicular',
#                           'Myeloid', 'Pericytes', 'Plasma', 'T_cells'))
# )
# print(p)
