# 03_cluster_stability_usage.R
#
# runs subsampling-based cluster stability assessment (scclusteval) across
# one or more analysis schemas, looping over a single shared annotation
# table instead of repeating a manual block of code per schema
#
# changes from original:
# - assess_stability() now takes `resolution` as a required arg + passes
#   it through explicitly, instead of run_stability_analysis() being
#   called w/ a hardcoded resolution = 0.5 regardless of what resolution
#   input clusters were actually generated at
# - per-schema cluster labels/numbers read from a single annotation table
#   (same format as SWANS's CLUSTER_ANNOTATION_FILE: cluster, celltype,
#   optionally schema + resolution) rather than duplicated as inline
#   vectors for every schema. same class of fix as resolution threading
#   above -- one source of truth instead of several hand-maintained
#   copies that can silently drift out of sync
# - schemas processed in a loop rather than copy-pasted blocks
#
# Requires: 02_cluster_stability_functions.R, Seurat, qs2, tidyverse

source('02_cluster_stability_functions.R')
library(tidyverse)


#' run stability assessment for a single schema
#'
#' @param seurat_obj Seurat object
#' @param output_dir base output dir
#' @param out_name subdir name for this schema's output (typically schema
#'   identifier, e.g. 'standard.rpca_snn_res.0.2')
#' @param cluster_labels character vector of cell type labels, in same
#'   order as cluster_numbers
#' @param cluster_numbers numeric/character vector of cluster IDs
#' @param resolution REQUIRED. must match resolution used to generate
#'   cluster identities currently set via Idents(seurat_obj)
#' @param jaccard_object optional path to a previously-saved jaccard_results
#'   .qs2 file, to skip rerunning the (expensive) subsampling step
#' @param stability_threshold Jaccard threshold for STABLE/UNSTABLE (default 0.75)
#' @return list w/ scclusteval_output, cluster_stability, celltype_stability,
#'   plots, manuscript_summary
assess_stability <- function(seurat_obj, output_dir, out_name, cluster_labels, cluster_numbers, resolution, jaccard_object = NULL, stability_threshold = 0.75)
{
  if (missing(resolution) || is.null(resolution))
  {
    stop('`resolution` required for assess_stability() + must match ',
         'resolution used to generate cluster identities currently ',
         'set via Idents(seurat_obj). this is the exact parameter that was ',
         'previously hardcoded + caused a resolution mismatch bug -- do ',
         'not reintroduce a default here.')
  }

  scclusteval_dir <- file.path(output_dir, out_name, 'scclusteval')
  stability_dir <- file.path(scclusteval_dir, 'stability')
  stability_plots <- file.path(stability_dir, 'figures')
  stability_tables <- file.path(stability_dir, 'tables')
  dir.create(stability_plots, recursive = TRUE, showWarnings = FALSE)
  dir.create(stability_tables, recursive = TRUE, showWarnings = FALSE)

  cat(sprintf('loaded: %d cells, %d clusters (resolution = %s)\n', ncol(seurat_obj), length(unique(Idents(seurat_obj))), resolution))
  print(table(Idents(seurat_obj)))

  if (!is.null(jaccard_object))
  {
    cat('loading existing jaccard results...\n')
    jaccard_results <- qs2::qs_read(jaccard_object)
    scclusteval_output <- list(jaccard_results = jaccard_results, output_dir = stability_dir, resolution = resolution)
  }
  else
  {
    cat('running scclusteval...\n')
    scclusteval_output <- run_stability_analysis(
      seurat_obj = seurat_obj,
      output_dir = scclusteval_dir,
      resolution = resolution,
      sampling_fractions = c(0.8, 0.9),
      n_reps = 10,
      seed = 42,
      verbose = TRUE
    )
  }

  cat('calculating cluster-level stability...\n')
  cluster_stability <- calculate_cluster_stability(scclusteval_output, stability_threshold = stability_threshold)

  # defensive checks -- this exact class of bug (celltype column silently
  # ending up numeric instead of the real label) has happened before + went
  # undetected across every schema until someone inspected the raw CSV.
  # fail loudly here instead of silently exporting wrong labels again.
  if (length(cluster_labels) != length(cluster_numbers))
  {
    stop('cluster_labels (n=', length(cluster_labels), ') + cluster_numbers (n=',
         length(cluster_numbers), ') have different lengths -- cannot map reliably.')
  }
  if (!is.character(cluster_labels) || all(grepl('^[0-9.]+$', cluster_labels)))
  {
    stop('cluster_labels looks numeric, not a real celltype label vector: ',
         paste(utils::head(cluster_labels), collapse = ', '),
         ' -- check the annotation file/config upstream of assess_stability().')
  }

  # direct match()-based mapping instead of a join -- avoids silent
  # type/format mismatches a left_join can mask
  # preserve raw cluster number as its own column BEFORE remapping celltype
  # -- plot_stability()'s heatmap needs both (celltype for x-axis, cluster
  # for y-axis). the old join-based approach added this as a side effect;
  # match()-based remapping below doesn't, so it's done explicitly here.
  cluster_stability$cluster <- cluster_stability$celltype

  celltype_lookup <- as.character(cluster_labels)
  names(celltype_lookup) <- as.character(cluster_numbers)
  cluster_stability$celltype <- celltype_lookup[as.character(cluster_stability$celltype)]

  n_unmapped <- sum(is.na(cluster_stability$celltype))
  if (n_unmapped > 0)
  {
    stop(n_unmapped, ' cluster(s) in cluster_stability had no matching entry in ',
         'cluster_numbers -- raw cluster IDs present: ',
         paste(names(celltype_lookup), collapse = ', '))
  }

  cluster_to_celltype <- data.frame(cluster = cluster_numbers, celltype = cluster_labels)

  cat('calculating cell-type-level stability...\n')
  celltype_stability <- calculate_celltype_stability(cluster_stability, cluster_to_celltype, stability_threshold = stability_threshold)

  cat('generating plots...\n')
  plots <- plot_stability(cluster_stability, stability_plots, stability_threshold = stability_threshold, save_plots = TRUE)

  cat('exporting results...\n')
  export_stability_results(cluster_stability, celltype_stability, stability_tables)

  manuscript_summary <- summarize_stability(celltype_stability, stability_threshold = stability_threshold)
  cat(manuscript_summary, '\n')

  return(list(
    scclusteval_output = scclusteval_output,
    cluster_stability = cluster_stability,
    celltype_stability = celltype_stability,
    plots = plots,
    manuscript_summary = manuscript_summary
  ))
}


#' run assess_stability() across every schema listed in an annotation table
#'
#' @param seurat_obj Seurat object
#' @param output_dir base output dir
#' @param annotation_table data frame w/ cols: schema, resolution,
#'   cluster, celltype. one row per (schema, cluster) pair. `schema` should
#'   match a meta.data col name (e.g. 'standard.rpca_snn_res.0.2'). optional
#'   `out_name` col overrides the output folder name -- useful when a
#'   differently-filtered Seurat object reuses the same schema/column name
#'   (e.g. a stringent-QC variant), which would otherwise collide w/ + overwrite
#'   the main run's output. defaults to `schema` if not provided.
#' @param jaccard_object optional path to a previously-saved raw jaccard
#'   results .qs2 (e.g. stability/<out_name>/scclusteval/jaccard_results.08_09.qs2)
#'   -- skips rerunning the expensive subsampling step, only redoes labeling/
#'   aggregation. only meaningful when annotation_table has exactly one schema.
#' @return named list of assess_stability() results, one per schema
run_stability_for_all_schemas <- function(seurat_obj, output_dir, annotation_table, jaccard_object = NULL)
{
  required_cols <- c('schema', 'resolution', 'cluster', 'celltype')
  missing_cols <- setdiff(required_cols, colnames(annotation_table))
  if (length(missing_cols) > 0)
  {
    stop('annotation_table missing required col(s): ', paste(missing_cols, collapse = ', '))
  }
  if (all(grepl('^[0-9.]+$', as.character(annotation_table$celltype))) &&
      !identical(as.character(annotation_table$celltype), as.character(annotation_table$cluster)))
  {
    # celltype column looks numeric but ISN'T just an intentional
    # --no-annotation run (where celltype == cluster by design) -- likely
    # the annotation file's celltypes col got swapped/misread upstream
    stop('annotation_table$celltype looks numeric, not real celltype labels: ',
         paste(utils::head(unique(annotation_table$celltype)), collapse = ', '),
         ' -- check the annotation file/config passed to run_validation.R.')
  }
  if (!('out_name' %in% colnames(annotation_table)))
  {
    annotation_table$out_name <- annotation_table$schema
  }

  results <- list()

  for (schema in unique(annotation_table$schema))
  {
    schema_ann <- annotation_table %>% filter(schema == !!schema)
    schema_resolution <- unique(schema_ann$resolution)
    schema_out_name <- unique(schema_ann$out_name)

    if (length(schema_resolution) != 1)
    {
      stop('schema "', schema, '" has inconsistent or missing resolution values ',
           'in annotation_table: ', paste(schema_resolution, collapse = ', '))
    }
    if (length(schema_out_name) != 1)
    {
      stop('schema "', schema, '" has inconsistent out_name values ',
           'in annotation_table: ', paste(schema_out_name, collapse = ', '))
    }

    if (!(schema %in% colnames(seurat_obj@meta.data)))
    {
      warning('schema "', schema, '" not found in seurat_obj meta.data, skipping.')
      next
    }

    cat('\n==== schema:', schema, '(out_name:', schema_out_name, ') ====\n')
    Idents(seurat_obj) <- seurat_obj@meta.data[[schema]]

    results[[schema_out_name]] <- assess_stability(
      seurat_obj = seurat_obj,
      output_dir = output_dir,
      out_name = schema_out_name,
      cluster_labels = schema_ann$celltype,
      cluster_numbers = schema_ann$cluster,
      resolution = schema_resolution,
      jaccard_object = jaccard_object
    )
  }

  return(results)
}


# ============================================================================
# example usage
# ============================================================================
# seurat_obj <- qs2::qs_read('path/to/analyzed_seurat_object.qs2')
# DefaultAssay(seurat_obj) <- 'RNA'
#
# # one row per (schema, cluster) -- same info previously duplicated as
# # separate _labels/_numbers vectors per schema
# annotation_table <- read.table('schema_annotations.txt', header = TRUE, sep = '\t')
# # expected cols: schema, resolution, cluster, celltype
# # e.g.:
# # schema                          resolution  cluster  celltype
# # standard.rpca_snn_res.0.2       0.2         0        Follicular
# # standard.rpca_snn_res.0.2       0.2         1        T_cells
# # ...
#
# all_results <- run_stability_for_all_schemas(
#   seurat_obj = seurat_obj,
#   output_dir = 'CLUSTER_STABILITY/PROJECT/',
#   annotation_table = annotation_table
# )
#
# # to rerun a single schema (e.g. after fixing an annotation), just filter
# # annotation_table to that schema + call run_stability_for_all_schemas()
# # again, or call assess_stability() directly w/ resolution set explicitly
