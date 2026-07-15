# 00_match_clusters_across_schemas.R
#
# Automated cross-schema cluster matching via cell-membership overlap and
# graph-based community detection. This removes the need for manual annotation
# of every analysis schema before running the gene-robustness (Jaccard)
# comparison in 01_gene_robustness.R.
#
# Two clusters from different schemas are linked if they share a large
# fraction of the same cells (Jaccard index on cell barcodes). Linked clusters
# are grouped into "meta-clusters" via community detection, which can recover
# many-to-one relationships (e.g., three over-clustered fragments in one
# schema all corresponding to a single cluster in another).
#
# NOTE: this produces ANONYMOUS meta-cluster IDs (meta_cluster_1, meta_cluster_2,
# ...), not biologically-named cell types. Matching cell membership tells you
# which clusters correspond to each other, not what they biologically are. If
# biologically-named labels are needed, annotate each meta-cluster once
# (using marker genes from whichever schema is most trusted) after matching,
# rather than annotating every schema separately.
#
# Requires: Seurat, igraph, tidyverse

library(Seurat)
library(igraph)
library(tidyverse)


#' detect every clustering resolution in a Seurat object's schema cols
#' (SWANS naming pattern), no need to know resolutions up front. used to
#' scan full norm x integration x resolution grid for gene-signature
#' robustness instead of requiring one resolution up front.
#'
#' @param seurat_obj Seurat object
#' @return sorted numeric vector of unique resolutions found across all
#'   schema cols in meta.data
detect_all_resolutions <- function(seurat_obj)
{
  all_cols <- colnames(seurat_obj@meta.data)
  pattern <- '^(sct|standard)\\.(cca|harmony|rpca)_snn_res\\.[0-9.]+$'
  all_schema_cols <- all_cols[str_detect(all_cols, pattern)]

  if (length(all_schema_cols) == 0)
  {
    stop('no schema cols matching expected SWANS naming pattern ',
         '({normalization}.{integration}_snn_res.{resolution}) found ',
         'in meta.data. found cols: ', paste(all_cols, collapse = ', '))
  }

  detected_res <- str_extract(all_schema_cols, '(?<=res\\.)[0-9.]+$')
  sort(unique(as.double(detected_res)))
}


#' Auto-detect schema columns in a SWANS-analyzed Seurat object's meta.data,
#' matching the pattern {normalization}.{integration}_snn_res.{resolution}
#'
#' @param seurat_obj Seurat object
#' @param resolution REQUIRED. The clustering resolution to match on. Cross-schema
#'   matching assumes all input clusters come from the same clustering
#'   resolution -- mixing resolutions silently produces meaningless matches.
#' @return Character vector of matching meta.data column names
detect_schema_columns <- function(seurat_obj, resolution)
{
  if (missing(resolution) || is.null(resolution))
  {
    stop('`resolution` is required. Cross-schema cluster matching compares ',
         'clusters that must come from the same clustering resolution -- ',
         'mixing resolutions silently produces meaningless matches. ',
         'Specify the resolution used for the schemas you want to compare ',
         '(e.g., resolution = 0.2).')
  }

  all_cols <- colnames(seurat_obj@meta.data)

  # SWANS naming pattern: {normalization}.{integration}_snn_res.{resolution}
  # normalization in {sct, standard}; integration in {cca, harmony, rpca}
  pattern <- '^(sct|standard)\\.(cca|harmony|rpca)_snn_res\\.[0-9.]+$'
  all_schema_cols <- all_cols[str_detect(all_cols, pattern)]

  if (length(all_schema_cols) == 0)
  {
    stop('No schema columns matching the expected SWANS naming pattern ',
         '({normalization}.{integration}_snn_res.{resolution}) were found ',
         'in meta.data. Found columns: ', paste(all_cols, collapse = ', '))
  }

  detected_res <- str_extract(all_schema_cols, '(?<=res\\.)[0-9.]+$')
  schema_cols <- all_schema_cols[detected_res == as.character(resolution)]

  if (length(schema_cols) == 0)
  {
    stop('No schema columns found at resolution ', resolution,
         '. Detected resolutions in this object: ',
         paste(unique(detected_res), collapse = ', '))
  }

  if (length(schema_cols) < 2)
  {
    stop('Only one schema found at resolution ', resolution,
         '. Cross-schema matching requires at least two schemas to compare.')
  }

  return(schema_cols)
}


#' Match clusters across multiple schemas using cell-membership Jaccard overlap.
#'
#' @param seurat_obj Seurat object with cluster identities stored in meta.data,
#'   one column per schema (e.g., 'standard.rpca_snn_res.0.2').
#' @param schema_cols Character vector of meta.data column names to match across.
#'   Typically the output of detect_schema_columns().
#' @param min_jaccard Minimum cell-overlap Jaccard index for two clusters to be
#'   considered linked (default 0.3). Lower = more permissive matching, more
#'   clusters grouped together. This is dataset-dependent -- calibrate against
#'   a dataset with existing manual annotation before trusting the default on
#'   new data (see calibrate_against_manual_annotation() below).
#' @param community_resolution Resolution parameter passed to
#'   igraph::cluster_louvain (default 1.0; higher = more, smaller communities).
#'   Distinct from Seurat clustering `resolution` -- this controls community
#'   detection granularity on the cluster-overlap graph, not cell clustering.
#' @return A data frame: schema | cluster | meta_cluster
match_clusters_across_schemas <- function(seurat_obj, schema_cols, min_jaccard = 0.3, community_resolution = 1.0)
{
  meta <- seurat_obj@meta.data

  # 1. build cell-membership list: one entry per (schema, cluster) -> cell barcodes
  cluster_cells <- list()
  for (schema in schema_cols)
  {
    clusters <- unique(meta[[schema]])
    for (cl in clusters)
    {
      node_id <- paste0(schema, '::', cl)
      cluster_cells[[node_id]] <- rownames(meta)[meta[[schema]] == cl]
    }
  }

  node_ids <- names(cluster_cells)
  n <- length(node_ids)

  # 2. pairwise cell-overlap Jaccard between every cluster pair ACROSS schemas
  #    (skip same-schema pairs -- a cluster never needs to match itself)
  edges <- list()
  for (i in seq_len(n - 1))
  {
    schema_i <- str_split(node_ids[i], '::', simplify = TRUE)[1]
    for (j in seq((i + 1), n))
    {
      schema_j <- str_split(node_ids[j], '::', simplify = TRUE)[1]
      if (schema_i == schema_j) next  # only cross-schema links matter

      inter <- length(intersect(cluster_cells[[i]], cluster_cells[[j]]))
      uni <- length(union(cluster_cells[[i]], cluster_cells[[j]]))
      jac <- if (uni > 0) inter / uni else 0

      if (jac >= min_jaccard)
      {
        edges[[length(edges) + 1]] <- c(node_ids[i], node_ids[j], jac)
      }
    }
  }

  if (length(edges) == 0)
  {
    stop('No cluster pairs met the min_jaccard threshold (', min_jaccard,
         '). Try lowering min_jaccard, or check that schema_cols reference ',
         'the same underlying cells.')
  }

  edge_df <- as.data.frame(do.call(rbind, edges), stringsAsFactors = FALSE)
  colnames(edge_df) <- c('from', 'to', 'weight')
  edge_df$weight <- as.numeric(edge_df$weight)

  # 3. build graph, run community detection to find meta-clusters
  #    include ALL nodes as vertices (even unmatched ones) so no cluster
  #    silently disappears from the output
  g <- graph_from_data_frame(edge_df, directed = FALSE, vertices = node_ids)
  communities <- cluster_louvain(g, weights = E(g)$weight, resolution = community_resolution)

  # 4. any (schema, cluster) node with no edges above threshold becomes its own
  #    singleton meta-cluster -- flagged for review, not dropped
  membership_df <- data.frame(
    node_id = node_ids,
    meta_cluster_num = membership(communities)[node_ids],
    stringsAsFactors = FALSE
  )

  # 5. split node_id back into schema + cluster, format output, flag singletons
  result <- membership_df %>%
    separate(node_id, into = c('schema', 'cluster'), sep = '::') %>%
    mutate(meta_cluster = paste0('meta_cluster_', meta_cluster_num)) %>%
    group_by(meta_cluster) %>%
    mutate(is_singleton = n() == 1) %>%
    ungroup() %>%
    select(-meta_cluster_num) %>%
    arrange(meta_cluster, schema)

  n_singletons <- sum(result$is_singleton)
  if (n_singletons > 0)
  {
    warning(n_singletons, ' cluster(s) had no cross-schema match above ',
            'min_jaccard = ', min_jaccard, ' and were assigned singleton ',
            'meta-clusters. Review these before including them in downstream ',
            'comparisons -- filter(result, is_singleton) to inspect.')
  }

  return(result)
}


#' Sanity-check automated matching against an existing manual annotation.
#'
#' Compares the meta-clusters produced by match_clusters_across_schemas()
#' against a known-good manual annotation (schema, cluster, celltype), to help
#' calibrate min_jaccard and community_resolution before trusting the
#' automated output on unannotated data.
#'
#' @param matched_result Output of match_clusters_across_schemas()
#' @param manual_annotation Data frame with columns: schema, cluster, celltype
#' @return A data frame showing, for each meta-cluster, the distribution of
#'   manual cell type labels it contains. A well-calibrated matching should
#'   show each meta-cluster dominated by a single manual celltype.
calibrate_against_manual_annotation <- function(matched_result, manual_annotation)
{
  manual_annotation$cluster <- as.character(manual_annotation$cluster)

  comparison <- matched_result %>%
    left_join(manual_annotation, by = c('schema', 'cluster')) %>%
    count(meta_cluster, celltype, name = 'n_schema_cluster_pairs') %>%
    group_by(meta_cluster) %>%
    mutate(pct_of_meta_cluster = round(100 * n_schema_cluster_pairs / sum(n_schema_cluster_pairs), 1)) %>%
    arrange(meta_cluster, desc(n_schema_cluster_pairs)) %>%
    ungroup()

  n_mixed <- comparison %>%
    group_by(meta_cluster) %>%
    summarise(n_celltypes = n_distinct(celltype), .groups = 'drop') %>%
    filter(n_celltypes > 1) %>%
    nrow()

  if (n_mixed > 0)
  {
    warning(n_mixed, ' meta-cluster(s) contain more than one manually-annotated ',
            'cell type. Consider adjusting min_jaccard or community_resolution. ',
            'Inspect the returned data frame for details.')
  }

  return(comparison)
}


# ============================================================================
# Example usage
# ============================================================================
# seurat_obj <- qs2::qs_read('path/to/analyzed_seurat_object.qs2')
#
# schema_cols <- detect_schema_columns(seurat_obj, resolution = 0.2)
#
# matched_clusters <- match_clusters_across_schemas(
#   seurat_obj = seurat_obj,
#   schema_cols = schema_cols,
#   min_jaccard = 0.3,
#   community_resolution = 1.0
# )
#
# # optional: calibrate against existing manual annotation before trusting
# # the automated output on a new, unannotated dataset
# # manual_annotation <- data.frame(
# #   schema = ..., cluster = ..., celltype = ...
# # )
# # calibration <- calibrate_against_manual_annotation(matched_clusters, manual_annotation)
# # View(calibration)
#
# write.table(matched_clusters, 'matched_clusters.txt', sep = '\t', quote = FALSE, row.names = FALSE)