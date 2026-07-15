#!/usr/bin/env Rscript
# run_validation.R
#
# driver script for SWANS multi-schema validation helper scripts. parses
# user params (seurat object path, output dir, etc), runs requested
# step(s): cross-schema cluster matching, gene-signature robustness,
# cluster stability, and/or combining stability results across schemas.
#
# not part of core SWANS pipeline. run after SWANS produces analyzed
# Seurat object + z-score tables for dataset to validate.
#
# Requires: optparse, Seurat, qs2, igraph, tidyverse, ggplot2, patchwork, scclusteval
# (pdftools optional -- only needed for combine_gene_robustness_across_resolutions())
#
# config file (--config-file, or auto-passed by run_validation.sh if
# validation_config.txt exists in cwd): key=value lines, # for comments.
# recognized keys: project, resolutions, steps, seurat_object, zscore_dir,
# calibration_annotation, resolution, normalization, integration,
# annotation_file, no_annotation (true/false), out_name. output_dir NOT a config key -- run_validation.sh owns
# + creates that dir directly, always passes via --output-dir. any value
# ALSO given as CLI flag takes precedence over config file; config file
# takes precedence over built-in defaults. project lowercased on read.
#
# match/robustness do NOT require single resolution: scan every
# resolution present in Seurat object automatically (or specific
# comma-separated list via --resolutions/resolutions=), producing
# <output-dir>/gene_robustness/res_<resolution>/ (matched_clusters.txt +
# jaccard panels together) per resolution. after loop, results
# auto-collate into <output-dir>/gene_robustness/all_resolutions_gene_jaccard.txt
# (one long-format table, all resolutions) +
# all_resolutions_gene_jaccard_summary.pdf (boxplot, cross-schema Jaccard
# by resolution). assesses full norm x integration x resolution grid,
# useful for narrowing down schema before annotating + running stability
# on just the one chosen.
#
# stability runs on exactly ONE schema, built from --resolution +
# --normalization + --integration (or config equivalents), paired w/
# --annotation-file -- same format as SWANS CLUSTER_ANNOTATION_FILE
# (cluster, celltypes -- no schema/resolution cols needed, unlike
# --calibration-annotation for match, which is a separate bespoke format).
#
# usage (see also run_validation.sh for bash entry point):
#   Rscript run_validation.R --config-file validation_config.txt
#   Rscript run_validation.R \
#     --seurat-object path/to/analyzed_seurat_object.qs2 \
#     --project prjna790856 \
#     --zscore-dir path/to/report/tables \
#     --output-dir path/to/validation_output \
#     --steps match,robustness \
#     --resolutions 0.1,0.2,0.3

suppressPackageStartupMessages({
  library(optparse)
})

script_dir <- dirname(sub('--file=', '', grep('--file=', commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == '') script_dir <- '.'

option_list <- list(
  make_option('--config-file', type = 'character', default = NULL,
              help = 'path to key=value config file. CLI flags below override anything set here.'),
  make_option('--seurat-object', type = 'character', default = NULL,
              help = 'path to analyzed Seurat object (.qs2). required for match/stability steps.'),
  make_option('--resolutions', type = 'character', default = NULL,
              help = 'comma-separated resolutions to scan for match/robustness (e.g. 0.1,0.2,0.3). default: auto-detect every resolution present in Seurat object.'),
  make_option('--project', type = 'character', default = NULL,
              help = 'SWANS project name, used for z-score file naming. required for robustness step.'),
  make_option('--zscore-dir', type = 'character', default = NULL,
              help = 'dir containing SWANS z-score tables. required for robustness step.'),
  make_option('--output-dir', type = 'character', default = NULL,
              help = 'dir to write all outputs to. default: validation_output'),
  make_option('--steps', type = 'character', default = NULL,
              help = 'comma-separated steps to run: match, robustness, stability, combine. default: match,robustness,stability,combine'),
  make_option('--min-jaccard', type = 'double', default = 0.3,
              help = 'min cell-overlap Jaccard for cross-schema cluster matching. default: %default'),
  make_option('--community-resolution', type = 'double', default = 1.0,
              help = 'community detection resolution for cluster matching (distinct from clustering resolution). default: %default'),
  make_option('--top-n-genes', type = 'integer', default = 100,
              help = 'num top genes per meta-cluster for robustness comparison. default: %default'),
  make_option('--calibration-annotation', type = 'character', default = NULL,
              help = 'optional: path to existing manual annotation file (cols: schema, cluster, celltype) to sanity-check automated matching against.'),
  make_option('--resolution', type = 'double', default = NULL,
              help = 'resolution of the ONE chosen schema to run stability on (distinct from --resolutions, the robustness grid scan). required for stability step.'),
  make_option('--normalization', type = 'character', default = NULL,
              help = 'normalization method of chosen schema (sct or standard). required for stability step.'),
  make_option('--integration', type = 'character', default = NULL,
              help = 'integration method of chosen schema (cca, harmony, or rpca). required for stability step.'),
  make_option('--annotation-file', type = 'character', default = NULL,
              help = 'plain cluster/celltypes annotation file, same format as SWANS CLUSTER_ANNOTATION_FILE (cols: cluster, celltypes), for the chosen schema. required for stability step, unless --no-annotation is set.'),
  make_option('--no-annotation', action = 'store_true', default = FALSE,
              help = 'run stability per raw cluster number instead of per celltype -- no annotation file needed. each cluster is its own group regardless of whether it would otherwise be combined w/ another cluster sharing a celltype label. useful for vignettes/demos where cluster identity doesn\'t need to be resolved.'),
  make_option('--out-name', type = 'character', default = NULL,
              help = 'override output folder name for stability results (default: same as schema, i.e. {normalization}.{integration}_snn_res.{resolution}). needed when a differently-filtered Seurat object reuses the same schema name (e.g. a stringent-QC variant) -- otherwise its output silently overwrites the main run\'s.'),
  make_option('--jaccard-object', type = 'character', default = NULL,
              help = 'path to a previously-saved raw jaccard results .qs2 (e.g. stability/<out_name>/scclusteval/jaccard_results.08_09.qs2) -- skips rerunning the expensive subsampling step, only redoes labeling/aggregation. use after fixing a labeling-only bug w/o repeating the actual scclusteval computation.'),
  make_option('--stability-threshold', type = 'double', default = 0.75,
              help = 'Jaccard threshold above which a cluster is considered stable. default: %default'),
  make_option('--force-rerun', action = 'store_true', default = FALSE,
              help = 'skip the re-run prompt when existing collated robustness results are found -- always re-run match/robustness.')
)

opt_parser <- OptionParser(option_list = option_list,
                            description = 'run SWANS multi-schema validation steps (cluster matching, gene-signature robustness, cluster stability).')
opt <- parse_args(opt_parser)

# ---- read config file, if given (key=value, # for comments) ----
read_config_file <- function(path)
{
  lines <- readLines(path, warn = FALSE)
  config <- list()
  for (line in lines)
  {
    line <- trimws(line)
    if (line == '' || startsWith(line, '#')) next
    eq_pos <- regexpr('=', line, fixed = TRUE)
    if (eq_pos < 0) next
    key <- trimws(substr(line, 1, eq_pos - 1))
    value <- trimws(substr(line, eq_pos + 1, nchar(line)))
    if (!nzchar(value)) next
    if (!is.null(config[[key]]))
    {
      warning('config key "', key, '" set more than once in ', path,
              ' -- using last value ("', value, '"), ignoring earlier value ("', config[[key]], '").')
    }
    config[[key]] <- value
  }
  config
}

config <- list()
if (!is.null(opt$`config-file`))
{
  if (!file.exists(opt$`config-file`))
  {
    stop('config file not found: ', opt$`config-file`)
  }
  cat('reading config file:', opt$`config-file`, '\n\n')
  config <- read_config_file(opt$`config-file`)
  if (!is.null(config$project)) config$project <- tolower(config$project)
}

# CLI flag > config file > built-in default. only applied when CLI value
# still NULL (i.e. user didn't pass that flag explicitly)
apply_config <- function(opt_value, config_key, config, default = NULL)
{
  if (!is.null(opt_value)) return(opt_value)
  if (!is.null(config[[config_key]])) return(config[[config_key]])
  default
}

opt$project                    <- apply_config(opt$project, 'project', config)
opt$resolutions                <- apply_config(opt$resolutions, 'resolutions', config)
if (is.null(opt$`output-dir`)) opt$`output-dir` <- 'validation_output'  # only when running standalone, w/o run_validation.sh
opt$steps                      <- apply_config(opt$steps, 'steps', config, default = 'match,robustness,stability,combine')
opt$`seurat-object`            <- apply_config(opt$`seurat-object`, 'seurat_object', config)
opt$`zscore-dir`               <- apply_config(opt$`zscore-dir`, 'zscore_dir', config)
opt$`calibration-annotation`   <- apply_config(opt$`calibration-annotation`, 'calibration_annotation', config)
opt$resolution                 <- apply_config(opt$resolution, 'resolution', config)
if (!is.null(opt$resolution) && !is.numeric(opt$resolution)) opt$resolution <- as.double(opt$resolution)
opt$normalization               <- apply_config(opt$normalization, 'normalization', config)
opt$integration                 <- apply_config(opt$integration, 'integration', config)
opt$`annotation-file`           <- apply_config(opt$`annotation-file`, 'annotation_file', config)
opt$`out-name`                  <- apply_config(opt$`out-name`, 'out_name', config)
if (!isTRUE(opt$`no-annotation`) && !is.null(config$no_annotation))
{
  opt$`no-annotation` <- tolower(trimws(config$no_annotation)) %in% c('true', 'yes', '1')
}

# derive seurat-object / zscore-dir from project if still unset (assumes
# cwd = SWANS home dir: data/endpoints/<project>/analysis/...)
if (is.null(opt$`seurat-object`) && !is.null(opt$project))
{
  opt$`seurat-object` <- file.path('data/endpoints', opt$project, 'analysis/RDS',
                                    paste0(opt$project, '_analyzed_seurat_object.qs2'))
}
if (is.null(opt$`zscore-dir`) && !is.null(opt$project))
{
  opt$`zscore-dir` <- file.path('data/endpoints', opt$project, 'analysis/report/tables')
}

steps <- trimws(strsplit(opt$steps, ',')[[1]])
valid_steps <- c('match', 'robustness', 'stability', 'combine')
bad_steps <- setdiff(steps, valid_steps)
if (length(bad_steps) > 0)
{
  stop('unrecognized step(s): ', paste(bad_steps, collapse = ', '),
       '. valid steps: ', paste(valid_steps, collapse = ', '))
}

dir.create(opt$`output-dir`, showWarnings = FALSE, recursive = TRUE)

cat('==== SWANS multi-schema validation ====\n')
cat('steps to run:', paste(steps, collapse = ', '), '\n')
cat('output dir:', opt$`output-dir`, '\n\n')

source(file.path(script_dir, '00_match_clusters_across_schemas.R'))
source(file.path(script_dir, '01_gene_robustness.R'))
source(file.path(script_dir, '02_cluster_stability_functions.R'))
source(file.path(script_dir, '03_cluster_stability_usage.R'))
source(file.path(script_dir, '04_combine_and_plot.R'))

seurat_obj <- NULL

need_seurat_obj <- any(c('match', 'robustness', 'stability') %in% steps)
if (need_seurat_obj)
{
  if (is.null(opt$`seurat-object`))
  {
    stop('--seurat-object required for "match" + "stability" steps.')
  }
  cat('loading Seurat object from', opt$`seurat-object`, '...\n')
  if (!requireNamespace('qs2', quietly = TRUE))
  {
    stop('package "qs2" required to load Seurat object. if object was ',
         'saved w/ older "qs" package, convert first (see README.md).')
  }
  seurat_obj <- qs2::qs_read(opt$`seurat-object`)
  cat('loaded:', ncol(seurat_obj), 'cells\n\n')
}

#' checks whether a step's output already exists + prompts to re-run if so.
#' shared by match/robustness + stability, so both steps skip cleanly when
#' existing results are present + user declines to overwrite them.
#'
#' @param existing_path file/dir whose existence signals "already ran"
#' @param step_label human-readable step name for prompt text
#' @param force_rerun if TRUE, skip prompt + always return TRUE (re-run)
#' @return TRUE if step should run, FALSE if existing results should be kept
should_run_step <- function(existing_path, step_label, force_rerun)
{
  if (!file.exists(existing_path)) return(TRUE)

  cat('existing', step_label, 'results found at\n  ', existing_path, '\n')

  if (isTRUE(force_rerun))
  {
    cat('--force-rerun set -- re-running', step_label, ', existing results will be overwritten.\n\n')
    return(TRUE)
  }

  cat('re-run', step_label, '? existing results will be overwritten. [y/N]: ')
  answer <- tryCatch(tolower(trimws(readLines(con = 'stdin', n = 1))), error = function(e) '')
  if (identical(answer, 'y') || identical(answer, 'yes'))
  {
    cat('\n')
    return(TRUE)
  }

  cat('skipping', step_label, ', keeping existing results. (pass --force-rerun to skip this prompt)\n\n')
  return(FALSE)
}


# ---- steps: match + robustness (looped across resolutions) ----
run_match_robustness_now <- any(c('match', 'robustness') %in% steps)

if (run_match_robustness_now)
{
  existing_combined_txt <- file.path(opt$`output-dir`, 'gene_robustness', 'all_resolutions_gene_jaccard.txt')
  run_match_robustness_now <- should_run_step(existing_combined_txt, 'match/robustness', opt$`force-rerun`)
}

if (run_match_robustness_now)
{
  if (!is.null(opt$resolutions))
  {
    resolutions_to_run <- as.double(trimws(strsplit(opt$resolutions, ',')[[1]]))
    cat('resolutions to scan (from --resolutions/config):', paste(resolutions_to_run, collapse = ', '), '\n\n')
  } else
  {
    resolutions_to_run <- detect_all_resolutions(seurat_obj)
    cat('no --resolutions specified -- auto-detected', length(resolutions_to_run),
        'resolution(s) in Seurat object:', paste(resolutions_to_run, collapse = ', '), '\n\n')
  }

  for (res in resolutions_to_run)
  {
    cat('==== resolution', res, '====\n')
    res_dir <- file.path(opt$`output-dir`, 'gene_robustness', paste0('res_', res))
    dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

    matched_clusters <- NULL

    if ('match' %in% steps)
    {
      cat('---- step: cross-schema cluster matching ----\n')
      schema_cols <- detect_schema_columns(seurat_obj, resolution = res)
      cat('detected schema cols:', paste(schema_cols, collapse = ', '), '\n')

      matched_clusters <- match_clusters_across_schemas(
        seurat_obj = seurat_obj,
        schema_cols = schema_cols,
        min_jaccard = opt$`min-jaccard`,
        community_resolution = opt$`community-resolution`
      )

      match_out <- file.path(res_dir, 'matched_clusters.txt')
      write.table(matched_clusters, match_out, sep = '\t', quote = FALSE, row.names = FALSE)
      cat('matched clusters written to', match_out, '\n')

      if (!is.null(opt$`calibration-annotation`))
      {
        cat('calibrating against manual annotation at', opt$`calibration-annotation`, '...\n')
        manual_annotation <- read.table(opt$`calibration-annotation`, header = TRUE, sep = '\t')
        calibration <- calibrate_against_manual_annotation(matched_clusters, manual_annotation)
        calib_out <- file.path(res_dir, 'calibration_check.txt')
        write.table(calibration, calib_out, sep = '\t', quote = FALSE, row.names = FALSE)
        cat('calibration results written to', calib_out, '-- review before trusting automated matching.\n')
      }
      cat('\n')
    }

    if ('robustness' %in% steps)
    {
      cat('---- step: gene-signature robustness ----\n')
      if (is.null(opt$project) || is.null(opt$`zscore-dir`))
      {
        stop('--project + --zscore-dir required for "robustness" step.')
      }
      if (is.null(matched_clusters))
      {
        stop('"robustness" step requires "match" step to have run first ',
             '(include "match" in --steps for this resolution).')
      }

      schema_short <- unique(sub('_snn_res\\.[0-9.]+$', '', matched_clusters$schema))
      cat('schemas for z-score loading:', paste(schema_short, collapse = ', '), '\n')

      zscore_data <- load_zscore_tables(
        z_score_dir = opt$`zscore-dir`,
        project = opt$project,
        schema_cols = schema_short,
        resolution = res
      )

      aggregated <- aggregate_genes_by_metacluster(zscore_data, matched_clusters, top_n = opt$`top-n-genes`)

      robustness_out <- res_dir
      results <- compute_gene_jaccard(aggregated, out_dir = robustness_out)
      cat('gene robustness results written to', robustness_out, '\n')
    }
    cat('\n')
  }

  if ('robustness' %in% steps)
  {
    cat('---- step: collate gene robustness across resolutions ----\n')
    gene_robustness_dir <- file.path(opt$`output-dir`, 'gene_robustness')
    collated <- combine_gene_robustness_across_resolutions(gene_robustness_dir)
    cat('combined table written to', collated$combined_txt, '\n')
    if (!is.null(collated$combined_pdf)) cat('combined panel PDF written to', collated$combined_pdf, '\n')
    cat('\n')
  }
}

# ---- step: stability ----
if ('stability' %in% steps)
{
  cat('---- step: cluster stability ----\n')
  missing_stability_opts <- c()
  if (is.null(opt$normalization)) missing_stability_opts <- c(missing_stability_opts, '--normalization')
  if (is.null(opt$integration)) missing_stability_opts <- c(missing_stability_opts, '--integration')
  if (is.null(opt$resolution)) missing_stability_opts <- c(missing_stability_opts, '--resolution')
  if (!isTRUE(opt$`no-annotation`) && is.null(opt$`annotation-file`))
  {
    missing_stability_opts <- c(missing_stability_opts, '--annotation-file (or set --no-annotation)')
  }
  if (length(missing_stability_opts) > 0)
  {
    stop('"stability" step requires: ', paste(missing_stability_opts, collapse = ', '))
  }
  if (isTRUE(opt$`no-annotation`) && !is.null(opt$`annotation-file`))
  {
    # both set at once is almost always a stale config leftover (e.g.
    # no_annotation=true left over from an earlier vignette/demo run)
    # silently overriding a real --annotation-file -- this exact silent
    # contradiction previously caused every schema to run unannotated
    # w/o anyone noticing until the raw CSV was inspected. fail loudly.
    stop('both --no-annotation and --annotation-file are set -- this is ',
         'contradictory + usually means a stale no_annotation=true is ',
         'still in your config file, silently overriding --annotation-file. ',
         'fix validation_config.txt (remove/flip no_annotation) or drop ',
         'one of the two flags before rerunning.')
  }

  schema <- paste0(opt$normalization, '.', opt$integration, '_snn_res.', opt$resolution)
  out_name <- if (!is.null(opt$`out-name`)) opt$`out-name` else schema
  cat('stability schema:', schema, ' out_name:', out_name, '\n')

  if (isTRUE(opt$`no-annotation`))
  {
    # each cluster is its own group, labeled by its own number -- no
    # celltype annotation needed, no combining of clusters that would
    # otherwise share a celltype label. useful when cluster identity
    # itself isn't the point (e.g. demos, vignettes, or sanity-checking
    # raw cluster reproducibility independent of annotation).
    if (!(schema %in% colnames(seurat_obj@meta.data)))
    {
      stop('schema "', schema, '" not found in seurat_obj meta.data.')
    }
    cluster_ids <- sort(unique(seurat_obj@meta.data[[schema]]))
    cat('--no-annotation set -- running stability per raw cluster number (',
        length(cluster_ids), 'clusters), skipping celltype annotation.\n')
    annotation_table <- data.frame(
      schema = schema,
      resolution = opt$resolution,
      cluster = as.character(cluster_ids),
      celltype = as.character(cluster_ids),
      out_name = out_name
    )
  } else
  {
    ann <- read.table(opt$`annotation-file`, header = TRUE, sep = '\t')
    required_ann_cols <- c('cluster', 'celltypes')
    missing_ann_cols <- setdiff(required_ann_cols, colnames(ann))
    if (length(missing_ann_cols) > 0)
    {
      stop('--annotation-file missing required col(s): ', paste(missing_ann_cols, collapse = ', '),
           '. expects same format as SWANS CLUSTER_ANNOTATION_FILE (cluster, celltypes) -- schema + resolution come from config, not the file.')
    }

    annotation_table <- data.frame(
      schema = schema,
      resolution = opt$resolution,
      cluster = as.character(ann$cluster),
      celltype = ann$celltypes,
      out_name = out_name
    )
  }

  stability_out <- file.path(opt$`output-dir`, 'stability')
  existing_celltype_csv <- file.path(stability_out, out_name, 'scclusteval', 'stability', 'tables', 'celltype_stability.csv')

  # providing --jaccard-object means "redo labeling only, reusing cached
  # subsampling results" -- that's an explicit rerun request, skip the prompt
  force_rerun_effective <- isTRUE(opt$`force-rerun`) || !is.null(opt$`jaccard-object`)

  if (should_run_step(existing_celltype_csv, 'stability', force_rerun_effective))
  {
    all_results <- run_stability_for_all_schemas(
      seurat_obj = seurat_obj,
      output_dir = stability_out,
      annotation_table = annotation_table,
      jaccard_object = opt$`jaccard-object`
    )
    cat('stability results written to', stability_out, '\n\n')
  }
}

# ---- step: combine ----
if ('combine' %in% steps)
{
  cat('---- step: combine stability results ----\n')
  stability_out <- file.path(opt$`output-dir`, 'stability')
  if (!dir.exists(stability_out))
  {
    stop('no stability output found at ', stability_out, '. run "stability" step first.')
  }

  combined <- combine_stability_results(stability_out)
  combined_csv <- file.path(stability_out, 'combined_celltype_stability.csv')
  write.csv(combined, combined_csv, row.names = FALSE)

  plot_out <- file.path(stability_out, 'combined_stability_facets.pdf')
  plot_combined_stability(combined, out_path = plot_out, stability_threshold = opt$`stability-threshold`)
  cat('combined results written to', combined_csv, 'and', plot_out, '\n\n')
}

cat('==== done ====\n')
