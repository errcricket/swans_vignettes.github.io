#!/usr/bin/env Rscript
# install_dependencies.R
#
# installs R packages required by SWANS multi-schema validation helper
# scripts (00-05) into user-specified lib dir, skips anything already
# installed + satisfying min version requirement. only Seurat has version
# floor enforced: validation scripts assume Seurat v5 object/metadata
# structure (e.g. schema-suffixed metadata cols from IntegrateLayers), NOT
# compatible w/ Seurat v4 or earlier, even if v4 already installed.
#
# standalone bootstrap script -- deliberately base R only (no optparse, no
# tidyverse), since optparse itself may not be installed yet when this
# needs to run.
#
# usage:
#   Rscript install_dependencies.R --lib-path /path/to/output_dir/r_libs
#
# all installs go into --lib-path only, nothing installed into system or
# default user lib. dir should also be exported as R_LIBS_USER by caller
# (run_validation.sh does this automatically) so 00-05 pick up packages
# installed here w/o any changes to their own library() calls.

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, args, default = NULL)
{
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop('missing value for ', flag)
  args[idx + 1]
}

lib_path <- get_arg('--lib-path', args)
if (is.null(lib_path))
{
  stop('--lib-path required, e.g.: Rscript install_dependencies.R --lib-path /path/to/output_dir/r_libs')
}

dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
lib_path <- normalizePath(lib_path, mustWork = TRUE)
cat('installing R package deps into:', lib_path, '\n\n')

# prepend so requireNamespace()/install.packages() see this location first,
# w/o disturbing rest of existing lib search path
.libPaths(c(lib_path, .libPaths()))

CRAN_REPO <- 'https://cloud.r-project.org'
MIN_SEURAT_VERSION <- '5.0.0'

is_installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

install_cran <- function(pkg)
{
  cat('installing', pkg, 'from CRAN into', lib_path, '...\n')
  install.packages(pkg, lib = lib_path, repos = CRAN_REPO)
}

report <- list()

# ---- Seurat: version-gated, v4 objects/metadata not compatible ----
cat('---- Seurat ----\n')
if (is_installed('Seurat'))
{
  current_version <- as.character(utils::packageVersion('Seurat'))
  if (utils::packageVersion('Seurat') >= MIN_SEURAT_VERSION)
  {
    cat('Seurat', current_version, 'already installed + satisfies >=', MIN_SEURAT_VERSION,
        '-- skipping.\n\n')
    report[['Seurat']] <- paste0('ok (', current_version, ', pre-existing)')
  } else
  {
    cat('found Seurat', current_version, '-- these scripts require >=', MIN_SEURAT_VERSION,
        '. installing compatible version into', lib_path,
        'so it takes precedence over existing (incompatible) install.\n')
    install_cran('Seurat')
    new_version <- tryCatch(as.character(utils::packageVersion('Seurat', lib.loc = lib_path)),
                             error = function(e) NA)
    report[['Seurat']] <- paste0('installed ', new_version, ' (shadowing existing ', current_version, ')')
  }
} else
{
  cat('Seurat not found -- installing.\n')
  install_cran('Seurat')
  new_version <- tryCatch(as.character(utils::packageVersion('Seurat', lib.loc = lib_path)),
                           error = function(e) NA)
  report[['Seurat']] <- paste0('installed ', new_version)
}
cat('\n')

# ---- standard CRAN deps (no version floor) ----
cran_packages <- c('optparse', 'tidyverse', 'igraph', 'qs2', 'ggplot2', 'patchwork')
for (pkg in cran_packages)
{
  cat('----', pkg, '----\n')
  if (is_installed(pkg))
  {
    cat(pkg, 'already installed -- skipping.\n\n')
    report[[pkg]] <- 'ok (pre-existing)'
  } else
  {
    install_cran(pkg)
    report[[pkg]] <- if (is_installed(pkg)) 'installed' else 'FAILED'
    cat('\n')
  }
}

# ---- scclusteval: github-only, needs "remotes" as bootstrap ----
cat('---- remotes (bootstrap for scclusteval) ----\n')
if (!is_installed('remotes'))
{
  install_cran('remotes')
}
cat('\n')

cat('---- scclusteval ----\n')
if (is_installed('scclusteval'))
{
  cat('scclusteval already installed -- skipping.\n\n')
  report[['scclusteval']] <- 'ok (pre-existing)'
} else
{
  cat('installing scclusteval from github (crazyhottommy/scclusteval)...\n')
  remotes::install_github('crazyhottommy/scclusteval', lib = lib_path, upgrade = 'never')
  report[['scclusteval']] <- if (is_installed('scclusteval')) 'installed' else 'FAILED'
  cat('\n')
}

# ---- summary ----
cat('==== dependency install summary ====\n')
for (pkg in names(report))
{
  cat(sprintf('%-15s %s\n', pkg, report[[pkg]]))
}

failed <- names(report)[grepl('FAILED', unlist(report))]
if (length(failed) > 0)
{
  cat('\nwarning: following pkg(s) failed to install:', paste(failed, collapse = ', '), '\n')
  quit(status = 1)
}

cat('\nall deps satisfied.\n')