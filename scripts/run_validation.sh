#!/usr/bin/env bash
# run_validation.sh -- wrapper for SWANS validation scripts.
#
# only handles what must happen BEFORE R starts (can't depend on any R
# package, incl optparse): output-dir + lib-path creation, R_LIBS_USER
# export, dependency install. everything else (project, resolutions,
# steps, seurat-object, zscore-dir, annotations) read directly by
# 05_run_validation.R from config file -- bash doesn't duplicate that
# parsing logic.
#
#   bash run_validation.sh                     (uses validation_config.txt in cwd, if present)
#   bash run_validation.sh --config other.txt  (different config file)
#   bash run_validation.sh --install-only      (deps only, no pipeline run)
#   bash run_validation.sh --resolutions 0.1,0.2  (CLI flags override config, forwarded to R as-is)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v Rscript &> /dev/null; then
  echo "error: Rscript not found on PATH. install R + ensure Rscript available." >&2
  exit 1
fi

# ---- pre-scan args: bash-only concerns handled here ----
CONFIG_FILE="validation_config.txt"   # default location, no flag needed
CONFIG_FILE_EXPLICIT="false"
OUTPUT_DIR=""
LIB_PATH=""
INSTALL_ONLY="false"
FORWARD_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      if [[ $# -lt 2 ]]; then echo "error: --config requires value." >&2; exit 1; fi
      CONFIG_FILE="$2"; CONFIG_FILE_EXPLICIT="true"; shift 2 ;;
    --output-dir)
      if [[ $# -lt 2 ]]; then echo "error: --output-dir requires value." >&2; exit 1; fi
      OUTPUT_DIR="$2"; shift 2 ;;
    --lib-path)
      if [[ $# -lt 2 ]]; then echo "error: --lib-path requires value." >&2; exit 1; fi
      LIB_PATH="$2"; shift 2 ;;
    --install-only)
      INSTALL_ONLY="true"; shift 1 ;;
    *)
      FORWARD_ARGS+=("$1"); shift 1 ;;
  esac
done

if [[ "${CONFIG_FILE_EXPLICIT}" == "true" && ! -f "${CONFIG_FILE}" ]]; then
  echo "error: config file not found: ${CONFIG_FILE}" >&2
  exit 1
fi

# output-dir / lib-path: CLI flag > default. NOT config-file keys --
# bash owns + creates these dirs directly, always passes resolved value
# straight to R instead of letting it be set two different ways.
[[ -z "${OUTPUT_DIR}" ]] && OUTPUT_DIR="validation_output"
[[ -z "${LIB_PATH}" ]] && LIB_PATH="${OUTPUT_DIR}/r_libs"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LIB_PATH}"
LIB_PATH="$(cd "${LIB_PATH}" && pwd)"
export R_LIBS_USER="${LIB_PATH}${R_LIBS_USER:+:${R_LIBS_USER}}"

if [[ -f "${CONFIG_FILE}" ]]; then
  echo "config file: ${CONFIG_FILE}"
else
  echo "config file: none found (using CLI args / R defaults only)"
fi
echo "output dir:  ${OUTPUT_DIR}"
echo "R lib path:  ${LIB_PATH}"
echo ""

# ---- install missing deps ----
Rscript "${SCRIPT_DIR}/install_dependencies.R" --lib-path "${LIB_PATH}"

if [[ "${INSTALL_ONLY}" == "true" ]]; then
  echo ""
  echo "--install-only specified -- dependency install complete, skipping pipeline run."
  exit 0
fi

# bash owns output-dir -- always pass resolved value straight to R,
# regardless of --output-dir vs default.
FORWARD_ARGS+=("--output-dir" "${OUTPUT_DIR}")

# hand config file straight to R -- R parses project/resolutions/steps/
# seurat-object/zscore-dir/annotations itself, not duplicated here.
if [[ -f "${CONFIG_FILE}" ]]; then
  FORWARD_ARGS+=("--config-file" "${CONFIG_FILE}")
fi

# ---- run actual validation pipeline ----
Rscript "${SCRIPT_DIR}/05_run_validation.R" "${FORWARD_ARGS[@]}"
