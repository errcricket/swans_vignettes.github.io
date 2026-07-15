---
title: "Evaluating Schema Robustness and Cluster Stability in SWANS"
layout: default
---

# Evaluating Schema Robustness and Cluster Stability in SWANS

This vignette walks through a set of standalone helper scripts used to validate a SWANS
multi-schema scRNA-seq/snRNA-seq analysis. They are **not** part of the core Snakemake
pipeline and are not run automatically — you run them yourself, after SWANS has produced
an analyzed Seurat object and its accompanying marker z-score tables, to ask two separate
questions:

1. **Gene signature robustness** — across different normalization × integration schemas
   (e.g. SCTransform vs. standard log-normalization, combined with CCA/Harmony/RPCA
   integration), do independent analyses converge on the *same marker genes* for
   corresponding cell populations?
2. **Cluster stability** — for the *one* schema you intend to use going forward, how
   reproducible are its cluster boundaries under subsampling?

These are complementary but different validation questions, and the scripts keep them
separate: robustness scans the *full schema grid*, stability evaluates *one chosen schema*
at a time.

---

## 1. What's in the toolkit

| File (in `scripts` folder) | Role |
|---|---|
| `run_validation.sh` | Entry point. Handles everything that must happen *before* R starts: creates the output directory, sets up `R_LIBS_USER`, installs dependencies. Then hands off to `05_run_validation.R`. |
| `install_dependencies.R` | Installs required R packages (including `scclusteval` from GitHub) into a user-writable library path. |
| `validation_config.txt` | Key=value config file specifying dataset, resolutions, chosen schema, annotation (or lack of it), and which steps to run. Auto-detected in the current working directory. |
| `00_match_clusters_across_schemas.R` | Automated, annotation-free cross-schema cluster matching by cell-membership overlap. |
| `01_gene_robustness.R` | Gene-signature (marker) Jaccard concordance across matched clusters. |
| `02_cluster_stability_functions.R` | Core `scclusteval` subsampling machinery for a single schema. |
| `03_cluster_stability_usage.R` | Loops the stability analysis over one or more schemas listed in an annotation table. |
| `04_combine_and_plot.R` | Combines per-schema stability results and per-resolution robustness results into single tables/plots. |
| `05_run_validation.R` | Main driver. Reads the config file, applies CLI overrides, dispatches to the four steps. |

**Precedence:** CLI flag > config file value > built-in default.

---

## 2. Installing dependencies

```bash
# from the directory containing run_validation.sh
bash run_validation.sh --install-only
```

This installs (into a local library under your output directory, not system-wide):
`Seurat` (≥ 5.0.0 — the validation scripts assume Seurat v5 metadata structure and are
**not** compatible with v4 objects), `optparse`, `tidyverse`, `igraph`, `qs2`, `ggplot2`,
`patchwork`, and `scclusteval` (installed from `crazyhottommy/scclusteval` on GitHub).

`pdftools` is *not* installed by default — it's only needed if you later merge
per-resolution robustness PDF panels, and is loaded lazily at that point. If you hit a
missing-package error there, that's the one to install manually; it doesn't indicate the
base install failed.

---

## 3. The config file

```txt
# validation_config.txt
project=prjna790856

zscore_dir=../SWANS/data/endpoints/prjna790856/analysis/report/tables
seurat_object=../SWANS/data/endpoints/prjna790856/analysis/RDS/prjna790856_analyzed_seurat_object.qs2

resolution=0.2
normalization=sct
integration=harmony
#no_annotation=true
#annotation_file=stability_annotation.txt

steps=match,robustness,stability,combine
```

A few things worth noting about the keys:

- `resolutions` (**plural**) controls the **match/robustness** grid scan — it can be a
  comma-separated list, or omitted entirely to auto-detect every resolution present in
  the Seurat object's metadata.
- `resolution` / `normalization` / `integration` (**singular**) specify the **one**
  schema the **stability** step runs on — this is the schema you actually intend to
  publish or carry forward, not the full grid.
- `output_dir` is **not** a config key. `run_validation.sh` always owns and creates it,
  and passes it explicitly on the command line.
- `zscore_dir` and `seurat_object` can be omitted if you're running from the SWANS home
  directory for `project` — they'll be derived automatically as
  `data/endpoints/<project>/analysis/...`.

Run it with:

```bash
bash run_validation.sh
```

Or override anything at the command line, e.g.:

```bash
bash run_validation.sh --resolutions 0.1,0.2,0.3,0.4,0.5
```

---

## 4. Step 1 — cross-schema cluster matching (`match`)

Before you can ask "do these schemas agree on this cell type's markers," you first need
to know *which cluster in schema A corresponds to which cluster in schema B* — and
cluster numbers are not comparable across schemas on their own.

`00_match_clusters_across_schemas.R` solves this **without requiring manual annotation of
every schema**. For each pair of clusters from two different schemas, it computes the
Jaccard index of their *cell membership* (shared cell barcodes, not marker genes), links
cluster pairs above a minimum overlap (`--min-jaccard`, default 0.3), and runs Louvain
community detection on the resulting graph to group linked clusters into **meta-clusters**.

A few properties of this step worth knowing:

- Meta-cluster IDs (`meta_cluster_1`, `meta_cluster_2`, ...) are **anonymous** — matching
  cell membership tells you clusters correspond to each other, not what they biologically
  are. If you want biologically-named labels, annotate each meta-cluster once after
  matching (using the schema you trust most), rather than annotating every schema
  separately.
- Any cluster with no cross-schema link above `--min-jaccard` becomes its own singleton
  meta-cluster and is flagged (not silently dropped) — inspect these before trusting
  downstream comparisons.
- `--min-jaccard` and `--community-resolution` (Louvain granularity) are dataset-dependent.
  If you already have a manual annotation for this dataset, you can sanity-check the
  automated matching against it with `--calibration-annotation <file>` (columns: `schema`,
  `cluster`, `celltype`) — a well-calibrated matching should show each meta-cluster
  dominated by a single manual cell type.
- Matching requires clusters to come from the **same clustering resolution** — mixing
  resolutions produces meaningless matches, so `--resolutions`/auto-detection loops the
  whole match+robustness process one resolution at a time.

---

## 5. Step 2 — gene-signature robustness (`robustness`): how marker-based identity is established

This is where cluster *identity* — as opposed to cluster *correspondence* — enters the
picture, via marker genes.

For each schema, SWANS's own per-cluster marker output (Wilcoxon rank-sum z-scores,
`FindAllMarkers`-derived, one table per schema at `{project}_z_scores.<schema>.txt`) is
loaded and joined to the meta-cluster assignments from Step 1. Within each
(schema, meta-cluster) group:

1. Genes are ranked by mean z-score, descending (if more than one raw cluster within a
   schema maps to the same meta-cluster, their z-scores are averaged first).
2. The top *N* genes are retained (`--top-n-genes`, default 100) — this gene list **is**
   the marker signature standing in for that meta-cluster's identity in that schema.
3. Meta-clusters with no cross-schema match (singletons) are excluded by default, since a
   Jaccard comparison needs at least two schemas to be meaningful.

Then, for every meta-cluster, pairwise Jaccard similarity of these top-*N* gene lists is
computed across every pair of schemas. **This Jaccard index is the robustness metric**:
if two schemas independently arrive at largely the same top marker genes for cells they
already agree belong together (from Step 1's cell-membership matching), that cell
population's identity is being called consistently regardless of normalization/integration
choice. Low gene-signature Jaccard for a given meta-cluster — even when cell-membership
matching was clean — is a sign that *what the cluster is* is schema-sensitive, distinct
from whether the cluster boundary itself is stable (that's Step 3).

Outputs per resolution (`<output-dir>/gene_robustness/res_<resolution>/`):
- `matched_clusters.txt` — Step 1's schema/cluster/meta-cluster table
- `<meta_cluster>_gene_jaccard.txt` / `.pdf` — per-meta-cluster Jaccard matrix + heatmap
- `all_metacluster_gene_jaccard.txt` — combined long-format table
- `all_metaclusters_gene_jaccard_panel.pdf` — all meta-clusters in one panel figure

Because match/robustness scan every resolution present (or an explicit `--resolutions`
list), results are also auto-collated across resolutions into
`gene_robustness/all_resolutions_gene_jaccard.txt` plus a boxplot summarizing
within-normalization vs. cross-normalization concordance by resolution — useful for
narrowing down a resolution/schema before committing to one for the stability step.

---

## 6. Step 3 — cluster stability (`stability`): with or without annotation

Stability asks a different question from robustness: for **one schema**, if you
repeatedly subsample the data and recluster, how consistently do the same cells end up
grouped together? This uses `scclusteval`: cells are subsampled at two rates (default 0.8,
0.9), 10 replicates each, reclustered at the *same resolution* as the original clustering,
and compared back to the full-data clustering via Jaccard index.

The schema is specified via `--resolution` + `--normalization` + `--integration` (or the
equivalent config keys) — exactly one combination, not a grid.

### If you have an annotation file

Supply `--annotation-file path/to/file.txt` (or `annotation_file=` in the config). This is
the **same format as SWANS's `CLUSTER_ANNOTATION_FILE`** — two columns, `cluster` and
`celltypes`, no schema or resolution columns needed (those come from `--resolution`/
`--normalization`/`--integration`):

```txt
cluster	celltypes
0	Follicular
1	T_cells
2	Myeloid
```

Raw per-cluster Jaccard stability values are then aggregated up to the cell-type level.
You only need to annotate the **one** schema you're actually running stability on — not
every schema in the comparison grid.

### If you don't have an annotation file

Pass `--no-annotation` instead. Stability is then computed **per raw cluster number**,
with no cell-type labels at all — each cluster is its own group, and no clusters are
combined under a shared label. This is the right choice for:
- demos/vignettes where you don't want to require a completed annotation,
- sanity-checking raw cluster reproducibility independent of any annotation call.

**`--no-annotation` and `--annotation-file` must never both be set.** The scripts throw a
hard error if both are supplied — this guards against a real failure mode where a stale
`no_annotation=true` left in a config file silently overrode an explicit
`--annotation-file`, silently producing raw-cluster-number output (uncombined, unlabeled)
across every schema with no warning. If you see this error, check
`validation_config.txt` for a leftover `no_annotation=true` line before assuming a bug in
your command.

### Other flags worth knowing

- `--out-name` overrides the output subfolder name for this run. Needed when a
  differently-filtered Seurat object (e.g. a stringent-QC variant) reuses the same
  normalization/integration/resolution schema name — without it, its output would silently
  overwrite the main run's.
- `--jaccard-object path/to/jaccard_results.08_09.qs2` reuses a previously computed (and
  expensive) subsampling result, redoing only the labeling/aggregation step — useful after
  fixing an annotation-only bug, without repeating the actual `scclusteval` computation.
- Existing results trigger an overwrite prompt; `--force-rerun` skips it (needed for
  non-interactive/scripted runs). Supplying `--jaccard-object` also skips the prompt
  automatically, since reusing cached results is itself an explicit rerun request.

Outputs land at `<output-dir>/stability/<out_name>/scclusteval/stability/`, with
`tables/celltype_stability.csv`, `tables/cluster_stability.csv`, and
`figures/stability_heatmap.pdf` / `figures/celltype_stability_barplot.pdf`.

---

## 7. Step 4 — combining across schemas (`combine`)

Once stability has been run for each schema you want to compare (each with its own
`--out-name` if applicable), the `combine` step:

- discovers every schema subfolder under `<output-dir>/stability/`,
- reads each `celltype_stability.csv`,
- aggregates duplicate cell types within a schema (in case an annotation collapsed
  multiple raw clusters into one label),
- writes `stability/combined_celltype_stability.csv`,
- and produces a faceted comparison plot, `stability/combined_stability_facets.pdf`, with
  one panel per schema and reference lines at `--stability-threshold` (default 0.75,
  "stable") and a secondary "moderately stable" line at 0.60.

The equivalent collation for gene-signature robustness across resolutions happens
automatically at the end of the `robustness` step (Section 5) — there's no separate
`combine` call needed for that half.

---

## 8. Putting it together: a full run

```bash
# 1. one-time dependency install
bash run_validation.sh --install-only

# 2. scan the full schema x resolution grid for cross-schema cluster matching
#    and gene-signature robustness
bash run_validation.sh --steps match,robustness --resolutions 0.1,0.2,0.3

# 3. run stability on each of six schemas at the chosen resolution, one
#    annotation file per schema
bash run_validation.sh --steps stability --resolution 0.2 --normalization sct      --integration harmony  --annotation-file stability_annotations/ann_sct_harmony.txt
bash run_validation.sh --steps stability --resolution 0.2 --normalization sct      --integration cca       --annotation-file stability_annotations/ann_sct_cca.txt
bash run_validation.sh --steps stability --resolution 0.2 --normalization sct      --integration rpca      --annotation-file stability_annotations/ann_sct_rpca.txt
bash run_validation.sh --steps stability --resolution 0.2 --normalization standard --integration cca       --annotation-file stability_annotations/ann_standard_cca.txt
bash run_validation.sh --steps stability --resolution 0.2 --normalization standard --integration harmony   --annotation-file stability_annotations/ann_standard_harmony.txt
bash run_validation.sh --steps stability --resolution 0.2 --normalization standard --integration rpca      --annotation-file stability_annotations/ann_standard_rpca.txt

# 3b. alternative: no annotation available yet -- run per raw cluster number
bash run_validation.sh --steps stability --resolution 0.2 --normalization sct --integration harmony --no-annotation

# 4. combine all six schemas' stability results into one comparison plot
bash run_validation.sh --steps combine
```

---

## 9. A caveat worth knowing before re-running published numbers

Re-running the stability analysis on the same schema, even with matched software and
package versions, can produce **numerically different** results from a previous run. This
traces to an upstream, unresolved Seurat issue in which `FindClusters()` can produce
non-reproducible cluster numbering between runs, which propagates into the raw
`scclusteval` subsampling matrices. This is not a bug in these validation scripts — if a
re-run doesn't exactly match a previously published figure, that discrepancy is the
expected symptom of this upstream issue, not necessarily a sign of a broken setup. For
internal consistency, all schemas being compared should be (re-)run together, in the same
session, when this matters.

---

## 10. Output directory structure

```
validation_output/
├── gene_robustness/
│   ├── res_0.1/
│   │   ├── matched_clusters.txt
│   │   ├── meta_cluster_1_gene_jaccard.txt
│   │   ├── meta_cluster_1_gene_jaccard.pdf
│   │   ├── all_metacluster_gene_jaccard.txt
│   │   └── all_metaclusters_gene_jaccard_panel.pdf
│   ├── res_0.2/
│   │   └── ...
│   ├── all_resolutions_gene_jaccard.txt
│   └── all_resolutions_gene_jaccard_summary.pdf
└── stability/
    ├── sct.harmony_snn_res.0.2/
    │   └── scclusteval/
    │       ├── jaccard_results.08_09.qs2
    │       └── stability/
    │           ├── figures/
    │           │   ├── stability_heatmap.pdf
    │           │   └── celltype_stability_barplot.pdf
    │           └── tables/
    │               ├── cluster_stability.csv
    │               └── celltype_stability.csv
    ├── standard.rpca_snn_res.0.2
    │   └── ...
    ├── combined_celltype_stability.csv
    └── combined_stability_facets.pdf
```
