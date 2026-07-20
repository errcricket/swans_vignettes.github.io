# SWANS Tutorial: Post-Annotation Analysis
### Worked example — GSE191288 / PRJNA790856 (Wang et al. 2022, papillary thyroid carcinoma)

> **Scope of this document.** This vignette continues directly from **SWANS
> Vignette 1: Preliminary Analysis**. It assumes you have already:
> - run the preliminary pipeline pass
> - reviewed the interactive report
> - chosen one final normalization method, integration method, and resolution
>
> This document covers annotating that chosen schema, running the
> post-annotation analysis pass, and interpreting the resulting report —
> including differential expression across experimental conditions, pathway
> analysis, z-scores, and optional trajectory analysis.
>
> **A naming note.** SWANS's output files, rule names, and the
> `FinalSnakefile` itself use `final` in their naming (`final_report.html`,
> `final_analyzed_seurat_object`, etc.). This document and the SWANS
> manuscript refer to this analysis phase as **post-annotation** analysis.
> The two terms refer to the same thing — don't be thrown by the mismatch
> when you see `final_` in actual output paths below.

---

## 1. What you need before starting

- Everything from Vignette 1, completed: `samples.sample_list`,
  `configs/prelim_configs.yaml`, and a completed prelim run producing
  `PROJECT_analyzed_seurat_object.qs2`
- `configs/post_annotation_configs.yaml`
- A completed **`CLUSTER_ANNOTATION_FILE`** (built in Section 2 below)

---

## 2. Building the `CLUSTER_ANNOTATION_FILE`

This file assigns a cell type label to every cluster number in your chosen
schema. It can have almost any filename, with one restriction: **the first
character of the filename must be a letter, not a number.**

Format: tab-separated, with a header row of `cluster`, `celltypes`, and
(optionally) `partition`.

```text
cluster	celltypes	partition
0	Follicular	2
1	Pericytes	1
2	T_cells	3
3	Follicular	2
4	Endothelial	4
5	T_cells	3
6	T_cells	3
7	Myeloid	5
8	Plasma	2
9	Fibroblasts	1
10	B_cells	6
11	Endothelial	4
```

This is the actual annotation file used for the PRJNA790856 run, at the
schema chosen in Vignette 1 (`sct` normalization, `harmony` integration,
resolution `0.2`) — saved as
`gene_files/cell_annotation_prjna790856_res.0.2.txt`. Twelve clusters, eight
distinct cell type labels. A few things worth noting about how this file is
structured:

- **Repeated labels across multiple clusters are expected and fine.**
  Clusters 0 and 3 are both `Follicular`; clusters 2, 5, and 6 are all
  `T_cells`; clusters 4 and 11 are both `Endothelial`. SWANS aggregates
  clusters sharing the same label when computing z-scores and downstream
  DGE/pathway comparisons by cell type, so there's no need to force distinct
  labels onto biologically identical clusters just because they were
  assigned separate cluster numbers during clustering.
- **Partition numbers group clusters for trajectory analysis**, not by cell
  type label. Note that partition 2 contains both `Follicular` clusters (0,
  3) *and* the `Plasma` cluster (8) — partitioning is about which clusters
  Monocle3 should treat as connected along a shared trajectory, which is a
  separate question from cell-type identity. Every one of the 12 clusters
  has a partition assigned, since `PARTITION_TRAJECTORY: y` in this run's
  config (Section 3).

**The `partition` column is optional** and only needed if you're running
Monocle3 trajectory analysis and want clusters grouped into shared
trajectory branches/roots (e.g., grouping several closely related immune
subtypes into one partition so Monocle3 builds a single connected trajectory
across them, rather than isolated per-cluster paths). If you include this
column, **every** cluster must have a partition assigned — none can be left
blank. If you omit the column, set `PARTITION_TRAJECTORY: n` in the config
(Section 3) and Monocle3 will run without partitioning.

**Important:** the path to this file in `configs/post_annotation_configs.yaml`
must be a **relative path**, not an absolute one.

---

## 3. Configuring `configs/post_annotation_configs.yaml`

Copy `configs/example_post_annotation_configs.yaml`, rename it to
`configs/post_annotation_configs.yaml`, and customize.

**The first time you run this**, `RUN_FINAL_ANALYSIS` must be `n`. Do not set
it to `y` until the preliminary analysis has fully completed and you have a
finished `CLUSTER_ANNOTATION_FILE` ready.

```yaml
# run post-annotation (final) analysis
RUN_FINAL_ANALYSIS: y

# minimum % of cells in one of the compared groups, for a gene to be tested (DGE)
MIN_PCT: 0.1

# filtering DEGs by avg_log2FC magnitude before pathway analysis
# 0 keeps all results (appropriate for pathway/GSEA, which wants the full ranked list);
# a non-zero value pre-filters to a shorter, higher-effect-size DEG list
AVG_LOG2FC_THRESHOLD: 1

# adj.p.value threshold for pathway/GSEA results
FINAL_FILTERING_THRESHOLD: 0.10

# ---- final schema selection (matches what was chosen in Vignette 1) ----
FINAL_SEURAT_NORMALIZATION_METHOD: sct
FINAL_SEURAT_INTEGRATION_METHOD: harmony
FINAL_RESOLUTION: 0.2

# relative path to the annotation file built in Section 2
CLUSTER_ANNOTATION_FILE: cell_annotation_prjna790856_res.0.2.txt

# ---- trajectory analysis ----
RUN_TRAJECTORY_ANALYSIS: y
PARTITION_TRAJECTORY: y

# ---- supplying a pre-analyzed Seurat object (not used in this tutorial) ----
PROVIDE_ANALYZED_SEURAT_OBJECT: n
USER_ANALYZED_SEURAT_OBJECT:
USER_ANALYZED_SEURAT_OBJECT_META_SAMPLE:
USER_ANALYZED_SEURAT_OBJECT_META_EXPERIMENT:
USER_ANALYZED_SEURAT_OBJECT_META_ANNOTATION:
USER_UMAP_REDUCTION:
USER_TNSE_REDUCTION:
ANNOTATE_PROVIDED_FINAL_SEURAT_OBJECT:

# ---- output ----
# cellchat: saves the final object in a format ready for downstream CellChat
# ligand-receptor / cell-cell communication analysis, beyond what this tutorial covers
FINAL_STORAGE: cellchat

# 10x end user license agreement — required if 'cloupe' is included in FINAL_STORAGE;
# set here regardless as a matter of course
EULA: y

# ---- optional gene visualization (not used for this run) ----
FINAL_USER_GENE_FILE:
FINAL_VISUALIZATION:

# multiple experimental conditions present (PTC vs. Normal), so conserved genes is meaningful here
FINAL_CONSERVED_GENES: y

FINAL_THREADS: 30
```

A note on `FINAL_STORAGE: cellchat` specifically: this tells SWANS to save
the final annotated object in a format ready to feed into
[CellChat](https://github.com/sqjin/CellChat) for cell-cell communication /
ligand-receptor analysis. That analysis itself is outside the scope of this
tutorial, but it's worth knowing this option exists if your next step after
annotation is exploring signaling between the cell types SWANS identifies —
Follicular-to-Myeloid signaling, for instance, would be a natural next
question for this dataset.

---

## 4. Running the pipeline (second pass)

From the same working directory as before:

```bash
bash run_snakemake.sh
```

Snakemake determines what actually needs to run based on which inputs have
changed since the last run — it does not rerun the entire preliminary
analysis just because you're now running the post-annotation pass. In
practice this means: if you need to tweak something small in the
post-annotation config (say, `FINAL_FILTERING_THRESHOLD`) and rerun, only
the rules downstream of that change are re-executed, not the full 18-schema
comparison from Vignette 1.

As before, an email notification is sent to your `CONTACT` address on
success or failure, including the exact failed command and log path if
something goes wrong.

---

## 5. Reading the post-annotation ("final") report

**Location:**
`data/endpoints/PROJECT/analysis/final_analysis/PROJECT_final_report.html`

The report is organized into the following sections, in order:

### Sample & parameter overview
A summary of samples analyzed and the full set of YAML parameters used for
this run — useful as a record of exactly how the analysis was configured
when interpreting or sharing results later.

### Annotated UMAPs
Multiple views of the same annotated embedding: overall, split by
experimental condition, split by sample, and colored by cell cycle phase
(G1/S/G2M). Comparing the "by experiment" and "by sample" views is a quick
way to check whether a cluster is driven by biology shared across the
cohort or is an artifact of one or two samples.

### Cell counts / proportions
Bar plots and tables of cluster composition — overall, by experimental
condition, and by sample. This is the section that most directly supports
comparison against a reference dataset: for PRJNA790856, this is where you'd
compare SWANS's recovered cell-type proportions against those reported by
Wang et al. 2022 (see Section 8 below).

### Differentially expressed genes (DEGs) by cluster
Unlike the preliminary report's DEGs (each cluster vs. all others), this
table shows within-cluster DEGs **across experimental conditions** — e.g.,
tumor vs. normal within the Follicular cell cluster specifically. This is
the core comparison most bench scientists will care about: not "what
defines this cluster" but "what changes in this cell type between
conditions."

### Pathway analysis (GSEA) results
Two tables of the same underlying GSEA results, ranked two different ways:
by adjusted p-value, and by `avg_log2FC` (NES). Both orderings are provided
because they surface different things — p-value ranking highlights
statistically confident pathways regardless of effect size, while NES/fold-
change ranking highlights the strongest effect sizes regardless of
significance. For a bench scientist without a stats background: p-value
ranking answers "what am I most sure is real," NES ranking answers "what
changed the most."

### Conserved genes (optional)
Genes with consistent expression in a given cluster across experimental 
conditions. 

### Z-score transformations
Three levels of granularity: by cluster, by cluster + experimental
condition, and by cluster + experimental condition + sample. These are
useful as a second lens on expression that's independent of the
`FindAllMarkers`-style pairwise comparison structure.

### Variable gene heatmap
A heatmap of the top 75 most variable genes across the dataset, plus a full
table of all variable genes.

### Trajectory analysis (if run)
A Monocle3 trajectory plot overlaid on the annotated UMAP, showing inferred
pseudotime paths between clusters — partitioned according to your
`CLUSTER_ANNOTATION_FILE`'s `partition` column if you supplied one.

---

## 6. Interpreting results as a bench scientist

A practical reading order if you're not a bioinformatician: start with the
annotated UMAP and cell-count tables to confirm the biology looks sensible
at a glance (right populations, reasonable proportions). Then go to the DEG
table filtered to your cell type of interest and sort by adjusted p-value.
Cross-check the top hits against the pathway table's NES ranking — genes
that show up consistently in both a significant DEG list and a
high-magnitude pathway hit are the strongest candidates worth following up
experimentally. The z-score tables are most useful when a DEG result looks
ambiguous or when you want to compare expression of a specific gene across
every level (cluster, condition, sample) at once.

---

## 7. Benchmarking report (brief)

**Location:**
`data/endpoints/PROJECT/analysis/report/benchmarks/PROJECT_benchmark_report.html`

This report is informational, not scientific — it logs runtime, memory
(`max_rss`, `max_vms`), I/O, and CPU time for every rule in the pipeline.
It's most useful for planning resource allocation on future runs (e.g.,
seeing that DoubletFinder and the multi-schema Seurat analysis are the most
resource-intensive steps) rather than for interpreting biological results.

---

## 8. Validating against Wang et al. 2022 — worked example

This section demonstrates the kind of outcome-based validation the
PRJNA790856 tutorial is designed to support: does SWANS recover the biology
already published for this dataset?

Wang et al. 2022 reported 8 cell types from 29,561 cells after their own QC
filtering (feature range 201–5,000, mito <25% — the same thresholds this
tutorial's prelim config was deliberately set to match, see Vignette 1
Section 5): Follicular cells, T cells, Endothelial, Pericytes, B cells,
Fibroblasts, Myeloid, Mast cells. A final note on analysis differences, 
we used doubletFinder as an additional QC step that was not employed in
the Wang et al. analysis. 

The annotation file in Section 2 above reflects what SWANS actually
recovered at the chosen schema (`sct`/`harmony`/res `0.2`):

**Concordant — recovered cleanly, matching Wang et al.:**
- Follicular cells, Fibroblasts, Endothelial cells, Myeloid, Pericytes

**More granular in SWANS (defensible, not contradictory):**
- Plasma cells — resolved as distinct from B cells (cluster 8 vs. 10),
  where Wang et al. grouped them together

**Same population, less subdivided than one might expect:**
- T cells — all three T-cell clusters (2, 5, 6) share the single label
  `T_cells` rather than being split into named subtypes at this resolution.
  This is a deliberate choice at res 0.2, not a limitation: the
  `USER_GENE_FILE` marker panel used in the prelim pass (NK markers KLRB1/
  NCAM1/FCGR3A, Treg markers IL2RA/FOXP3, tissue-resident memory markers
  ZNF683/ITGAE/CXCR6) is available for characterizing these three clusters
  further via dot plot without needing a higher-resolution reclustering.

**Not split into its own cluster, but not missing either:**
- Mast cells — no cluster in the res 0.2 annotation is labeled `Mast`; TPSAB1-expressing
  cells co-cluster with the broader myeloid population (cluster 7). This is
  visible directly in the standard dot plot SWANS generates for every
  schema from the `USER_GENE_FILE` marker panel (found under
  `analysis/normalization/.../PROJECT_DimPlot_Proportions_sct.harmony_snn_res.0.2.pdf`):
  TPSAB1 shows a small but distinct signal concentrated in cluster 7,
  alongside the other myeloid markers CD68/CD14. So the marker signal is
  there in this tutorial's own res 0.2 output — it just isn't broken out as
  its own cluster label, consistent with mast cells being a low-abundance
  population (~1.4% in Wang et al.) sharing myeloid lineage.

  Note: the manuscript's own validation of this (Table 9, Figure 8) further
  confirms mast cell signal at resolution 0.5 using `Seurat::AddModuleScore()`
  with a TPSAB1-based gene set. That module-scoring step is a supplementary
  analysis performed for the manuscript's validation narrative, not a
  standard SWANS report output — the dot plot referenced above, generated
  automatically from `USER_GENE_FILE`, is the tool this tutorial's pipeline
  run actually produces for this purpose.

If you're using this tutorial as a template for validating your own
dataset against a reference, this comparison structure — concordant /
more granular / resolution-dependent — is a useful framework for organizing
that discussion.

---
