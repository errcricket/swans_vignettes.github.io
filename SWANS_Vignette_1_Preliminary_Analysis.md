# SWANS Tutorial: Preliminary Analysis
### Worked example — GSE191288 / PRJNA790856 (Wang et al. 2022, papillary thyroid carcinoma)

> **Scope of this document.** This vignette walks a first-time user through the
> **preliminary analysis** phase of SWANS — from raw input data through schema
> selection in the interactive report. It stops at the point where you have chosen
> one clustering schema and are ready to annotate it. Annotation, differential
> expression across conditions, GSEA, and trajectory analysis are covered in the
> companion document, **SWANS Vignette 2: Post-Annotation Analysis**.
>
> This tutorial uses GSE191288 (SRA project PRJNA790856), a bilateral papillary
> thyroid carcinoma single-cell dataset from Wang et al. 2022, as a public, fully
> worked example. It is the same dataset used for validation in the SWANS
> methods manuscript.

---

## 1. Requirements

SWANS uses Snakemake as a workflow manager, with Singularity executing each rule
inside a container. You do not need to install Seurat, R packages, or Cell Ranger
yourself — the required Docker images (`POND` for all analysis rules, `cellranger`
for the `cellranger_counts` rule only) are pulled automatically by Singularity.

For this tutorial specifically, you will **not** need the `cellranger` image at
all, because PRJNA790856 is run starting from feature-barcode matrix files
(see Section 3).

Minimum environment:
- Snakemake (tested on 7.32.4)
- Singularity / singularity-ce (tested on 4.3.2-1.el8)

---

## 2. What you need before starting

You must have the following in your working directory before running the
pipeline:

- `samples.sample_list`
- `configs/prelim_configs.yaml`
- `configs/post_annotation_configs.yaml`

And your starting single-cell data — see Section 3 for the three supported
formats.

---

## 3. Starting data: three supported input types

SWANS supports **three different starting points**, controlled by the
`STARTING_DATA` field in `configs/prelim_configs.yaml`. Which one you use
depends on what you have available.

| `STARTING_DATA` value | What you need | `path_to_starting_data` should point to | Notes |
|---|---|---|---|
| `fastq` | Raw FASTQ files | Folder containing the FASTQ files | Triggers `RUN_CELLRANGER: y`. Requires exact naming convention (see README, "Naming conventions for running Cell Ranger"). SWANS runs Cell Ranger's `count` pipeline for you. |
| `cellranger` | Completed Cell Ranger output | Folder **containing** the `outs` directory (not `outs` itself) | Use this if Cell Ranger has already been run outside SWANS. |
| `matrix` | Feature-barcode matrix files | Folder containing `features.tsv.gz`, `barcodes.tsv.gz`, `matrix.mtx.gz` (named exactly this) | No Cell Ranger step is run. MultiQC cannot be generated from this starting point. |

**This tutorial uses `matrix`.** GSE191288's public deposit provides
filtered feature-barcode matrices per sample, so PRJNA790856 is run with:

```yaml
STARTING_DATA: matrix
RUN_CELLRANGER: n
RUN_MULTIQC: n
```

> **A note on why SoupX is not used for this dataset.** SoupX (ambient RNA
> correction) is optional but recommended by SWANS, and normally you should
> run it if you have the data to support it. SoupX needs either Cell Ranger
> `outs` (recommended), raw feature-barcode matrix files, or an `.h5` file —
> critically, it needs access to the **unfiltered (raw)** matrix alongside the
> filtered one to estimate the ambient RNA profile. GEO's deposit of
> GSE191288 provides only the **filtered** matrices, not raw. Because the raw
> matrix isn't available for this dataset, SoupX is **not run** here:
>
> ```yaml
> RUN_SOUPX: n
> ```
>
> If you are starting your own project from `fastq` or `cellranger` output,
> you will typically have both raw and filtered data available, and running
> SoupX is recommended in that case.

---

## 4. Setting up `samples.sample_list`

This file lives at the top of your working directory and defines which
samples will be analyzed, what experimental condition each belongs to, and
where its starting data lives. Fields are **tab-separated**, with one sample
per line and no blank lines.

```text
samples condition   path_to_starting_data
T1L PTC PRJNA790856/GSM5743021_T1L/
T1R PTC PRJNA790856/GSM5743022_T1R/
T2L PTC PRJNA790856/GSM5743023_T2L/
T2R PTC PRJNA790856/GSM5743024_T2R/
T3L PTC PRJNA790856/GSM5743025_T3L/
T3R PTC PRJNA790856/GSM5743026_T3R/
NT  Normal  PRJNA790856/GSM5743027_NT/
```

Seven samples total: three patients, each with paired **L**eft and **R**ight
papillary thyroid carcinoma (PTC) tumor samples, plus one normal thyroid
(NT) sample used as the comparison condition. This is the bilateral-tumor
design referenced throughout the SWANS manuscript's validation section.

Each `path_to_starting_data` entry points to the folder containing that
sample's `features.tsv.gz`, `barcodes.tsv.gz`, and `matrix.mtx.gz`, per
Section 3. 

---

## 5. Configuring `configs/prelim_configs.yaml`

Copy `configs/example_prelim_configs.yaml`, rename it to
`configs/prelim_configs.yaml`, and customize. Below is the annotated,
complete configuration used for this tutorial's PRJNA790856 run.

```yaml
# contact (email will be sent when jobs complete)
contact: user_name@email.com

# project name (IN LOWER CASE) — determines output directory under data/endpoints/
PROJECT: prjna790856

# organism
ORGANISM: human

# sequencing type — 10X Genomics Standard/GEM or Flex
SEQUENCING_TYPE: standard

# ---- starting data ----
STARTING_DATA: matrix
RUN_CELLRANGER: n
RUN_MULTIQC: n
OUTPUT_BAM: n

# ---- SoupX / DoubletFinder ----
# SoupX not run — raw (unfiltered) matrices unavailable for this dataset (see Section 3)
RUN_SOUPX: n
SOUPX_START: no_clusters
RUN_DOUBLETFINDER: y

# ---- filtering thresholds ----
# these intentionally match the thresholds used by Wang et al. 2022 for this same
# dataset (<25% mito, 201-5,000 genes) — matching the original publication's QC
# criteria is a deliberate validation choice, so that any differences in downstream
# results reflect the pipeline/clustering, not a different starting cell population.
MITO: 25
RIBO: 100

MIN_FEATURE_THRESHOLD: 200
MAX_FEATURE_THRESHOLD: 5000

# ---- normalization / scaling ----
SPLIT_LAYERS_BY: Sample
COMPONENTS: 50
NUM_VARIABLE_FEATURES: 3000
SCALE_DATA_FEATURES: variable

# ---- regression (none applied) ----
MITO_REGRESSION: n
RIBO_REGRESSION: n
REGRESSION_FILE:
CELL_CYCLE_REGRESSION: n
CELL_CYCLE_METHOD: standard

# ---- schema comparison: six normalization x integration combinations, seven resolutions ----
# {standard, sct} normalization x {cca, harmony, rpca} integration
SEURAT_NORMALIZATION_METHOD: sct,standard
SEURAT_INTEGRATION_METHOD: cca,harmony,rpca
RESOLUTION: 0.1,0.2,0.3,0.4,0.5,0.6,0.7

# which FindClusters algorithm (blank = default Louvain; 2 = Louvain w/ multilevel
# refinement; 3 = SLM; 4 = Leiden)
FIND_CLUSTERS_ALGORITHM:

# ---- integration / annotation options not used in this tutorial ----
REFERENCE_BASED_INTEGRATION: n
REFERENCE_SAMPLES:
RUN_AZIMUTH: n
AZIMUTH_REFERENCE:
RUN_TRANSFERDATA: n
TRANSFERDATA_REDUCTION: standard.pca
TRANSFERDATA_ANNOCOL: celltypes

# ---- misc ----
TSNE: n
CONSERVED_GENES: n
STORAGE:
THREADS: 30
MEMORY: 1389064417920

# candidate marker genes for visualization across all schemas — see Section 8a below
USER_GENE_FILE: gene_files/prjna1185392.txt
VISUALIZATION: dot

# ---- environment-specific paths — leave as configured for your institution's environment ----
CELLRANGER: /usr/local/bin/cellranger
MULTIQC: /usr/local/bin/multiqc
RPATH: /usr/local/lib/R/site-library/
CELLRANGER_REFERENCE:
```

A few fields worth understanding, not just copying:

- **`RESOLUTION: 0.1,0.2,0.3,0.4,0.5,0.6,0.7`** — SWANS clusters at every
  resolution listed, in a single pass, for every normalization × integration
  combination. With 2 normalization methods, 3 integration methods, and 7
  resolutions, this produces **42 distinct clustering schemas**, each with
  its own DEG table, z-score table, and cluster proportion table. This is
  the core design feature that lets you compare schemas side-by-side later
  in the interactive report rather than committing to one combination up
  front — and running a wide resolution sweep like this one gives a much
  clearer Clustree picture of exactly where over-clustering starts.
- **`MITO` / `RIBO`** — set the cutoff to `0` to exclude all cells expressing
  any amount of that gene class, or `100` for no filtering at all.
- **`COMPONENTS: 50`** — this is a starting number of PCs to *investigate*,
  not the number ultimately used. SWANS quantifies the number of PCs needed
  to explain 90% of variance and uses that value downstream.
- **`TRANSFERDATA_REDUCTION` / `TRANSFERDATA_ANNOCOL`** are populated here
  even though `RUN_TRANSFERDATA: n` — they're inert while that flag is off,
  left filled in simply because this config was adapted from a template
  where TransferData was used. Harmless to leave populated, but they have
  no effect unless `RUN_TRANSFERDATA` is set to `y`.
- **`USER_GENE_FILE`** points to `gene_files/prjna1185392.txt` — a curated
  list of thyroid, immune, stromal, and vascular marker genes (TG, TSHR,
  TPO for follicular cells; CD3D, CD8B for T cells; IL2RA/FOXP3 for Tregs;
  KLRB1/NCAM1/FCGR3A for NK cells; TPSAB1 for mast cells; PDGFRA for
  fibroblasts; VWF for endothelial; HIGD1B for pericytes; CD68/CD14/S100A8/
  S100A9 for myeloid, among others). This is a text file with no header,
  and one gene per line. With `VISUALIZATION: dot`, SWANS
  generates a dot plot of this gene panel across every schema, which is
  genuinely useful for cross-checking whether a cluster's identity holds up
  under known markers — independent of, and a good complement to, the
  automatically generated DEG/z-score tables.

---

## 6. Running the pipeline

From the directory containing `samples.sample_list` and
`configs/prelim_configs.yaml`:

```bash
bash run_snakemake.sh
```

Snakemake manages parallel execution of rules according to your `THREADS`
setting. You do not need to babysit the run — an email is sent to the
`CONTACT` address on completion, whether successful or not.

**For reference**, on a comparably-sized 9-sample test run (a different
project, shown here only to illustrate what the benchmarking report looks
like), the heaviest single rule — the multi-schema Seurat analysis step —
took roughly 1h 53m and peaked around 16.6 GB RSS on 30 threads. Runtime 
and memory will scale with your own dataset size, number of samples, and 
— notably for this tutorial's 42-schema sweep — how many 
normalization/integration/resolution combinations you're comparing 
in a single pass. 

If a rule fails, both the terminal output and the failure email will include
the exact shell command that was run and the path to that rule's log file
(under `logs/PROJECT/rule_name/`). Read that log first — it almost always
identifies the specific error.

---

## 7. Reading the QC report

**Location:** `data/endpoints/PROJECT/analysis/report/PROJECT_qc_report.html`

This report is generated per sample and shows:

- DoubletFinder classification (UMAP of doublets vs. singlets, plus a
  violin plot of `nFeature_RNA` split by classification) — for each sample,
  you'll see total cells, doublet rate, doublet count, and cells remaining
  after doublet removal
- `nFeature_RNA`, `nCount_RNA`, `percent.mito`, `percent.ribo` — both
  unfiltered and filtered, across all samples together, so you can spot
  outlier samples at a glance
- Scatter plots correlating mito%/ribo% against RNA counts/features, before
  and after filtering, which is a useful sanity check that your thresholds
  are removing what you intend (low-quality/dying cells) without removing
  legitimate biology

**What to look for:** a healthy sample generally shows a tight,
near-linear relationship between `nCount_RNA` and `nFeature_RNA` after
filtering, and a clear separation between the doublet and singlet violin
distributions. If one sample's doublet rate or mito% distribution looks very
different from the others, that's worth a closer look before proceeding —
it may indicate a technical issue specific to that sample rather than
biological signal.

---

## 8. Reading the Interactive (Shiny) report

**Location:** `data/endpoints/PROJECT/analysis/report/Interactive_report.Rmd`

This is the report you'll spend the most time in, since it's what lets you
compare all 42 schemas and choose one.

**Launching it:** set your R working directory to the `report/` folder, open
`Interactive_report.Rmd` in RStudio, install the required packages if you
haven't already —

```
yaml_2.3.10, magick_2.8.4, pdftools_3.4.0, rsconnect_1.3.1,
dplyr_1.1.4, DT_0.33, shinyfullscreen_1.1.0, shinyjs_2.1.0, shiny_1.9.0
```

— then click **Run Document**.

The report has three independent panels, each with its own dropdown
selections. You can compare up to three different schemas side-by-side, or
set panels to the same normalization/integration and vary only
resolution.

**Plots.** Each panel: choose normalization (standard/sct), integration
(cca/rpca/harmony), and resolution (0.1/0.2/0.3/0.4/0.5/06/0.7), and the matching UMAP
renders. Changing any dropdown updates that panel's UMAP; the other two
panels are unaffected.

**Tables.** Each panel independently lets you choose between:
- **Top 100 DEGs by cluster** (`avg_log2FC`) — upregulated markers per
  cluster from `FindAllMarkers`
- **Z-scores by cluster** — genes ranked by scaled expression rather than
  fold-change; useful when `avg_log2FC` results are dominated by comparisons
  between one dominant cluster and its nearest neighbor (see the
  "Suggestions" note in Section 9 below)
- **Cell counts/proportions by cluster and experimental condition**
- **Cell cycle phase distribution.** For the selected schema, G1/S/G2M phase
proportions are shown per experimental condition, broken down by cluster —
useful for confirming that a cluster split isn't simply tracking cell cycle
phase rather than a distinct cell type.

Enter `-1` in the cluster number field to see all clusters; enter a specific
cluster number to restrict the table to that cluster; an invalid positive
number shows nothing.

**Clustree.** For the panel's selected normalization + integration, a
Clustree figure shows how cells move between clusters as resolution
increases — this is where over-clustering becomes visually obvious (a
cluster splitting into two children that are barely distinguishable in
composition, for instance).

---

## 9. A caution on `FindAllMarkers` and large, closely related clusters

Worth knowing before you interpret DEG tables: `FindAllMarkers` compares
each cluster against *all other cells combined*. If your UMAP has one large
body split into several sub-clusters (a very common pattern), DGE results
will tend to highlight differences between the largest two sub-clusters
rather than markers distinguishing the smaller ones — because the smaller
clusters are outnumbered in the "all remaining cells" comparison. This is
one reason SWANS also computes z-scores (scaled expression by cluster,
independent of this comparison structure) as a second lens on the same
question. If you suspect this is happening, the z-score table or a
user-supplied candidate gene list (via `USER_GENE_FILE`) are more reliable
than the DEG table alone.

---

## 10. Choosing a schema

This is the practical decision point the rest of the pipeline depends on.
Using PRJNA790856 as a worked example, comparing SWANS's clustering against
the eight cell types reported by Wang et al. 2022 (Follicular cells, T
cells, Endothelial, Pericytes, B cells, Fibroblasts, Myeloid, Mast cells),
here is what the schema ultimately chosen for this dataset — `sct`
normalization, `harmony` integration, resolution `0.2` — actually recovered
(12 clusters, annotated in Vignette 2):

- **Follicular cells, Fibroblasts, Endothelial cells, Myeloid, Pericytes** —
  all recovered cleanly and consistently
- **T cells** — resolved into three separate clusters (2, 5, 6), all
  labeled `T_cells` at this resolution rather than split further into named
  subtypes. The `USER_GENE_FILE` marker panel (Section 5) includes NK
  (KLRB1, NCAM1, FCGR3A), Treg (IL2RA, FOXP3), and tissue-resident memory
  (ZNF683, ITGAE, CXCR6) markers specifically so these three clusters can be
  characterized further via the dot plot without needing to formally split
  them in the annotation file — a lighter-weight alternative to increasing
  resolution further just to separate T cell subtypes
- **Plasma cells** — resolved as a cluster distinct from B cells (cluster 8
  vs. cluster 10) — more granular than Wang et al., who grouped these
  together, and a defensible split rather than a contradiction
- **Mast cells** — **not** split out into their own cluster at this
  resolution; TPSAB1-expressing cells co-cluster with the broader myeloid
  population (cluster 7). The signal is still visible, though: the standard
  dot plot SWANS generates from the `USER_GENE_FILE` marker panel for this
  schema shows a small but distinct TPSAB1 signal in cluster 7 alongside
  CD68/CD14, so mast cell identity here comes from the marker panel rather
  than from a dedicated cluster label — expected for a low-abundance
  population (~1.4% in Wang et al.) sharing myeloid lineage.

In general: prefer the schema that recovers biologically expected
populations cleanly, is stable across nearby resolutions in the Clustree
view, and doesn't require excessive resolution to separate rare populations
you don't specifically need resolved. Chapter/Section 2.4 of the SWANS
manuscript ("Guidance for Schema Selection") covers the full reasoning
behind this in more depth.

---

## 11. Where to go from here

- **Robustness / stability analysis.** SWANS ships (or will ship — check the
  repository for current status) a separate suite of generic robustness and
  stability scripts, built to assess how consistent gene signatures and
  cluster assignments are across schemas and resolutions. These are a more
  advanced, standalone resource and are intentionally **not** covered in
  this tutorial — see the scripts' own documentation if you want to go
  beyond schema selection into quantitative schema robustness.
- **Full statistical and design rationale.** The SWANS methods manuscript
  covers the reasoning behind default parameter choices, the six-schema
  comparison design, and the outcome-based validation approach in more
  depth than this tutorial.

## 12. Next steps

Once you've settled on one normalization × integration × resolution
combination, you're ready to annotate clusters and move into the
post-annotation phase — covered in **SWANS Vignette 2: Post-Annotation
Analysis**, which requires a completed `CLUSTER_ANNOTATION_FILE` as its
starting point.
