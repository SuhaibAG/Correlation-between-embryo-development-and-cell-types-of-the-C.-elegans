# Correlation-between-embryo-development-and-cell-types-of-the-C.-elegans

R workflows for analyzing *C. elegans* embryonic gene expression with two complementary datasets:

1. **Bulk RNA-seq:** a developmental time course analyzed with Salmon + DESeq2 to find genes that change across stages.
2. **Single-cell RNA-seq:** STARsolo counts for ~70,000 cells from three time points (300, 400, 500), clustered and annotated by tissue, then linked back to the bulk stage-specific genes.

## Repository structure

```
├── Bulk_Scripts/
│   ├── import_data_project.R
│   ├── Differential_Expression_project.R
│   ├── Annotation_project.R
│   ├── PCA_project.R
│   ├── clustering_project.R
│   └── GO_term.R
├── SC_Scripts/Single-cell-project/
│   ├── Single-cell-project.Rproj
│   ├── data_loading.R
│   ├── pre_processing.R
│   ├── normalisation.R
│   ├── batching.R
│   ├── dimensionality_reduction.R
│   ├── clustering.R
│   ├── Marker_genes.R
│   ├── heat-mapping.R
│   ├── my_graphs.R
│   └── *.jpeg                  # output figures
└── extra-files/
    ├── sample_info.csv
    ├── cluster_names.csv
    ├── sig_genes_table.csv
    └── stages.csv
```

---

## 1. Bulk RNA-seq (`Bulk_Scripts/`)

Run the scripts in order, in the same R session (later scripts reuse `txi`, `dds`, `res`, `res_shrink_annot`, `plot_data`).

| # | Script | What it does |
|---|--------|--------------|
| 1 | `import_data_project.R` | Imports Salmon `quant.sf` files with `tximport`, summarizes to genes, removes sample `mpfc_700`, filters low-count genes (row sums ≤ 5), and compares raw, log2 and rlog counts. |
| 2 | `Differential_Expression_project.R` | `DESeq2` with design `~ developmental_stage` and a likelihood ratio test (LRT) against `~ 1`. Checks dispersion estimates. |
| 3 | `Annotation_project.R` | Annotates results with gene names and Entrez IDs from the *C. elegans* EnsDb, applies `ashr` LFC shrinkage, and draws MA plots and a Venn diagram of up- vs. down-regulated genes (padj < 0.05). |
| 4 | `PCA_project.R` | VST transformation and PCA colored by developmental stage. |
| 5 | `clustering_project.R` | Heatmap of significant genes (padj < 0.01) with hierarchical clustering of samples and genes. |
| 6 | `GO_term.R` | Per stage, takes the top 51 genes by mean VST difference vs. all other samples, then runs `enrichGO` (BP, CC, MF) and `compareCluster` across stages. |

## 2. Single-cell RNA-seq (`SC_Scripts/Single-cell-project/`)

Run the scripts in order, in the same R session (`sce` is passed between them). Counts come from STARsolo (`Solo.out/Gene/filtered`). Open `Single-cell-project.Rproj` in RStudio to work in this folder.

| # | Script | What it does |
|---|--------|--------------|
| 1 | `data_loading.R` | Reads all samples with `DropletUtils::read10xCounts` into one `SingleCellExperiment` (46,926 genes × 72,591 cells). |
| 2 | `pre_processing.R` | Adds metadata (sample, SRA run, time point), removes undetected genes (21,547 remain), annotates chromosomes via EnsDb, flags mitochondrial genes, and removes low-quality cells with `quickPerCellQC` per sample (70,055 cells remain). |
| 3 | `normalisation.R` | Pooled size factors (`quickCluster` + `computePooledFactors`) and log-normalized counts with `scran`/`scater`. |
| 4 | `batching.R` | Batch correction across samples with `batchelor::quickCorrect`, stored as a `corrected` reduced dimension. |
| 5 | `dimensionality_reduction.R` | Models gene variance, selects the top 10% highly variable genes, runs PCA (elbow point and denoised PCs) and UMAP. |
| 6 | `clustering.R` | Louvain clustering over a sweep of k = 5–30 (nearest-neighbor graph), picks k = 25 using cluster counts and silhouette width (21 clusters). |
| 7 | `Marker_genes.R` | Finds cluster markers with `scoreMarkers` (blocked by sample), checks known tissue markers (e.g. `unc-119`, `myo-3`, `myo-2`, `elt-2`, `dpy-7`, `pie-1`, `skn-1`), and plots bulk stage genes across clusters. |
| 8 | `heat-mapping.R` | Assigns each bulk stage gene to its peak stage (Early → Late), builds cluster-level heatmaps with time point, tissue and stage annotations, and finds each cluster's dominant time point. |
| 9 | `my_graphs.R` | Final figure styling and color palettes (UMAPs and heatmap). |

### Figures

| UMAP by time point | UMAP by tissue group |
|:--:|:--:|
| ![UMAP colored by time point](SC_Scripts/Single-cell-project/timepoint.jpeg) | ![UMAP colored by tissue group](SC_Scripts/Single-cell-project/Major_group.jpeg) |

| Time point and tissue side by side | Cluster heatmap of bulk stage genes |
|:--:|:--:|
| ![UMAPs by time point and tissue](SC_Scripts/Single-cell-project/stage_vs_time.jpeg) | ![Heatmap of bulk stage genes across clusters](SC_Scripts/Single-cell-project/Heatmap.jpeg) |

## Requirements

R with these packages:

- **Bulk:** `tximport`, `DESeq2`, `tidyverse`, `ashr`, `AnnotationHub`, `ensembldb`, `pheatmap`, `dendextend`, `ggvenn`, `RColorBrewer`, `ggrepel`, `clusterProfiler`, `org.Ce.eg.db`, `enrichplot`, `DOSE`
- **Single-cell:** `DropletUtils`, `scater`, `scran`, `batchelor`, `bluster`, `PCAtools`, `BiocParallel`, `igraph`, `patchwork`, `tidyverse`, `AnnotationHub`, `ensembldb`, `RColorBrewer`

## Input data

Supporting tables are in `extra-files/`:

| File | Used by | Contents |
|------|---------|----------|
| `sample_info.csv` | Single-cell | 32 samples with `sample`, `sra_run`, `timepoint` (12 at 300, 12 at 400, 8 at 500). |
| `sig_genes_table.csv` | Single-cell | Mean VST expression of ~7,100 significant genes across five stage groups (Early, Intermediate1–3, Late). |
| `stages.csv` | Single-cell | The same format for a ~100-gene subset. |
| `cluster_names.csv` | Single-cell | Tissue labels for the 21 single-cell clusters (e.g. body wall, cuticle, intestine, pharynx, somatic gonad) and a major group for each. |

Large inputs are not included in the repo:

- **Bulk:** a sample metadata file with a `developmental_stage` column, `references/tx2gene.tsv`, and Salmon output in `preprocessed/rnaseq/star_salmon/`.
- **Single-cell:** STARsolo output in `preprocessed/scrnaseq/star/`.
