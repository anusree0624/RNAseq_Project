# RNA-Seq Differential Gene Expression Analysis of T Cells During *Toxoplasma gondii* Infection

## Project Overview

This project analyzes RNA-Seq data from mouse T-cell populations isolated from the spleen and brain during *Toxoplasma gondii* ME49 infection. The objective is to identify transcriptional changes associated with infection progression and tissue-specific immune responses. Publicly available RNA-Seq datasets were processed to generate gene expression count matrices, perform differential expression analysis, and visualize gene expression patterns.

The analysis was conducted using the DESeq2 framework in R, including normalization, differential expression testing, and generation of multiple quality assessment and exploratory visualizations.
## Objectives

* Analyze RNA-Seq datasets from infected mouse T-cell populations.
* Perform quality assessment of sequencing data.
* Generate and process gene count matrices.
* Normalize gene expression counts.
* Identify differentially expressed genes (DEGs) between experimental conditions.
* Visualize sample relationships and expression patterns.
* Generate publication-ready figures for downstream biological interpretation.

## Dataset Information

### Study Details

| Parameter           | Description                                               |
| ------------------- | --------------------------------------------------------- |
| Project ID          | PRJNA1465895                                              |
| Organism            | Mus musculus (Mouse)                                      |
| Study Title         | Spleen and Brain T Cells in *Toxoplasma gondii* Infection |
| Sequencing Platform | Illumina NovaSeq 6000                                     |
| Library Type        | Paired-End RNA-Seq                                        |
| Tissue Types        | Brain and Spleen                                          |
| Cell Types          | CD4+ and CD8+ T Cells                                     |
| Infection Model     | *Toxoplasma gondii* ME49                                  |
| Mouse Strain        | C57BL/6NCrl                                               |
| Sex                 | Female                                                    |
| Time Points         | 9 dpi and 23 dpi                                          |

### Sample Accessions

* SRR38588400
* SRR38588403
* SRR38588406
* SRR38588409
* SRR38588412
* SRR38588418
* SRR38588421
* SRR38588424
* SRR38588427
* SRR38588415

## Tools and Packages Used

### Bioinformatics Tools

* SRA Toolkit
* FastQC
* MultiQC
* HISAT2
* SAMtools
* featureCounts

### R Packages

* DESeq2
* ggplot2
* pheatmap


## Project Directory Structure

```text
RNAseq_Project/
│
├── accession_list/
├── data/
├── counts/
├── qc/
├── results/
├── fig/
├── scripts/
└── README.md
```

---

## Workflow

```text
Public RNA-Seq Data (SRA)
            │
            ▼
      Data Download
            │
            ▼
      Quality Control
     (FastQC/MultiQC)
            │
            ▼
       Read Alignment
          (HISAT2)
            │
            ▼
      BAM Processing
         (SAMtools)
            │
            ▼
     Gene Quantification
      (featureCounts)
            │
            ▼
      Count Matrix
            │
            ▼
    Differential Expression
          (DESeq2)
            │
            ▼
       Significant DEGs
            │
            ▼
      Data Visualization
            │
            ▼
    Biological Interpretation
```


## Differential Expression Analysis

Gene expression analysis was performed using DESeq2.

### Metadata Processing

* Imported sample metadata.
* Assigned experimental conditions.
* Matched metadata with count matrix sample names.
* Converted conditions into factors for statistical analysis.

### DESeq2 Workflow

1. Constructed DESeq2 dataset object.
2. Filtered genes with low read counts.
3. Performed normalization and dispersion estimation.
4. Conducted differential expression testing.
5. Exported complete and significant DEG results.

Filtering criteria:

* Adjusted p-value (padj) < 0.05
* |log2 Fold Change| > 1



## Commands Used

### Load Required Libraries

```r
library(DESeq2)
library(ggplot2)
library(pheatmap)
```

### Create DESeq2 Dataset

```r
dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ condition
)
```

### Filter Low-Count Genes

```r
dds <- dds[rowSums(counts(dds)) > 10, ]
```

### Run Differential Expression Analysis

```r
dds <- DESeq(dds)
res <- results(dds)
```

### Extract Significant Genes

```r
sig <- subset(
    as.data.frame(res),
    padj < 0.05 &
    abs(log2FoldChange) > 1
)
```

### Save Results

```r
write.csv(as.data.frame(res),
          "results/DEG_results.csv")

write.csv(sig,
          "results/significant_genes.csv")
```

---

## Visualizations Generated

### Principal Component Analysis (PCA)

Used to evaluate clustering of samples based on global gene expression patterns.

Output:

```text
figures/PCA_plot.png
```

### Volcano Plot

Displays significantly upregulated and downregulated genes.

Output:

```text
figures/volcano_plot.png
```

### Heatmap

Visualization of expression profiles for the top differentially expressed genes.

Output:

```text
figures/heatmap_top20.png
```

### MA Plot

Shows expression changes relative to mean gene expression levels.

Output:

```text
figures/MA_plot.png
```

### Sample Distance Heatmap

Assesses similarity between samples.

Output:

```text
figures/sample_distance_heatmap.png
```

### Dispersion Plot

Evaluates variance estimates used by DESeq2.

Output:

```text
figures/dispersion_plot.png
```

### Normalized Counts Boxplot

Shows distribution of normalized counts across samples.

Output:

```text
figures/boxplot_normalized_counts.png
```

### Top DEG Expression Plot

Visualizes expression of the most significant gene across conditions.

Output:

```text
figures/top_gene_expression.png
```

### Adjusted P-value Histogram

Examines statistical significance distribution.

Output:

```text
figures/padj_histogram.png
```

---

## Outputs

| Folder  | Description                                                  |
| ------- | ------------------------------------------------------------ |
| counts  | Raw count matrices                                           |
| qc      | Quality control reports                                      |
| results | DEG results and normalized counts                            |
| fig     | PCA, volcano plots, heatmaps, and statistical visualizations |
| scripts | Analysis scripts and workflow code                           |

---

## Key Deliverables

* Differentially expressed gene list
* Normalized expression matrix
* PCA analysis
* Volcano plot visualization
* Heatmap of top DEGs
* Sample clustering analysis
* Statistical summaries of gene expression changes

This workflow provides a reproducible framework for identifying transcriptional signatures associated with T-cell responses during *Toxoplasma gondii* infection in mouse spleen and brain tissues.

This version is ready to paste directly into your `README.md` and matches the DESeq2 code you actually used.
