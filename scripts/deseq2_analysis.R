# Load libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)

#Check dimensions
length(colnames(counts))
length(rownames(meta))
#View sample names clearly
colnames(counts)
rownames(meta)
#Clean count sample names properly
colnames(counts) <- gsub(
    ".*alignment[./]",
    "",
    colnames(counts)
)
colnames(counts) <- sub("_sorted.bam", "", colnames(counts))
colnames(counts)
#Rebuild metadata correctly
meta <- read.table(
    "data/metadata.txt",
    header = FALSE,
    stringsAsFactors = FALSE
)
Meta
#assign metadata columns
colnames(meta) <- c("sample", "condition")
#set row names 
rownames(meta) <- meta$sample
#keep only condition column 
meta <- meta["condition"]
#convert to factor 
meta$condition <- as.factor(meta$condition)
#record metadata 
meta <- meta[colnames(counts), , drop = FALSE]
#final verification 
all(colnames(counts) == rownames(meta))
#run deg 
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ condition
)
#filter low counts 
dds <- dds[rowSums(counts(dds)) > 10, ]
#run deseq2
dds <- DESeq(dds)
#results 
res <- results(dds)
write.csv(
    as.data.frame(res),
    "results/DEG_results.csv"
)
#significant genes 
sig <- subset(
    as.data.frame(res),
    padj < 0.05 &
    abs(log2FoldChange) > 1
)

write.csv(
    sig,
    "results/significant_genes.csv"
)
#Full deg result
write.csv(
    as.data.frame(res),
    file = "results/DEG_results.csv",
    quote = FALSE
)
#save normalized count 
vsd <- vst(dds, blind = FALSE)

norm_counts <- assay(vsd)

write.csv(
    norm_counts,
    file = "results/normalized_counts.csv",
    quote = FALSE
)
#Save pca 
png(
    "figures/PCA_plot.png",
    width = 1200,
    height = 900,
    res = 150
)

plotPCA(vsd, intgroup = "condition")

dev.off()
#save volcano 
library(ggplot2)

res_df <- as.data.frame(res)

res_df$significant <- ifelse(
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) > 1,
    "Significant",
    "Not Significant"
)
#createplot 
volcano <- ggplot(
    res_df,
    aes(
        x = log2FoldChange,
        y = -log10(padj),
        color = significant
    )
) +
geom_point(size = 2) +
theme_minimal() +
labs(
    title = "Volcano Plot",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
)
#save 
ggsave(
    filename = "figures/volcano_plot.png",
    plot = volcano,
    width = 8,
    height = 6
)
#save pheatmap
library(pheatmap)

topgenes <- head(order(res$padj), 20)

mat <- assay(vsd)[topgenes, ]
png(
    "figures/heatmap_top20.png",
    width = 1200,
    height = 1200,
    res = 150
)

pheatmap(
    mat,
    scale = "row"
)

dev.off()
#ma plot
png(
    "figures/MA_plot.png",
    width = 1200,
    height = 900,
    res = 150
)

plotMA(res, ylim = c(-5, 5))

dev.off()
#sample distance pheatmap
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- colnames(vsd)

colnames(sampleDistMatrix) <- colnames(vsd)
png(
    "figures/sample_distance_heatmap.png",
    width = 1200,
    height = 1200,
    res = 150
)

pheatmap(
    sampleDistMatrix
)

dev.off()
#dispersion plot 
png(
    "figures/dispersion_plot.png",
    width = 1200,
    height = 900,
    res = 150
)

plotDispEsts(dds)

dev.off()
#boxplot of normalized counts 
norm_counts <- counts(
    dds,
    normalized = TRUE
)
png(
    "figures/boxplot_normalized_counts.png",
    width = 1400,
    height = 900,
    res = 150
)

boxplot(
    log2(norm_counts + 1),
    las = 2,
    main = "Normalized Counts",
    ylab = "Log2 Counts"
)

dev.off()
#gene expression plot for deg 
topgene <- rownames(res)[which.min(res$padj)]

topgene
plotCounts(
    dds,
    gene = topgene,
    intgroup = "condition"
)
png(
    "figures/top_gene_expression.png",
    width = 1200,
    height = 900,
    res = 150
)

plotCounts(
    dds,
    gene = topgene,
    intgroup = "condition"
)

dev.off()
#histogram of adjusted p value 
png(
    "figures/padj_histogram.png",
    width = 1200,
    height = 900,
    res = 150
)

hist(
    res$padj,
    breaks = 50,
    col = "grey",
    main = "Adjusted P-value Distribution",
    xlab = "Adjusted P-value"
)

dev.off()

