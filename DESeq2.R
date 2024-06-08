# 设置个人库路径
.libPaths("~/R/library")

# 安装必要的R包到个人库路径中
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib="~/R/library")
BiocManager::install("DESeq2", lib="~/R/library")
BiocManager::install("EnhancedVolcano", lib="~/R/library")
BiocManager::install("pheatmap", lib="~/R/library")

# 加载R包
library("DESeq2", lib.loc="~/R/library")
library("EnhancedVolcano", lib.loc="~/R/library")
library("pheatmap", lib.loc="~/R/library")

# 读取数据
counts_matrix <- read.csv("path/to/counts_matrix.csv", row.names=1)
sample_info <- read.csv("path/to/sample_info.csv")

# 检查数据
head(counts_matrix)
head(sample_info)

# 创建DESeq2对象并进行差异表达分析
dds <- DESeqDataSetFromMatrix(countData = counts_matrix, 
                              colData = sample_info, 
                              design = ~ stage + group)
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
res <- results(dds)
head(res)

# 筛选结果并保存
res_filtered <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(as.data.frame(res_filtered), file="differential_expression_results.csv")

# 可视化结果
EnhancedVolcano(res_filtered,
                lab = rownames(res_filtered),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0)
EnhancedVolcano(res_filtered,
                lab = rownames(res_filtered),
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 1,
                legendLabels = c('NS','Log2 FC','p-value','p-value & Log2 FC'),
                legendPosition = 'top',
                title = "Volcano Plot")
sig_genes <- rownames(res_filtered)
norm_counts <- counts(dds, normalized=TRUE)
sig_norm_counts <- norm_counts[sig_genes,]
pheatmap(sig_norm_counts, scale = "row", clustering_distance_rows = "correlation")
