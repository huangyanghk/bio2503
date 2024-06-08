data <- read.csv("expression_data.csv", row.names = 1, header = TRUE)
datExpr <- t(data)

# 样本聚类和数据清理（与之前步骤相同）
sampleTree <- hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
cutHeight <- 15
clust <- cutreeStatic(sampleTree, cutHeight = cutHeight)
keepSamples <- (clust == 1)
datExpr <- datExpr[keepSamples, ]

# 过滤低表达基因
geneVar <- apply(datExpr, 2, var)
nGenes <- 5000
topGenes <- order(geneVar, decreasing = TRUE)[1:nGenes]
datExprFiltered <- datExpr[, topGenes]

# 标准化数据
datExprFiltered <- scale(datExprFiltered)


# 选择软阈值功率
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)


# 绘制软阈值选择图
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

# 根据软阈值选择图，选择softpower=20
softPower <- 20

# 计算邻接矩阵
adjacency <- adjacency(datExprFiltered, power = softPower)

# 计算拓扑重叠矩阵（TOM）及其距离矩阵
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# 聚类基因
geneTree <- hclust(as.dist(dissTOM), method = "average")

# 模块检测
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)

# 将动态剪切得到的基因模块颜色化
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# 绘制基因树及模块颜色
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

# 读取表型数据
traitData <- read.csv("trait_data.csv", row.names = 1, header = TRUE)

# 计算模块特征基因（MEs）
MEs <- moduleEigengenes(datExpr, colors = dynamicColors)$eigengenes
MEs <- orderMEs(MEs)

# 计算模块与表型的相关性
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 可视化相关性矩阵
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(traitData), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = ""),
               setStdMargins = FALSE, cex.text = 0.5, zlim = c(-1,1),
               main = "Module-trait relationships")
