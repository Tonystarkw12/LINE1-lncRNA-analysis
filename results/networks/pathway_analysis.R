
# 通路富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# 基因列表
genes <- c("TP53", "MYC", "EGFR", "KRAS", "BRCA1")

# 转换为Entrez ID
entrez_ids <- bitr(genes, fromType="SYMBOL", 
                   toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO富集分析
go_results <- enrichGO(gene = entrez_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2,
                       readable = TRUE)

# KEGG通路分析
kegg_results <- enrichKEGG(gene = entrez_ids$ENTREZID,
                           organism = 'hsa',
                           pvalueCutoff = 0.05)

# 保存结果
write.csv(as.data.frame(go_results), file="results/networks/pathways_go.csv")
write.csv(as.data.frame(kegg_results), file="results/networks/pathways_kegg.csv")

# 绘图
png("results/networks/pathways_go.png", width=1200, height=800, res=300)
dotplot(go_results, showCategory=20) + ggtitle("GO Biological Process Enrichment")
dev.off()

png("results/networks/pathways_kegg.png", width=1200, height=800, res=300)
dotplot(kegg_results, showCategory=20) + ggtitle("KEGG Pathway Enrichment")
dev.off()

cat("通路分析完成\n")
