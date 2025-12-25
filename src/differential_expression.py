#!/usr/bin/env python3
"""
差异表达分析模块
使用DESeq2进行差异表达分析
"""

import os
import subprocess
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
import seaborn as sns

class DifferentialExpression:
    """差异表达分析器"""
    
    def __init__(self, config_path: str = "config/config.yaml"):
        """初始化差异表达分析器"""
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
    
    def count_reads(self, bam_files: List[str], gtf_file: str, output_file: str):
        """
        使用featureCounts计算基因表达量
        
        Args:
            bam_files: BAM文件列表
            gtf_file: GTF注释文件
            output_file: 输出文件路径
        """
        print("计算基因表达量...")
        
        # 构建featureCounts命令
        bam_str = " ".join(bam_files)
        cmd = f"""featureCounts -T 8 -p -t exon -g gene_id \\
            -a {gtf_file} -o {output_file} {bam_str}"""
        
        subprocess.run(cmd, shell=True, check=True)
        print(f"✓ 表达量计算完成: {output_file}")
    
    def create_deseq2_script(self, counts_file: str, output_dir: str = "results/tables/"):
        """
        创建DESeq2分析脚本
        
        Args:
            counts_file: 计数文件路径
            output_dir: 输出目录
        """
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        padj_threshold = self.config['analysis']['differential_expression']['padj_threshold']
        log2fc_threshold = self.config['analysis']['differential_expression']['log2fc_threshold']

        script_content = f'''
# DESeq2差异表达分析脚本
suppressPackageStartupMessages({{
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(EnhancedVolcano)
}})

padj_threshold <- {padj_threshold}
log2fc_threshold <- {log2fc_threshold}

# 读取计数数据
# 兼容两种输入：
# 1) featureCounts 输出（包含 Geneid, Chr, Start, End, Strand, Length + 样本列）
# 2) 纯计数矩阵（包含 Geneid + 样本列，或第一列为基因ID）
counts_raw <- read.delim("{counts_file}", header=TRUE, comment.char="#", check.names=FALSE)

if (nrow(counts_raw) == 0) {{
  stop("counts 文件为空或无法解析: {counts_file}")
}}

featurecounts_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")

if (all(featurecounts_cols %in% colnames(counts_raw))) {{
  gene_ids <- counts_raw$Geneid
  count_df <- counts_raw[, setdiff(colnames(counts_raw), featurecounts_cols), drop=FALSE]
}} else {{
  if ("Geneid" %in% colnames(counts_raw)) {{
    gene_ids <- counts_raw$Geneid
    count_df <- counts_raw[, setdiff(colnames(counts_raw), "Geneid"), drop=FALSE]
  }} else {{
    gene_ids <- counts_raw[[1]]
    count_df <- counts_raw[, -1, drop=FALSE]
  }}
}}

count_matrix <- as.matrix(count_df)
mode(count_matrix) <- "numeric"
count_matrix[is.na(count_matrix)] <- 0
count_matrix <- round(count_matrix)

rownames(count_matrix) <- make.unique(as.character(gene_ids))

# 创建样本信息：默认前一半 control，后一半 L1OE
n_samples <- ncol(count_matrix)
if (n_samples < 2) {{
  stop("样本列数不足，无法运行 DESeq2。至少需要 2 个样本列。")
}}

n_control <- floor(n_samples / 2)
n_treat <- n_samples - n_control
condition <- factor(c(rep("control", n_control), rep("L1OE", n_treat)))

coldata <- data.frame(condition = condition, row.names = colnames(count_matrix))

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)

# 预过滤低表达基因
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]

dds <- DESeq(dds)

res <- results(dds, contrast=c("condition", "L1OE", "control"))

# lfcShrink 的 apeglm/ashr 需要额外依赖；这里使用 "normal" 以保证默认可运行
res <- lfcShrink(dds, coef="condition_L1OE_vs_control", type="normal")

# 保存结果
write.csv(as.data.frame(res), file="{output_dir}/all_genes_deseq2.csv")

res_sig <- res[which(!is.na(res$padj) & res$padj < padj_threshold & abs(res$log2FoldChange) > log2fc_threshold), ]
write.csv(as.data.frame(res_sig), file="{output_dir}/significant_genes_deseq2.csv")

# 绘制火山图
png("{output_dir}/volcano_plot.png", width=1200, height=800, res=300)
EnhancedVolcano(as.data.frame(res),
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'padj',
    xlab = bquote(~Log[2]~ 'fold change'),
    ylab = bquote(~-Log[10]~ italic(P)),
    pCutoff = padj_threshold,
    FCcutoff = log2fc_threshold,
    pointSize = 2.0,
    labSize = 3.0,
    title = 'L1OE vs Control',
    subtitle = 'Differential Expression Analysis')
dev.off()

# 绘制热图（按 padj 取前 50）
res_order <- order(res$padj, na.last=NA)
top_genes <- head(res_order, 50)

if (length(top_genes) > 2) {{
  rld <- rlog(dds, blind=FALSE)
  mat <- assay(rld)[top_genes, , drop=FALSE]
  mat <- mat - rowMeans(mat)
  png("{output_dir}/heatmap_top50.png", width=1200, height=1000, res=300)
  pheatmap(mat,
           annotation_col=coldata,
           show_rownames=TRUE,
           show_colnames=TRUE,
           fontsize_row=8,
           fontsize_col=12,
           main="Top 50 Differential Expressed Genes")
  dev.off()
}}

# PCA分析
png("{output_dir}/pca_plot.png", width=1000, height=800, res=300)
vsd <- vst(dds, blind=FALSE)
print(plotPCA(vsd, intgroup="condition") +
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal())
dev.off()

cat("DESeq2分析完成!\\n")
cat("显著差异基因数量:", nrow(res_sig), "\\n")
'''
        
        script_path = f"{output_dir}/deseq2_analysis.R"
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(script_content)
        
        return script_path
    
    def run_deseq2_analysis(self, script_path: str):
        """
        运行DESeq2分析
        
        Args:
            script_path: R脚本路径
        """
        print("运行DESeq2分析...")
        
        cmd = f"Rscript {script_path}"
        subprocess.run(cmd, shell=True, check=True)
        
        print("✓ DESeq2分析完成")
    
    def filter_lncrnas(self, diff_file: str, gtf_file: str, output_file: str):
        """
        筛选差异表达的lncRNA
        
        Args:
            diff_file: 差异分析结果文件
            gtf_file: GTF注释文件
            output_file: 输出文件路径
        """
        print("筛选差异表达lncRNA...")
        
        # 读取差异分析结果
        diff_df = pd.read_csv(diff_file, index_col=0)
        
        # 从GTF文件提取lncRNA信息
        lncrnas = set()
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    attributes = fields[8]
                    if 'gene_biotype "lncRNA"' in attributes:
                        # 提取gene_id
                        for attr in attributes.split(';'):
                            if 'gene_id' in attr:
                                gene_id = attr.split('"')[1]
                                lncrnas.add(gene_id)
                                break
        
        # 筛选lncRNA
        lncrna_diff = diff_df[diff_df.index.isin(lncrnas)]
        
        # 保存结果
        lncrna_diff.to_csv(output_file)
        
        print(f"✓ 筛选完成，发现 {len(lncrna_diff)} 个差异表达lncRNA")
        
        return lncrna_diff
    
    def analyze_line1_association(self, lncrna_file: str, line1_bed: str, gtf_file: str, 
                                output_file: str, distance_threshold: int = 10000):
        """
        分析lncRNA与LINE-1的位置关联
        
        Args:
            lncrna_file: lncRNA差异表达文件
            line1_bed: LINE-1 BED文件
            gtf_file: GTF注释文件
            output_file: 输出文件
            distance_threshold: 距离阈值(bp)
        """
        print("分析lncRNA与LINE-1位置关联...")
        
        # 创建lncRNA BED文件
        lncrna_bed = "data/processed/lncrnas.bed"
        self._create_lncrna_bed(lncrna_file, gtf_file, lncrna_bed)
        
        # 使用bedtools分析位置关系
        cmd = f"""bedtools window -a {lncrna_bed} -b {line1_bed} -w {distance_threshold} > data/processed/lncrna_line1_overlap.txt"""
        subprocess.run(cmd, shell=True, check=True)
        
        # 分析结果
        overlap_df = pd.read_csv("data/processed/lncrna_line1_overlap.txt", 
                                sep='\t', header=None,
                                names=['lncRNA_chr', 'lncRNA_start', 'lncRNA_end', 'lncRNA_id',
                                       'LINE1_chr', 'LINE1_start', 'LINE1_end', 'LINE1_id'])
        
        # 计算距离和关联类型
        overlap_df['distance'] = overlap_df.apply(self._calculate_distance, axis=1)
        overlap_df['association_type'] = overlap_df.apply(self._determine_association_type, axis=1)
        
        # 保存结果
        overlap_df.to_csv(output_file, index=False)
        
        print(f"✓ 关联分析完成，发现 {len(overlap_df)} 个关联关系")
        
        return overlap_df
    
    def _create_lncrna_bed(self, lncrna_file: str, gtf_file: str, output_bed: str):
        """创建lncRNA BED文件"""
        lncrnas = set(pd.read_csv(lncrna_file, index_col=0).index)
        
        with open(output_bed, 'w') as out_f:
            with open(gtf_file, 'r') as in_f:
                for line in in_f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) >= 9:
                        attributes = fields[8]
                        if 'gene_biotype "lncRNA"' in attributes:
                            # 提取gene_id
                            gene_id = None
                            for attr in attributes.split(';'):
                                if 'gene_id' in attr:
                                    gene_id = attr.split('"')[1]
                                    break
                            
                            if gene_id and gene_id in lncrnas:
                                chrom = fields[0]
                                start = int(fields[3]) - 1  # BED格式是0-based
                                end = fields[4]
                                out_f.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t.\t{fields[6]}\n")
    
    def _calculate_distance(self, row):
        """计算lncRNA与LINE-1之间的距离"""
        lncrna_center = (row['lncRNA_start'] + row['lncRNA_end']) / 2
        line1_center = (row['LINE1_start'] + row['LINE1_end']) / 2
        return abs(lncrna_center - line1_center)
    
    def _determine_association_type(self, row):
        """确定关联类型"""
        # 检查是否重叠
        if (row['lncRNA_start'] <= row['LINE1_end'] and 
            row['lncRNA_end'] >= row['LINE1_start']):
            return "overlap"
        
        # 检查是否在LINE-1内部
        if (row['lncRNA_start'] >= row['LINE1_start'] and 
            row['lncRNA_end'] <= row['LINE1_end']):
            return "within"
        
        # 检查LINE-1是否在lncRNA内部
        if (row['LINE1_start'] >= row['lncRNA_start'] and 
            row['LINE1_end'] <= row['lncRNA_end']):
            return "contains"
        
        # 否则为邻近
        return "proximity"


def main():
    """主函数"""
    de_analyzer = DifferentialExpression()
    
    # 示例流程
    bam_files = ["data/processed/sample1.bam", "data/processed/sample2.bam"]
    gtf_file = "data/annotation/gencode.v45.annotation.gtf"
    
    # 计算表达量
    de_analyzer.count_reads(bam_files, gtf_file, "data/processed/gene_counts.txt")
    
    # 创建并运行DESeq2分析
    script_path = de_analyzer.create_deseq2_script("data/processed/gene_counts.txt")
    de_analyzer.run_deseq2_analysis(script_path)
    
    # 筛选lncRNA
    de_analyzer.filter_lncrnas("results/tables/significant_genes_deseq2.csv", 
                              gtf_file, "results/tables/lncrna_differential.csv")
    
    print("差异表达分析完成!")


if __name__ == "__main__":
    main()