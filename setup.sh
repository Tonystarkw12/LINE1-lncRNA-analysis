#!/bin/bash

# LINE-1 lncRNA项目环境设置脚本
# 用于快速配置项目环境和依赖

echo "=========================================="
echo "LINE-1 lncRNA项目环境设置"
echo "=========================================="

# 检查conda是否安装
if ! command -v conda &> /dev/null; then
    echo "错误: 未找到conda，请先安装Anaconda或Miniconda"
    exit 1
fi

# 创建conda环境
echo "正在创建conda环境..."
conda env create -f environment.yml

# 激活环境
echo "激活环境..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate line1_lncrna

# 检查R环境
echo "检查R环境..."
R --version

# 安装R包（如果需要）
echo "检查R包..."
Rscript -e "
if (!require('BiocManager', quietly = TRUE)) {
    install.packages('BiocManager')
}
BiocManager::install(c('DESeq2', 'GenomicFeatures', 'GenomicRanges', 
                      'AnnotationDbi', 'org.Hs.eg.db', 'clusterProfiler', 
                      'enrichplot'), ask = FALSE)
"

# 创建必要目录
echo "创建项目目录..."
mkdir -p data/{raw,processed,reference,annotation,line1}
mkdir -p results/{figures,tables,reports,networks}
mkdir -p logs

# 设置权限
chmod +x src/*.py

echo "=========================================="
echo "环境设置完成!"
echo "=========================================="
echo "使用方法:"
echo "1. 激活环境: conda activate line1_lncrna"
echo "2. 运行分析: jupyter notebook notebooks/LINE1_lncRNA_analysis.ipynb"
echo "3. 或直接运行: python src/data_download.py"
echo "=========================================="