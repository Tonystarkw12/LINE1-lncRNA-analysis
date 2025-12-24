#!/usr/bin/env python3
"""
工具函数模块
提供项目中常用的辅助函数
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import yaml
import logging
from typing import List, Dict, Tuple, Optional
import warnings

def setup_logging(log_file: str = "logs/project.log", level: str = "INFO"):
    """设置日志记录"""
    Path(log_file).parent.mkdir(parents=True, exist_ok=True)
    
    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    
    return logging.getLogger(__name__)

def load_config(config_path: str = "config/config.yaml") -> Dict:
    """加载配置文件"""
    with open(config_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)

def check_file_exists(file_path: str, description: str = "文件") -> bool:
    """检查文件是否存在"""
    if os.path.exists(file_path):
        print(f"✓ {description}存在: {file_path}")
        return True
    else:
        print(f"✗ {description}不存在: {file_path}")
        return False

def create_directory(dir_path: str) -> bool:
    """创建目录"""
    try:
        Path(dir_path).mkdir(parents=True, exist_ok=True)
        print(f"✓ 目录创建成功: {dir_path}")
        return True
    except Exception as e:
        print(f"✗ 目录创建失败: {e}")
        return False

def validate_dataframe(df: pd.DataFrame, required_columns: List[str], 
                       df_name: str = "DataFrame") -> bool:
    """验证DataFrame是否包含必需的列"""
    missing_columns = set(required_columns) - set(df.columns)
    if missing_columns:
        print(f"✗ {df_name}缺少必需列: {missing_columns}")
        return False
    else:
        print(f"✓ {df_name}验证通过，包含所有必需列")
        return True

def calculate_gene_statistics(df: pd.DataFrame, 
                            gene_col: str = "gene_id",
                            expression_cols: List[str] = None) -> Dict:
    """计算基因表达统计信息"""
    if expression_cols is None:
        expression_cols = [col for col in df.columns if col != gene_col]
    
    stats = {}
    for col in expression_cols:
        if col in df.columns:
            stats[col] = {
                'mean': df[col].mean(),
                'median': df[col].median(),
                'std': df[col].std(),
                'min': df[col].min(),
                'max': df[col].max(),
                'zeros': (df[col] == 0).sum(),
                'non_zeros': (df[col] > 0).sum()
            }
    
    return stats

def filter_genes_by_expression(df: pd.DataFrame, 
                              min_count: int = 10,
                              min_samples: int = 2,
                              gene_col: str = "gene_id") -> pd.DataFrame:
    """根据表达量过滤基因"""
    expression_cols = [col for col in df.columns if col != gene_col]
    
    # 计算每个基因在至少min_samples个样本中表达量大于min_count
    mask = (df[expression_cols] >= min_count).sum(axis=1) >= min_samples
    
    filtered_df = df[mask].copy()
    print(f"✓ 基因过滤完成: {len(df)} -> {len(filtered_df)}")
    
    return filtered_df

def normalize_expression_data(df: pd.DataFrame, 
                            method: str = "cpm",
                            gene_col: str = "gene_id") -> pd.DataFrame:
    """标准化表达数据"""
    expression_cols = [col for col in df.columns if col != gene_col]
    df_normalized = df.copy()
    
    if method == "cpm":
        # Counts per million
        lib_sizes = df[expression_cols].sum(axis=0)
        df_normalized[expression_cols] = df[expression_cols].div(lib_sizes, axis=1) * 1e6
    elif method == "tpm":
        # Transcripts per million (需要基因长度信息)
        warnings.warn("TPM标准化需要基因长度信息，当前使用CPM方法")
        lib_sizes = df[expression_cols].sum(axis=0)
        df_normalized[expression_cols] = df[expression_cols].div(lib_sizes, axis=1) * 1e6
    elif method == "quantile":
        # Quantile normalization
        from sklearn.preprocessing import QuantileTransformer
        qt = QuantileTransformer(output_distribution='normal')
        df_normalized[expression_cols] = qt.fit_transform(df[expression_cols])
    
    print(f"✓ 表达数据标准化完成，方法: {method}")
    return df_normalized

def plot_data_summary(df: pd.DataFrame, 
                     gene_col: str = "gene_id",
                     output_file: str = None,
                     figsize: Tuple[int, int] = (15, 10)):
    """绘制数据汇总图"""
    expression_cols = [col for col in df.columns if col != gene_col]
    
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    fig.suptitle('Expression Data Summary', fontsize=16)
    
    # 1. 样本表达量分布
    df[expression_cols].boxplot(ax=axes[0, 0])
    axes[0, 0].set_title('Sample Expression Distribution')
    axes[0, 0].tick_params(axis='x', rotation=45)
    
    # 2. 基因表达量分布
    gene_means = df[expression_cols].mean(axis=1)
    axes[0, 1].hist(gene_means, bins=50, alpha=0.7)
    axes[0, 1].set_title('Gene Mean Expression Distribution')
    axes[0, 1].set_xlabel('Mean Expression')
    axes[0, 1].set_ylabel('Gene Count')
    
    # 3. 零表达基因比例
    zero_ratios = (df[expression_cols] == 0).sum(axis=1) / len(expression_cols)
    axes[0, 2].hist(zero_ratios, bins=20, alpha=0.7)
    axes[0, 2].set_title('Zero Expression Ratio Distribution')
    axes[0, 2].set_xlabel('Zero Expression Ratio')
    axes[0, 2].set_ylabel('Gene Count')
    
    # 4. 样本间相关性热图
    corr_matrix = df[expression_cols].corr()
    sns.heatmap(corr_matrix, annot=True, cmap='coolwarm', center=0, 
                ax=axes[1, 0], cbar_kws={'shrink': 0.8})
    axes[1, 0].set_title('Sample Correlation Matrix')
    
    # 5. 主成分分析
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    # 标准化数据
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(df[expression_cols].T)
    
    # PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(scaled_data)
    
    axes[1, 1].scatter(pca_result[:, 0], pca_result[:, 1], alpha=0.7)
    axes[1, 1].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    axes[1, 1].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    axes[1, 1].set_title('PCA of Samples')
    
    # 添加样本标签
    for i, sample in enumerate(expression_cols):
        axes[1, 1].annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                           alpha=0.7, fontsize=8)
    
    # 6. 表达量统计摘要
    summary_stats = df[expression_cols].describe()
    axes[1, 2].axis('off')
    axes[1, 2].table(cellText=summary_stats.round(2).values,
                     rowLabels=summary_stats.index,
                     colLabels=summary_stats.columns,
                     cellLoc='center',
                     loc='center')
    axes[1, 2].set_title('Expression Statistics Summary')
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"✓ 数据汇总图已保存: {output_file}")
    
    plt.show()

def save_analysis_summary(results: Dict, output_file: str):
    """保存分析结果摘要"""
    summary_lines = ["LINE-1 lncRNA项目分析摘要", "=" * 40, ""]
    
    for key, value in results.items():
        if isinstance(value, dict):
            summary_lines.append(f"{key}:")
            for sub_key, sub_value in value.items():
                summary_lines.append(f"  {sub_key}: {sub_value}")
        else:
            summary_lines.append(f"{key}: {value}")
        summary_lines.append("")
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write("\\n".join(summary_lines))
    
    print(f"✓ 分析摘要已保存: {output_file}")

def check_memory_usage():
    """检查内存使用情况"""
    import psutil
    memory = psutil.virtual_memory()
    print(f"内存使用情况: {memory.percent}% ({memory.used/1024/1024/1024:.1f}GB / {memory.total/1024/1024/1024:.1f}GB)")

def format_file_size(size_bytes: int) -> str:
    """格式化文件大小"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f}{unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f}PB"

def get_file_info(file_path: str) -> Dict:
    """获取文件信息"""
    if not os.path.exists(file_path):
        return {"exists": False}
    
    stat = os.stat(file_path)
    return {
        "exists": True,
        "size": stat.st_size,
        "size_formatted": format_file_size(stat.st_size),
        "modified": pd.Timestamp.fromtimestamp(stat.st_mtime),
        "is_file": os.path.isfile(file_path),
        "is_dir": os.path.isdir(file_path)
    }

def create_progress_bar(total: int, description: str = "Processing"):
    """创建进度条"""
    from tqdm import tqdm
    return tqdm(total=total, desc=description)

def validate_sequence_format(sequence: str, seq_type: str = "DNA") -> bool:
    """验证序列格式"""
    if seq_type.upper() == "DNA":
        valid_bases = set("ATCGN")
    elif seq_type.upper() == "RNA":
        valid_bases = set("AUCGN")
    elif seq_type.upper() == "PROTEIN":
        valid_bases = set("ACDEFGHIKLMNPQRSTVWY*")
    else:
        return False
    
    sequence = sequence.upper()
    return all(base in valid_bases for base in sequence)

def reverse_complement(sequence: str) -> str:
    """计算反向互补序列"""
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(complement_dict.get(base, base) for base in reversed(sequence))

def calculate_gc_content(sequence: str) -> float:
    """计算GC含量"""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence) * 100) if len(sequence) > 0 else 0

def batch_process_files(file_paths: List[str], 
                       process_function, 
                       batch_size: int = 10,
                       **kwargs):
    """批量处理文件"""
    results = []
    
    for i in range(0, len(file_paths), batch_size):
        batch = file_paths[i:i + batch_size]
        print(f"处理批次 {i//batch_size + 1}/{(len(file_paths)-1)//batch_size + 1}")
        
        for file_path in batch:
            try:
                result = process_function(file_path, **kwargs)
                results.append(result)
            except Exception as e:
                print(f"处理文件失败 {file_path}: {e}")
                results.append(None)
    
    return results

# 常用常量
DNA_BASES = "ATCG"
RNA_BASES = "AUCG"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

if __name__ == "__main__":
    # 测试工具函数
    print("工具函数模块测试")
    
    # 测试配置加载
    config = load_config()
    print(f"配置加载成功: {list(config.keys())}")
    
    # 测试目录创建
    create_directory("test_output")
    
    # 测试文件信息
    if os.path.exists("config/config.yaml"):
        info = get_file_info("config/config.yaml")
        print(f"配置文件信息: {info}")
    
    print("✓ 工具函数测试完成")