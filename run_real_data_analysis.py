#!/usr/bin/env python3
"""
LINE-1 lncRNA 项目 - 真实RNA-seq数据分析流程
使用真实的GEO数据(GSE159217)进行分析
"""

import os
import sys
import subprocess
import warnings
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
from pathlib import Path
import logging
from datetime import datetime

# 添加src目录到路径
sys.path.append('src')

from differential_expression import DifferentialExpression
from network_analysis import NetworkAnalyzer
from visualization import Visualizer

# 设置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('logs/real_data_analysis.log'),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 配置参数
CONFIG = {
    'threads': 8,  # 多线程加速
    'data_dir': 'data/raw/fastq',
    'processed_dir': 'data/processed',
    'reference_dir': 'data/reference/hg38_hisat2',
    'genome_fasta': 'data/raw/hg38.fa',
    'gtf_file': 'data/raw/gencode.v45.annotation.gtf',
    'rmsk_file': 'data/raw/rmsk.txt',
    'sample_name': 'SRR12793631',
    'results_dir': 'results_real_data',
}

def run_command(cmd, description=""):
    """运行shell命令并处理错误"""
    logger.info(f"运行: {description}")
    logger.info(f"命令: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True,
                              capture_output=True, text=True)
        logger.info(f"✓ {description}完成")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"✗ {description}失败")
        logger.error(f"错误信息: {e.stderr}")
        return False

def step1_quality_control():
    """步骤1: 质量控制 (fastp)"""
    logger.info("="*60)
    logger.info("步骤 1: 质量控制 (fastp)")
    logger.info("="*60)

    # 检查是否已经运行过fastp
    clean_1 = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}_1_clean.fastq.gz"
    clean_2 = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}_2_clean.fastq.gz"

    if os.path.exists(clean_1) and os.path.exists(clean_2):
        logger.info(f"✓ 质控数据已存在: {clean_1}")
        return True

    # 运行fastp
    cmd = f"""fastp \
        -i {CONFIG['data_dir']}/{CONFIG['sample_name']}_1.fastq.gz \
        -I {CONFIG['data_dir']}/{CONFIG['sample_name']}_2.fastq.gz \
        -o {clean_1} \
        -O {clean_2} \
        -h fastp.html \
        -j fastp.json \
        -w {CONFIG['threads']} \
        --detect_adapter_for_pe"""

    return run_command(cmd, "fastp质量控制")

def step2_build_hisat2_index():
    """步骤2: 构建HISAT2索引"""
    logger.info("="*60)
    logger.info("步骤 2: 构建HISAT2基因组索引")
    logger.info("="*60)

    index_prefix = f"{CONFIG['reference_dir']}/hg38"
    index_files = [f"{index_prefix}.{i}.ht2" for i in range(1, 9)]

    if all(os.path.exists(f) for f in index_files):
        logger.info(f"✓ HISAT2索引已存在: {index_prefix}")
        return True

    os.makedirs(CONFIG['reference_dir'], exist_ok=True)

    cmd = f"hisat2-build -p {CONFIG['threads']} {CONFIG['genome_fasta']} {index_prefix}"
    return run_command(cmd, "构建HISAT2索引")

def step3_alignment():
    """步骤3: 序列比对 (HISAT2)"""
    logger.info("="*60)
    logger.info("步骤 3: 序列比对 (HISAT2)")
    logger.info("="*60)

    sam_file = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}.sam"
    bam_file = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}.bam"

    # 如果已经存在BAM文件,跳过
    if os.path.exists(bam_file):
        logger.info(f"✓ 比对结果已存在: {bam_file}")
        return True

    # HISAT2比对
    index_prefix = f"{CONFIG['reference_dir']}/hg38"
    clean_1 = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}_1_clean.fastq.gz"
    clean_2 = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}_2_clean.fastq.gz"

    cmd = f"""hisat2 -p {CONFIG['threads']} \
        -x {index_prefix} \
        -1 {clean_1} \
        -2 {clean_2} \
        -S {sam_file} \
        --rna-strandness FR \
        --ignore-quals"""

    if not run_command(cmd, "HISAT2序列比对"):
        return False

    # SAM转BAM并排序
    logger.info("将SAM转换为BAM并排序...")
    cmd = f"samtools view -@ {CONFIG['threads']} -bS {sam_file} | \
            samtools sort -@ {CONFIG['threads']} -o {bam_file}"
    if not run_command(cmd, "SAM转BAM并排序"):
        return False

    # 建立索引
    cmd = f"samtools index {bam_file}"
    if not run_command(cmd, "建立BAM索引"):
        return False

    # 删除SAM文件节省空间
    if os.path.exists(sam_file):
        os.remove(sam_file)
        logger.info(f"删除中间文件: {sam_file}")

    # 生成比对统计
    cmd = f"samtools flagstat {bam_file}"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    stats_file = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}_alignment_stats.txt"
    with open(stats_file, 'w') as f:
        f.write(result.stdout)
    logger.info(f"比对统计已保存: {stats_file}")
    logger.info(result.stdout)

    return True

def step4_quantification():
    """步骤4: 基因定量 (featureCounts)"""
    logger.info("="*60)
    logger.info("步骤 4: 基因定量 (featureCounts)")
    logger.info("="*60)

    bam_file = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}.bam"
    counts_file = f"{CONFIG['processed_dir']}/gene_counts_real.txt"

    if os.path.exists(counts_file):
        logger.info(f"✓ 定量结果已存在: {counts_file}")
        return True

    cmd = f"""featureCounts -T {CONFIG['threads']} -p -B -C \
        -a {CONFIG['gtf_file']} \
        -o {counts_file} \
        {bam_file}"""

    return run_command(cmd, "featureCounts基因定量")

def step5_differential_expression():
    """步骤5: 差异表达分析 (Python模拟DESeq2)"""
    logger.info("="*60)
    logger.info("步骤 5: 差异表达分析")
    logger.info("="*60)

    # 读取featureCounts结果
    counts_file = f"{CONFIG['processed_dir']}/gene_counts_real.txt"

    if not os.path.exists(counts_file):
        logger.error(f"定量文件不存在: {counts_file}")
        return False

    # 读取counts
    logger.info("读取基因表达数据...")
    counts_df = pd.read_csv(counts_file, sep='\t', comment='#',
                          skiprows=1, index_col=0)

    # 提取GeneID和Counts列
    # featureCounts输出格式: Geneid Chr Start End Strand Length (Sample)
    count_col = [col for col in counts_df.columns if 'SRR' in col or 'bam' in col.lower()][0]

    # 创建表达矩阵
    expression_df = pd.DataFrame({
        'Geneid': counts_df.index,
        'Length': counts_df['Length'],
        'counts': counts_df[count_col]
    })

    # 过滤低表达基因 (counts >= 10)
    expression_df = expression_df[expression_df['counts'] >= 10]
    logger.info(f"过滤后基因数: {len(expression_df)}")

    # 计算RPKM/FPKM
    total_reads = expression_df['counts'].sum()
    gene_lengths_kb = expression_df['Length'] / 1000
    library_size_millions = total_reads / 1e6

    expression_df['fpkm'] = expression_df['counts'] / (gene_lengths_kb * library_size_millions)

    # 识别高表达lncRNA
    # 从GTF中提取lncRNA基因
    logger.info("识别lncRNA基因...")

    lncrna_genes = set()
    with open(CONFIG['gtf_file'], 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'gene':
                attributes = fields[8]
                if 'gene_biotype "lncRNA"' in attributes or 'gene_type "lncRNA"' in attributes:
                    gene_id = None
                    for attr in attributes.split(';'):
                        if 'gene_id' in attr:
                            gene_id = attr.split('"')[1]
                            break
                    if gene_id:
                        lncrna_genes.add(gene_id)

    logger.info(f"发现 {len(lncrna_genes)} 个lncRNA基因")

    # 筛选lncRNA
    lncrna_df = expression_df[expression_df['Geneid'].isin(lncrna_genes)].copy()
    logger.info(f"表达数据中包含 {len(lncrna_df)} 个lncRNA")

    # 保存结果
    os.makedirs(f"{CONFIG['results_dir']}/tables", exist_ok=True)

    # 保存所有基因表达数据
    expression_df.to_csv(f"{CONFIG['results_dir']}/tables/all_genes_expression.csv", index=False)

    # 保存lncRNA表达数据
    lncrna_df.to_csv(f"{CONFIG['results_dir']}/tables/lncrna_expression.csv", index=False)

    # 识别高表达lncRNA (FPKM > 1)
    highly_expressed_lncrna = lncrna_df[lncrna_df['fpkm'] > 1].sort_values('fpkm', ascending=False)
    highly_expressed_lncrna.to_csv(f"{CONFIG['results_dir']}/tables/highly_expressed_lncrna.csv", index=False)

    logger.info(f"高表达lncRNA (FPKM>1): {len(highly_expressed_lncrna)}")
    logger.info(f"\nTop 10 高表达lncRNA:")
    for idx, row in highly_expressed_lncrna.head(10).iterrows():
        logger.info(f"  {row['Geneid']}: FPKM={row['fpkm']:.2f}, Counts={row['counts']}")

    return True, highly_expressed_lncrna

def step6_line1_association_analysis(lncrna_df):
    """步骤6: LINE-1关联分析"""
    logger.info("="*60)
    logger.info("步骤 6: LINE-1关联分析")
    logger.info("="*60)

    # 读取rmsk.txt获取LINE-1元件位置
    logger.info("读取LINE-1元件注释...")

    line1_elements = []
    with open(CONFIG['rmsk_file'], 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                # rmsk实际格式: (前几列是统计信息) ... genoName genoStart genoEnd strand repName repClass repFamily ...
                geno_name = fields[5]
                geno_start = int(fields[6])
                geno_end = int(fields[7])
                strand = fields[8]
                rep_name = fields[9]
                rep_class = fields[10]
                rep_family = fields[11]

                if rep_class == 'LINE' and rep_family.startswith('L1'):
                    line1_elements.append({
                        'chr': geno_name,
                        'start': geno_start,
                        'end': geno_end,
                        'strand': strand,
                        'name': rep_name,
                        'family': rep_family
                    })

    logger.info(f"发现 {len(line1_elements)} 个LINE-1元件")

    # 读取lncRNA位置信息
    logger.info("读取lncRNA位置信息...")

    lncrna_coords = {}
    with open(CONFIG['gtf_file'], 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'gene':
                attributes = fields[8]
                if 'gene_biotype "lncRNA"' in attributes or 'gene_type "lncRNA"' in attributes:
                    gene_id = None
                    for attr in attributes.split(';'):
                        if 'gene_id' in attr:
                            gene_id = attr.split('"')[1]
                            break
                    if gene_id and gene_id in lncrna_df['Geneid'].values:
                        lncrna_coords[gene_id] = {
                            'chr': fields[0],
                            'start': int(fields[3]),
                            'end': int(fields[4]),
                            'strand': fields[6]
                        }

    # 分析lncRNA与LINE-1的关联
    logger.info("分析lncRNA与LINE-1的位置关联...")

    associations = []
    proximity_threshold = 50000  # 50kb范围内认为有关联

    for lncrna_id in lncrna_df['Geneid']:
        if lncrna_id not in lncrna_coords:
            continue

        lncrna_coord = lncrna_coords[lncrna_id]

        # 检查每个LINE-1元件
        for line1 in line1_elements:
            if line1['chr'] != lncrna_coord['chr']:
                continue

            # 检查重叠
            if not (line1['end'] < lncrna_coord['start'] or line1['start'] > lncrna_coord['end']):
                associations.append({
                    'lncRNA_id': lncrna_id,
                    'line1_name': line1['name'],
                    'line1_family': line1['family'],
                    'chr': line1['chr'],
                    'lncrna_start': lncrna_coord['start'],
                    'lncrna_end': lncrna_coord['end'],
                    'line1_start': line1['start'],
                    'line1_end': line1['end'],
                    'association_type': 'overlap',
                    'distance': 0
                })
            # 检查临近
            else:
                distance = min(abs(line1['end'] - lncrna_coord['start']),
                             abs(line1['start'] - lncrna_coord['end']))
                if distance <= proximity_threshold:
                    associations.append({
                        'lncRNA_id': lncrna_id,
                        'line1_name': line1['name'],
                        'line1_family': line1['family'],
                        'chr': line1['chr'],
                        'lncrna_start': lncrna_coord['start'],
                        'lncrna_end': lncrna_coord['end'],
                        'line1_start': line1['start'],
                        'line1_end': line1['end'],
                        'association_type': 'proximal',
                        'distance': distance
                    })

    # 保存结果
    associations_df = pd.DataFrame(associations)
    if len(associations_df) > 0:
        associations_df.to_csv(f"{CONFIG['results_dir']}/tables/lncrna_line1_associations.csv", index=False)
        logger.info(f"✓ 发现 {len(associations_df)} 个lncRNA-LINE1关联")
        logger.info(f"  - 关联lncRNA数: {associations_df['lncRNA_id'].nunique()}")
        logger.info(f"  - 重叠关联: {len(associations_df[associations_df['association_type']=='overlap'])}")
        logger.info(f"  - 临近关联: {len(associations_df[associations_df['association_type']=='proximal'])}")
    else:
        logger.warning("未发现lncRNA-LINE1关联")
        # 创建空DataFrame
        associations_df = pd.DataFrame(columns=['lncRNA_id', 'line1_name', 'line1_family',
                                               'chr', 'lncrna_start', 'lncrna_end',
                                               'line1_start', 'line1_end',
                                               'association_type', 'distance'])
        associations_df.to_csv(f"{CONFIG['results_dir']}/tables/lncrna_line1_associations.csv", index=False)

    return associations_df

def step7_visualization(lncrna_df, associations_df):
    """步骤7: 可视化"""
    logger.info("="*60)
    logger.info("步骤 7: 生成可视化结果")
    logger.info("="*60)

    os.makedirs(f"{CONFIG['results_dir']}/figures", exist_ok=True)

    # 1. lncRNA表达分布图
    fig, axes = plt.subplots(1, 2, figsize=(15, 5))

    # FPKM分布
    axes[0].hist(np.log10(lncrna_df['fpkm'] + 0.01), bins=50, color='skyblue', edgecolor='black')
    axes[0].set_xlabel('Log10(FPKM+0.01)')
    axes[0].set_ylabel('Number of lncRNAs')
    axes[0].set_title('lncRNA Expression Distribution')
    axes[0].axvline(np.log10(1.01), color='red', linestyle='--', label='FPKM=1')
    axes[0].legend()

    # Top 20 高表达lncRNA
    top20 = lncrna_df.nlargest(20, 'fpkm')
    axes[1].barh(range(len(top20)), top20['fpkm'], color='coral')
    axes[1].set_yticks(range(len(top20)))
    axes[1].set_yticklabels([gid.split('.')[0] for gid in top20['Geneid']], fontsize=8)
    axes[1].set_xlabel('FPKM')
    axes[1].set_title('Top 20 Expressed lncRNAs')
    axes[1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(f"{CONFIG['results_dir']}/figures/lncrna_expression.png", dpi=300, bbox_inches='tight')
    plt.close()
    logger.info("✓ lncRNA表达分布图已生成")

    # 2. LINE-1关联统计
    if len(associations_df) > 0:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))

        # 关联类型分布
        assoc_counts = associations_df['association_type'].value_counts()
        axes[0].pie(assoc_counts, labels=assoc_counts.index, autopct='%1.1f%%', startangle=90)
        axes[0].set_title('LINE-1 Association Types')

        # LINE-1家族分布
        family_counts = associations_df['line1_family'].value_counts().head(10)
        axes[1].bar(range(len(family_counts)), family_counts.values, color='lightgreen')
        axes[1].set_xticks(range(len(family_counts)))
        axes[1].set_xticklabels(family_counts.index, rotation=45, ha='right')
        axes[1].set_ylabel('Count')
        axes[1].set_title('Top 10 LINE-1 Families')
        plt.tight_layout()
        plt.savefig(f"{CONFIG['results_dir']}/figures/line1_associations.png", dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("✓ LINE-1关联统计图已生成")

    # 3. 基因组分布图
    if len(associations_df) > 0:
        # 按染色体统计
        chr_counts = associations_df['chr'].value_counts().sort_index()

        # 过滤掉非标准染色体
        chr_counts = chr_counts[chr_counts.index.str.match(r'chr\d+|chrX|chrY')]

        plt.figure(figsize=(12, 5))
        plt.bar(range(len(chr_counts)), chr_counts.values, color='steelblue')
        plt.xticks(range(len(chr_counts)), [c.replace('chr', '') for c in chr_counts.index], rotation=45)
        plt.xlabel('Chromosome')
        plt.ylabel('Number of Associations')
        plt.title('lncRNA-LINE1 Associations by Chromosome')
        plt.tight_layout()
        plt.savefig(f"{CONFIG['results_dir']}/figures/chromosome_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        logger.info("✓ 染色体分布图已生成")

    return True

def step8_generate_report(lncrna_df, associations_df):
    """步骤8: 生成分析报告"""
    logger.info("="*60)
    logger.info("步骤 8: 生成分析报告")
    logger.info("="*60)

    os.makedirs(f"{CONFIG['results_dir']}/reports", exist_ok=True)

    # 读取比对统计
    stats_file = f"{CONFIG['processed_dir']}/{CONFIG['sample_name']}_alignment_stats.txt"
    alignment_stats = ""
    if os.path.exists(stats_file):
        with open(stats_file, 'r') as f:
            alignment_stats = f.read()

    # 生成报告
    report = f"""
LINE-1 lncRNA 项目 - 真实数据分析报告
{'='*60}

分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
数据来源: GSE159217 (GEO数据库)
样本: {CONFIG['sample_name']}

分析流程:
{'-'*60}
1. 质量控制: fastp
2. 序列比对: HISAT2 (hg38)
3. 基因定量: featureCounts (GENCODE v45)
4. 差异表达分析: 基于FPKM
5. LINE-1关联分析: 基于基因组位置

比对统计:
{'-'*60}
{alignment_stats}

表达分析结果:
{'-'*60}
总lncRNA基因数: {len(lncrna_df)}
高表达lncRNA (FPKM>1): {len(lncrna_df[lncrna_df['fpkm'] > 1])}
平均FPKM: {lncrna_df['fpkm'].mean():.2f}
中位数FPKM: {lncrna_df['fpkm'].median():.2f}

LINE-1关联分析结果:
{'-'*60}
"""

    if len(associations_df) > 0:
        report += f"""
总关联数: {len(associations_df)}
关联lncRNA数: {associations_df['lncRNA_id'].nunique()}
重叠关联: {len(associations_df[associations_df['association_type']=='overlap'])}
临近关联: {len(associations_df[associations_df['association_type']=='proximal'])}

Top 10 LINE-1家族:
{'-'*60}
"""
        family_counts = associations_df['line1_family'].value_counts().head(10)
        for family, count in family_counts.items():
            report += f"  {family}: {count}\n"

        report += f"""
高优先级lncRNA (高表达且与LINE-1关联):
{'-'*60}
"""
        # 找出高表达且与LINE-1关联的lncRNA
        associated_lncrnas = set(associations_df['lncRNA_id'].unique())
        high_expr_lncrna = lncrna_df[lncrna_df['fpkm'] > 1]
        priority_lncrna = high_expr_lncrna[high_expr_lncrna['Geneid'].isin(associated_lncrnas)].sort_values('fpkm', ascending=False)

        for idx, row in priority_lncrna.head(20).iterrows():
            assoc_count = len(associations_df[associations_df['lncRNA_id'] == row['Geneid']])
            report += f"  {row['Geneid']}: FPKM={row['fpkm']:.2f}, 关联数={assoc_count}\n"
    else:
        report += "未发现lncRNA-LINE1关联\n"

    report += f"""

科学发现:
{'-'*60}
"""

    if len(associations_df) > 0:
        # 计算关联lncRNA的比例
        assoc_ratio = associations_df['lncRNA_id'].nunique() / len(lncrna_df) * 100
        report += f"1. {assoc_ratio:.1f}% 的lncRNA与LINE-1元件存在位置关联\n"
        report += f"2. {len(associations_df[associations_df['association_type']=='overlap'])} 个lncRNA直接与LINE-1重叠\n"
        report += f"3. 提示这些lncRNA可能通过LINE-1元件发挥调控作用\n"

    report += f"""

后续分析建议:
{'-'*60}
1. 如果有多个样本,进行差异表达分析(DESeq2)
2. 对高优先级lncRNA进行功能富集分析
3. 构建lncRNA-mRNA共表达网络
4. 进行保守性分析(PhyloP)
5. 预测lncRNA二级结构(RNAfold)

{'='*60}
生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
"""

    # 保存报告
    report_file = f"{CONFIG['results_dir']}/reports/analysis_report.txt"
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write(report)

    logger.info(f"✓ 分析报告已生成: {report_file}")

    # 打印报告到控制台
    logger.info("\n" + report)

    return True

def main():
    """主函数"""
    logger.info("="*60)
    logger.info("LINE-1 lncRNA 项目 - 真实数据分析开始")
    logger.info("="*60)

    # 创建必要的目录
    os.makedirs(CONFIG['processed_dir'], exist_ok=True)
    os.makedirs(f"{CONFIG['results_dir']}/tables", exist_ok=True)
    os.makedirs(f"{CONFIG['results_dir']}/figures", exist_ok=True)
    os.makedirs(f"{CONFIG['results_dir']}/reports", exist_ok=True)
    os.makedirs('logs', exist_ok=True)

    # 步骤1: 质量控制
    if not step1_quality_control():
        logger.error("质量控制失败,退出")
        return False

    # 步骤2: 构建HISAT2索引
    if not step2_build_hisat2_index():
        logger.error("构建索引失败,退出")
        return False

    # 步骤3: 序列比对
    if not step3_alignment():
        logger.error("序列比对失败,退出")
        return False

    # 步骤4: 基因定量
    if not step4_quantification():
        logger.error("基因定量失败,退出")
        return False

    # 步骤5: 差异表达分析(单样本情况下只做表达分析)
    success, lncrna_df = step5_differential_expression()
    if not success:
        logger.error("差异表达分析失败,退出")
        return False

    # 步骤6: LINE-1关联分析
    associations_df = step6_line1_association_analysis(lncrna_df)

    # 步骤7: 可视化
    step7_visualization(lncrna_df, associations_df)

    # 步骤8: 生成报告
    step8_generate_report(lncrna_df, associations_df)

    logger.info("="*60)
    logger.info("分析完成!")
    logger.info("="*60)
    logger.info(f"\n主要输出目录: {CONFIG['results_dir']}/")
    logger.info(f"  - tables/: 表格数据")
    logger.info(f"  - figures/: 图片结果")
    logger.info(f"  - reports/: 分析报告")

    return True

if __name__ == "__main__":
    try:
        success = main()
        sys.exit(0 if success else 1)
    except Exception as e:
        logger.error(f"分析过程中发生错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
