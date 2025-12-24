#!/usr/bin/env python3
"""
可视化模块
生成各种分析结果的可视化图表
"""

import os
import pandas as pd
import numpy as np
import yaml
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from pathlib import Path
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Arial Unicode MS']
plt.rcParams['axes.unicode_minus'] = False

class Visualizer:
    """可视化器"""
    
    def __init__(self, config_path: str = "config/config.yaml"):
        """初始化可视化器"""
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        # 设置绘图样式
        self.setup_style()
        
        # 创建输出目录
        Path(self.config['output']['figures_dir']).mkdir(parents=True, exist_ok=True)
    
    def setup_style(self):
        """设置绘图样式"""
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
        # 设置图形大小和DPI
        self.figure_size = self.config['visualization']['figure_size']
        self.figure_dpi = self.config['visualization']['figure_dpi']
        self.figure_format = self.config['visualization']['figure_format']
    
    def plot_volcano(self, diff_data: str, output_file: str = None):
        """
        绘制火山图
        
        Args:
            diff_data: 差异分析结果文件
            output_file: 输出文件路径
        """
        print("绘制火山图...")
        
        # 读取数据
        df = pd.read_csv(diff_data, index_col=0)
        
        # 设置阈值
        padj_threshold = self.config['analysis']['differential_expression']['padj_threshold']
        log2fc_threshold = self.config['analysis']['differential_expression']['log2fc_threshold']
        
        # 分类基因
        df['significance'] = 'not_significant'
        df.loc[(df['padj'] < padj_threshold) & (df['log2FoldChange'] > log2fc_threshold), 'significance'] = 'up'
        df.loc[(df['padj'] < padj_threshold) & (df['log2FoldChange'] < -log2fc_threshold), 'significance'] = 'down'
        
        # 创建图形
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        # 绘制散点
        colors = {'up': '#ff6b6b', 'down': '#4dabf7', 'not_significant': '#868e96'}
        for sig_type in colors:
            subset = df[df['significance'] == sig_type]
            ax.scatter(subset['log2FoldChange'], -np.log10(subset['padj']), 
                      c=colors[sig_type], alpha=0.6, s=20, label=sig_type)
        
        # 添加阈值线
        ax.axhline(y=-np.log10(padj_threshold), color='red', linestyle='--', alpha=0.7)
        ax.axvline(x=log2fc_threshold, color='red', linestyle='--', alpha=0.7)
        ax.axvline(x=-log2fc_threshold, color='red', linestyle='--', alpha=0.7)
        
        # 设置标签和标题
        ax.set_xlabel('Log2 Fold Change', fontsize=12)
        ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
        ax.set_title('Volcano Plot: Differential Expression Analysis', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 添加统计信息
        up_count = len(df[df['significance'] == 'up'])
        down_count = len(df[df['significance'] == 'down'])
        ax.text(0.02, 0.98, f'Up: {up_count}\\nDown: {down_count}', 
                transform=ax.transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        plt.tight_layout()
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/volcano_plot.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 火山图已保存: {output_file}")
    
    def plot_heatmap(self, expression_data: str, top_genes: int = 50, 
                    output_file: str = None):
        """
        绘制热图
        
        Args:
            expression_data: 表达量数据文件
            top_genes: 显示的基因数量
            output_file: 输出文件路径
        """
        print("绘制热图...")
        
        # 读取数据
        df = pd.read_csv(expression_data, index_col=0)
        
        # 选择变异最大的基因
        gene_var = df.var(axis=1)
        top_var_genes = gene_var.nlargest(top_genes).index
        df_subset = df.loc[top_var_genes]
        
        # 标准化数据
        df_scaled = (df_subset - df_subset.mean(axis=1).values.reshape(-1, 1)) / df_subset.std(axis=1).values.reshape(-1, 1)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(15, 10))
        
        # 绘制热图
        sns.heatmap(df_scaled, cmap='RdBu_r', center=0, 
                   cbar_kws={'label': 'Z-score'}, ax=ax)
        
        ax.set_title(f'Expression Heatmap - Top {top_genes} Variable Genes', 
                    fontsize=14, fontweight='bold')
        ax.set_xlabel('Samples', fontsize=12)
        ax.set_ylabel('Genes', fontsize=12)
        
        plt.tight_layout()
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/heatmap.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 热图已保存: {output_file}")
    
    def plot_network(self, network_file: str, output_file: str = None, 
                    layout: str = 'spring'):
        """
        绘制网络图
        
        Args:
            network_file: 网络文件路径
            output_file: 输出文件路径
            layout: 布局算法
        """
        print("绘制网络图...")
        
        # 读取网络
        G = nx.read_edgelist(network_file)
        
        # 设置布局
        if layout == 'spring':
            pos = nx.spring_layout(G, k=0.15, iterations=50)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        else:
            pos = nx.random_layout(G)
        
        # 创建图形
        fig, ax = plt.subplots(figsize=self.figure_size)
        
        # 获取节点类型
        node_colors = []
        for node in G.nodes():
            node_type = G.nodes[node].get('type', 'unknown')
            if node_type == 'lncRNA':
                node_colors.append('#ff6b6b')
            elif node_type == 'mRNA':
                node_colors.append('#4dabf7')
            elif node_type == 'LINE1':
                node_colors.append('#51cf66')
            else:
                node_colors.append('#868e96')
        
        # 绘制网络
        nx.draw(G, pos, ax=ax, 
               node_color=node_colors, 
               node_size=50,
               edge_color='gray',
               alpha=0.7,
               with_labels=False,
               width=0.5)
        
        # 添加图例
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#ff6b6b', label='lncRNA'),
            Patch(facecolor='#4dabf7', label='mRNA'),
            Patch(facecolor='#51cf66', label='LINE-1'),
            Patch(facecolor='#868e96', label='Other')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        ax.set_title('Regulatory Network: lncRNA-mRNA-LINE-1', 
                    fontsize=14, fontweight='bold')
        ax.axis('off')
        
        plt.tight_layout()
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/network.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 网络图已保存: {output_file}")
    
    def plot_interactive_network(self, network_file: str, output_file: str = None):
        """
        绘制交互式网络图
        
        Args:
            network_file: 网络文件路径
            output_file: 输出文件路径
        """
        print("绘制交互式网络图...")
        
        # 读取网络
        G = nx.read_edgelist(network_file)
        
        # 设置布局
        pos = nx.spring_layout(G, k=0.15)
        
        # 提取边和节点信息
        edge_x = []
        edge_y = []
        for edge in G.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
        
        # 创建边的轨迹
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines')
        
        # 提取节点信息
        node_x = []
        node_y = []
        node_text = []
        node_colors = []
        
        for node in G.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
            
            # 节点信息
            node_type = G.nodes[node].get('type', 'unknown')
            degree = G.degree(node)
            node_text.append(f'{node}<br>Type: {node_type}<br>Degree: {degree}')
            
            # 节点颜色
            if node_type == 'lncRNA':
                node_colors.append('#ff6b6b')
            elif node_type == 'mRNA':
                node_colors.append('#4dabf7')
            elif node_type == 'LINE1':
                node_colors.append('#51cf66')
            else:
                node_colors.append('#868e96')
        
        # 创建节点的轨迹
        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers',
            hoverinfo='text',
            text=node_text,
            marker=dict(
                showscale=True,
                colorscale='YlGnBu',
                reversescale=True,
                color=[],
                size=10,
                colorbar=dict(
                    thickness=15,
                    len=0.5,
                    x=1.02,
                    title='Node Connections'
                ),
                line_width=2))
        
        # 设置节点颜色
        node_adjacencies = []
        for node, adjacencies in enumerate(G.adjacency()):
            node_adjacencies.append(len(adjacencies[1]))
        
        node_trace.marker.color = node_adjacencies
        node_trace.marker.colorscale = 'Viridis'
        
        # 创建图形
        fig = go.Figure(data=[edge_trace, node_trace],
                       layout=go.Layout(
                           title='Interactive Regulatory Network',
                           titlefont_size=16,
                           showlegend=False,
                           hovermode='closest',
                           margin=dict(b=20,l=5,r=5,t=40),
                           annotations=[ dict(
                               text="Network visualization showing lncRNA-mRNA-LINE-1 interactions",
                               showarrow=False,
                               xref="paper", yref="paper",
                               x=0.005, y=-0.002,
                               xanchor='left', yanchor='bottom',
                               font=dict(color="#888", size=12)
                           )],
                           xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                           yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                       )
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/interactive_network.html"
        fig.write_html(output_file)
        
        print(f"✓ 交互式网络图已保存: {output_file}")
    
    def plot_conservation_scores(self, conservation_data: str, output_file: str = None):
        """
        绘制保守性得分分布
        
        Args:
            conservation_data: 保守性数据文件
            output_file: 输出文件路径
        """
        print("绘制保守性得分分布...")
        
        # 读取数据
        df = pd.read_csv(conservation_data)
        
        # 创建子图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. 平均保守性得分分布
        axes[0, 0].hist(df['avg_phylop'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
        axes[0, 0].set_xlabel('Average PhyloP Score')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Distribution of Average Conservation Scores')
        
        # 2. 保守性vs基因长度
        axes[0, 1].scatter(df['length'], df['avg_phylop'], alpha=0.6, s=30)
        axes[0, 1].set_xlabel('Gene Length (bp)')
        axes[0, 1].set_ylabel('Average PhyloP Score')
        axes[0, 1].set_title('Conservation vs Gene Length')
        
        # 3. 高保守区域比例
        axes[1, 0].hist(df['high_cons_ratio'], bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
        axes[1, 0].set_xlabel('High Conservation Ratio')
        axes[1, 0].set_ylabel('Frequency')
        axes[1, 0].set_title('Distribution of High Conservation Regions')
        
        # 4. 保守性统计箱线图
        cons_metrics = ['avg_phylop', 'max_phylop', 'min_phylop', 'std_phylop']
        df_box = df[cons_metrics].melt(var_name='Metric', value_name='Score')
        sns.boxplot(data=df_box, x='Metric', y='Score', ax=axes[1, 1])
        axes[1, 1].set_title('Conservation Metrics Distribution')
        axes[1, 1].tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/conservation_analysis.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 保守性分析图已保存: {output_file}")
    
    def plot_sequence_features(self, sequence_data: str, output_file: str = None):
        """
        绘制序列特征分析图
        
        Args:
            sequence_data: 序列特征数据文件
            output_file: 输出文件路径
        """
        print("绘制序列特征分析图...")
        
        # 读取数据
        df = pd.read_csv(sequence_data)
        
        # 创建子图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. GC含量分布
        axes[0, 0].hist(df['gc_content'], bins=30, alpha=0.7, color='lightgreen', edgecolor='black')
        axes[0, 0].set_xlabel('GC Content (%)')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('GC Content Distribution')
        axes[0, 0].axvline(df['gc_content'].mean(), color='red', linestyle='--', 
                          label=f'Mean: {df["gc_content"].mean():.1f}%')
        axes[0, 0].legend()
        
        # 2. 基因长度分布
        axes[0, 1].hist(df['length'], bins=30, alpha=0.7, color='lightblue', edgecolor='black')
        axes[0, 1].set_xlabel('Gene Length (bp)')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Gene Length Distribution')
        
        # 3. GC含量vs基因长度
        axes[1, 0].scatter(df['length'], df['gc_content'], alpha=0.6, s=30)
        axes[1, 0].set_xlabel('Gene Length (bp)')
        axes[1, 0].set_ylabel('GC Content (%)')
        axes[1, 0].set_title('GC Content vs Gene Length')
        
        # 4. 重复序列比例
        if 'repeat_ratio' in df.columns:
            axes[1, 1].hist(df['repeat_ratio'], bins=30, alpha=0.7, color='lightcoral', edgecolor='black')
            axes[1, 1].set_xlabel('Repeat Sequence Ratio')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title('Repeat Sequence Distribution')
        else:
            axes[1, 1].text(0.5, 0.5, 'No repeat data available', 
                           ha='center', va='center', transform=axes[1, 1].transAxes)
            axes[1, 1].set_title('Repeat Sequence Analysis')
        
        plt.tight_layout()
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/sequence_features.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 序列特征分析图已保存: {output_file}")
    
    def plot_secondary_structure_summary(self, structure_data: str, output_file: str = None):
        """
        绘制二级结构分析汇总图
        
        Args:
            structure_data: 二级结构数据文件
            output_file: 输出文件路径
        """
        print("绘制二级结构分析汇总图...")
        
        # 读取数据
        df = pd.read_csv(structure_data)
        
        # 创建子图
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # 1. 自由能分布
        axes[0, 0].hist(df['free_energy'], bins=30, alpha=0.7, color='orange', edgecolor='black')
        axes[0, 0].set_xlabel('Free Energy (kcal/mol)')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title('Distribution of Free Energy')
        axes[0, 0].axvline(df['free_energy'].mean(), color='red', linestyle='--',
                          label=f'Mean: {df["free_energy"].mean():.2f}')
        axes[0, 0].legend()
        
        # 2. 配对比例分布
        axes[0, 1].hist(df['pairing_ratio'], bins=30, alpha=0.7, color='purple', edgecolor='black')
        axes[0, 1].set_xlabel('Base Pairing Ratio')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title('Distribution of Base Pairing Ratio')
        
        # 3. 自由能vs配对比例
        axes[1, 0].scatter(df['pairing_ratio'], df['free_energy'], alpha=0.6, s=30)
        axes[1, 0].set_xlabel('Base Pairing Ratio')
        axes[1, 0].set_ylabel('Free Energy (kcal/mol)')
        axes[1, 0].set_title('Free Energy vs Base Pairing Ratio')
        
        # 4. 茎环结构数量分布
        if 'stem_loops' in df.columns:
            axes[1, 1].hist(df['stem_loops'], bins=20, alpha=0.7, color='brown', edgecolor='black')
            axes[1, 1].set_xlabel('Number of Stem-Loops')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title('Distribution of Stem-Loop Structures')
        else:
            axes[1, 1].text(0.5, 0.5, 'No stem-loop data available',
                           ha='center', va='center', transform=axes[1, 1].transAxes)
            axes[1, 1].set_title('Stem-Loop Analysis')
        
        plt.tight_layout()
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/secondary_structure.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 二级结构分析图已保存: {output_file}")
    
    def create_summary_dashboard(self, output_file: str = None):
        """
        创建项目汇总仪表板
        
        Args:
            output_file: 输出文件路径
        """
        print("创建项目汇总仪表板...")
        
        # 创建子图
        fig = plt.figure(figsize=(20, 15))
        
        # 定义网格布局
        gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
        
        # 添加标题
        fig.suptitle('LINE-1 lncRNA Project Summary Dashboard', 
                    fontsize=20, fontweight='bold', y=0.95)
        
        # 各个子图的标题（示例）
        subplot_titles = [
            'Differential Expression Analysis',
            'Network Topology',
            'Conservation Scores',
            'Sequence Features',
            'Secondary Structures',
            'Pathway Enrichment',
            'LINE-1 Associations',
            'Sample Clustering',
            'Key Regulators'
        ]
        
        # 创建占位子图
        for i, title in enumerate(subplot_titles):
            row = i // 3
            col = i % 3
            ax = fig.add_subplot(gs[row, col])
            
            # 创建示例数据或占位内容
            if i == 0:  # 差异分析
                x = np.random.randn(100)
                y = np.random.randn(100)
                ax.scatter(x, y, alpha=0.6)
                ax.set_title(title)
                ax.set_xlabel('Log2FC')
                ax.set_ylabel('-Log10(padj)')
                
            elif i == 1:  # 网络拓扑
                sizes = np.random.randint(10, 100, 50)
                ax.pie(sizes[:5], labels=['Hub', 'Bottleneck', 'Bridge', 'Peripheral', 'Other'], autopct='%1.1f%%')
                ax.set_title(title)
                
            else:  # 其他占位图
                ax.text(0.5, 0.5, f'{title}\\n(Analysis Results)', 
                       ha='center', va='center', transform=ax.transAxes, fontsize=12)
                ax.set_title(title)
        
        # 保存图形
        if output_file is None:
            output_file = f"{self.config['output']['figures_dir']}/summary_dashboard.{self.figure_format}"
        plt.savefig(output_file, dpi=self.figure_dpi, bbox_inches='tight')
        plt.close()
        
        print(f"✓ 汇总仪表板已保存: {output_file}")


def main():
    """主函数"""
    visualizer = Visualizer()
    
    # 示例可视化流程
    visualizer.plot_volcano("results/tables/significant_genes_deseq2.csv")
    visualizer.plot_heatmap("data/processed/gene_counts.txt")
    visualizer.plot_network("results/networks/integrated_network.edgelist")
    visualizer.plot_conservation_scores("results/tables/conservation_scores.csv")
    visualizer.plot_sequence_features("results/tables/sequence_features.csv")
    visualizer.create_summary_dashboard()
    
    print("可视化完成!")


if __name__ == "__main__":
    main()