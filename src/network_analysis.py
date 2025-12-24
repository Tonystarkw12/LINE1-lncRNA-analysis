#!/usr/bin/env python3
"""
调控网络构建模块
构建lncRNA-mRNA-LINE-1调控网络，进行网络分析
"""

import os
import subprocess
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple, Set
import requests
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

class NetworkAnalyzer:
    """网络分析器"""
    
    def __init__(self, config_path: str = "config/config.yaml"):
        """初始化网络分析器"""
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
    
    def build_ppi_network(self, gene_list: List[str], output_file: str):
        """
        构建蛋白质-蛋白质相互作用网络
        
        Args:
            gene_list: 基因列表
            output_file: 输出文件
        """
        print("构建PPI网络...")
        
        # 使用STRING数据库获取PPI信息
        ppi_edges = self._get_string_interactions(gene_list)
        
        # 创建网络图
        G = nx.Graph()
        G.add_edges_from(ppi_edges)
        
        # 计算网络拓扑特征
        self._calculate_network_metrics(G)
        
        # 保存网络
        nx.write_edgelist(G, output_file)
        
        print(f"✓ PPI网络构建完成，节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}")
        
        return G
    
    def _get_string_interactions(self, gene_list: List[str], score_threshold: float = 0.7) -> List[Tuple[str, str]]:
        """
        从STRING数据库获取蛋白质相互作用
        
        Args:
            gene_list: 基因列表
            score_threshold: 相互作用得分阈值
            
        Returns:
            相互作用边列表
        """
        ppi_edges = []
        
        # STRING API (示例，需要实际API密钥)
        string_url = "https://string-db.org/api/json/network"
        
        # 分批查询以避免API限制
        batch_size = 100
        for i in range(0, len(gene_list), batch_size):
            batch_genes = gene_list[i:i+batch_size]
            genes_str = "%0d".join(ord(c) for c in ",".join(batch_genes))
            
            params = {
                "identifiers": ",".join(batch_genes),
                "species": 9606,  # 人类
                "required_score": int(score_threshold * 1000),
                "limit": 1000
            }
            
            try:
                response = requests.get(string_url, params=params)
                if response.status_code == 200:
                    interactions = response.json()
                    for interaction in interactions:
                        gene1 = interaction['preferredName_A']
                        gene2 = interaction['preferredName_B']
                        score = interaction['score']
                        
                        if score >= score_threshold * 1000:
                            ppi_edges.append((gene1, gene2))
                
            except Exception as e:
                print(f"STRING API查询失败: {e}")
                # 使用模拟数据
                ppi_edges = self._generate_mock_ppi(gene_list)
                break
        
        return ppi_edges
    
    def _generate_mock_ppi(self, gene_list: List[str]) -> List[Tuple[str, str]]:
        """生成模拟PPI数据用于测试"""
        import random
        
        ppi_edges = []
        # 随机生成一些相互作用
        for i in range(min(100, len(gene_list) * 2)):
            gene1 = random.choice(gene_list)
            gene2 = random.choice(gene_list)
            if gene1 != gene2:
                ppi_edges.append((gene1, gene2))
        
        return list(set(ppi_edges))  # 去重
    
    def build_lncrna_mrna_network(self, lncrna_file: str, mrna_file: str, 
                                expression_data: str, output_file: str):
        """
        构建lncRNA-mRNA共表达网络
        
        Args:
            lncrna_file: lncRNA表达文件
            mrna_file: mRNA表达文件
            expression_data: 表达量数据
            output_file: 输出文件
        """
        print("构建lncRNA-mRNA共表达网络...")
        
        # 读取表达数据
        lncrna_expr = pd.read_csv(lncrna_file, index_col=0)
        mrna_expr = pd.read_csv(mrna_file, index_col=0)
        
        # 计算相关性
        correlations = []
        for lncrna in lncrna_expr.index:
            for mrna in mrna_expr.index:
                corr, pvalue = stats.pearsonr(lncrna_expr.loc[lncrna], mrna_expr.loc[mrna])
                if abs(corr) > 0.7 and pvalue < 0.05:  # 相关性阈值
                    correlations.append({
                        'lncRNA': lncrna,
                        'mRNA': mrna,
                        'correlation': corr,
                        'pvalue': pvalue
                    })
        
        # 创建网络
        G = nx.Graph()
        
        # 添加节点和边
        for corr in correlations:
            lncrna = corr['lncRNA']
            mrna = corr['mRNA']
            weight = abs(corr['correlation'])
            
            G.add_node(lncrna, type='lncRNA')
            G.add_node(mrna, type='mRNA')
            G.add_edge(lncrna, mrna, weight=weight, correlation=corr['correlation'])
        
        # 保存网络
        nx.write_edgelist(G, output_file)
        
        # 保存相关性数据
        corr_df = pd.DataFrame(correlations)
        corr_df.to_csv(output_file.replace('.edgelist', '_correlations.csv'), index=False)
        
        print(f"✓ lncRNA-mRNA网络构建完成，节点数: {G.number_of_nodes()}, 边数: {G.number_of_edges()}")
        
        return G, corr_df
    
    def integrate_line1_network(self, lncrna_mrna_network: nx.Graph, 
                              line1_associations: pd.DataFrame, output_file: str):
        """
        整合LINE-1关联信息到调控网络
        
        Args:
            lncrna_mrna_network: lncRNA-mRNA网络
            line1_associations: LINE-1关联数据
            output_file: 输出文件
        """
        print("整合LINE-1关联...")
        
        # 复制网络
        G = lncrna_mrna_network.copy()
        
        # 添加LINE-1节点和边
        line1_nodes = set()
        for _, row in line1_associations.iterrows():
            lncrna = row['lncRNA_id']
            line1_id = f"LINE1_{row['LINE1_id']}"
            
            # 添加LINE-1节点
            G.add_node(line1_id, type='LINE1', 
                      association_type=row['association_type'],
                      distance=row['distance'])
            
            # 添加lncRNA-LINE-1边
            weight = 1.0 / (1.0 + row['distance']/1000)  # 距离越近权重越高
            G.add_edge(lncrna, line1_id, weight=weight, 
                      association_type=row['association_type'],
                      distance=row['distance'])
            
            line1_nodes.add(line1_id)
        
        # 保存整合网络
        nx.write_edgelist(G, output_file)
        
        print(f"✓ LINE-1整合完成，新增LINE-1节点数: {len(line1_nodes)}")
        
        return G
    
    def _calculate_network_metrics(self, G: nx.Graph):
        """计算网络拓扑特征"""
        # 度中心性
        degree_centrality = nx.degree_centrality(G)
        nx.set_node_attributes(G, degree_centrality, 'degree_centrality')
        
        # 介数中心性
        betweenness_centrality = nx.betweenness_centrality(G)
        nx.set_node_attributes(G, betweenness_centrality, 'betweenness_centrality')
        
        # 紧密中心性
        closeness_centrality = nx.closeness_centrality(G)
        nx.set_node_attributes(G, closeness_centrality, 'closeness_centrality')
        
        # 特征向量中心性
        try:
            eigenvector_centrality = nx.eigenvector_centrality(G)
            nx.set_node_attributes(G, eigenvector_centrality, 'eigenvector_centrality')
        except:
            print("无法计算特征向量中心性")
    
    def identify_key_regulators(self, network: nx.Graph, output_file: str):
        """
        识别关键调控因子
        
        Args:
            network: 网络图
            output_file: 输出文件
        """
        print("识别关键调控因子...")
        
        # 提取网络指标
        node_metrics = []
        for node, data in network.nodes(data=True):
            metrics = {
                'node': node,
                'type': data.get('type', 'unknown'),
                'degree': network.degree(node),
                'degree_centrality': data.get('degree_centrality', 0),
                'betweenness_centrality': data.get('betweenness_centrality', 0),
                'closeness_centrality': data.get('closeness_centrality', 0),
                'eigenvector_centrality': data.get('eigenvector_centrality', 0)
            }
            node_metrics.append(metrics)
        
        # 转换为DataFrame
        metrics_df = pd.DataFrame(node_metrics)
        
        # 计算综合得分
        scaler = StandardScaler()
        features = ['degree_centrality', 'betweenness_centrality', 
                   'closeness_centrality', 'eigenvector_centrality']
        
        # 处理缺失值
        for feature in features:
            metrics_df[feature] = metrics_df[feature].fillna(0)
        
        # 标准化并计算综合得分
        scaled_features = scaler.fit_transform(metrics_df[features])
        metrics_df['hub_score'] = np.mean(scaled_features, axis=1)
        
        # 排序
        metrics_df = metrics_df.sort_values('hub_score', ascending=False)
        
        # 分类关键节点
        key_regulators = {
            'hubs': metrics_df[metrics_df['degree_centrality'] > 0.8]['node'].tolist(),
            'bottlenecks': metrics_df[metrics_df['betweenness_centrality'] > 0.8]['node'].tolist(),
            'top_hubs': metrics_df.head(20)['node'].tolist()
        }
        
        # 保存结果
        metrics_df.to_csv(output_file, index=False)
        
        with open(output_file.replace('.csv', '_key_regulators.txt'), 'w') as f:
            f.write("关键调控因子分析结果\\n\\n")
            f.write(f"Top 20 关键节点:\\n")
            for node in key_regulators['top_hubs']:
                f.write(f"- {node}\\n")
            
            f.write(f"\\nHub节点 (度中心性 > 0.8):\\n")
            for node in key_regulators['hubs']:
                f.write(f"- {node}\\n")
            
            f.write(f"\\nBottleneck节点 (介数中心性 > 0.8):\\n")
            for node in key_regulators['bottlenecks']:
                f.write(f"- {node}\\n")
        
        print(f"✓ 关键调控因子识别完成")
        
        return metrics_df, key_regulators
    
    def perform_network_clustering(self, network: nx.Graph, output_file: str, 
                                 n_clusters: int = 5):
        """
        网络聚类分析
        
        Args:
            network: 网络图
            output_file: 输出文件
            n_clusters: 聚类数量
        """
        print("进行网络聚类分析...")
        
        # 使用社区发现算法
        try:
            import community as community_louvain
            communities = community_louvain.best_partition(network)
        except ImportError:
            print("警告: python-louvain未安装，使用简单的聚类方法")
            communities = self._simple_clustering(network, n_clusters)
        
        # 添加聚类信息到网络
        nx.set_node_attributes(network, communities, 'community')
        
        # 分析聚类特征
        cluster_stats = []
        for cluster_id in set(communities.values()):
            nodes_in_cluster = [node for node, comm in communities.items() 
                              if comm == cluster_id]
            
            # 计算聚类内连接密度
            subgraph = network.subgraph(nodes_in_cluster)
            density = nx.density(subgraph)
            
            # 节点类型统计
            node_types = {}
            for node in nodes_in_cluster:
                node_type = network.nodes[node].get('type', 'unknown')
                node_types[node_type] = node_types.get(node_type, 0) + 1
            
            cluster_stats.append({
                'cluster_id': cluster_id,
                'size': len(nodes_in_cluster),
                'density': density,
                'lncRNA_count': node_types.get('lncRNA', 0),
                'mRNA_count': node_types.get('mRNA', 0),
                'LINE1_count': node_types.get('LINE1', 0),
                'nodes': nodes_in_cluster
            })
        
        # 保存聚类结果
        cluster_df = pd.DataFrame(cluster_stats)
        cluster_df.to_csv(output_file, index=False)
        
        # 保存节点聚类标签
        node_clusters = pd.DataFrame([
            {'node': node, 'cluster': communities[node]} 
            for node in network.nodes()
        ])
        node_clusters.to_csv(output_file.replace('.csv', '_nodes.csv'), index=False)
        
        print(f"✓ 网络聚类完成，共识别 {len(cluster_stats)} 个聚类")
        
        return cluster_df, communities
    
    def _simple_clustering(self, network: nx.Graph, n_clusters: int) -> Dict:
        """简单的聚类方法"""
        # 使用节点度数进行K-means聚类
        degrees = np.array([network.degree(node) for node in network.nodes()]).reshape(-1, 1)
        
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(degrees)
        
        communities = {}
        for i, node in enumerate(network.nodes()):
            communities[node] = int(cluster_labels[i])
        
        return communities
    
    def analyze_pathways(self, gene_list: List[str], output_file: str):
        """
        通路富集分析
        
        Args:
            gene_list: 基因列表
            output_file: 输出文件
        """
        print("进行通路富集分析...")
        
        # 创建R脚本进行通路分析
        r_script = f'''
# 通路富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# 基因列表
genes <- c({", ".join([f'"{gene}"' for gene in gene_list[:100]])})

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
write.csv(as.data.frame(go_results), file="{output_file.replace('.csv', '_go.csv')}")
write.csv(as.data.frame(kegg_results), file="{output_file.replace('.csv', '_kegg.csv')}")

# 绘图
png("{output_file.replace('.csv', '_go.png')}", width=1200, height=800, res=300)
dotplot(go_results, showCategory=20) + ggtitle("GO Biological Process Enrichment")
dev.off()

png("{output_file.replace('.csv', '_kegg.png')}", width=1200, height=800, res=300)
dotplot(kegg_results, showCategory=20) + ggtitle("KEGG Pathway Enrichment")
dev.off()

cat("通路分析完成\\n")
'''
        
        script_path = "results/networks/pathway_analysis.R"
        with open(script_path, 'w', encoding='utf-8') as f:
            f.write(r_script)
        
        # 运行R脚本
        try:
            subprocess.run(f"Rscript {script_path}", shell=True, check=True)
            print("✓ 通路富集分析完成")
        except subprocess.CalledProcessError as e:
            print(f"通路分析失败: {e}")


def main():
    """主函数"""
    analyzer = NetworkAnalyzer()
    
    # 示例基因列表
    gene_list = ["TP53", "MYC", "EGFR", "KRAS", "BRCA1"]
    
    # 构建PPI网络
    ppi_network = analyzer.build_ppi_network(gene_list, "results/networks/ppi_network.edgelist")
    
    # 识别关键调控因子
    metrics_df, key_regulators = analyzer.identify_key_regulators(ppi_network, 
                                                                "results/networks/network_metrics.csv")
    
    # 网络聚类
    cluster_df, communities = analyzer.perform_network_clustering(ppi_network, 
                                                                "results/networks/clusters.csv")
    
    # 通路分析
    analyzer.analyze_pathways(gene_list, "results/networks/pathways.csv")
    
    print("网络分析完成!")


if __name__ == "__main__":
    main()