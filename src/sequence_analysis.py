#!/usr/bin/env python3
"""
序列分析与保守性分析模块
包括保守性分析、二级结构预测、序列特征提取
"""

import os
import subprocess
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple
import requests

class SequenceAnalyzer:
    """序列分析器"""
    
    def __init__(self, config_path: str = "config/config.yaml"):
        """初始化序列分析器"""
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
    
    def extract_lncrna_sequences(self, lncrna_file: str, gtf_file: str, 
                               genome_fasta: str, output_file: str):
        """
        提取lncRNA序列
        
        Args:
            lncrna_file: lncRNA列表文件
            gtf_file: GTF注释文件
            genome_fasta: 参考基因组FASTA文件
            output_file: 输出FASTA文件
        """
        print("提取lncRNA序列...")
        
        # 读取lncRNA列表
        lncrnas = set(pd.read_csv(lncrna_file, index_col=0).index)
        
        # 创建lncRNA坐标字典
        lncrna_coords = {}
        with open(gtf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\\t')
                if len(fields) >= 9 and fields[2] == 'gene':
                    attributes = fields[8]
                    if 'gene_biotype "lncRNA"' in attributes:
                        gene_id = None
                        gene_name = None
                        for attr in attributes.split(';'):
                            if 'gene_id' in attr:
                                gene_id = attr.split('"')[1]
                            elif 'gene_name' in attr:
                                gene_name = attr.split('"')[1]
                        
                        if gene_id and gene_id in lncrnas:
                            lncrna_coords[gene_id] = {
                                'chr': fields[0],
                                'start': int(fields[3]),
                                'end': int(fields[4]),
                                'strand': fields[6],
                                'name': gene_name if gene_name else gene_id
                            }
        
        # 提取序列
        with open(output_file, 'w') as out_f:
            for gene_id, coords in lncrna_coords.items():
                # 使用bedtools提取序列
                cmd = f"""bedtools getfasta -fi {genome_fasta} \\
                    -chr {coords['chr']} -start {coords['start']} -end {coords['end']} \\
                    -name -s > temp_seq.fa"""
                subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL)
                
                # 读取临时序列文件
                with open("temp_seq.fa", 'r') as temp_f:
                    seq_record = next(SeqIO.parse(temp_f, "fasta"))
                    seq_record.id = f"{gene_id}|{coords['name']}"
                    seq_record.description = f"{coords['chr']}:{coords['start']}-{coords['end']}({coords['strand']})"
                    SeqIO.write(seq_record, out_f, "fasta")
                
                # 清理临时文件
                os.remove("temp_seq.fa")
        
        print(f"✓ 序列提取完成，共提取 {len(lncrna_coords)} 个lncRNA序列")
    
    def calculate_conservation_scores(self, lncrna_fasta: str, output_file: str):
        """
        计算保守性得分
        
        Args:
            lncrna_fasta: lncRNA FASTA文件
            output_file: 输出文件
        """
        print("计算保守性得分...")
        
        # 下载PhyloP数据 (示例：使用UCSC的bigWig文件)
        phylop_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw"
        phylop_file = "data/annotation/hg38.phyloP100way.bw"
        
        if not os.path.exists(phylop_file):
            print("下载PhyloP保守性数据...")
            self._download_file(phylop_url, phylop_file)
        
        # 使用pyBigWig计算保守性得分
        conservation_scores = []
        
        try:
            import pyBigWig
            
            bw = pyBigWig.open(phylop_file)
            
            for record in SeqIO.parse(lncrna_fasta, "fasta"):
                # 解析染色体和位置信息
                header_parts = record.description.split(':')
                if len(header_parts) >= 2:
                    chrom = header_parts[0].split()[0]  # 去掉染色体前缀
                    positions = header_parts[1].split('(')[0].split('-')
                    start = int(positions[0]) - 1  # 转换为0-based
                    end = int(positions[1])
                    
                    # 获取保守性得分
                    scores = bw.values(chrom, start, end)
                    scores = [s if s is not None else 0 for s in scores]
                    
                    # 计算统计量
                    avg_score = np.mean(scores)
                    max_score = np.max(scores)
                    min_score = np.min(scores)
                    std_score = np.std(scores)
                    
                    # 计算高保守区域比例
                    high_cons_threshold = self.config['analysis']['conservation']['phylop_threshold']
                    high_cons_ratio = sum(1 for s in scores if s > high_cons_threshold) / len(scores)
                    
                    conservation_scores.append({
                        'gene_id': record.id.split('|')[0],
                        'gene_name': record.id.split('|')[1],
                        'avg_phylop': avg_score,
                        'max_phylop': max_score,
                        'min_phylop': min_score,
                        'std_phylop': std_score,
                        'high_cons_ratio': high_cons_ratio,
                        'length': len(record.seq)
                    })
            
            bw.close()
            
        except ImportError:
            print("警告: pyBigWig未安装，跳过保守性分析")
            return None
        
        # 保存结果
        conservation_df = pd.DataFrame(conservation_scores)
        conservation_df.to_csv(output_file, index=False)
        
        print(f"✓ 保守性分析完成")
        
        return conservation_df
    
    def predict_secondary_structure(self, lncrna_fasta: str, output_dir: str = "results/structures/"):
        """
        预测RNA二级结构
        
        Args:
            lncrna_fasta: lncRNA FASTA文件
            output_dir: 输出目录
        """
        print("预测RNA二级结构...")
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        structure_results = []
        
        for record in SeqIO.parse(lncrna_fasta, "fasta"):
            gene_id = record.id.split('|')[0]
            
            # 使用RNAfold预测二级结构
            temp_file = f"temp_{gene_id}.fa"
            with open(temp_file, 'w') as f:
                SeqIO.write(record, f, "fasta")
            
            try:
                # 运行RNAfold
                cmd = f"RNAfold < {temp_file}"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                
                if result.returncode == 0:
                    output_lines = result.stdout.strip().split('\\n')
                    if len(output_lines) >= 2:
                        structure = output_lines[1].split(' ')[0]  # 二级结构
                        energy = float(output_lines[1].split('(')[1].split(')')[0])  # 自由能
                        
                        # 计算结构特征
                        gc_content = GC(record.seq)
                        length = len(record.seq)
                        
                        # 计算碱基配对比例
                        paired_bases = structure.count('(') + structure.count(')')
                        pairing_ratio = paired_bases / length if length > 0 else 0
                        
                        # 计算茎环结构数量
                        stem_loops = self._count_stem_loops(structure)
                        
                        structure_results.append({
                            'gene_id': gene_id,
                            'gene_name': record.id.split('|')[1],
                            'length': length,
                            'gc_content': gc_content,
                            'structure': structure,
                            'free_energy': energy,
                            'pairing_ratio': pairing_ratio,
                            'stem_loops': stem_loops
                        })
                        
                        # 保存结构文件
                        structure_file = f"{output_dir}/{gene_id}_structure.txt"
                        with open(structure_file, 'w') as f:
                            f.write(f">{record.description}\\n")
                            f.write(f"{str(record.seq)}\\n")
                            f.write(f"{structure}\\n")
                            f.write(f"Free Energy: {energy} kcal/mol\\n")
                
            except subprocess.CalledProcessError as e:
                print(f"RNAfold预测失败 {gene_id}: {e}")
            
            finally:
                # 清理临时文件
                if os.path.exists(temp_file):
                    os.remove(temp_file)
        
        # 保存结果汇总
        structure_df = pd.DataFrame(structure_results)
        structure_df.to_csv(f"{output_dir}/secondary_structures.csv", index=False)
        
        print(f"✓ 二级结构预测完成，共预测 {len(structure_results)} 个lncRNA")
        
        return structure_df
    
    def _count_stem_loops(self, structure: str) -> int:
        """计算茎环结构数量"""
        stem_loops = 0
        stack = []
        
        for char in structure:
            if char == '(':
                stack.append(char)
            elif char == ')':
                if stack:
                    stack.pop()
                    if len(stack) == 0:  # 完成一个茎环
                        stem_loops += 1
        
        return stem_loops
    
    def analyze_sequence_features(self, lncrna_fasta: str, output_file: str):
        """
        分析序列特征
        
        Args:
            lncrna_fasta: lncRNA FASTA文件
            output_file: 输出文件
        """
        print("分析序列特征...")
        
        sequence_features = []
        
        for record in SeqIO.parse(lncrna_fasta, "fasta"):
            seq = str(record.seq)
            gene_id = record.id.split('|')[0]
            
            # 基本特征
            length = len(seq)
            gc_content = GC(seq)
            at_content = 100 - gc_content
            
            # k-mer分析
            kmers = self._analyze_kmers(seq, k=4)
            
            # 重复序列分析
            repeats = self._find_repeats(seq)
            
            # 基序搜索 (示例：poly-A信号)
            poly_a_signals = self._find_motifs(seq, "AATAAA")
            
            # CpG岛分析
            cpg_islands = self._find_cpg_islands(seq)
            
            sequence_features.append({
                'gene_id': gene_id,
                'gene_name': record.id.split('|')[1],
                'length': length,
                'gc_content': gc_content,
                'at_content': at_content,
                'top_kmer': kmers[0] if kmers else None,
                'repeat_count': repeats['count'],
                'repeat_ratio': repeats['ratio'],
                'poly_a_signals': poly_a_signals,
                'cpg_islands': cpg_islands['count'],
                'cpg_length': cpg_islands['total_length']
            })
        
        # 保存结果
        features_df = pd.DataFrame(sequence_features)
        features_df.to_csv(output_file, index=False)
        
        print(f"✓ 序列特征分析完成")
        
        return features_df
    
    def _analyze_kmers(self, sequence: str, k: int = 4) -> List[Tuple[str, int]]:
        """分析k-mer频率"""
        kmers = {}
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            kmers[kmer] = kmers.get(kmer, 0) + 1
        
        # 返回频率最高的k-mers
        sorted_kmers = sorted(kmers.items(), key=lambda x: x[1], reverse=True)
        return sorted_kmers[:10]
    
    def _find_repeats(self, sequence: str, min_length: int = 10) -> Dict:
        """查找重复序列"""
        repeats = {}
        
        for length in range(min_length, min(50, len(sequence)//2)):
            for i in range(len(sequence) - length + 1):
                repeat = sequence[i:i+length]
                if sequence.count(repeat) > 1:
                    repeats[repeat] = sequence.count(repeat)
        
        total_repeats = sum(repeats.values())
        repeat_ratio = total_repeats / len(sequence) if len(sequence) > 0 else 0
        
        return {
            'count': len(repeats),
            'ratio': repeat_ratio,
            'repeats': repeats
        }
    
    def _find_motifs(self, sequence: str, motif: str) -> int:
        """查找特定基序"""
        count = 0
        start = 0
        while True:
            pos = sequence.find(motif, start)
            if pos == -1:
                break
            count += 1
            start = pos + 1
        return count
    
    def _find_cpg_islands(self, sequence: str, min_length: int = 200, 
                         gc_content_threshold: float = 0.5, 
                         cpg_ratio_threshold: float = 0.6) -> Dict:
        """查找CpG岛"""
        cpg_islands = []
        i = 0
        
        while i < len(sequence):
            # 寻找CpG富集区域
            if sequence[i:i+2] == "CG":
                start = i
                gc_count = 0
                cpg_count = 0
                length = 0
                
                # 扩展区域
                while i < len(sequence) and length < min_length:
                    if sequence[i] in ['G', 'C']:
                        gc_count += 1
                    if i < len(sequence) - 1 and sequence[i:i+2] == "CG":
                        cpg_count += 1
                    length += 1
                    i += 1
                
                # 检查是否满足CpG岛标准
                if length >= min_length:
                    gc_content = gc_count / length
                    observed_cpg = cpg_count
                    expected_cpg = (gc_count/2) * (gc_count/2) / length
                    cpg_ratio = observed_cpg / expected_cpg if expected_cpg > 0 else 0
                    
                    if gc_content >= gc_content_threshold and cpg_ratio >= cpg_ratio_threshold:
                        cpg_islands.append({
                            'start': start,
                            'end': i,
                            'length': length,
                            'gc_content': gc_content,
                            'cpg_ratio': cpg_ratio
                        })
            
            i += 1
        
        total_length = sum(island['length'] for island in cpg_islands)
        
        return {
            'count': len(cpg_islands),
            'total_length': total_length,
            'islands': cpg_islands
        }
    
    def _download_file(self, url: str, output_path: str):
        """下载文件"""
        response = requests.get(url, stream=True)
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)


def main():
    """主函数"""
    analyzer = SequenceAnalyzer()
    
    # 示例流程
    lncrna_file = "results/tables/lncrna_differential.csv"
    gtf_file = "data/annotation/gencode.v45.annotation.gtf"
    genome_fasta = "data/reference/hg38.fa"
    
    # 提取lncRNA序列
    analyzer.extract_lncrna_sequences(lncrna_file, gtf_file, genome_fasta, 
                                     "data/processed/lncrnas.fa")
    
    # 保守性分析
    analyzer.calculate_conservation_scores("data/processed/lncrnas.fa", 
                                         "results/tables/conservation_scores.csv")
    
    # 二级结构预测
    analyzer.predict_secondary_structure("data/processed/lncrnas.fa")
    
    # 序列特征分析
    analyzer.analyze_sequence_features("data/processed/lncrnas.fa", 
                                      "results/tables/sequence_features.csv")
    
    print("序列分析完成!")


if __name__ == "__main__":
    main()