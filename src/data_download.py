#!/usr/bin/env python3
"""
数据下载与预处理模块
用于下载公共数据集、参考基因组和注释文件
"""

import os
import subprocess
import requests
from pathlib import Path
import yaml
from typing import List, Dict
from tqdm import tqdm

class DataDownloader:
    """数据下载器"""
    
    def __init__(self, config_path: str = "config/config.yaml"):
        """初始化数据下载器"""
        with open(config_path, 'r', encoding='utf-8') as f:
            self.config = yaml.safe_load(f)
        
        # 创建必要的目录
        self._create_directories()
    
    def _create_directories(self):
        """创建数据目录结构"""
        dirs = [
            "data/raw",
            "data/processed", 
            "data/reference",
            "data/annotation",
            "data/line1",
            "results/figures",
            "results/tables",
            "results/reports",
            "results/networks"
        ]
        
        for dir_path in dirs:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
            print(f"✓ 创建目录: {dir_path}")
    
    def download_sra_data(self, sra_ids: List[str], output_dir: str = "data/raw/"):
        """
        下载SRA数据
        
        Args:
            sra_ids: SRA登录号列表
            output_dir: 输出目录
        """
        print("开始下载SRA数据...")
        
        for sra_id in tqdm(sra_ids, desc="下载进度"):
            try:
                # 使用prefetch下载
                cmd = f"prefetch {sra_id} --output-directory {output_dir}"
                subprocess.run(cmd, shell=True, check=True)
                
                # 转换为fastq格式
                fastq_dir = os.path.join(output_dir, "fastq")
                Path(fastq_dir).mkdir(exist_ok=True)
                
                cmd = f"fastq-dump {output_dir}/{sra_id}/{sra_id}.sra --gzip --outdir {fastq_dir}"
                subprocess.run(cmd, shell=True, check=True)
                
                print(f"✓ 成功下载并转换: {sra_id}")
                
            except subprocess.CalledProcessError as e:
                print(f"✗ 下载失败 {sra_id}: {e}")
    
    def download_reference_genome(self, genome: str = "hg38"):
        """
        下载参考基因组
        
        Args:
            genome: 基因组版本
        """
        print(f"下载参考基因组: {genome}")
        
        if genome == "hg38":
            fasta_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
            fasta_path = "data/reference/hg38.fa.gz"
            
            # 下载基因组FASTA文件
            self._download_file(fasta_url, fasta_path)
            
            # 解压缩
            cmd = "gunzip -c data/reference/hg38.fa.gz > data/reference/hg38.fa"
            subprocess.run(cmd, shell=True, check=True)
            
            # 创建BWA索引
            cmd = "bwa index data/reference/hg38.fa"
            subprocess.run(cmd, shell=True, check=True)
            
            print("✓ 参考基因组下载并索引完成")
    
    def download_gencode_annotation(self, version: str = "v45"):
        """
        下载GENCODE注释文件
        
        Args:
            version: GENCODE版本
        """
        print(f"下载GENCODE注释文件: {version}")
        
        gtf_url = f"https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"
        gtf_path = f"data/annotation/gencode.{version}.annotation.gtf.gz"
        
        # 下载GTF文件
        self._download_file(gtf_url, gtf_path)
        
        # 解压缩
        cmd = f"gunzip -c {gtf_path} > data/annotation/gencode.{version}.annotation.gtf"
        subprocess.run(cmd, shell=True, check=True)
        
        print("✓ GENCODE注释文件下载完成")
    
    def download_line1_sequences(self):
        """下载LINE-1序列信息"""
        print("下载LINE-1序列信息...")
        
        # 从UCSC下载RepeatMasker注释
        repeat_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz"
        repeat_path = "data/line1/rmsk.txt.gz"
        
        self._download_file(repeat_url, repeat_path)
        
        # 解压缩并过滤LINE-1
        cmd = "gunzip -c data/line1/rmsk.txt.gz | grep 'LINE1' > data/line1/LINE1_raw.txt"
        subprocess.run(cmd, shell=True, check=True)
        
        # 创建BED文件
        self._create_line1_bed()
        
        # 提取LINE-1序列
        self._extract_line1_sequences()
        
        print("✓ LINE-1序列信息下载完成")
    
    def _create_line1_bed(self):
        """创建LINE-1 BED文件"""
        cmd = """awk 'BEGIN{OFS="\\t"} {
            if($12 == "LINE1") {
                print $6, $7, $8, $11, $10, $3
            }
        }' data/line1/LINE1_raw.txt > data/line1/LINE1.bed"""
        
        subprocess.run(cmd, shell=True, check=True)
    
    def _extract_line1_sequences(self):
        """提取LINE-1序列"""
        cmd = """bedtools getfasta -fi data/reference/hg38.fa -bed data/line1/LINE1.bed -fo data/line1/LINE1_sequences.fa"""
        subprocess.run(cmd, shell=True, check=True)
    
    def _download_file(self, url: str, output_path: str):
        """下载文件"""
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        
        with open(output_path, 'wb') as f, tqdm(
            desc=os.path.basename(output_path),
            total=total_size,
            unit='B',
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    bar.update(len(chunk))
    
    def quality_control(self, fastq_file: str, output_dir: str = "data/processed/"):
        """
        质量控制
        
        Args:
            fastq_file: FASTQ文件路径
            output_dir: 输出目录
        """
        print(f"质控处理: {fastq_file}")
        
        filename = os.path.basename(fastq_file).replace('.fastq.gz', '')
        output_file = os.path.join(output_dir, f"{filename}_clean.fastq.gz")
        
        cmd = f"fastp -i {fastq_file} -o {output_file} --detect_adapter_for_pe --qualified_quality_phred 20 --length_required 30 --thread 4"
        subprocess.run(cmd, shell=True, check=True)
        
        print(f"✓ 质控完成: {output_file}")
    
    def align_reads(self, fastq_file: str, output_bam: str):
        """
        序列比对
        
        Args:
            fastq_file: FASTQ文件路径
            output_bam: 输出BAM文件路径
        """
        print(f"序列比对: {fastq_file}")
        
        # BWA比对
        sam_file = output_bam.replace('.bam', '.sam')
        cmd = f"bwa mem -t 8 data/reference/hg38.fa {fastq_file} > {sam_file}"
        subprocess.run(cmd, shell=True, check=True)
        
        # 转换为BAM并排序
        cmd = f"samtools view -bS {sam_file} | samtools sort -o {output_bam}"
        subprocess.run(cmd, shell=True, check=True)
        
        # 建立索引
        cmd = f"samtools index {output_bam}"
        subprocess.run(cmd, shell=True, check=True)
        
        # 清理临时文件
        os.remove(sam_file)
        
        print(f"✓ 比对完成: {output_bam}")


def main():
    """主函数"""
    downloader = DataDownloader()
    
    # 下载参考基因组和注释文件
    downloader.download_reference_genome()
    downloader.download_gencode_annotation()
    downloader.download_line1_sequences()
    
    # 下载SRA数据 (示例)
    sra_ids = ["SRR1234567", "SRR1234568", "SRR1234569"]
    downloader.download_sra_data(sra_ids)
    
    print("数据下载完成!")


if __name__ == "__main__":
    main()