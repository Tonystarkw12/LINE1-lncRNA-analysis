#!/usr/bin/env python3
"""
从GEO下载已处理的表达数据
"""
import requests
import os
import pandas as pd
import sys

def download_geodata(gse_id):
    """
    下载GEO已处理的数据

    Args:
        gse_id: GEO编号，如GSE159217
    """
    print(f"正在下载 {gse_id} 的处理数据...")

    # GEO系列矩阵文件URL
    matrix_url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&targ=gsm&form=text&view=data=full"

    try:
        print(f"访问: {matrix_url}")
        response = requests.get(matrix_url, timeout=60)
        response.raise_for_status()

        # 保存原始数据
        output_dir = "data/processed"
        os.makedirs(output_dir, exist_ok=True)

        raw_file = f"{output_dir}/{gse_id}_raw.txt"
        with open(raw_file, 'w', encoding='utf-8') as f:
            f.write(response.text)

        print(f"✓ 原始数据已下载: {raw_file}")
        print(f"  文件大小: {len(response.text):,} 字符")

        # 尝试解析为表格
        try:
            # GEO数据通常是制表符分隔的
            from io import StringIO
            df = pd.read_csv(StringIO(response.text), sep='\t', index_col=0)

            # 保存为CSV
            csv_file = f"{output_dir}/{gse_id}_expression.csv"
            df.to_csv(csv_file)

            print(f"✓ 表格数据已保存: {csv_file}")
            print(f"  基因数: {len(df)}")
            print(f"  样本数: {len(df.columns)}")
            print(f"\n前5个样本: {list(df.columns[:5])}")

            return df

        except Exception as e:
            print(f"⚠️ 解析表格失败: {e}")
            print("原始数据已保存，您可以手动处理")
            return None

    except Exception as e:
        print(f"❌ 下载失败: {e}")
        return None

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python download_geodata.py <GSE_ID>")
        print("示例: python download_geodata.py GSE159217")
        sys.exit(1)

    gse_id = sys.argv[1]
    download_geodata(gse_id)
