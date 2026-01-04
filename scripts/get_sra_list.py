#!/usr/bin/env python3
"""
从SRA项目页面提取SRA编号列表
"""
import requests
import re
import sys

def get_sra_runs(srp_id):
    """
    从SRA项目获取所有run编号

    Args:
        srp_id: SRA项目编号，如SRP286734
    """
    url = f"https://www.ncbi.nlm.nih.gov/sra/{srp_id}"
    print(f"正在访问: {url}")

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # 查找所有SRR编号
        srr_pattern = r'SRR\d+'
        srr_ids = re.findall(srr_pattern, response.text)

        # 去重并排序
        srr_ids = sorted(list(set(srr_ids)))

        print(f"\n找到 {len(srr_ids)} 个SRA样本:")
        for i, srr in enumerate(srr_ids, 1):
            print(f"  {i}. {srr}")

        # 保存到文件
        output_file = f"data/raw/{srp_id}_sra_list.txt"
        with open(output_file, 'w') as f:
            for srr in srr_ids:
                f.write(f"{srr}\n")
        print(f"\n✓ SRA列表已保存到: {output_file}")

        return srr_ids

    except Exception as e:
        print(f"错误: {e}")
        return []

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python get_sra_list.py <SRP_ID>")
        print("示例: python get_sra_list.py SRP286734")
        sys.exit(1)

    srp_id = sys.argv[1]
    get_sra_runs(srp_id)
