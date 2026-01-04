#!/usr/bin/env python3
"""
使用ENA API获取SRA样本列表
ENA (European Nucleotide Archive) 通常比NCBI SRA更快
"""
import requests
import sys

def get_sra_from_ena(srp_id):
    """
    使用ENA Web API获取SRA样本列表

    Args:
        srp_id: SRA项目编号，如SRP286734
    """
    # ENA API endpoint
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport"
    params = {
        'accession': srp_id,
        'result': 'read_run',
        'fields': 'run_accession,experiment_title',
        'format': 'tsv'
    }

    print(f"正在从ENA查询: {srp_id}")
    print(f"URL: {url}")

    try:
        response = requests.get(url, params=params, timeout=30)
        response.raise_for_status()

        lines = response.text.strip().split('\n')
        if len(lines) < 2:
            print("未找到样本")
            return []

        # 解析TSV
        srr_list = []
        for line in lines[1:]:  # 跳过标题行
            fields = line.split('\t')
            if len(fields) > 0:
                srr_id = fields[0]
                title = fields[1] if len(fields) > 1 else "N/A"
                srr_list.append((srr_id, title))

        print(f"\n找到 {len(srr_list)} 个SRA样本:")
        for i, (srr, title) in enumerate(srr_list, 1):
            print(f"  {i}. {srr} - {title[:60]}...")

        # 只保存SRR编号
        srr_ids = [srr for srr, _ in srr_list]

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
        print("用法: python get_sra_from_ena.py <SRP_ID>")
        print("示例: python get_sra_from_ena.py SRP286734")
        sys.exit(1)

    srp_id = sys.argv[1]
    get_sra_from_ena(srp_id)
