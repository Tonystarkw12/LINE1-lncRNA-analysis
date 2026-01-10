#!/usr/bin/env python3
"""
使用Python下载GEO数据并解压
"""
import requests
import tarfile
import os
import sys

def download_and_extract(gse_id="GSE159217"):
    """下载并解压GEO数据"""

    # 构建下载URL
    gse_num = gse_id[3:]  # 提取数字部分
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{gse_id[:-3]}nnn/{gse_id}/suppl/{gse_id}_RAW.tar"

    output_dir = "/home/tony/line1/LINE1_lncRNA_project/data/processed"
    tar_path = os.path.join(output_dir, f"{gse_id}_RAW.tar")

    print(f"正在下载: {url}")
    print(f"保存到: {tar_path}")

    try:
        # 流式下载
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))
        print(f"文件大小: {total_size / 1024 / 1024:.2f} MB")

        downloaded = 0
        with open(tar_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        progress = (downloaded / total_size) * 100
                        print(f"\r下载进度: {progress:.1f}%", end='')

        print(f"\n✓ 下载完成: {tar_path}")

        # 解压
        print("\n正在解压...")
        with tarfile.open(tar_path, 'r') as tar:
            tar.extractall(path=output_dir)

        print(f"✓ 解压完成")

        # 列出文件
        print("\n解压后的文件:")
        files = os.listdir(output_dir)
        gz_files = [f for f in files if f.endswith('.gz')]
        print(f"找到 {len(gz_files)} 个gzip文件")

        return gz_files

    except Exception as e:
        print(f"\n✗ 错误: {e}")
        import traceback
        traceback.print_exc()
        return []

if __name__ == "__main__":
    download_and_extract()
