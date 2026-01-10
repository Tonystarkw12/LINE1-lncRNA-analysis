#!/usr/bin/env python3
"""
从GEO FTP下载supplementary files
"""
import requests
import os
import re
import sys

def get_supplementary_files(gse_id):
    """
    获取并下载GEO supplementary files

    Args:
        gse_id: GEO编号，如GSE159217
    """
    print(f"正在查找 {gse_id} 的supplementary files...")

    # GEO查询页面
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}"

    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()

        # 查找supplementary file链接
        # GEO的supplementary files通常在ftp://ftp.ncbi.nlm.nih.gov/geo/series/
        pattern = r'ftp://ftp\.ncbi\.nlm\.nih\.gov/geo/series/[^\s"]+\.tar\.gz'
        files = re.findall(pattern, response.text)

        if not files:
            print("未找到supplementary files")
            # 尝试查找其他格式的文件
            pattern2 = r'(GSE\d+_RAW\.tar|GSE\d+_RAW\.tar\.gz)'
            files2 = re.findall(pattern2, response.text)
            if files2:
                print(f"找到: {files2}")
                # 构建完整FTP路径
                gse_num = gse_id[3:]  # 提取数字部分
                file_path = f"GSE{gse_num}_RAW/GSE{gse_num}_RAW.tar"
                ftp_url = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE{gse_num[:-3]}nn/{gse_id}/raw/{file_path}"
                print(f"FTP URL: {ftp_url}")
                files = [ftp_url]

        if not files:
            return []

        print(f"找到 {len(files)} 个文件:")
        for i, f in enumerate(files, 1):
            print(f"  {i}. {f}")

        # 下载文件
        output_dir = "data/processed"
        os.makedirs(output_dir, exist_ok=True)

        downloaded = []
        for file_url in files:
            filename = os.path.basename(file_url)
            output_path = os.path.join(output_dir, filename)

            print(f"\n正在下载: {filename}")
            try:
                # 对于FTP，使用stream下载
                response = requests.get(file_url, stream=True, timeout=60)
                response.raise_for_status()

                total_size = int(response.headers.get('content-length', 0))

                with open(output_path, 'wb') as f:
                    if total_size > 0:
                        downloaded_size = 0
                        for chunk in response.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                                downloaded_size += len(chunk)
                                progress = (downloaded_size / total_size) * 100
                                print(f"\r进度: {progress:.1f}%", end='')
                    else:
                        f.write(response.content)

                print(f"\n✓ 下载完成: {output_path}")
                downloaded.append(output_path)

            except Exception as e:
                print(f"\n✗ 下载失败: {e}")

        return downloaded

    except Exception as e:
        print(f"错误: {e}")
        return []

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("用法: python download_geo_supplementary.py <GSE_ID>")
        print("示例: python download_geo_supplementary.py GSE159217")
        sys.exit(1)

    gse_id = sys.argv[1]
    get_supplementary_files(gse_id)
