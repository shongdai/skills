#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
check_expMatrix_all.py — DSH-GEO 下载结果统一检查器（CLI 版）
=============================================================

适配本技能 4 个下载脚本的所有产物，自动识别每个 GSE 文件夹的下载类型：

  - probe      : GEO_download_probe.R       (Series Matrix 已归一化芯片)
  - cel        : GEO_download_cel.R         (CEL RAW + RMA)
  - ncbicount  : GEO_download_ncbicount.R   (NCBI 整理 raw_counts / FPKM / TPM)
  - count      : GEO_download_count.R       (作者 supplementary count)
  - raw        : 仅有原始下载文件，未生成 expMatrix
  - unknown    : 既无 series_matrix 又无 expMatrix

检查项：
  - 文件齐全度：clinical / expMatrix / probe2gene / group / diff
  - 数据质量：形状 / 数值范围 / 是否含 NA / 是否含负值 / symbol 含 '---'
  - 多平台：series_matrix 数量 vs expMatrix 数量
  - 未处理完：series_matrix 数 > expMatrix 数

CLI 用法：
  python GEO_expMatrix_check.py --dir ./
  python GEO_expMatrix_check.py --dir ./data --output ./check.txt
  python GEO_expMatrix_check.py                  # 默认扫描当前目录

输出：
  控制台打印 + 默认在扫描目录下生成 check_result.txt
"""

import argparse
import os
import re
import sys
import pandas as pd


# ============================================================
# 1. 单文件质量检查
# ============================================================
def check_single_matrix(file_path, folder_name):
    """检查单个 expMatrix*.csv 的数据质量。"""
    fname = os.path.basename(file_path)
    type_name = (fname.replace('expMatrix_', '')
                      .replace('expMatrix', '')
                      .replace('.csv', '')) or 'Single'

    try:
        df = pd.read_csv(file_path, index_col=0)

        numeric_df = df.select_dtypes(include=['number'])
        has_na = bool(numeric_df.isna().any().any())
        na_count = int(numeric_df.isna().sum().sum())

        has_negative = bool((numeric_df < 0).any().any())
        negative_count = int((numeric_df < 0).sum().sum())

        min_value = float(numeric_df.min().min()) if not numeric_df.empty else None
        max_value = float(numeric_df.max().max()) if not numeric_df.empty else None

        return {
            'folder': folder_name,
            'type': type_name,
            'file_name': fname,
            'exists': True,
            'error': None,
            'data_shape': df.shape,
            'has_na': has_na,
            'na_count': na_count,
            'has_negative': has_negative,
            'negative_count': negative_count,
            'min_value': min_value,
            'max_value': max_value,
        }
    except Exception as e:
        return {
            'folder': folder_name,
            'type': type_name,
            'file_name': fname,
            'exists': True,
            'error': str(e),
            'data_shape': None,
            'has_na': None,
            'na_count': None,
            'has_negative': None,
            'negative_count': None,
            'min_value': None,
            'max_value': None,
        }


# ============================================================
# 2. 文件夹类型自动识别
# ============================================================
NCBI_TYPES = {'Count', 'FPKM', 'TPM'}
EXP_RE = re.compile(r'^expMatrix(?:_(.+))?\.csv$')


def classify_folder(files, exp_types):
    """
    根据文件夹内容推断下载类型。
    参数：
      files     : set，文件夹内全部文件名（不含子目录）
      exp_types : set，所有 expMatrix_<type>.csv 中 <type> 部分的集合
                  （type 为空表示 expMatrix.csv 单平台）
    返回：probe / cel / ncbicount / count / raw / unknown
    """
    fnames_lc = ' '.join(f.lower() for f in files)
    has_series_matrix = '_series_matrix.txt.gz' in fnames_lc
    has_cel = '.cel.gz' in fnames_lc or '.cel"' in fnames_lc or '_raw.tar' in fnames_lc
    has_probe2gene = any(f.startswith('probe2gene') for f in files)
    has_ncbi_count = bool(exp_types & NCBI_TYPES)

    if has_ncbi_count:
        return 'ncbicount'

    if has_cel and has_probe2gene:
        return 'cel'

    if has_series_matrix and has_probe2gene:
        return 'probe'

    if exp_types and not has_series_matrix and not has_probe2gene:
        return 'count'

    if has_series_matrix and not exp_types:
        return 'raw'

    return 'unknown'


# ============================================================
# 3. 主扫描
# ============================================================
def check_expMatrix(base_dir, gse_filter=None):
    """遍历 base_dir 下所有 GSE 子文件夹，返回 (results, folder_info)。

    gse_filter : 可选的 GSE 白名单（大小写不敏感）。传入后只检查匹配的子文件夹。
    """
    results = []
    folder_info = {}

    if not os.path.isdir(base_dir):
        raise FileNotFoundError(f'目录不存在: {base_dir}')

    whitelist = {g.upper() for g in gse_filter} if gse_filter else None

    for folder_name in sorted(os.listdir(base_dir)):
        folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue
        if not folder_name.upper().startswith('GSE'):
            continue
        if whitelist is not None and folder_name.upper() not in whitelist:
            continue

        files = os.listdir(folder_path)
        files_set = set(files)

        series_matrix_files = sorted([f for f in files if f.endswith('_series_matrix.txt.gz')])
        exp_matrix_files = sorted([f for f in files if EXP_RE.match(f)])
        probe2gene_files = sorted([f for f in files if f.startswith('probe2gene')])

        exp_types = set()
        for f in exp_matrix_files:
            m = EXP_RE.match(f)
            suffix = m.group(1) if m and m.group(1) else ''
            exp_types.add(suffix)

        clinical_files = sorted([f for f in files if f.startswith('clinical') and f.endswith('.csv')])
        group_files = sorted([f for f in files if f.startswith('group') and f.endswith('.txt')])
        diff_files = sorted([f for f in files if f.startswith('diff') and f.endswith('.txt')])

        download_type = classify_folder(files_set, exp_types)

        folder_info[folder_name] = {
            'download_type': download_type,
            'series_matrix_count': len(series_matrix_files),
            'series_matrix_files': series_matrix_files,
            'exp_matrix_count': len(exp_matrix_files),
            'clinical_count': len(clinical_files),
            'probe2gene_count': len(probe2gene_files),
            'group_count': len(group_files),
            'diff_count': len(diff_files),
        }

        if not exp_matrix_files:
            results.append({
                'folder': folder_name,
                'type': None,
                'file_name': None,
                'exists': False,
                'error': None,
                'data_shape': None,
                'has_na': None, 'na_count': None,
                'has_negative': None, 'negative_count': None,
                'min_value': None, 'max_value': None,
            })
            continue

        for exp_file in exp_matrix_files:
            exp_path = os.path.join(folder_path, exp_file)
            results.append(check_single_matrix(exp_path, folder_name))

    return results, folder_info


# ============================================================
# 4. 报告生成
# ============================================================
def _num(v, fmt='{}'):
    if v is None:
        return 'N/A'
    return fmt.format(v)


def render_lines(results, folder_info):
    """生成报告行列表，供打印 / 写文件复用。"""
    lines = []

    def pr(s=''):
        lines.append(s)

    pr('=' * 160)
    pr('GEO下载结果统计')
    pr('=' * 160)

    folders_with_file = [r for r in results if r['exists']]
    folders_without_file = [r for r in results if not r['exists']]
    unique_missing = sorted({r['folder'] for r in folders_without_file})

    # ---------- 4.1 类型分布 ----------
    type_counter = {}
    for info in folder_info.values():
        type_counter[info['download_type']] = type_counter.get(info['download_type'], 0) + 1

    pr('\n【下载类型分布】')
    pr('-' * 60)
    pr(f"{'类型':<15} {'GSE 数量':<10}")
    pr('-' * 60)
    type_label = {
        'probe':     'probe (芯片)',
        'cel':       'cel (CEL RAW)',
        'ncbicount': 'ncbicount (NCBI RNA-seq)',
        'count':     'count (作者 supplementary)',
        'raw':       'raw (已下未处理)',
        'unknown':   'unknown',
    }
    for k in ('probe', 'cel', 'ncbicount', 'count', 'raw', 'unknown'):
        if k in type_counter:
            pr(f"  {type_label[k]:<25} {type_counter[k]:<10}")

    # ---------- 4.2 有 expMatrix 的文件夹质量检查 ----------
    unique_done = sorted({r['folder'] for r in folders_with_file})
    pr(f'\n【存在 expMatrix*.csv 的文件夹】({len(unique_done)} 个，共 {len(folders_with_file)} 个矩阵）')
    pr('-' * 140)
    header = (f"{'文件夹':<14} {'下载类型':<11} {'矩阵类型':<10} "
              f"{'基因×样本':<16} {'NA数量':<9} {'负值数量':<9} {'范围':<20}")
    pr(header)
    pr('-' * 140)

    def row_key(r):
        return (r['folder'], r['type'] or '')

    def _range_str(r):
        if r['min_value'] is None or r['max_value'] is None:
            return 'N/A'
        return f"{r['min_value']:.2f}-{r['max_value']:.2f}"

    for r in sorted(folders_with_file, key=row_key):
        if r.get('error'):
            pr(f"  {r['folder']:<14} [读取失败] {r['file_name']}: {r['error']}")
            continue
        dtype = folder_info[r['folder']]['download_type']
        pr(f"{r['folder']:<14} {dtype:<11} {r['type']:<10} "
           f"{str(r['data_shape']):<16} "
           f"{_num(r['na_count']):<9} {_num(r['negative_count']):<9} "
           f"{_range_str(r):<20}")

    pr('\n' + '=' * 160)

    # ---------- 4.4 汇总 ----------
    folders_with_na = [r for r in folders_with_file if r['has_na']]
    folders_with_neg = [r for r in folders_with_file if r['has_negative']]
    folders_with_err = [r for r in folders_with_file if r.get('error')]

    multi_platform = [(f, info) for f, info in folder_info.items()
                      if info['series_matrix_count'] > 1]

    pr('汇总统计:')
    pr(f"  总文件夹数:              {len(folder_info)}")
    pr(f"  有 expMatrix*.csv:      {len(unique_done)} (共 {len(folders_with_file)} 个矩阵)")
    missing_detail = f"({','.join(unique_missing)})" if unique_missing else ""
    pr(f"  无 expMatrix*.csv:      {len(unique_missing)}{missing_detail}")
    pr(f"  读取失败的矩阵:           {len(folders_with_err)}")
    pr(f"  含 NA 的矩阵:             {len(folders_with_na)}")
    pr(f"  含负值的矩阵:             {len(folders_with_neg)}")

    def _gpl_tag(filename):
        m = re.search(r'-GPL\d+', filename, re.IGNORECASE)
        return m.group(0)[1:].upper() if m else ''

    if multi_platform:
        parts = []
        for folder, info in sorted(multi_platform):
            gpls = '/'.join(_gpl_tag(f) for f in info['series_matrix_files'] if _gpl_tag(f))
            parts.append(f"{folder}:{gpls}" if gpls else folder)
        multi_detail = ';'.join(parts)
        pr(f"  有个多平台:               {len(multi_platform)}({multi_detail})")
    else:
        pr(f"  有个多平台:               0")

    # ---------- 4.5 未处理完 ----------
    incomplete = [(f, info) for f, info in folder_info.items()
                  if info['series_matrix_count'] > 0
                  and info['exp_matrix_count'] < info['series_matrix_count']]
    if incomplete:
        pr(f'\n【未处理完的文件夹】({len(incomplete)} 个)')
        pr('-' * 90)
        pr(f"{'文件夹':<14} {'series_matrix':<14} {'clinical':<10} {'expMatrix':<10}")
        pr('-' * 90)
        for folder, info in sorted(incomplete):
            pr(f"  {folder:<12} {info['series_matrix_count']:<14} "
               f"{info['clinical_count']:<10} {info['exp_matrix_count']:<10}")

    # ---------- 4.6 质量异常详情 ----------
    if folders_with_err:
        pr(f'\n【读取失败的矩阵】{len(folders_with_err)} 个')
        pr('-' * 100)
        for r in folders_with_err:
            pr(f"  {r['folder']}/{r['file_name']}: {r['error']}")

    if folders_with_neg:
        pr(f'\n【含负值的矩阵详情】{len(folders_with_neg)} 个')
        pr('-' * 80)
        pr(f"{'文件夹':<14} {'矩阵类型':<10} {'负值数量':<10} {'范围':<20}")
        pr('-' * 80)
        for r in sorted(folders_with_neg, key=row_key):
            pr(f"  {r['folder']:<12} {r['type']:<10} {r['negative_count']:<10} "
               f"{_range_str(r):<20}")

    if folders_with_na:
        pr(f'\n【含 NA 的矩阵详情】{len(folders_with_na)} 个')
        pr('-' * 70)
        pr(f"{'文件夹':<14} {'矩阵类型':<10} {'NA 数量':<10}")
        pr('-' * 70)
        for r in sorted(folders_with_na, key=row_key):
            pr(f"  {r['folder']:<12} {r['type']:<10} {r['na_count']:<10}")

    return lines


def print_results(results, folder_info):
    for line in render_lines(results, folder_info):
        print(line)


def save_results(results, folder_info, output_path):
    lines = render_lines(results, folder_info)
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


# ============================================================
# 5. CLI 入口
# ============================================================
def build_parser():
    p = argparse.ArgumentParser(
        prog="GEO_expMatrix_check.py",
        description="DSH-GEO 下载结果统一检查器 — 扫描指定目录下所有 GSE 子文件夹，"
                    "检查表达矩阵文件齐全度与数据质量。",
    )
    p.add_argument("--dir", default=".",
                   help="扫描的根目录（默认当前目录）")
    p.add_argument("--gse", default=None,
                   help="只检查指定的 GSE（逗号/空格分隔），不传则扫描全部 GSE 子文件夹")
    p.add_argument("--output", default=None,
                   help="结果报告输出路径（默认 <dir>/check_result.txt）")
    p.add_argument("--no-save", action="store_true",
                   help="仅打印到控制台，不写文件")
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()

    base_dir = os.path.abspath(args.dir)
    gse_filter = None
    if args.gse:
        gse_filter = [g.strip() for g in re.split(r"[,\s]+", args.gse) if g.strip()]
        print(f'扫描目录: {base_dir}  仅检查: {",".join(gse_filter)}')
    else:
        print(f'扫描目录: {base_dir}')

    results, folder_info = check_expMatrix(base_dir, gse_filter=gse_filter)

    if not folder_info:
        print('[警告] 当前目录下未找到任何 GSE 开头的子文件夹')
        return 1

    print_results(results, folder_info)

    if not args.no_save:
        output_file = args.output or os.path.join(base_dir, 'check_result.txt')
        save_results(results, folder_info, output_file)
        print(f'\n结果已保存至: {output_file}')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[中断]")
        sys.exit(130)
