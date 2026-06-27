#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
GEO_smart_sniff.py - GSE 元数据嗅探器（v2.0 轻量版）
====================================================

仅查询 NCBI 元数据（KB 级），不下载 SOFT/Matrix。秒级响应。

输出：
  - data_type       : microarray / rna_seq / non_coding_rna / methylation / unknown
  - species         : human / mouse / rat / ...
  - org_db          : org.Hs.eg.db / org.Mm.eg.db / ...
  - id_type         : ENSEMBL / ENTREZID / SYMBOL（仅 rna_seq 且无 NCBI 标准化）
  - script + args   : 推荐的 R 脚本及参数
  - command         : 完整命令行

用法：
    python GEO_smart_sniff.py --gse GSE266899
    python GEO_smart_sniff.py --gse GSE266899,GSE5281 --json

依赖：
    pip install biopython   (GEOparse 不再需要)
"""

import argparse
import json
import os
import re
import sys
import time
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import xml.etree.ElementTree as ET

try:
    from Bio import Entrez
except ImportError:
    print("[错误] 需要安装 biopython: pip install biopython", file=sys.stderr)
    sys.exit(1)

Entrez.email = "geo_smart_sniff@example.com"  # NCBI 要求

# ---------------------------------------------------------------------------
# 默认代理（可在 CLI 覆盖）
# ---------------------------------------------------------------------------
DEFAULT_PROXY = "http://127.0.0.1:7897"


# ---------------------------------------------------------------------------
# 平台 / 物种映射
# ---------------------------------------------------------------------------
MICROARRAY_PLATFORMS = {
    "GPL570", "GPL571", "GPL96", "GPL1261", "GPL1355", "GPL6244",
    "GPL6883", "GPL6884", "GPL8321", "GPL11532", "GPL14668",
    "GPL6097", "GPL6102", "GPL6104", "GPL29849",
    "GPL6480", "GPL4133", "GPL4134", "GPL7202", "GPL15331",
    "GPL6947", "GPL10558", "GPL13376", "GPL8490",
}

RNA_SEQ_PLATFORMS = {
    "GPL11154", "GPL13112", "GPL16791", "GPL18573", "GPL20301",
    "GPL21273", "GPL21697", "GPL23227", "GPL24676", "GPL24943",
    "GPL25248", "GPL25410",
}

SPECIES_MAP = {
    "Homo sapiens":              ("human",     "org.Hs.eg.db"),
    "Mus musculus":              ("mouse",     "org.Mm.eg.db"),
    "Rattus norvegicus":         ("rat",       "org.Rn.eg.db"),
    "Danio rerio":               ("zebrafish", "org.Dr.eg.db"),
    "Saccharomyces cerevisiae":  ("yeast",     "org.Sc.sgd.db"),
    "Caenorhabditis elegans":    ("worm",      "org.Ce.eg.db"),
    "Drosophila melanogaster":   ("fly",       "org.Dm.eg.db"),
    "Macaca mulatta":            ("macaque",   "org.Mmu.eg.db"),
    "Sus scrofa":                ("pig",       "org.Ss.eg.db"),
    "Gallus gallus":             ("chicken",   "org.Gg.eg.db"),
}


# ---------------------------------------------------------------------------
# 1. Entrez esearch: GSE -> GDS UID
# ---------------------------------------------------------------------------
def gse_to_gds_uid(gse_id, timeout=10):
    """轻量级：把 GSE 号转成 GDS 数据库的 UID。"""
    h = Entrez.esearch(db="gds", term=f"{gse_id}[Accession]", retmode="xml", timeout=timeout)
    rec = Entrez.read(h); h.close()
    if not rec.get("IdList"):
        return None
    return rec["IdList"][0]


# ---------------------------------------------------------------------------
# 2. Entrez esummary: UID -> 元数据 dict
# ---------------------------------------------------------------------------
def gds_summary(uid, timeout=10):
    """拉 GDS summary（含 platform、taxon、sample count、title、gdsType）。"""
    h = Entrez.esummary(db="gds", id=uid, retmode="xml", timeout=timeout)
    rec = Entrez.read(h); h.close()
    if not rec:
        return None
    s = rec[0] if isinstance(rec, list) else rec
    gpl_raw = s.get("GPL", "")
    if isinstance(gpl_raw, list):
        gpl_raw = gpl_raw[0] if gpl_raw else ""
    gpl = str(gpl_raw)
    if gpl and not gpl.upper().startswith("GPL"):
        gpl = "GPL" + gpl
    taxon = s.get("taxon", "")
    if isinstance(taxon, list):
        taxon = taxon[0] if taxon else ""
    return {
        "accession":  str(s.get("Accession", "")),
        "title":      str(s.get("title", "")),
        "gpl":        gpl,
        "gpl_title":  str(s.get("platform_title", "") or s.get("PlatformTitle", "")),
        "n_samples":  int(s.get("n_samples", 0) or 0),
        "gds_type":   str(s.get("gdsType", "")),
        "entry_type": str(s.get("entryType", "")),
        "taxon":      str(taxon),
    }


# ---------------------------------------------------------------------------
# 3. Entrez efetch: UID -> full XML（含 supplementary file / entry type）
# ---------------------------------------------------------------------------
def gds_full_xml(uid, timeout=15):
    """拉 GDS 完整 XML，用于准确判断 entryType 和 sample 物种。"""
    h = Entrez.efetch(db="gds", id=uid, retmode="xml", timeout=timeout)
    xml_text = h.read(); h.close()
    if isinstance(xml_text, bytes):
        xml_text = xml_text.decode("utf-8", errors="ignore")
    return xml_text


def parse_organism_from_xml(xml_text):
    """从 GDS XML 中提取第一个 sample 的 organism。"""
    try:
        root = ET.fromstring(xml_text)
        # 路径：./DocumentSummary/Records/Record/Sample/Channel/Organism
        for elem in root.iter():
            tag = elem.tag.split("}")[-1]
            if tag == "Organism" and elem.text:
                return elem.text.strip()
    except ET.ParseError:
        pass
    return ""


# ---------------------------------------------------------------------------
# 4. 嗅探 supplementary 文件（HTTP GET GEO FTP 目录列表）
# ---------------------------------------------------------------------------
def _http_get(url, timeout=10):
    """带代理的 HTTP GET（urllib 自动读 HTTP_PROXY/HTTPS_PROXY 环境变量）。"""
    req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as r:
        return r.read()


def _http_head(url, timeout=8):
    """带代理的 HTTP HEAD。"""
    req = Request(url, method="HEAD", headers={"User-Agent": "Mozilla/5.0"})
    with urlopen(req, timeout=timeout) as r:
        return r.status


def _check_proxy_connectivity(proxy_url, test_url="https://www.ncbi.nlm.nih.gov/",
                              timeout=4):
    """
    快速检测代理是否可用。返回 (ok, latency_ms, reason)。
    - ok=True  : 代理可用
    - ok=False : 代理不通（reason 给出原因）

    检测策略（两步）：
      1. TCP 握手到代理 host:port（3s 超时）→ 代理进程是否在跑
      2. 实际 HTTP GET 一个稳定 URL（4s 超时）→ 代理能否转发 HTTPS
         使用 https://www.ncbi.nlm.nih.gov/ 因为它是稳定的首页，
         而 eutils endpoint 会 404 误导检测。
    """
    import socket
    import time as _t

    # 1) 解析代理 host:port
    try:
        if "://" in proxy_url:
            proxy_url = proxy_url.split("://", 1)[1]
        host, port = proxy_url.rsplit(":", 1)
        port = int(port)
    except Exception as e:
        return False, 0, f"代理 URL 解析失败: {e}"

    # 2) TCP 握手测活（代理进程是否在跑）
    t0 = _t.time()
    try:
        with socket.create_connection((host, port), timeout=timeout):
            pass
    except Exception as e:
        return False, 0, f"代理 {host}:{port} 不可达: {e.__class__.__name__}"

    tcp_ms = round((_t.time() - t0) * 1000, 1)

    # 3) 实际 HTTP 测试（确认代理能正确转发 HTTPS 到目标）
    t0 = _t.time()
    try:
        req = Request(test_url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=timeout) as r:
            if 200 <= r.status < 400:
                return True, round((_t.time() - t0) * 1000, 1), "OK"
            return False, tcp_ms, f"代理返回 {r.status}"
    except HTTPError as e:
        # HTTPError 也算"代理工作"（它能解析目标主机）
        if e.code in (403, 404):
            return True, round((_t.time() - t0) * 1000, 1), f"OK (代理可达，目标 {e.code})"
        return False, tcp_ms, f"代理 HTTP 错误 {e.code}"
    except Exception as e:
        return False, tcp_ms, f"代理 HTTP 失败: {e.__class__.__name__}"


# ---------------------------------------------------------------------------
# 4. 嗅探 supplementary 文件（HTTP GET GEO FTP 目录列表）
# ---------------------------------------------------------------------------
def fetch_suppl_files(gse_id, timeout=10):
    """通过 HTTP GET FTP 目录的 HTML 列出所有 supplementary 文件名。"""
    digits = gse_id.replace("GSE", "")
    if len(digits) < 3:
        return []
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{digits[:-3]}nnn/{gse_id}/suppl/"
    try:
        html = _http_get(url, timeout=timeout).decode("utf-8", errors="ignore")
    except Exception:
        return []

    # 提取 href="..." 中的文件名（排除父目录与外链）
    files = re.findall(r'href="([^"]+)"', html)
    files = [f for f in files
             if not f.startswith("?") and not f.startswith("/")
             and not f.startswith("http")
             and f not in ("../", "/")]
    return files


# ---------------------------------------------------------------------------
# 4b. 探测 NCBI 标准化文件（HEAD 请求）
#    URL 格式参考 GEO_download_ncbicount.R：
#    https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc={GSE}&format=file&file={GSE}_{type}_{genome}_NCBI.tsv.gz
# ---------------------------------------------------------------------------
NCBI_GENOME_TAG = {
    "human": "GRCh38.p13",
    "mouse": "GRCm39",
    "rat":   "mRatBN7.2",
}


def check_ncbi_norm_exists(gse_id, species, timeout=8):
    """用 HEAD 请求探测 NCBI 标准化文件是否存在。返回 bool。"""
    genome = NCBI_GENOME_TAG.get(species, "")
    if not genome:
        return False
    url = (f"https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts"
           f"&acc={gse_id}&format=file"
           f"&file={gse_id}_raw_counts_{genome}_NCBI.tsv.gz")
    try:
        return _http_head(url, timeout=timeout) == 200
    except Exception:
        return False


# ---------------------------------------------------------------------------
# 5. 嗅探 series matrix 第一列前几行（用于推断 ID 类型）
# ---------------------------------------------------------------------------
def fetch_series_matrix_first_col(gse_id, timeout=15):
    """从 series matrix 文件流式读取第一列前 5 行（不全量下载）。"""
    import gzip
    import io
    digits = gse_id.replace("GSE", "")
    if len(digits) < 3:
        return []
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE{digits[:-3]}nnn/{gse_id}/matrix/{gse_id}_series_matrix.txt.gz"
    try:
        req = Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(req, timeout=timeout) as resp:
            # 流式解压：逐行读取，读到 5 条数据行即停止
            with gzip.GzipFile(fileobj=resp) as gz:
                lines = []
                for line in gz:
                    if line.startswith(b"!") or line.startswith(b'"'):
                        continue
                    parts = line.decode("utf-8", errors="ignore").split("\t")
                    if parts:
                        lines.append(parts[0])
                    if len(lines) >= 5:
                        break
                return lines
    except Exception:
        return []


# ---------------------------------------------------------------------------
# 6. 主嗅探函数
# ---------------------------------------------------------------------------
def sniff_gse(gse_id, verbose=False):
    gse_id = gse_id.strip().upper()
    if not gse_id.startswith("GSE"):
        gse_id = "GSE" + gse_id

    t0 = time.time()
    info = {
        "gse_id":            gse_id,
        "ok":                False,
        "error":             None,
        "elapsed_sec":       None,
        "data_type":         None,
        "species":           None,
        "org_db":            None,
        "id_type":           None,
        "platform_id":       None,
        "n_samples":         0,
        "suppl_files":       [],
        "has_ncbi_norm":     False,
        "has_author_counts": False,
        "has_cel":           False,
        "script":            None,
        "args":              [],
        "command":           None,
        "title":             None,
    }

    # ---- Step 1: GSE -> GDS UID ----
    try:
        uid = gse_to_gds_uid(gse_id)
    except Exception as e:
        info["error"] = f"esearch 失败: {e}"
        info["elapsed_sec"] = round(time.time() - t0, 2)
        return info
    if not uid:
        info["error"] = f"未找到 {gse_id}"
        info["elapsed_sec"] = round(time.time() - t0, 2)
        return info

    # ---- Step 2: UID -> summary ----
    try:
        summ = gds_summary(uid)
    except Exception as e:
        info["error"] = f"esummary 失败: {e}"
        info["elapsed_sec"] = round(time.time() - t0, 2)
        return info
    if not summ:
        info["error"] = "esummary 返回空"
        info["elapsed_sec"] = round(time.time() - t0, 2)
        return info

    info["title"]       = summ["title"]
    info["platform_id"] = summ["gpl"]
    info["n_samples"]   = summ["n_samples"]

    # ---- Step 3: data_type 推断（基于 entryType + gdsType + platform） ----
    et = (summ.get("entry_type") or "").lower()
    gt = (summ.get("gds_type") or "").lower()
    platform = info["platform_id"] or ""

    if "expression profiling by array" in et or "expression profiling by array" in gt:
        info["data_type"] = "microarray"
    elif "expression profiling by high throughput sequencing" in et or \
         "expression profiling by high throughput sequencing" in gt:
        info["data_type"] = "rna_seq"
    elif "non-coding rna profiling" in et or "non-coding rna profiling" in gt:
        info["data_type"] = "non_coding_rna"
    elif "methylation profiling" in et or "methylation profiling" in gt:
        info["data_type"] = "methylation"
    elif platform in MICROARRAY_PLATFORMS:
        info["data_type"] = "microarray"
    elif platform in RNA_SEQ_PLATFORMS:
        info["data_type"] = "rna_seq"
    else:
        info["data_type"] = "unknown"

    # ---- Step 4: organism (优先用 esummary.taxon，回退到 XML) ----
    org = summ.get("taxon", "")
    if not org:
        try:
            xml = gds_full_xml(uid)
            org = parse_organism_from_xml(xml)
        except Exception:
            org = ""
    if org in SPECIES_MAP:
        info["species"], info["org_db"] = SPECIES_MAP[org]
    elif org:
        info["species"] = org.lower().replace(" ", "_")

    # ---- Step 5: supplementary files (FTP suppl/) ----
    try:
        info["suppl_files"] = fetch_suppl_files(gse_id)
    except Exception:
        info["suppl_files"] = []

    fnames_lc = " ".join(info["suppl_files"]).lower()
    info["has_author_counts"] = bool(re.search(r"(raw_?count|gene[_ ]?count|read[_ ]?count|expression)", fnames_lc))
    info["has_cel"]           = any(f.lower().endswith((".cel", ".cel.gz")) for f in info["suppl_files"])

    # ---- Step 5b: NCBI 标准化文件探测（HEAD 请求，仅 rna_seq 时） ----
    if info["data_type"] == "rna_seq" and info["species"] in ("human", "mouse", "rat"):
        try:
            info["has_ncbi_norm"] = check_ncbi_norm_exists(gse_id, info["species"])
        except Exception:
            info["has_ncbi_norm"] = False

    # ---- Step 6: ID 类型（rna_seq 且无 NCBI 标准化时） ----
    if info["data_type"] == "rna_seq" and not info["has_ncbi_norm"]:
        try:
            first_ids = fetch_series_matrix_first_col(gse_id)
            sample_id = first_ids[0] if first_ids else ""
        except Exception:
            sample_id = ""
        if sample_id.startswith("ENSMUSG") or sample_id.startswith("ENSG"):
            info["id_type"] = "ENSEMBL"
        elif re.match(r"^\d+(\.\d+)?$", sample_id):
            info["id_type"] = "ENTREZID"
        elif re.match(r"^[A-Za-z][A-Za-z0-9_.\-]+$", sample_id) and not sample_id.upper().startswith("ENS"):
            info["id_type"] = "SYMBOL"
        else:
            info["id_type"] = "ENSEMBL"

    # ---- Step 7: 决策 ----
    if info["data_type"] == "microarray":
        if info["has_cel"]:
            info["script"] = "GEO_download_cel.R"
            info["args"]  = ["--gse", gse_id]
        else:
            info["script"] = "GEO_download_probe.R"
            info["args"]  = ["--gse", gse_id, "--diff", "TRUE"]
    elif info["data_type"] == "rna_seq":
        # ncbicount.R v2.1 仅支持 human；非 human 强制走 count.R
        use_ncbicount = info["has_ncbi_norm"] and info["species"] == "human"
        if use_ncbicount:
            info["script"] = "GEO_download_ncbicount.R"
            info["args"]  = ["--gse", gse_id]
            # GEO_download_ncbicount.R v2.1 已 human-only，不再接受 --species
        else:
            info["script"] = "GEO_download_count.R"
            info["args"]  = ["--gse", gse_id, "--diff", "TRUE"]
            if info["species"]:
                info["args"] += ["--species", info["species"]]
            if info["id_type"]:
                info["args"] += ["--from-type", info["id_type"]]
    else:
        info["script"] = "GEO_download_probe.R"
        info["args"]  = ["--gse", gse_id]

    info["command"] = f"Rscript scripts/{info['script']} " + " ".join(info["args"])
    info["ok"] = True
    info["elapsed_sec"] = round(time.time() - t0, 2)
    return info


# ---------------------------------------------------------------------------
# 输出
# ---------------------------------------------------------------------------
def print_text(info):
    print("=" * 64)
    print(f"GSE        : {info['gse_id']}    [{info['data_type']}]")
    print(f"标题       : {(info['title'] or '')[:80]}")
    print(f"平台       : {info['platform_id']}    样本数: {info['n_samples']}")
    print(f"物种       : {info['species']}    -> OrgDb = {info['org_db']}")
    print(f"ID 类型    : {info['id_type']}")
    print("-" * 64)
    print(f"特征       : NCBI_norm={info['has_ncbi_norm']}, "
          f"author_cnt={info['has_author_counts']}, "
          f"CEL={info['has_cel']}")
    if info["suppl_files"]:
        print(f"补充文件   : {len(info['suppl_files'])} 个")
        for f in info["suppl_files"][:6]:
            print(f"   - {f}")
        if len(info["suppl_files"]) > 6:
            print(f"   ... 还有 {len(info['suppl_files']) - 6} 个")
    print("-" * 64)
    if info["ok"]:
        print(f"[推荐脚本]  {info['script']}")
        print(f"[推荐命令]  {info['command']}")
    else:
        print(f"[错误]  {info['error']}")
    print(f"[耗时]    {info['elapsed_sec']}s")
    print("=" * 64)


def main():
    ap = argparse.ArgumentParser(description="GSE 元数据嗅探器（轻量版）")
    ap.add_argument("--gse", required=True, help="GSE 号（多个用逗号分隔）")
    ap.add_argument("--json", action="store_true", help="输出 JSON 格式")
    ap.add_argument("--proxy", default=DEFAULT_PROXY,
                    help=f"HTTP/HTTPS 代理 URL（默认 {DEFAULT_PROXY}）")
    ap.add_argument("--no-proxy", dest="no_proxy", action="store_true",
                    help="禁用代理（直连）")
    ap.add_argument("--skip-proxy-check", dest="skip_check", action="store_true",
                    help="跳过代理可用性检查（节省 ~3s）")
    args = ap.parse_args()

    # 决定最终是否使用代理
    use_proxy = False
    proxy_url = None
    proxy_info = ("(未配置)", 0, "")

    if args.no_proxy:
        # 用户强制禁用
        use_proxy = False
    elif args.proxy:
        proxy_url = args.proxy
        if args.skip_check:
            use_proxy = True
            proxy_info = (proxy_url, 0, "跳过检查")
        else:
            # 默认先检查代理是否可用
            ok, latency, reason = _check_proxy_connectivity(proxy_url)
            proxy_info = (proxy_url, latency, reason)
            if ok:
                use_proxy = True
            else:
                # 代理不可用 → 自动回退到直连
                use_proxy = False
                print(f"[代理检测] {proxy_url} 不可用 ({reason})，自动回退到直连\n",
                      file=sys.stderr)

    # 配置代理环境变量
    if use_proxy and proxy_url:
        os.environ["HTTP_PROXY"]  = proxy_url
        os.environ["HTTPS_PROXY"] = proxy_url
        os.environ["http_proxy"]  = proxy_url
        os.environ["https_proxy"] = proxy_url
    else:
        for k in ("HTTP_PROXY", "HTTPS_PROXY", "http_proxy", "https_proxy"):
            os.environ.pop(k, None)

    gse_list = [g.strip() for g in re.split(r"[,\s]+", args.gse) if g.strip()]

    results = []
    for i, g in enumerate(gse_list):
        if i > 0:
            time.sleep(0.34)  # NCBI 速率限制：3 req/s
        results.append(sniff_gse(g))

    # 把代理信息注入到每条结果
    for r in results:
        r["proxy_used"] = proxy_url if use_proxy else None
        r["proxy_check_ms"] = proxy_info[1]
        r["proxy_check_status"] = proxy_info[2]

    if args.json:
        print(json.dumps(results, ensure_ascii=False, indent=2))
    else:
        # 顶部打印代理状态
        if use_proxy:
            print(f"[代理] 启用 {proxy_url}（{proxy_info[1]}ms）\n",
                  file=sys.stderr)
        else:
            print(f"[代理] 直连模式\n", file=sys.stderr)
        for info in results:
            print_text(info)
            print()

    if any(not r["ok"] for r in results):
        sys.exit(1)


if __name__ == "__main__":
    main()
