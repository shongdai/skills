#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEO 预下载器（CLI 版）
---------------------
按"预设 / 类别"批量下载 GEO 公共数据，作为后续 dsh-geo 技能 R 脚本的输入素材。

类别（共 5 个）：
  soft    : SOFT 完整记录（{GSE}_family.soft.gz）
  miniml  : MINiML XML 完整记录（{GSE}_family.xml.tgz）
  matrix  : Series Matrix（{GSE}_series_matrix.txt.gz）  ← R 芯片脚本必需
  suppl   : suppl/ 目录所有补充文件（CEL/作者上传 count/...）  ← R cel/count 脚本必需
  gpl     : GPL 平台注释（NCBI 文本）  ← R 芯片脚本的本地 GPL 回退

预设（对应后续要跑哪个 R 脚本）：
  all       → soft + miniml + matrix + suppl + gpl   （全量备份）
  probe     → matrix + gpl                            （配合 GEO_download_probe.R）
  cel       → matrix + suppl + gpl                    （配合 GEO_download_cel.R）
  ncbicount → matrix                                  （配合 GEO_download_ncbicount.R；NCBI 整理矩阵由 R 脚本另行下载）
  count     → matrix + suppl                          （配合 GEO_download_count.R）
  custom    → 由 --types 指定

CLI 示例：
  python download_GEO.py --gse GSE76262                                  # 默认 all
  python download_GEO.py --gse GSE76262 --preset probe
  python download_GEO.py --gse GSE76262,GSE12345 --preset cel --out ./data
  python download_GEO.py --gse GSE76262 --preset custom --types matrix,gpl
  python download_GEO.py --gse GSE76262 --no-proxy
  python download_GEO.py --gse GSE76262 --proxy http://127.0.0.1:7897 --workers 16

依赖: pip install requests
设计：python 失败不影响 R，R 脚本本身能独立下载（这里只是预下载加速 + 离线兜底）。
"""

import argparse
import re
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

try:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
except ImportError:
    sys.exit("[ERROR] 缺少 requests, 请执行: pip install requests")


# ==================== 默认配置 ====================
DEFAULT_PROXY = "http://127.0.0.1:7897"

ALL_TYPES = ["soft", "miniml", "matrix", "suppl", "gpl"]

PRESETS = {
    "all":       ["soft", "miniml", "matrix", "suppl", "gpl"],
    "probe":     ["matrix", "gpl"],
    "cel":       ["matrix", "suppl", "gpl"],
    "ncbicount": ["matrix"],
    "count":     ["matrix", "suppl"],
}
# ==================================================


# ==================== 运行时配置 ====================
# 由 main() 从 CLI 参数填充；下载函数都读这里
CFG = {
    "gses": [],
    "email": "your.email@example.com",
    "workers": 8,
    "timeout": 30,
    "retries": 3,
    "api_gap": 0.4,
    "proxy": DEFAULT_PROXY,
    "types": set(ALL_TYPES),
    "suppl_block": ["filelist"],
    "out": ".",
}
# ====================================================


_print_lock = Lock()


# ── 工具 ───────────────────────────────────────
def _prefix(gse_id: str) -> str:
    m = re.match(r"^GSE(\d+)$", gse_id)
    if not m or len(m[1]) < 4:
        raise ValueError(f"非法 GSE: {gse_id}")
    return f"GSE{m[1][:-3]}nnn"


def _hsize(n) -> str:
    n = float(n)
    for u in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} TB"


def _norm_gses(val) -> list:
    if isinstance(val, str):
        val = re.split(r"[,\s]+", val)
    return [x.strip().upper() for x in val
            if x.strip() and re.match(r"^GSE\d+$", x.strip().upper())]


# ── URL 构造 ───────────────────────────────────
def _fixed_urls(gse_id: str) -> dict:
    b = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{_prefix(gse_id)}/{gse_id}"
    return {
        "soft":   f"{b}/soft/{gse_id}_family.soft.gz",
        "miniml": f"{b}/miniml/{gse_id}_family.xml.tgz",
        "matrix": f"{b}/matrix/{gse_id}_series_matrix.txt.gz",
    }


def _gpl_url(gpl: str) -> str:
    return (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
            f"?acc={gpl}&targ=self&form=text&view=data")


# ── HTTP Session ───────────────────────────────
def _session() -> requests.Session:
    s = requests.Session()
    if CFG["proxy"]:
        s.proxies = {"http": CFG["proxy"], "https": CFG["proxy"]}
    retry = Retry(total=CFG["retries"], backoff_factor=1.0,
                  status_forcelist=(429, 500, 502, 503, 504),
                  allowed_methods=frozenset(["GET", "HEAD"]))
    adapter = HTTPAdapter(max_retries=retry,
                          pool_connections=CFG["workers"],
                          pool_maxsize=CFG["workers"] * 2)
    s.mount("https://", adapter)
    s.mount("http://", adapter)
    s.headers["User-Agent"] = "Mozilla/5.0 (GEO-Downloader)"
    return s


# ── 下载（带进度条） ────────────────────────────
def _download(session, url: str, dest: Path) -> tuple:
    """返回 (ok, size_bytes, fname, status_msg)。"""
    fname = dest.name
    try:
        head = session.head(url, allow_redirects=True, timeout=CFG["timeout"])
        if head.status_code == 404:
            return False, 0, fname, "下载失败: 404 文件不存在"
        if head.status_code >= 400:
            return False, 0, fname, f"下载失败: HTTP {head.status_code}"

        total = int(head.headers.get("Content-Length", 0))
        already = dest.stat().st_size if dest.exists() else 0
        if total and already >= total:
            return True, total, fname, "已下载,跳过"

        headers = {}
        mode = "wb"
        if already and "bytes" in head.headers.get("Accept-Ranges", "").lower():
            headers["Range"] = f"bytes={already}-"
            mode = "ab"

        with session.get(url, stream=True, timeout=CFG["timeout"], headers=headers) as r:
            r.raise_for_status()
            downloaded = already
            t0 = time.time()
            last_tick = t0
            with open(dest, mode) as f:
                for chunk in r.iter_content(chunk_size=65536):
                    if not chunk:
                        continue
                    f.write(chunk)
                    downloaded += len(chunk)

                    if total:
                        now = time.time()
                        if now - last_tick >= 0.5:
                            elapsed = now - t0
                            speed = (downloaded - already) / elapsed if elapsed > 0 else 0
                            pct = downloaded * 100 / total
                            bar_w = 20
                            filled = int(bar_w * downloaded / total)
                            bar = "#" * filled + "-" * (bar_w - filled)
                            with _print_lock:
                                print(f"\r    {fname}  [{bar}] {pct:5.1f}%  "
                                      f"{_hsize(downloaded)}/{_hsize(total)}  "
                                      f"{_hsize(int(speed))}/s",
                                      end="", flush=True)
                            last_tick = now

        return True, dest.stat().st_size, fname, "下载完成"
    except Exception as e:
        return False, 0, fname, f"下载失败: {e}"


def _batch(session, urls: list, out_dir: Path,
           names: Optional[list] = None) -> list:
    tasks = [(u, out_dir / (names[i] if names else Path(u).name))
             for i, u in enumerate(urls)]
    results = []
    with ThreadPoolExecutor(max_workers=CFG["workers"]) as pool:
        futures = {pool.submit(_download, session, u, d): (u, d)
                   for u, d in tasks}
        for fut in as_completed(futures):
            results.append(fut.result())
    return results


# ── Suppls 抓取 ────────────────────────────────
def _list_suppl(gse_id: str, session) -> list:
    url = (f"https://ftp.ncbi.nlm.nih.gov/geo/series/"
           f"{_prefix(gse_id)}/{gse_id}/suppl/")
    try:
        txt = session.get(url, timeout=CFG["timeout"]).text
        pre = re.search(r"<pre>(.*?)</pre>", txt, re.DOTALL)
        if not pre:
            return []
        files = []
        for m in re.finditer(r'<a href="([^"]+)"', pre.group(1)):
            href = m.group(1).strip()
            if not href or href.startswith(("/", "http")):
                continue
            files.append((href, url + href))
        return files
    except Exception as e:
        print(f"    suppl 抓取失败: {e}")
        return []


# ── GPL 查询 ───────────────────────────────────
def _query_gpls(gse_id: str, session) -> list:
    eu = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    try:
        r = session.get(f"{eu}/esearch.fcgi", params={
            "db": "gds", "term": f"{gse_id}[Accession]",
            "retmax": 1, "email": CFG["email"]
        }, timeout=CFG["timeout"])
        gds = ET.fromstring(r.text).findtext(".//Id")
        if not gds:
            return []
        time.sleep(CFG["api_gap"])

        r = session.get(f"{eu}/esummary.fcgi", params={
            "db": "gds", "id": gds, "email": CFG["email"]
        }, timeout=CFG["timeout"])
        gpls = []
        for item in ET.fromstring(r.text).iter("Item"):
            if item.get("Name") == "GPL" and item.text:
                for p in re.split(r"[;,/\s]+", item.text):
                    p = p.strip().upper()
                    if p.isdigit():
                        p = "GPL" + p
                    if re.match(r"^GPL\d+$", p) and p not in gpls:
                        gpls.append(p)
        return gpls
    except Exception as e:
        print(f"    GPL 查询失败: {e}")
        return []


# ── 单 GSE 处理 ────────────────────────────────
def _process_one(gse_id: str, session) -> dict:
    print(f"\n{'='*60}")
    print(f"  {gse_id}")
    print(f"{'='*60}")

    out_dir = Path(CFG["out"]) / gse_id
    out_dir.mkdir(parents=True, exist_ok=True)

    records = []

    def _record(ftype: str, r: list):
        seen = set()
        for ok, b, name, status in r:
            if name in seen:
                continue
            seen.add(name)
            records.append((ftype, name, b, status))

    # 1) 固定文件
    urls = _fixed_urls(gse_id)
    for t in ["soft", "miniml", "matrix"]:
        if t not in CFG["types"]:
            continue
        r = _batch(session, [urls[t]], out_dir)
        _record(t, r)

    # 2) suppl
    if "suppl" in CFG["types"]:
        suppl = _list_suppl(gse_id, session)
        blocks = [b.lower() for b in CFG["suppl_block"] if b]
        suppl = [(n, u) for n, u in suppl
                 if not any(b in n.lower() for b in blocks)]
        if suppl:
            names, sul = zip(*suppl)
            r = _batch(session, list(sul), out_dir, names=list(names))
            _record("suppl", r)

    # 3) GPL
    gpls = []
    if "gpl" in CFG["types"]:
        gpls = _query_gpls(gse_id, session)
        if gpls:
            r = _batch(session, [_gpl_url(g) for g in gpls], out_dir,
                       names=[f"{g}.txt" for g in gpls])
            _record("GPL", r)

    ok_count = sum(1 for _, _, _, s in records if s in ("下载完成", "已下载,跳过"))
    fail_count = len(records) - ok_count
    total_bytes = sum(b for _, _, b, _ in records)

    print(f"\n  {gse_id} 下载 {len(records)} 个文件, "
          f"{ok_count} 个成功, {fail_count} 个失败, 总大小 {_hsize(total_bytes)}")
    print(f"  {'文件类型':<8s} {'文件名称':<50s} {'文件大小':>10s}  {'状态'}")
    print(f"  {'-'*8} {'-'*50} {'-'*10}  {'-'*12}")
    for ftype, name, b, status in records:
        sz = _hsize(b) if b else "-"
        st = "成功" if status in ("下载完成", "已下载,跳过") else "失败"
        print(f"  {ftype:<8s} {name:<50s} {sz:>10s}  {st}")

    return {gse_id: {"records": records, "gpls": gpls,
                     "ok": ok_count, "fail": fail_count,
                     "bytes": total_bytes}}


# ── CLI 入口 ───────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="download_GEO.py",
        description="GEO 预下载器 — 按预设/类别批量下载，为后续 R 脚本提供素材",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "预设映射（对应后续 R 脚本）:\n"
            "  all       = soft+miniml+matrix+suppl+gpl   全量备份\n"
            "  probe     = matrix+gpl                     → GEO_download_probe.R\n"
            "  cel       = matrix+suppl+gpl               → GEO_download_cel.R\n"
            "  ncbicount = matrix                         → GEO_download_ncbicount.R\n"
            "  count     = matrix+suppl                   → GEO_download_count.R\n"
            "  custom    = 由 --types 指定\n"
        ),
    )
    p.add_argument("--gse", required=True,
                   help="GSE 号，多个用逗号或空格分隔，如 GSE76262,GSE12345")
    p.add_argument("--preset", default="all",
                   choices=list(PRESETS.keys()) + ["custom"],
                   help="下载预设（默认 all 全 5 类）")
    p.add_argument("--types", default="",
                   help="自定义类别（preset=custom 时生效），如 matrix,gpl，"
                        f"可选: {','.join(ALL_TYPES)}")
    p.add_argument("--out", default=".", help="输出根目录（默认当前目录）")
    p.add_argument("--proxy", default=DEFAULT_PROXY,
                   help=f"HTTP/HTTPS 代理（默认 {DEFAULT_PROXY}；空串关闭）")
    p.add_argument("--no-proxy", action="store_true",
                   help="禁用代理（等价 --proxy \"\"）")
    p.add_argument("--workers", type=int, default=8, help="并行线程数（默认 8）")
    p.add_argument("--timeout", type=int, default=30, help="单请求超时秒数（默认 30）")
    p.add_argument("--retries", type=int, default=3, help="失败重试次数（默认 3）")
    p.add_argument("--api-gap", type=float, default=0.4,
                   help="E-utilities 调用间隔秒（默认 0.4）")
    p.add_argument("--email", default="your.email@example.com",
                   help="NCBI E-utilities 邮箱（GPL 查询用）")
    p.add_argument("--suppl-block", default="filelist",
                   help="suppl 屏蔽关键字，逗号分隔（默认 filelist）")
    return p


def main() -> int:
    args = build_parser().parse_args()

    # 解析类别
    if args.preset == "custom":
        types = [t.strip().lower() for t in re.split(r"[,\s]+", args.types) if t.strip()]
        unknown = [t for t in types if t not in ALL_TYPES]
        if unknown:
            sys.exit(f"[ERROR] 未知类别: {unknown}；可选: {ALL_TYPES}")
        if not types:
            sys.exit(f"[ERROR] --preset custom 需配合 --types，可选: {ALL_TYPES}")
    else:
        types = PRESETS[args.preset]

    # 填充运行时配置
    CFG["gses"]        = _norm_gses(args.gse)
    CFG["email"]       = args.email
    CFG["workers"]     = max(1, args.workers)
    CFG["timeout"]     = args.timeout
    CFG["retries"]     = args.retries
    CFG["api_gap"]     = args.api_gap
    CFG["proxy"]       = "" if args.no_proxy else (args.proxy or "")
    CFG["types"]       = set(types)
    CFG["suppl_block"] = [b.strip() for b in args.suppl_block.split(",") if b.strip()]
    CFG["out"]         = args.out

    if not CFG["gses"]:
        sys.exit("[ERROR] --gse 解析后为空或非法（需形如 GSE12345）")

    Path(CFG["out"]).mkdir(parents=True, exist_ok=True)

    print(f"{'='*60}")
    print(f"  GEO 预下载器（CLI）")
    print(f"{'='*60}")
    print(f"  预设      : {args.preset}")
    print(f"  类别      : {', '.join(sorted(CFG['types']))}")
    print(f"  GSE 列表  : {', '.join(CFG['gses'])}")
    print(f"  输出目录  : {Path(CFG['out']).resolve()}")
    print(f"  代理      : {CFG['proxy'] or '(直连)'}")
    print(f"  并行线程  : {CFG['workers']}")
    print(f"  屏蔽关键字: {', '.join(CFG['suppl_block']) or '(无)'}")

    session = _session()
    summary = {}
    for g in CFG["gses"]:
        try:
            summary.update(_process_one(g, session))
        except Exception as e:
            print(f"\n[ERROR] {g} 处理异常: {e}（跳过；R 脚本仍可独立下载）")
            summary[g] = {"records": [], "gpls": [],
                          "ok": 0, "fail": 1, "bytes": 0}

    print(f"\n{'='*60}")
    print(f"  最终统计")
    print(f"{'='*60}")
    print(f"  {'GSE':<15s} {'成功':>5s}/{'失败':<5s}  {'数据量':>10s}  {'GPL'}")
    print(f"  {'-'*15} {'-'*11}  {'-'*10}  {'-'*20}")

    grand_ok = grand_fail = grand_bytes = 0
    for g, v in summary.items():
        ok = v["ok"]; fail = v["fail"]
        gpl_str = ", ".join(v["gpls"]) if v["gpls"] else "-"
        grand_ok += ok; grand_fail += fail; grand_bytes += v["bytes"]
        print(f"  {g:<15s} {ok:>4d}/{fail:<4d}  {_hsize(v['bytes']):>10s}  {gpl_str}")

    total = grand_ok + grand_fail
    print(f"\n  共 {len(summary)} 个 GSE | 总文件 {total} 个 | "
          f"成功 {grand_ok} | 失败 {grand_fail} | 总数据 {_hsize(grand_bytes)}")

    if grand_fail:
        print("\n[提示] 部分文件失败，R 脚本可独立下载这些文件，不影响后续流程。")
    return 0 if grand_fail == 0 else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\n[中断]")
        sys.exit(130)