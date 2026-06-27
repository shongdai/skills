#!/usr/bin/env python3
"""GEO_run.py — 嗅探→下载→检查 一键编排（输出精简，适配 AI token 预算）

设计要点（v1.2）：
1. 内置 UTF-8 日志双写（--log），子进程 stdout/stderr 实时 tee 到日志；
2. stdout 设为行缓冲（reconfigure(line_buffering=True)），AI 中途读 log 能看进度；
3. 各阶段打印明确标记：[START] [SNIFF_DONE] [GSE_BEGIN] [GSE_END] [CHECK_DONE] [ALL_DONE]
   便于 AI 通过关键字定位当前进度；
4. 文件末尾追加 ==== ALL_DONE exit=<code> ====，AI 据此判断"完成"。
"""

import argparse, io, json, os, re, subprocess, sys, time
from datetime import datetime
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = None       # 由 main() 从 --root 参数或 cwd 填充
LOG_FH = None     # 日志文件句柄（main 中打开）

# ============================================================
# 日志双写：所有 print 经 _log() 同时落地 stdout + log 文件
# ============================================================
def _log(msg=""):
    """同步打印到 stdout 和日志文件（若开启）。"""
    print(msg, flush=True)
    if LOG_FH is not None:
        LOG_FH.write(msg + "\n")
        LOG_FH.flush()

# ============================================================
# 子进程：流式读取 stdout（合并 stderr），逐行 tee 到日志
# ============================================================
def _stream(cmd, cwd=None, timeout=1800, prefix=""):
    """运行子进程并实时 tee 输出。返回 (returncode, full_output_str)。

    - 编码强制 UTF-8 + errors=replace，避免 Windows GBK 解码崩溃；
    - PYTHONIOENCODING=utf-8 让 Python 子进程也用 UTF-8 输出；
    - 每行立即 flush 到 log，AI 中途读文件能看到实时进度；
    - prefix 用于在日志中标识来源（如 "[sniff]"）。
    """
    env = os.environ.copy()
    env.setdefault("PYTHONIOENCODING", "utf-8")

    try:
        proc = subprocess.Popen(
            cmd, cwd=cwd or str(ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            encoding="utf-8", errors="replace",
            env=env, bufsize=1,  # 行缓冲
        )
    except Exception as e:
        _log(f"  [启动失败] {e}")
        return 1, str(e)

    lines = []
    start = time.time()
    try:
        for line in proc.stdout:
            line = line.rstrip("\r\n")
            lines.append(line)
            # 仅把 prefix 写入 log（控制台不打印每行，避免 token 爆炸）
            if LOG_FH is not None:
                LOG_FH.write(f"{prefix}{line}\n")
                LOG_FH.flush()
            # 超时检查
            if time.time() - start > timeout:
                proc.kill()
                _log(f"  [超时] {timeout}s")
                return 124, "\n".join(lines)
        proc.wait()
    except KeyboardInterrupt:
        proc.kill()
        raise

    return proc.returncode, "\n".join(lines)

# ============================================================
def _sniff(gses, proxy_url):
    """批量嗅探，返回 {gse: info}"""
    gse_str = ",".join(gses)
    cmd = [sys.executable, str(SCRIPTS / "GEO_smart_sniff.py"),
           "--gse", gse_str, "--json"]
    cmd += (["--proxy", proxy_url] if proxy_url else ["--no-proxy"])

    rc, out = _stream(cmd, timeout=120, prefix="[sniff] ")
    if rc != 0:
        _log(f"[sniff失败 rc={rc}] {out[-200:]}")
        return {}
    try:
        items = json.loads(out)
    except json.JSONDecodeError:
        m = re.search(r'\[.*\]', out, re.DOTALL)
        if not m:
            _log(f"[sniff解析失败] {out[-200:]}")
            return {}
        items = json.loads(m.group())
    return {it["gse_id"]: it for it in items if it.get("ok")}

def _check_done(gse, out_dir):
    """检查 expMatrix*.csv 是否已生成"""
    d = Path(out_dir) / gse
    return bool(list(d.glob("expMatrix*.csv")))

def _cel_pre_download(gse, out_dir, proxy_url):
    """CEL 路径：预下载 RAW 文件"""
    cmd = [sys.executable, str(SCRIPTS / "GEO_pre_download.py"),
           "--gse", gse, "--preset", "cel", "--out", out_dir]
    cmd += (["--proxy", proxy_url] if proxy_url else ["--no-proxy"])
    _log(f"  {gse:<12} [预下载] 开始")
    rc, _ = _stream(cmd, cwd=str(ROOT), timeout=3600, prefix=f"[{gse}|pre] ")
    _log(f"  {gse:<12} [预下载] {'OK' if rc == 0 else 'FAIL'}")
    return rc == 0

def _run_r(gse, info, out_dir, diff, timeout, proxy_url):
    """运行 R 脚本处理一个 GSE"""
    script = info.get("script", "")
    args   = list(info.get("args") or [])
    if not script:
        _log(f"  {gse:<12} [跳过] 无推荐脚本")
        return False

    extra = ["--out", out_dir]
    if diff and "--diff" not in args:
        extra += ["--diff", "TRUE"]
    extra += ["--timeout", str(timeout)]
    extra += ["--proxy", proxy_url if proxy_url else ""]

    r_cmd = ["Rscript", str(SCRIPTS / script)] + args + extra

    _log(f"  {gse:<12} [R运行] {script.replace('.R','')} ...")
    rc, out = _stream(r_cmd, cwd=str(ROOT),
                      timeout=timeout + 300, prefix=f"[{gse}|R] ")
    done = _check_done(gse, out_dir)
    if rc == 0 and done:
        _log(f"  {gse:<12} [R运行] OK")
        return True
    else:
        if rc != 0:
            tail = out[-500:] if out else "(无输出)"
            _log(f"  {gse:<12} [R运行] FAIL rc={rc}\n    {tail}")
        elif not done:
            _log(f"  {gse:<12} [R运行] FAIL (退出码 0 但无 expMatrix 产物)")
        return False

def _report_one(gse, out_dir):
    """单 GSE 完成报告"""
    d = Path(out_dir) / gse
    exps = sorted(d.glob("expMatrix*.csv"))
    if not exps:
        return f"{gse:<12} [无产物]"
    try:
        import pandas as pd
        ep = exps[0]
        df = pd.read_csv(ep, index_col=0).select_dtypes("number")
        rng = f"{df.min().min():.2f}~{df.max().max():.2f}"
        return f"{gse:<12} [OK] {df.shape[0]}×{df.shape[1]} {rng}"
    except Exception:
        return f"{gse:<12} [OK] {exps[0].name}"

# ============================================================
def main():
    p = argparse.ArgumentParser(
        description="GEO 一键编排：嗅探→下载→检查",
        epilog="日志：默认写 <out>/_GEO_run.log（UTF-8），可用 --log 自定义。"
    )
    p.add_argument("--gse", required=True)
    p.add_argument("--out", default=".")
    p.add_argument("--diff", action="store_true")
    p.add_argument("--timeout", type=int, default=1800)
    p.add_argument("--skip-check", action="store_true")
    p.add_argument("--root", default=None, help="项目根目录（默认当前工作目录）")
    p.add_argument("--proxy", default="http://127.0.0.1:7897",
                   help="HTTP 代理（默认 http://127.0.0.1:7897）")
    p.add_argument("--no-proxy", action="store_true", help="禁用代理")
    p.add_argument("--log", default=None,
                   help="日志路径（默认 <out>/_GEO_run.log，UTF-8 编码）")
    p.add_argument("--no-log", action="store_true",
                   help="禁用日志文件（仅控制台输出）")
    args = p.parse_args()

    global ROOT, LOG_FH
    ROOT = Path(args.root) if args.root else Path.cwd()
    proxy_url = "" if args.no_proxy else args.proxy

    gses = [g.strip().upper() for g in re.split(r"[, ]+", args.gse) if g.strip()]
    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)

    # —— 让 stdout 也行缓冲，便于实时 tail——
    try:
        sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
    except Exception:
        pass

    # —— 打开日志（UTF-8 强制）——
    if not args.no_log:
        log_path = args.log or os.path.join(out_dir, "_GEO_run.log")
        LOG_FH = open(log_path, "w", encoding="utf-8", buffering=1)  # 行缓冲
        _log(f"[START] {datetime.now():%Y-%m-%d %H:%M:%S}  pid={os.getpid()}")
        _log(f"[LOG] {log_path}")
        _log(f"[GSE] {','.join(gses)}  out={out_dir}  proxy={proxy_url or '(none)'}")

    rc_final = 0
    try:
        # 1. 嗅探
        _log("\n=== 嗅探 ===")
        info_map = _sniff(gses, proxy_url)
        _log(f"[SNIFF_DONE] {len(info_map)}/{len(gses)} ok")
        if not info_map:
            _log("[错误] 嗅探全部失败")
            rc_final = 1
            return rc_final

        # 2. 逐 GSE 下载
        _log("\n=== 下载 ===")
        results = {}
        for idx, gse in enumerate(gses, 1):
            _log(f"\n[GSE_BEGIN] {idx}/{len(gses)} {gse}")
            info = info_map.get(gse)
            if not info:
                _log(f"  {gse:<12} [跳过] 嗅探未返回")
                results[gse] = False
                _log(f"[GSE_END] {gse} FAIL")
                continue
            script = info.get("script", "")
            smp = info.get("n_samples", "?")
            plat = info.get("platform_id", "?")
            dtype = info.get("data_type", "?")
            _log(f"  {gse:<12} [>>>] {script.replace('.R','')}, {dtype}, {smp}smp, {plat}")

            if script == "GEO_download_cel.R":
                if not _cel_pre_download(gse, out_dir, proxy_url):
                    results[gse] = False
                    _log(f"[GSE_END] {gse} FAIL (pre-download)")
                    continue

            ok = _run_r(gse, info, out_dir, args.diff, args.timeout, proxy_url)
            results[gse] = ok
            _log(f"[GSE_END] {gse} {'OK' if ok else 'FAIL'}")

        # 3. 汇总
        _log("\n=== 结果 ===")
        for gse in gses:
            if results.get(gse):
                _log("  " + _report_one(gse, out_dir))
            else:
                _log(f"  {gse:<12} [FAIL]")
        ok_n = sum(1 for v in results.values() if v)
        fail_n = len(gses) - ok_n
        _log(f"\n==== {ok_n}/{len(gses)} done, {fail_n} fail ====")
        rc_final = 0 if fail_n == 0 else 1

        # 4. 检查（流式 tee 子进程输出）
        if not args.skip_check:
            _log("\n=== 检查 ===")
            check_cmd = [sys.executable,
                         str(SCRIPTS / "GEO_expMatrix_check.py"),
                         "--dir", out_dir,
                         "--gse", ",".join(gses),
                         "--no-save"]
            _stream(check_cmd, cwd=str(ROOT),
                    timeout=120, prefix="[check] ")
            _log("[CHECK_DONE]")

    except KeyboardInterrupt:
        _log("\n[中断]")
        rc_final = 130
    finally:
        _log(f"\n==== ALL_DONE exit={rc_final}  {datetime.now():%Y-%m-%d %H:%M:%S} ====")
        if LOG_FH is not None:
            LOG_FH.close()

    return rc_final

if __name__ == "__main__":
    sys.exit(main())
