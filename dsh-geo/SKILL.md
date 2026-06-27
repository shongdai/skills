---
name: "dsh-geo"
description: "GEO data skill. Query GSE metadata and download/preprocess expression matrices (probe/CEL/RNA-seq) via bundled CLI R scripts. Invoke when user mentions GSE/GEO or needs expression matrix download."
---

# DSH-GEO 技能（命令行调用版）

> **核心理念**：技能自带可执行 R 脚本（在 `scripts/` 目录下），通过命令行参数调用，**不需要修改源码**。

---

## 一句话定位

| 做什么 | 用什么 |
|---|---|
| **查询 / 检索 / 元数据** | Python（`Bio.Entrez` + `GEOparse`） |
| **下载 / 预处理 / ID映射 / log2** | 本技能自带 R 脚本（`scripts/*.R`） |

| 数据形态 | 调用脚本（CLI） |
|---|---|
| **一键自动（推荐）** | `python scripts/GEO_run.py --gse <GSE>`  ← 嗅探→下载→检查 全自动 |
| **不知道用哪个脚本？** | `python scripts/GEO_smart_sniff.py --gse <GSE>`  ← 先嗅探，5-7s |
| **想先把原始文件拉到本地**（离线/加速） | `python scripts/GEO_pre_download.py --gse <GSE> --preset <probe/cel/ncbicount/count/all>` |
| Series Matrix（已归一化矩阵）| `Rscript scripts/GEO_download_probe.R --gse <GSE>` |
| RAW（.CEL / .txt 原始文件）| `Rscript scripts/GEO_download_cel.R --gse <GSE>` |
| RNA-seq Count/FPKM/TPM（**NCBI 整理版**） | `Rscript scripts/GEO_download_ncbicount.R --gse <GSE> --species <human/mouse>` |
| RNA-seq Count（**作者上传 supplementary**） | `Rscript scripts/GEO_download_count.R --gse <GSE> --species <human/mouse> --from-type <ID>` |
| miRNA | _（未实现，见扩展计划）_ |

> ⚠️ Python 仅用于查询；下载与预处理一律走 `scripts/` 中的 R 脚本（命令行调用，无需改源码）。
>
> 🌐 **代理默认开启**：`http://127.0.0.1:7897`（4 个 R 脚本 + sniff 脚本）。不需代理时传 `--proxy ""` / `--no-proxy`。
>
> 🔍 **sniff 脚本启动时自动检测代理**：先 TCP 握手 + 1 次 HTTP 测活（~1.2s），可用就用，不可用自动回退直连并打印警告。脚本不卡死。

---

## 何时使用本技能

只要用户提到以下任意一项，**立即调用本技能**：

- 搜索 / 检索 GEO 数据集（不知道 GSE 号）
- 查看某个 GSE 的元数据 / 实验设计 / 样本信息
- 下载 GEO（GSE/GSM）表达数据
- 探针 ID → 基因 symbol
- 表达矩阵预处理 / log2 转换 / 基因去重
- 多平台 GSE 批量处理
- 临床表型 / 分组文件生成

---

## 决策树

```
用户给了 GSE 号了吗？
   │
   ├─ 没有 / 想"搜索 / 看元数据"
   │     └─→ 跳到「一、查询 GEO（Python）」
   │
   ├─ 有了，想"一键下载"（推荐）
   │     └─→ python scripts/GEO_run.py --gse <GSE> [--diff]    ← 全自动
   │           内部：嗅探 → (CEL预下载) → R脚本 → 检查报告
   │
   ├─ 有了，想"手动分步"（调试 / 单步控制）
   │     └─→ 看下方分步路径
   │
   │
   ├─ 有了，想"下载表达矩阵 + ID映射 + log2"（最常见）
   │     └─→ Rscript scripts/GEO_download_probe.R --gse <GSE_ID>
   │         → python scripts/GEO_expMatrix_check.py            [检查结果]
   │
   ├─ 有了，想"RAW / CEL / RMA 重新归一化"
   │     └─→ python scripts/GEO_pre_download.py --gse <GSE_ID> --preset cel   （必需！cel.R 默认不下载 suppl）
   │         → Rscript scripts/GEO_download_cel.R --gse <GSE_ID>
   │         → python scripts/GEO_expMatrix_check.py            [检查结果]
   │
   ├─ RNA-seq，NCBI 整理过的标准 raw_counts / FPKM / TPM
   │     └─→ Rscript scripts/GEO_download_ncbicount.R --gse <GSE> --species human/mouse
   │         → python scripts/GEO_expMatrix_check.py            [检查结果]
   │
   ├─ RNA-seq，NCBI 未整理，需要从作者上传的 supplementary 解析 *count* 文件
   │     └─→ Rscript scripts/GEO_download_count.R --gse <GSE> --species human/mouse --from-type ENSEMBL
   │         → python scripts/GEO_expMatrix_check.py            [检查结果]
   │
   └─ miRNA？
         └─ TODO：等待 scripts/GEO_download_mirna.R 加入
```

### 0. 自动嗅探（推荐先用）

```powershell
# 单个 GSE
python scripts/GEO_smart_sniff.py --gse GSE266899

# 多个 GSE 批量
python scripts/GEO_smart_sniff.py --gse GSE5281,GSE266899,GSE190451

# JSON 输出（机器可读，便于集成到 workflow）
python scripts/GEO_smart_sniff.py --gse GSE266899 --json
```

**嗅探流程**（5-7s/GSE，纯元数据，不下载大文件）：

| 步骤 | 调用 | 用途 |
|---|---|---|
| 1 | `Entrez.esearch(db="gds", term="GSE[Accession]")` | GSE → GDS UID |
| 2 | `Entrez.esummary(db="gds", id=UID)` | 拿 platform / taxon / n_samples / gdsType |
| 3 | FTP `suppl/` 目录 HTML 解析 | 是否有作者 count / CEL 文件 |
| 4 | HEAD `ncbi.nlm.nih.gov/geo/download/...` | 是否有 NCBI 标准化 raw_counts |
| 5 | `Entrez.efetch` + FTP `matrix/*.gz` 第一列 | 推断 ID 类型（ENSEMBL/ENTREZID/SYMBOL） |

**输出**：

```
GSE        : GSE266899    [rna_seq]
平台       : GPL13112    样本数: 18
物种       : mouse    -> OrgDb = org.Mm.eg.db
ID 类型    : ENSEMBL
特征       : NCBI_norm=False, author_cnt=True, CEL=False
[推荐脚本]  GEO_download_count.R
[推荐命令]  Rscript scripts/GEO_download_count.R --gse GSE266899 --diff TRUE --species mouse --from-type ENSEMBL
[耗时]    6.69s
```

**依赖**：`uv pip install biopython`（不需要 GEOparse）

---

### 0.5 预下载（仅 CEL 路径必需，其他可选加速）

`scripts/GEO_pre_download.py` 使用 Python `requests` 多线程下载，**支持 HTTP Range 断点续传**，比 R 脚本内置的 `getGEOSuppFiles`（不支持续传）更稳定。

| 路径 | 是否必需预下载 | 原因 |
|---|---|---|
| **CEL** (`GEO_download_cel.R`) | **必需** | RAW.tar 文件大（数百 MB），`cel.R` 默认不调用 `getGEOSuppFiles` |
| probe / ncbicount / count | 可选 | 数据量小，R 脚本内置下载即可 |

```powershell
# 配合 R CEL 路径（cel.R）：matrix + suppl + gpl  ← 必需！
python scripts/GEO_pre_download.py --gse GSE5281 --preset cel

#### 预设 ↔ R 脚本对应关系

| 预设 `--preset` | 下载类别 | 后续 R 脚本 | 用途 |
|---|---|---|---|
| `all`（默认） | soft + miniml + matrix + suppl + gpl | 任意 | 全量归档备份 |
| `probe` | matrix + gpl | `GEO_download_probe.R` | 已归一化芯片矩阵 + 本地 GPL 注释回退 |
| `cel` | matrix + suppl + gpl | `GEO_download_cel.R` | CEL 原始芯片 + RMA |
| `ncbicount` | matrix | `GEO_download_ncbicount.R` | RNA-seq NCBI 整理矩阵（raw_counts 由 R 单独下载） |
| `count` | matrix + suppl | `GEO_download_count.R` | RNA-seq 作者上传 count 文件 |
| `custom` | 由 `--types` 指定 | — | 任意组合（可选 `soft,miniml,matrix,suppl,gpl`）|

#### 主要参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--gse <ids>` | **必填** | GSE 号，逗号/空格分隔 |
| `--preset <name>` | `all` | `all` / `probe` / `cel` / `ncbicount` / `count` / `custom` |
| `--types <a,b>` | — | `custom` 时指定（`soft,miniml,matrix,suppl,gpl` 任选） |
| `--out <dir>` | `.` | 输出根目录（每个 GSE 一个子目录） |
| `--proxy <url>` | `http://127.0.0.1:7897` | HTTP 代理；`--no-proxy` 禁用 |
| `--workers <N>` | `8` | 并行线程数 |
| `--timeout <sec>` | `30` | 单请求超时 |
| `--retries <N>` | `3` | 失败重试 |
| `--suppl-block <a,b>` | `filelist` | 屏蔽 suppl 文件名包含的关键字 |

> 完整流程见 [§2.0 GEO_run.py](#20-georunpy-—-一键自动编排推荐)。CEL 路径需额外预下载，详见 [§2.4](#24-geodownload_celr-cel-raw-路径)。

#### R 脚本运行耗时参考

R 脚本执行分为两个阶段：**网络下载**（依赖代理和带宽）和**本地计算**（依赖 CPU/内存）。建议调用时设置充足的超时。

| 脚本 | 典型耗时 | 瓶颈 |
|---|---|---|
| `GEO_download_probe.R` | 1-5 分钟 | GPL 注释下载 + ID 映射（大 GPL 如 GPL570 > 10 万行注释） |
| `GEO_download_cel.R` | 5-30 分钟 | RAW.tar 下载 + RMA 归一化（样本数 × CEL 大小） |
| `GEO_download_ncbicount.R` | 1-5 分钟 | NCBI count 矩阵下载 + 基因注释 |
| `GEO_download_count.R` | 2-10 分钟 | suppl 文件下载 + 解压 + 识别 count 文件 |

> 建议命令调用时带 `--timeout 1800`（30 分钟），并确保运行环境不要提前终止子进程。多 GSE 批量时总耗时 = Σ(各 GSE 耗时)。

---

## 一、查询 GEO（Python）

> 本章节代码**仅用于查询、检索、查看元数据**，**不要用于下载数据**。

### 0. 准备

```python
from Bio import Entrez
Entrez.email = "your.email@example.com"
# 可选：Entrez.api_key = "<your_key>"  无 key 3 req/s，有 key 10 req/s
```

### 1.1 搜索 GEO DataSets

```python
def search_geo_datasets(query, retmax=20):
    handle = Entrez.esearch(db="gds", term=query, retmax=retmax, usehistory="y")
    results = Entrez.read(handle); handle.close()
    return results

# 例
search_geo_datasets("breast cancer[MeSH] AND Homo sapiens[Organism]")
search_geo_datasets("GPL570[Accession]")
search_geo_datasets("RNA-seq[DataSet Type] AND 2024[Publication Date]")
```

> 💡 **仅查询**。要下载请用本技能 R 脚本。

### 1.2 搜索 GEO Profiles（基因维度）

```python
def search_geo_profiles(gene, organism="Homo sapiens", retmax=100):
    h = Entrez.esearch(db="geoprofiles",
                       term=f"{gene}[Gene Name] AND {organism}[Organism]",
                       retmax=retmax)
    r = Entrez.read(h); h.close(); return r

search_geo_profiles("TP53")
```

> 💡 **仅查询**。

### 1.3 高级搜索（多条件组合）

```python
def advanced_geo_search(terms, op="AND"):
    return search_geo_datasets(f" {op} ".join(terms))

advanced_geo_search([
    "RNA-seq[DataSet Type]",
    "Homo sapiens[Organism]",
    "2024[Publication Date]",
])
```

### 1.4 查看 GSE 元数据

**方式 A：esummary（轻量）**

```python
def fetch_geo_summary(gse_id):
    sh = Entrez.esearch(db="gds", term=f"{gse_id}[Accession]")
    sr = Entrez.read(sh); sh.close()
    if not sr['IdList']: return None
    uh = Entrez.esummary(db="gds", id=sr['IdList'][0])
    s = Entrez.read(uh); uh.close(); return s[0]

info = fetch_geo_summary("GSE76262")
print(info.get('title'), info.get('n_samples'), info.get('GPL'))
```

**方式 B：GEOparse（完整元数据）**

```python
import GEOparse
gse = GEOparse.get_GEO(geo="GSE76262", destdir="./_meta_cache")
print(gse.metadata['title'], gse.metadata['summary'])
```

> 💡 仅用于查询；要下载分析，调用 R 脚本。

### 1.5 NCBI 速率限制

| 模式 | 速率 | sleep |
|---|---|---|
| 无 key | 3 req/s | `time.sleep(0.34)` |
| 有 key | 10 req/s | `time.sleep(0.1)` |

---

## 二、下载与预处理（R CLI 脚本）

### 2.0 GEO_run.py — 一键自动编排（推荐）

将嗅探→下载→检查串成一条命令，**输出精简**，适合批量处理。

```powershell
# 单 GSE
python scripts/GEO_run.py --gse GSE5281

# 多 GSE + 生成 diff 模板
python scripts/GEO_run.py --gse GSE5281,GSE266899 --diff

# 指定输出目录 + 自定义超时
python scripts/GEO_run.py --gse GSE5281 --out ./data --timeout 3600

# 跳过最终检查（只下载）
python scripts/GEO_run.py --gse GSE5281 --skip-check

# 自定义日志文件 / 禁用日志
python scripts/GEO_run.py --gse GSE5281 --log ./run.log
python scripts/GEO_run.py --gse GSE5281 --no-log
```

内部自动判断：
- CEL 路径 → 先 `pre_download --preset cel` → 再 `cel.R`
- 其他路径 → 直接调对应 R 脚本
- 全部完成后 → 自动跑 `GEO_expMatrix_check.py`

#### 🤖 AI 调用范式（重要）

`GEO_run.py` v1.2 起内置 **UTF-8 双写日志 + 阶段标记**，AI **不再依赖终端缓冲**，应通过日志文件判断进度：

| 项 | 约定 |
|---|---|
| 日志路径 | 默认 `<out>/_GEO_run.log`（UTF-8） |
| 编码 | 强制 UTF-8（避免 Windows GBK 乱码） |
| 实时性 | 行缓冲，可边跑边读 |
| 阶段标记 | `[START]` `[SNIFF_DONE]` `[GSE_BEGIN]` `[GSE_END]` `[CHECK_DONE]` `[ALL_DONE exit=N]` |

**推荐 AI 工作流**：

1. **启动**：用非阻塞方式跑 `GEO_run.py`，不要重定向（脚本自带日志双写）：
   ```powershell
   python scripts/GEO_run.py --gse GSE5281,GSE266899 --diff
   ```
2. **轮询日志**：通过 Read 工具读 `<out>/_GEO_run.log` 检查进度。
   - 找到 `==== ALL_DONE exit=0 ====` 即视为成功结束
   - 找到 `==== ALL_DONE exit=1 ====` 视为有失败 → 看 `[GSE_END] <GSE> FAIL` 上方的 `[R运行] FAIL rc=...` 块
3. **建议的检查时间节点**（按 GSE 数量调整）：

   | 步骤 | 等待时长 | 期望日志关键字 |
   |---|---|---|
   | 启动后 15 秒 | 短 | `[SNIFF_DONE] N/N ok` |
   | + 30 秒（每 GSE） | 中 | `[GSE_BEGIN] 1/N <GSE>` |
   | + 2-10 分钟（每 GSE，依数据量） | 长 | `[GSE_END] <GSE> OK/FAIL` |
   | 最终 | — | `[CHECK_DONE]` + `==== ALL_DONE exit=N ====` |

   未见 `ALL_DONE` 时勿提前判定失败 —— 继续等待或读最新一行 prefix（如 `[GSE266899|R]`）确认仍在运行。

4. **乱码处理**：日志由 Python 直接写 UTF-8，**不要**用 PowerShell `>` 重定向（会变 GBK）。

---

> 一律通过命令行调用，**绝不修改源码**。

### 2.1 GEO_download_probe.R（Series Matrix 路径，最常用）

```powershell
# 单个 GSE
Rscript scripts/GEO_download_probe.R --gse GSE76262

# 多个 GSE（逗号分隔，无空格）
Rscript scripts/GEO_download_probe.R --gse GSE76262,GSE12345

# 自定义输出目录 + 生成分组文件
Rscript scripts/GEO_download_probe.R --gse GSE76262 --out ./data --diff TRUE

# 启用 idmap() 在线注释（默认 FALSE，走本地 GPL）
Rscript scripts/GEO_download_probe.R --gse GSE76262 --use-idmap TRUE

# 本地 GPL 列名不同
Rscript scripts/GEO_download_probe.R --gse GSE76262 --probe-col probe_id --gene-col "Gene Symbol"

# 通过代理下载
Rscript scripts/GEO_download_probe.R --gse GSE76262 --proxy http://127.0.0.1:7897

# 查看帮助
Rscript scripts/GEO_download_probe.R --help
```

**参数表**：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--gse <ids>` | **必填** | GSE 号，多个用逗号分隔 |
| `--out <dir>` | `.` | 输出根目录 |
| `--use-idmap <T/F>` | `FALSE` | 是否调用 AnnoProbe::idmap() |
| `--probe-col <name>` | `ID` | 本地 GPL 探针 ID 列名（模糊匹配）|
| `--gene-col <name>` | `symbol` | 本地 GPL 基因列名（模糊匹配）|
| `--diff <T/F>` | `FALSE` | 生成空 group.txt / diff.txt |
| `--use-bitr <T/F>` | `FALSE` | 当 GPL 注释列不是 SYMBOL 时，调用 `clusterProfiler::bitr` 二次转换 |
| `--from-type <type>` | `ENTREZID` | bitr 输入 ID 类型（ENTREZID/REFSEQ/ENSEMBL/SYMBOL/...）；仅 `--use-bitr TRUE` 生效 |
| `--org-db <pkg>` | `org.Hs.eg.db` | 物种注释包；仅 `--use-bitr TRUE` 生效 |
| `--proxy <url>` | `http://127.0.0.1:7897` | HTTP 代理；不需代理时传空串 `--proxy ""` |
| `--timeout <sec>` | `600` | 下载超时 |
| `--help` | | 显示帮助 |

#### 🩹 当默认参数失败（GPL 没有 `Symbol/symbol` 列）

脚本会打印「列匹配警告」。打开 `{out}/{GSE}/GPLxxxx.txt` 看表头，按下表挑列：

| GPL 实际列名 | 内容样例 | 应加的参数 |
|---|---|---|
| `Gene Symbol` / `symbol` | `TP53` | （默认即可，无需 bitr） |
| `GB_ACC` | `NM_014332.1`（RefSeq） | `--use-bitr TRUE --gene-col GB_ACC --from-type REFSEQ` |
| `ENTREZ_GENE_ID` | `780`（Entrez 数字） | `--use-bitr TRUE --gene-col ENTREZ_GENE_ID --from-type ENTREZID` |
| `SPOT_ID`（纯数字） | `1`, `780`, `5982` | `--use-bitr TRUE --gene-col SPOT_ID --from-type ENTREZID` |
| `ENSEMBL` / 含 `ENSG`/`ENSMUSG` | `ENSG00000141510` | `--use-bitr TRUE --gene-col <列名> --from-type ENSEMBL` |
| 含 `UNIPROT` | `P04637` | `--use-bitr TRUE --gene-col <列名> --from-type UNIPROT` |

> 小鼠数据加 `--org-db org.Mm.eg.db`，大鼠加 `--org-db org.Rn.eg.db`。

**实战例子**：

```powershell
# GSE8479 (GPL2700)：只有 GB_ACC 列（RefSeq）
Rscript scripts/GEO_download_probe.R --gse GSE8479 --use-bitr TRUE --gene-col GB_ACC --from-type REFSEQ

# GSE117525 (GPL20880)：只有 SPOT_ID 列，内容是 Entrez 数字
Rscript scripts/GEO_download_probe.R --gse GSE117525 --use-bitr TRUE --gene-col SPOT_ID --from-type ENTREZID

# GSE14520 (GPL571)：有 ENTREZ_GENE_ID 列
Rscript scripts/GEO_download_probe.R --gse GSE14520 --use-bitr TRUE --gene-col ENTREZ_GENE_ID --from-type ENTREZID
```

### 2.2 GEO_download_ncbicount.R（NCBI 标准化 RNA-seq counts，**推荐用于 RNA-seq**）

NCBI 自 2023 年起为部分 GEO RNA-seq 数据集提供官方整理的 **raw_counts / FPKM / TPM** 三种矩阵，并配套基因注释文件，比作者上传的原始 supplementary 文件更规范。

```powershell
# 单 GSE（默认 human）
Rscript scripts/GEO_download_ncbicount.R --gse GSE56545

# 多 GSE 批量
Rscript scripts/GEO_download_ncbicount.R --gse GSE56545,GSE70089,GSE33294 --diff TRUE

# 小鼠数据（自动切换 GRCm39 + Mouse 注释文件 + org.Mm.eg.db）
Rscript scripts/GEO_download_ncbicount.R --gse GSE100001 --species mouse

# 自定义输出 + 代理
Rscript scripts/GEO_download_ncbicount.R --gse GSE56545 --out ./data --proxy http://127.0.0.1:7897
```

**参数表**：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--gse <ids>` | **必填** | GSE 号，多个用逗号分隔 |
| `--species <name>` | `human` | `human` / `mouse`（**自动切换 genome_tag + annot 文件 + OrgDb**）|
| `--out <dir>` | `.` | 输出根目录 |
| `--diff <T/F>` | `FALSE` | 生成 group/diff 模板 |
| `--proxy <url>` | `http://127.0.0.1:7897` | HTTP 代理；不需代理时传空串 `--proxy ""` |
| `--timeout <sec>` | `600` | 下载超时 |
| `--help` | | 显示帮助 |

**物种自动切换的映射**：

| species | genome_tag | annot 文件 | OrgDb |
|---|---|---|---|
| human | `GRCh38.p13` | `Human.GRCh38.p13.annot.tsv.gz` | `org.Hs.eg.db` |
| mouse | `GRCm39` | `Mouse.GRCm39.annot.tsv.gz` | `org.Mm.eg.db` |

### 2.3 GEO_download_count.R（作者上传的 supplementary count 文件）

当 NCBI 没整理 RNA-seq counts 时（`ncbicount` 脚本下载 404），从 FTP `suppl/` 目录拉所有附件，自动解压并识别 `*count*` 文件。

```powershell
# 单 GSE，默认 human + ENSEMBL ID
Rscript scripts/GEO_download_count.R --gse GSE152418

# 小鼠 + ENSEMBL
Rscript scripts/GEO_download_count.R --gse GSE99999 --species mouse --from-type ENSEMBL

# 输入是 SYMBOL 或 ENTREZID
Rscript scripts/GEO_download_count.R --gse GSE99999 --from-type ENTREZID

# 自定义物种（例如大鼠）
Rscript scripts/GEO_download_count.R --gse GSE99999 --species custom --org-db org.Rn.eg.db

# 文件名关键字不是 count 而是 expression
Rscript scripts/GEO_download_count.R --gse GSE99999 --keyword expression
```

**参数表**：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--gse <ids>` | **必填** | GSE 号，多个用逗号分隔 |
| `--species <name>` | `human` | `human` / `mouse` / `custom`（custom 时必须配 `--org-db`）|
| `--org-db <pkg>` | `(空)` | OrgDb 包名（custom 时必填，如 `org.Rn.eg.db`）|
| `--from-type <id>` | `ENSEMBL` | 输入基因 ID 类型（ENSEMBL/ENTREZID/SYMBOL/...）|
| `--keyword <str>` | `count` | 文件名关键字（不区分大小写）|
| `--out <dir>` | `.` | 输出根目录 |
| `--diff <T/F>` | `FALSE` | 生成 group/diff 模板 |
| `--proxy <url>` | `http://127.0.0.1:7897` | HTTP 代理；不需代理时传空串 `--proxy ""` |
| `--timeout <sec>` | `600` | 下载超时 |
| `--help` | | 显示帮助 |

### 2.4 GEO_download_cel.R（CEL RAW 路径）

> 必须先执行预下载！本脚本默认不调用 `getGEOSuppFiles`（避免下载中断），改用 `GEO_pre_download.py` 断点续传下载 RAW 文件。

```powershell
# 推荐流程
python scripts/GEO_pre_download.py --gse GSE5281 --preset cel     # 1. 先预下载
Rscript scripts/GEO_download_cel.R --gse GSE5281                  # 2. 再跑 R

# 完整选项
Rscript scripts/GEO_download_cel.R --gse GSE5281 --out ./data --diff TRUE --use-idmap TRUE --proxy http://127.0.0.1:7897

# 如果一定要在线下载（不推荐，不支持断点续传）
Rscript scripts/GEO_download_cel.R --gse GSE5281 --download-suppl TRUE

# 查看帮助
Rscript scripts/GEO_download_cel.R --help
```

**参数表**（与 probe 版基本一致；CEL 默认 `--gene-col gene_assignment`，更贴合 Affymetrix 平台）。

### 2.5 输出文件

每个 GSE 在 `{out}/{gse_id}/` 下生成：

| 文件 | 内容 |
|---|---|
| `clinical{_GPL}.csv` | 临床表型（pData）|
| `expMatrix{_GPL}.csv` | 处理后的表达矩阵（symbol × sample）|
| `probe2gene{_GPL}.csv` | 探针 → 基因 symbol 映射 |
| `group{_GPL}.txt` | 分组文件（`--diff TRUE` 时）|
| `diff{_GPL}.txt` | 差异比较定义（`--diff TRUE` 时）|

根目录还生成 `download_status.txt`（批量状态汇总）。

### 2.6 关键设计

- **断点续跑**：检测到 `expMatrix.csv` 已存在自动跳过
- **健壮性**：单 GSE / 单平台失败不影响其它
- **多平台**：单平台输出无后缀；多平台自动加 `_GPL570` 等
- **idmap 回退**：pipe → soft → bioc → 本地 GPL → NCBI 自动下载

### 2.7 下载结果检查

`scripts/GEO_expMatrix_check.py` 扫描输出目录下所有 GSE 子文件夹，自动识别下载类型并检查数据质量。

```powershell
# 扫描当前目录
python scripts/GEO_expMatrix_check.py

# 指定目录 + 自定义输出路径
python scripts/GEO_expMatrix_check.py --dir ./data --output ./check.txt

# 仅打印到控制台
python scripts/GEO_expMatrix_check.py --dir ./data --no-save
```

**参数**：

| 参数 | 默认 | 说明 |
|---|---|---|
| `--dir <path>` | `.` | 扫描的根目录 |
| `--output <path>` | `<dir>/check_result.txt` | 结果报告输出路径 |
| `--no-save` | — | 仅控制台打印，不写文件 |

**检查项**：

| 类别 | 检查内容 |
|---|---|
| 文件齐全度 | clinical / expMatrix / probe2gene / group / diff 是否存在 |
| 下载类型自动识别 | 根据产物推断 probe / cel / ncbicount / count / raw / unknown |
| 数据质量 | 矩阵形状、数值范围、是否含 NA、是否含负值 |
| 多平台检测 | series_matrix 数量 vs expMatrix 数量 |
| 未处理完 | series_matrix 数 > expMatrix 数 |

**输出示例**：

```
【下载类型分布】
  probe (芯片)               3
  ncbicount (NCBI RNA-seq)   2

【存在 expMatrix*.csv 的文件夹】(5 个，共 7 个矩阵)
文件夹         下载类型     矩阵类型   基因×样本         NA数量   负值数量   范围
GSE12345       probe        Single     (20531, 12)      0        0          2.10-14.55
GSE56545       ncbicount    Count      (18472, 8)       0        0          1.00-45678.00

汇总统计:
  总文件夹数:              5
  有 expMatrix*.csv:      5 (共 7 个矩阵)
  含 NA 的矩阵:             0
  含负值的矩阵:             0
```

### 2.8 下游对接

```r
exp     <- fread("data/GSE76262/expMatrix.csv") |> column_to_rownames("symbol")
clin    <- fread("data/GSE76262/clinical.csv")
group   <- fread("data/GSE76262/group.txt")
diff_df <- fread("data/GSE76262/diff.txt")
```

---

## 三、边界分工（Python vs R 脚本）

| 任务 | 用什么 | 命令 |
|---|---|---|
| 按关键词搜索 GEO | **Python** | `Entrez.esearch(db="gds", ...)` |
| 按基因找跨研究 Profiles | **Python** | `Entrez.esearch(db="geoprofiles", ...)` |
| 看 GSE 标题/设计/样本数 | **Python** | `Entrez.esummary` |
| 批量元数据汇总成表 | **Python** | `Entrez.esummary` + `pandas.DataFrame` |
| **下载表达矩阵（芯片）** | **R CLI** | `Rscript scripts/GEO_download_probe.R --gse ...` |
| **下载 RAW + RMA（芯片）** | **R CLI** | `python scripts/GEO_pre_download.py ... --preset cel` → `Rscript scripts/GEO_download_cel.R --gse ...` |
| **下载 NCBI RNA-seq counts** | **R CLI** | `Rscript scripts/GEO_download_ncbicount.R --gse ... --species human/mouse` |
| **下载作者 supplementary counts** | **R CLI** | `Rscript scripts/GEO_download_count.R --gse ... --species human/mouse` |
| **下载结果质量检查** | **Python CLI** | `python scripts/GEO_expMatrix_check.py --dir <dir>` |
| 物种切换（人/小鼠） | **R CLI** | `--species human` / `--species mouse` |
| 探针 → 基因 symbol | **R CLI** | 上述脚本内置 |
| log2 / 去重 / 过滤 | **R CLI** | 上述脚本内置 |
| 多 GSE 批量下载 | **R CLI** | `--gse GSE1,GSE2,GSE3` |
| 多平台 GSE 处理 | **R CLI** | 自动 `_GPLxxx` 后缀 |
| 临床表型 + 分组模板 | **R CLI** | `--diff TRUE` |
| 下游差异分析 | **R 项目脚本** | 兄弟目录 `DEG_limma.R` 等 |

---

## 四、典型组合工作流

### 场景：从 0 开始

**步骤 1：Python 检索 GSE**

```python
from Bio import Entrez
Entrez.email = "your.email@example.com"
h = Entrez.esearch(db="gds",
    term="breast cancer[MeSH] AND Homo sapiens[Organism] AND GPL570[Accession]",
    retmax=10)
ids = Entrez.read(h)['IdList']; h.close()
```

**步骤 2：Python esummary 看元数据**

```python
for uid in ids[:5]:
    sh = Entrez.esummary(db="gds", id=uid)
    s = Entrez.read(sh)[0]; sh.close()
    print(s.get('Accession'), s.get('n_samples'), s.get('title')[:80])
```

**步骤 3：直接调用 R 脚本下载**

```powershell
Rscript scripts/GEO_download_probe.R --gse GSE76262,GSE12345 --out ./data --diff TRUE
```

**步骤 4：检查结果**

```powershell
python scripts/GEO_expMatrix_check.py --dir ./data
```

### 场景：已知 GSE 号

跳过步骤 1-2，一键下载见 [§2.0](#20-georunpy-—-一键自动编排推荐)。或手动分步：
```powershell
python scripts/GEO_smart_sniff.py --gse GSE76262                        # 1. 确认类型
Rscript scripts/<推荐脚本> --gse GSE76262 [选项]                        # 2. R 处理
python scripts/GEO_expMatrix_check.py                                   # 3. 检查结果
```

### 场景：只看元数据

只跑步骤 1-2 Python 部分，**不调用 R 脚本**。

---

## 依赖

### R 依赖（下载与预处理）

通过技能内置的 [scripts/install_dependencies.R](./scripts/install_dependencies.R) 安装：

| 类型 | 包 |
|---|---|
| CRAN | tidyverse, data.table, archive, devtools, remotes |
| Bioconductor | limma, affy, GEOquery, clusterProfiler, org.Hs.eg.db / org.Mm.eg.db |
| GitHub | jmzeng1314/AnnoProbe |

首次运行前（CLI）：
```powershell
# 装人+小鼠（默认）
Rscript scripts/install_dependencies.R

# 只装人
Rscript scripts/install_dependencies.R --species human

# 只装小鼠
Rscript scripts/install_dependencies.R --species mouse

# 额外装大鼠 / 斑马鱼
Rscript scripts/install_dependencies.R --species all --org-db org.Rn.eg.db,org.Dr.eg.db
```

### Python 依赖（查询 + 预下载）

```bash
# 一次性安装全部
uv pip install biopython GEOparse pandas requests
```

| 包 | 用途 | 脚本 |
|---|---|---|
| `biopython` | NCBI E-utilities 查询 | `GEO_smart_sniff.py` / 1.x 节示例 |
| `requests` | HTTP 下载（带重试 / 代理） | `GEO_pre_download.py` |
| `GEOparse` | 完整 SOFT 元数据解析（可选） | 1.4 节方式 B |
| `pandas` | 元数据汇总表格（可选） | 三、边界分工表 |

---

## 扩展计划

后续新脚本按相同 CLI 风格放在 `scripts/` 目录，并同步更新本 SKILL.md：

- [x] `scripts/GEO_download_ncbicount.R` — RNA-seq Count/FPKM/TPM（NCBI 标准化）✅
- [x] `scripts/GEO_download_count.R` — RNA-seq Count（作者 supplementary）✅
- [ ] `scripts/GEO_download_mirna.R` — miRNA 表达
- [ ] `scripts/GEO_download_methylation.R` — DNA 甲基化
- [ ] `scripts/GEO_download_singlecell.R` — 单细胞表达
- [ ] `scripts/GEO_merge_combat.R` — 多 GSE ComBat 批次校正合并

---

## 版本表

| 组件 | 版本 | 备注 |
|---|---|---|
| `SKILL.md`                | 3.6 | 新增 AI 调用范式 + 日志检查时间节点 |
| `GEO_run.py`              | 1.2 | 内置 UTF-8 双写日志 + 流式 tee + 阶段标记 |
| `GEO_smart_sniff.py`      | 2.1 | ncbicount 路径仅 human，不再传 --species |
| `GEO_pre_download.py`     | 1.0 | 5 类预设 + 自动代理 + 断点续传 |
| `GEO_download_probe.R`    | 3.0 | bitr 二次转换 |
| `GEO_download_cel.R`      | 2.1 | 默认跳过 getGEOSuppFiles，需预下载 |
| `GEO_download_ncbicount.R`| 2.1 | 仅 human（移除 mouse）|
| `GEO_download_count.R`    | 2.0 | 物种切换 + custom OrgDb |
| `GEO_expMatrix_check.py`  | 1.0 | 下载结果质量检查 |
| `install_dependencies.R`  | 2.0 | 按物种条件安装 |

---

## FAQ

**Q1: 下载很慢或超时**
A: 加 `--proxy http://127.0.0.1:7897`；或加大 `--timeout 1800`。

**Q2: idmap() 全部失败**
A: 脚本自动回退到本地 GPL → NCBI 下载；NCBI 也连不上时手动放 `GPL.txt` 到 `{out}/{gse_id}/`。

**Q3: 列名找不到（"列匹配警告"）/ probe2gene 全空 / 基因符号大量为 NA**
A: 这是 GPL 没有 `Symbol/symbol` 列，需要走 bitr 二次转换。
  1. 打开 `{out}/{GSE}/GPLxxxx.txt`，找含基因 ID 的列（GB_ACC / ENTREZ_GENE_ID / SPOT_ID / ENSEMBL 等）。
  2. 重跑脚本时加 `--use-bitr TRUE --gene-col <列名> --from-type <ID类型>`。
  3. 详细对照表见 [2.1 节 🩹 故障恢复速查](#-当默认参数失败gpl-没有-symbolsymbol-列)。

**Q4: 怎么从 Python 检索结果无缝传给 R 脚本？**
A:
```python
gse_ids = ["GSE76262", "GSE12345", "GSE99999"]
import subprocess
subprocess.run(["Rscript", "scripts/GEO_download_probe.R",
                "--gse", ",".join(gse_ids), "--out", "./data"])
```

**Q5: 能否绕过 R 脚本用 Python 下载？**
A: 不推荐。ID 映射 / log2 / 过滤都在 R 里。预下载仅 CEL 路径必需（`GEO_pre_download.py --preset cel`），其他路径数据量小，R 脚本内置下载足够。

**Q6: 出错了如何调试？**
A: 加 `--help` 查看用法；脚本输出含详细 `cat()` 日志便于定位。

**Q7: R 脚本要跑多久？会不会卡死？**
A: 见 [§0.5 R 脚本运行耗时参考](#r-脚本运行耗时参考)。建议加 `--timeout 1800`（30 分钟上限）。如果卡在"加载依赖"阶段超过 2 分钟，通常是 Bioconductor 包未安装，先跑 `Rscript scripts/install_dependencies.R`。

---

## 参考资料

- [GEO 官网](https://www.ncbi.nlm.nih.gov/geo/)
- [GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)
- [AnnoProbe](https://github.com/jmzeng1314/AnnoProbe)
- [GEOparse](https://geoparse.readthedocs.io/)
- [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [references/geo_reference.md](./references/geo_reference.md) — GEO 数据格式、FTP 路径、SOFT/MINiML 解析

---

## 作者

科研木鱼（闲鱼/小红书：科研木鱼）
