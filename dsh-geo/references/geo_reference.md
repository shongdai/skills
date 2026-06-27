# GEO 深度知识参考

> **定位**：本文档专注于 GEO 数据库的**底层格式 / 命名规则 / 平台特性 / ID 体系**等"深度知识"。
> 日常操作命令、CLI 用法、决策流程见 [../SKILL.md](../SKILL.md)；本文不重复操作示例。

---

## 1. GEO 五大记录类型

| 类型 | 前缀 | 含义 | 一句话 |
|---|---|---|---|
| Platform | `GPL` | 平台（探针定义/参考基因组） | "用什么仪器/什么探针测的" |
| Sample | `GSM` | 单个样本 | "一份原始数据" |
| Series | `GSE` | 实验集合（n × GSM） | "一篇文章的所有样本" |
| DataSet | `GDS` | 编辑过的可比较 GSE 子集 | NCBI 人工整理（很少新增） |
| Profile | — | 单基因在多 Series 的表达模式 | 仅在 `geoprofiles` 库可查 |

> **重要差异**：
> - 同一个 GSE 可包含**多个 GPL**（混合平台实验），脚本要按 GPL 拆分输出（本技能用 `_GPL570` 后缀处理）。
> - `gds` 数据库同时包含 GSE 和 GDS，esearch 时 `[Accession]` 字段会自动消歧。

---

## 2. GEO 数据形态决策矩阵

| 形态 | Series Matrix 含表达值？ | suppl/ 含什么？ | NCBI raw_counts？ | 推荐脚本 |
|---|---|---|---|---|
| **已归一化芯片** | ✅（log2 完成） | 一般空或 RAW.tar | — | `GEO_download_probe.R` |
| **未处理芯片** | ❌（占位/空） | .CEL.gz × N | — | `GEO_download_cel.R` |
| **RNA-seq（NCBI 整理）** | ❌（部分有 normalized） | 各种 | ✅ | `GEO_download_ncbicount.R` |
| **RNA-seq（仅作者上传）** | ❌ | `*counts.tsv.gz` / `*FPKM.txt` | ❌ | `GEO_download_count.R` |
| **甲基化** | ✅（β-value 矩阵） | IDAT × 2N | — | TODO（minfi 管线） |
| **单细胞** | ❌ | barcodes/features/matrix.mtx | — | TODO（10x 三件套） |

**Series Matrix 是否含表达值的判定**：解压后看前几个非 `!` 行是否是数值矩阵；若 `series_matrix` 文件 < 1 MB（仅元数据），多半是 RNA-seq/单细胞，需要去 `suppl/` 找原始文件。

---

## 3. SOFT 格式深度解析

### 3.1 三种 SOFT 子类型

| 文件 | 大小 | 用途 |
|---|---|---|
| `GSExxx_family.soft.gz` | 数十~数百 MB | 完整 Series（含所有 GSM 的表达矩阵） |
| `GPLxxx.soft.gz` | 几 MB ~ 几十 MB | 单平台注释（探针 → 基因映射） |
| `GSExxx_series_matrix.txt.gz` | 几 MB | 紧凑表达矩阵 + 精简元数据（**GEOquery 默认读这个**） |

### 3.2 SOFT 行前缀语义

| 前缀 | 含义 | 例 |
|---|---|---|
| `^` | 实体起始（新对象） | `^SERIES = GSE76262` |
| `!` | 属性键值对 | `!Series_title = ...` |
| `#` | 数据表列定义 | `#ID = Probe identifier` |
| `_table_begin` / `_table_end` | 表数据边界 | `!platform_table_begin` |
| 其它 | 表数据行（tab 分隔） | `1007_s_at\tU48705\t...` |

### 3.3 跳过头注释的安全方式

R 读 GPL 文本时常因不知道注释行数报错。本技能 `GEO_download_probe.R` 用了一段循环：

```r
con <- file(local_gpl_file, "r"); skip_lines <- 0L
repeat {
  line <- readLines(con, n = 1)
  if (length(line) == 0) break
  if (!grepl("^[#!^]", line)) break  # 第一行非 SOFT 标记 = 表头
  skip_lines <- skip_lines + 1L
}
close(con)
fread(local_gpl_file, skip = skip_lines, ...)
```

逻辑：跳过所有以 `^` / `!` / `#` 开头的元数据行，第一行非注释即视为表头。

### 3.4 重要属性键

#### Series 级
- `!Series_type`：与 esummary 的 `gdsType` **不完全一致**，例如 `Expression profiling by high throughput sequencing` ≠ `entryType` 字段。本技能 sniff 脚本同时检查两者。
- `!Series_supplementary_file`：FTP suppl/ 的绝对路径（一般 = `ftp://.../suppl/GSE..._RAW.tar`）。

#### Platform 级
- `!Platform_distribution`：`commercial` / `custom` / `non-commercial`，custom 平台的注释完整度最差。
- `!Platform_data_row_count`：探针总数，超 50w 多半是甲基化（450K / EPIC）。

#### Sample 级
- `!Sample_characteristics_ch1`：可重复出现的"标签:值"对，本技能脚本未默认解析（直接落到 clinical.csv 的 `characteristics_ch1.x` 列）；下游差异分析需 grep。

---

## 4. MINiML 格式（XML 等价物）

### 4.1 何时优先用 MINiML

| 场景 | 推荐 |
|---|---|
| 程序化遍历元数据 | MINiML（结构化） |
| Excel/文本查看 | SOFT |
| 极大 Series（> 1000 GSM） | MINiML 流式解析（lxml.iterparse）省内存 |

### 4.2 关键 XPath

```
./Series/Title
./Series/Overall-Design
./Sample/Channel/Characteristics[@tag='cell type']
./Platform/Data-Table/Column[@position='1']/Name
```

> 注意：MINiML 用了默认命名空间 `http://www.ncbi.nlm.nih.gov/geo/info/MINiML`，所有 XPath 必须带 `{ns}` 前缀，或用 `.iter()` 取本地名（本技能 sniff 脚本采用后者：`tag.split("}")[-1]`）。

---

## 5. FTP 目录命名规则

### 5.1 通用模式

```
ftp.ncbi.nlm.nih.gov/geo/{kind}/{prefix}/{accession}/{subdir}/
```

| 变量 | 取值 | 示例 |
|---|---|---|
| `kind` | `series` / `samples` / `platforms` | series |
| `prefix` | `{ACC[:-3]}nnn`（去掉末三位 + `nnn`） | GSE123nnn |
| `subdir` | `matrix` / `soft` / `miniml` / `suppl` / `annot` | suppl |

### 5.2 边界情况

| 情况 | 处理 |
|---|---|
| GSE 编号 < 1000（如 GSE12） | `prefix = "GSEnnn"`（不裁切） |
| GSE 编号刚好 4 位（如 GSE1234） | `prefix = "GSE1nnn"` |
| 5+ 位 | 标准裁切：`GSE76262 → GSE76nnn` |

本技能 `GEO_pre_download.py` 用正则强制 4+ 位（`re.match(r"^GSE(\d+)$", ...)`，长度 < 4 直接报错），实战中 GEO 编号都 ≥ 4 位。

### 5.3 NCBI 整理 RNA-seq counts URL

非 FTP 路径，独立于 `/geo/series/...`：

```
https://www.ncbi.nlm.nih.gov/geo/download/?
    type=rnaseq_counts
    &acc={GSE}
    &format=file
    &file={GSE}_{kind}_{genome}_NCBI.tsv.gz
```

| `{kind}` | 含义 |
|---|---|
| `raw_counts` | 整数 count 矩阵 |
| `norm_counts_FPKM` | FPKM |
| `norm_counts_TPM` | TPM |

`{genome}` 见下表（与本技能 `GEO_download_ncbicount.R` 的 `species_config` 对齐）：

| 物种 | genome_tag | annot 文件 |
|---|---|---|
| human | `GRCh38.p13` | `Human.GRCh38.p13.annot.tsv.gz` |
| mouse | `GRCm39` | `Mouse.GRCm39.annot.tsv.gz` |
| rat | `mRatBN7.2` | `Rat.mRatBN7.2.annot.tsv.gz`（仅嗅探脚本支持） |

> 何时存在：NCBI 自 2023 起逐步回填，2018 年后投稿的 RNA-seq 多数有；甲基化 / ChIP-seq / 单细胞**不会有**。

---

## 6. E-utilities 字段完整参考

### 6.1 `db=gds` 常用搜索字段

| 字段 | 例 | 备注 |
|---|---|---|
| `[Accession]` | `GSE76262[Accession]` | 精确匹配；用 `GSE` 不需要 `GDS` |
| `[Title]` | `"breast cancer"[Title]` | 引号锁词组 |
| `[Organism]` | `Homo sapiens[Organism]` | 拉丁名 |
| `[DataSet Type]` | `"Expression profiling by high throughput sequencing"[DataSet Type]` | 注意是 entryType 不是 gdsType |
| `[Platform]` | `GPL570[Platform]` | 含子平台版本 |
| `[Publication Date]` | `2020:2024[Publication Date]` | YYYY 或 YYYY/MM/DD |
| `[Number of Samples]` | `10:30[Number of Samples]` | 范围语法 |
| `[Supplementary Files]` | `CEL[Supplementary Files]` | 按附件扩展名筛 |

### 6.2 `db=geoprofiles` 字段

| 字段 | 备注 |
|---|---|
| `[Gene Name]` | HGNC 符号 |
| `[Reporter Type]` | 探针类型 |
| `[ANOVA p-value]` | 数值范围 |

### 6.3 esummary 关键返回字段（v2.0）

| 字段 | 类型 | 用途 |
|---|---|---|
| `Accession` | str | "GSE76262" |
| `gdsType` | str | "Expression profiling by array" |
| `entryType` | str | 多数同 gdsType，但 RNA-seq 时更细 |
| `GPL` | str/list | 平台列表 |
| `taxon` | str | 物种拉丁名 |
| `n_samples` | int | GSM 数 |
| `PDAT` | str | 发表日期 |
| `FTPLink` | str | 直达 series 目录 |

> 本技能 sniff 脚本同时读 `gdsType` 与 `entryType` 后做 OR 判断，规避 NCBI 偶发的字段缺失。

---

## 7. 探针 ID 类型对照

GEO 平台注释列千差万别，本技能 `GEO_download_probe.R` 通过 `--use-bitr TRUE --from-type ...` 二次转换。完整对照如下：

| GPL 列名样例 | 内容样例 | bitr `--from-type` |
|---|---|---|
| `Gene Symbol` / `symbol` / `Symbol` | `TP53` | （无需 bitr） |
| `GB_ACC` | `NM_014332.1` | `REFSEQ` |
| `ENTREZ_GENE_ID` | `780` | `ENTREZID` |
| `SPOT_ID`（纯数字） | `1, 780, 5982` | `ENTREZID` |
| `ENSEMBL_ID` / `ENSG...` | `ENSG00000141510` | `ENSEMBL` |
| `ENSMUSG...` | 鼠 ENSEMBL | `ENSEMBL`（加 `--org-db org.Mm.eg.db`） |
| `UNIPROT` / `UniProt` | `P04637` | `UNIPROT` |
| `gene_assignment`（Affymetrix HuGene/MoGene） | `NM_xxx // SYMBOL // ...` | 默认解析（CEL 脚本走这条） |
| `Composite Element` | Illumina/MEX 自定义 ID | 多数需 GPL.txt 手映射 |

### 7.1 ID 类型 sniff 规则（matrix 首行）

| 形如 | 推断为 |
|---|---|
| `ENSG\d+` / `ENSMUSG\d+` | `ENSEMBL` |
| 纯数字 `^\d+(\.\d+)?$` | `ENTREZID` |
| `^[A-Za-z][A-Za-z0-9_.\-]+` 且非 ENS 开头 | `SYMBOL` |
| 其它 | 兜底 `ENSEMBL` |

本技能 sniff 脚本仅取 `series_matrix` 第一列前 5 行做正则匹配，**不读全文**（5–10 KB 流量）。

### 7.2 版本号清理

ENSEMBL/REFSEQ 在 GPL 表里常带版本号（`ENSG00000001.4` / `NM_014332.1`），bitr 转换前必须去掉：

```r
probe2gene %>% mutate(symbol = str_remove(symbol, "\\..*$"))
```

否则 bitr 会丢失 30%+ 映射。

---

## 8. 平台技术特性

### 8.1 Affymetrix（GPL96/570/1261/6244/13112...）

| 维度 | 特性 |
|---|---|
| 探针后缀 | `_at`（unique） / `_s_at`（cross-hybrid） / `_x_at`（low spec） / `_a_at` |
| 多探针/基因 | 常见，**本技能默认按 rowmean 降序取首条** |
| 推荐流程 | CEL → `affy::ReadAffy` → `rma()`（背景校正 + quantile + log2） |
| 注释优先级 | `Gene Symbol` > `gene_assignment` > 调 bitr 转 `Entrez_Gene_ID` |
| 数据范围 | RMA 后 log2 空间 ~ [2, 14] |

### 8.2 Illumina BeadChip（GPL6883/6884/10558...）

| 维度 | 特性 |
|---|---|
| 探针 | `ILMN_\d+` |
| 重复探针 | 同基因往往多 ILMN ID，去重前先看 detection p-value |
| 推荐流程 | 优先用 author normalized；若有 IDAT 走 `lumi` 包 |
| 注释列 | 常见 `Symbol` 直接可用 |

### 8.3 Agilent 两通道（GPL4133/4134/6480/7202）

| 维度 | 特性 |
|---|---|
| 数据形态 | 双通道（Cy3 vs Cy5），需配对 |
| 文件 | `*.txt`（Feature Extraction 输出），非 CEL |
| 推荐流程 | `limma::read.maimages` → `normexp` 背景校正 → `loess` 通道内归一化 → `Aquantile` 通道间 |
| 本技能 | `GEO_download_cel.R` 在检测到 `.txt` 时自动走 Agilent 路径 |

### 8.4 RNA-seq（HiSeq/NovaSeq → GPL11154/16791/18573/20301/24676 ...）

| 阶段 | 文件 | 处理 |
|---|---|---|
| Raw FASTQ | 大多不在 GEO（在 SRA） | SRA-Toolkit 拉取 |
| Aligned BAM | 偶有 | 一般不重复 align |
| **Counts**（本技能关注） | `*_raw_counts.tsv.gz` / `*_counts.txt` | `ncbicount.R` / `count.R` |
| 归一化矩阵 | `*FPKM.tsv` / `*TPM.tsv` | 仅可视化用；差异分析必须用 counts |

> **counts vs FPKM/TPM 红线**：DESeq2/edgeR/limma-voom 要求**整数原始 counts**，喂 FPKM 会导致 dispersion 估计错误。本技能产出的 `expMatrix_Count.csv` 是合法输入。

### 8.5 甲基化（GPL13534 / GPL21145 / GPL8490）

| 平台 | 名称 | 探针数 |
|---|---|---|
| GPL8490 | 27K | ~27,000 |
| GPL13534 | 450K | ~485,000 |
| GPL21145 | EPIC (850K) | ~865,000 |
| GPL33022 | EPIC v2 | ~935,000 |

> Series Matrix 提供的多是 β-value（0~1），原始 IDAT 才能跑 BMIQ/SWAN 校正。本技能未支持，规划见 [SKILL.md 扩展计划](../SKILL.md#扩展计划)。

---

## 9. NCBI 限速与稳定性

### 9.1 速率上限（2024 起）

| 模式 | 上限 | 推荐 sleep |
|---|---|---|
| 无 API key | 3 req/s | 0.34s |
| 有 API key | 10 req/s | 0.10s |
| `epost` 批量 | 算 1 次 | 优于多次 esearch |

> 申请 key：https://account.ncbi.nlm.nih.gov/ → Account settings → API Key Management

### 9.2 常见错误码

| 状态 | 含义 | 对策 |
|---|---|---|
| 429 | 限速 | sleep + retry |
| 500/502/503/504 | NCBI 后端波动 | 指数退避；本技能 `GEO_pre_download.py` 已配 `Retry(backoff_factor=1.0)` |
| 414 | URL 太长（id 太多） | 改 `epost` + history server |
| 403 / 404 | 资源不存在或被撤回 | 看 GSE 是否 superseded |

### 9.3 history server（大批量必用）

```python
h = Entrez.esearch(db="gds", term=q, usehistory="y", retmax=10000)
r = Entrez.read(h); h.close()
key, env = r["QueryKey"], r["WebEnv"]
# 后续 esummary 直接传 key+env，免传 id 列表
Entrez.esummary(db="gds", query_key=key, WebEnv=env, retstart=0, retmax=500)
```

---

## 10. 大数据集处理建议

| 场景 | 建议 |
|---|---|
| 单 GSE > 500 GSM | 用 `data.table::fread` 读 matrix；不要 `read.table` |
| 多 GSE 合并（ComBat） | float32 + 仅保留交集基因；2 万基因 × 1000 样本约 80 MB |
| Series Matrix 巨大（> 500 MB） | 改用 family.soft.gz + GEOparse 流式 |
| 平台注释频繁查询 | 把 `GPL.txt` 缓存到本地 `~/.cache/geo/` |

---

## 11. 引用规范

> 使用 GEO 数据应在论文中：
> 1. 注明 GSE 号（"... obtained from GEO under accession GSE76262 ..."）
> 2. 引用原始文献（每个 GSE 的 PubMed ID 在 esummary 的 `PubMedIds` 字段）
> 3. 引用 GEO 本身：Barrett T. et al. *NCBI GEO: archive for functional genomics data sets — update.* Nucleic Acids Res. 2013.

---

## 参考链接

- [GEO 主页](https://www.ncbi.nlm.nih.gov/geo/)
- [E-utilities 文档](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [SOFT 格式规范](https://www.ncbi.nlm.nih.gov/geo/info/soft.html)
- [MINiML 格式规范](https://www.ncbi.nlm.nih.gov/geo/info/MINiML.html)
- [GEOquery（Bioconductor）](https://bioconductor.org/packages/release/bioc/html/GEOquery.html)
- [GEOparse（Python）](https://geoparse.readthedocs.io/)
- [AnnoProbe](https://github.com/jmzeng1314/AnnoProbe)
