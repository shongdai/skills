#!/usr/bin/env Rscript
# =============================================================================
# GEO_download_cel.R - GEO 芯片 RAW 数据下载与预处理（CLI 版）
#
# 描述: 从 GEO 下载 supplementary RAW（.CEL/.txt），自动选择处理路径：
#         - .CEL → affy::ReadAffy + rma()（背景校正 + quantile + log2）
#         - .txt → limma::read.maimages + normexp（Agilent 路径）
#       然后做探针 → 基因 symbol 转换、基因过滤、log2 检测。
# 作者: 科研木鱼（闲鱼/小红书）
# 版本: 2.0  —— 命令行参数化
#
# 使用方法（CLI）:
#   Rscript GEO_download_cel.R --gse GSE5281
#   Rscript GEO_download_cel.R --gse GSE5281 --out ./data --diff TRUE
#   Rscript GEO_download_cel.R --gse GSE5281 --use-idmap TRUE --proxy http://127.0.0.1:7897
#   Rscript GEO_download_cel.R --help
#
# 参数:
#   --gse <id>           必填，单个 GSE 号（CEL 路径一次处理一个 GSE）
#   --out <dir>          输出根目录，默认 "."
#   --use-idmap <T/F>    是否优先 AnnoProbe::idmap() 在线注释，默认 FALSE
#   --probe-col <name>   本地 GPL 探针 ID 列名，默认 "ID"
#   --gene-col <name>    本地 GPL 基因注释列名，默认 "gene_assignment"
#   --diff <T/F>         是否生成 group/diff 模板，默认 FALSE
#   --proxy <url>        HTTP 代理 URL，默认 http://127.0.0.1:7897；不需代理时传空串 --proxy ""
#   --timeout <sec>      下载超时秒数，默认 600
#   --help               显示帮助
# =============================================================================

# ---------- CLI 参数解析 ----------
parse_cli_args <- function(args) {
  defaults <- list(
    gse        = NULL,
    out        = ".",
    use_idmap  = FALSE,
    probe_col  = "ID",
    gene_col   = "gene_assignment",
    diff       = FALSE,
    proxy      = "http://127.0.0.1:7897",
    timeout    = 600L,
    help       = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    switch(a,
      "--gse"        = { defaults$gse       <- val; i <- i + 2 },
      "--out"        = { defaults$out       <- val; i <- i + 2 },
      "--use-idmap"  = { defaults$use_idmap <- toupper(val) %in% c("T", "TRUE", "1", "YES"); i <- i + 2 },
      "--probe-col"  = { defaults$probe_col <- val; i <- i + 2 },
      "--gene-col"   = { defaults$gene_col  <- val; i <- i + 2 },
      "--diff"       = { defaults$diff      <- toupper(val) %in% c("T", "TRUE", "1", "YES"); i <- i + 2 },
      "--proxy"      = { defaults$proxy     <- val; i <- i + 2 },
      "--timeout"    = { defaults$timeout   <- as.integer(val); i <- i + 2 },
      "--help"       = { defaults$help      <- TRUE; i <- i + 1 },
      "-h"           = { defaults$help      <- TRUE; i <- i + 1 },
      {
        cat("[警告] 未识别参数: ", a, "（已跳过）\n", sep = "")
        i <- i + 1
      }
    )
  }
  defaults
}

ARGS <- parse_cli_args(commandArgs(trailingOnly = TRUE))

if (isTRUE(ARGS$help) || is.null(ARGS$gse)) {
  cat("用法: Rscript GEO_download_cel.R --gse GSE_ID [选项]\n\n")
  cat("必填:\n")
  cat("  --gse <id>           单个 GSE 号（RAW 路径一次处理一个）\n\n")
  cat("可选:\n")
  cat("  --out <dir>          输出根目录（默认 .）\n")
  cat("  --use-idmap <T/F>    在线注释 idmap()（默认 FALSE）\n")
  cat("  --probe-col <name>   本地 GPL 探针列名（默认 ID）\n")
  cat("  --gene-col <name>    本地 GPL 基因列名（默认 gene_assignment）\n")
  cat("  --diff <T/F>         生成 group/diff 模板（默认 FALSE）\n")
  cat("  --proxy <url>        HTTP 代理 URL（默认 http://127.0.0.1:7897）\n")
  cat("  --timeout <sec>      下载超时（默认 600）\n")
  cat("  --help               显示帮助\n\n")
  cat("示例:\n")
  cat("  Rscript GEO_download_cel.R --gse GSE5281\n")
  cat("  Rscript GEO_download_cel.R --gse GSE5281 --out ./data --diff TRUE\n")
  if (is.null(ARGS$gse) && !isTRUE(ARGS$help)) {
    cat("\n[错误] 缺少必填参数 --gse\n")
    quit(status = 1)
  }
  quit(status = 0)
}

gse_id          <- ARGS$gse
use_idmap       <- isTRUE(ARGS$use_idmap)
probe_id_col    <- ARGS$probe_col
gene_symbol_col <- ARGS$gene_col
diff <- if (isTRUE(ARGS$diff)) "TRUE" else "FALSE"
out_root <- ARGS$out
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
output_dir <- file.path(out_root, gse_id)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("========== CLI 参数 ==========\n")
cat("GSE      :", gse_id, "\n")
cat("输出目录 :", normalizePath(output_dir, mustWork = FALSE), "\n")
cat("use_idmap:", use_idmap, "\n")
cat("probe-col:", probe_id_col, "\n")
cat("gene-col :", gene_symbol_col, "\n")
cat("diff     :", diff, "\n")
cat("proxy    :", ifelse(nzchar(ARGS$proxy), ARGS$proxy, "(未设)"), "\n")
cat("timeout  :", ARGS$timeout, "\n\n")

# ---------- 加载依赖 ----------
# 如缺少包，请运行: Rscript scripts/install_dependencies.R
suppressPackageStartupMessages({
  library(limma)
  library(affy)
  library(GEOquery)
  library(AnnoProbe)
  library(devtools)
  library(data.table)
  library(archive)
  library(tidyverse)
})
options(timeout = ARGS$timeout)
options(download.file.method.GEOquery = "libcurl")
if (nzchar(ARGS$proxy)) {
  Sys.setenv(http_proxy = ARGS$proxy, https_proxy = ARGS$proxy)
  cat("[代理已设置]", ARGS$proxy, "\n\n")
}

# ============================================================
# 第一步：下载 Series Matrix 与临床数据
# ============================================================
cat("========== 第一步：下载 Series Matrix ==========\n")
gse <- getGEO(gse_id, destdir = output_dir, getGPL = FALSE)
gpl <- annotation(gse[[1]])
cat("自动检测到GPL平台:", gpl, "\n")

cat("\n========== 第二步：保存临床数据 ==========\n")
clinical_data <- pData(gse[[1]])
cat("临床数据维度:", dim(clinical_data), "\n")
write.csv(clinical_data, file = file.path(output_dir, "clinical.csv"), row.names = FALSE)
cat("临床数据已保存至", file.path(output_dir, "clinical.csv"), "\n")

if (diff == "TRUE") {
  group_file <- file.path(output_dir, "group.txt")
  diff_file  <- file.path(output_dir, "diff.txt")
  if (file.exists(group_file) && file.exists(diff_file)) {
    cat("分组文件与差异比较文件已存在，跳过生成\n")
  } else {
    write.table(data.frame(sample = clinical_data$geo_accession,
                           group  = clinical_data$title,
                           stringsAsFactors = FALSE),
                file = group_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
    cat("分组文件已生成:", group_file, "\n")
    write.table(data.frame(control = character(0), treat = character(0)),
                file = diff_file, sep = "\t", quote = FALSE,
                row.names = FALSE, col.names = TRUE)
    cat("差异比较文件已生成:", diff_file, "\n")
  }
} else {
  cat("diff =", diff, "，跳过生成分组文件\n")
}

# 下载 supplementary 包含 RAW 数据
getGEOSuppFiles(gse_id, makeDirectory = FALSE, baseDir = output_dir)

# 解压 _RAW.tar
raw_dir <- file.path(output_dir, paste0(gse_id, "_RAW"))
untar(file.path(output_dir, paste0(gse_id, "_RAW.tar")), exdir = raw_dir)

# 自动检测样本文件类型
cel_files <- list.files(raw_dir, pattern = "\\.CEL(\\.gz)?$", full.names = FALSE)
txt_files <- list.files(raw_dir, pattern = "\\.txt(\\.gz)?$", full.names = FALSE)
cat("检测到 .CEL 文件:", length(cel_files), "  .txt 文件:", length(txt_files), "\n")

if (length(cel_files) > 0) {
  cat("走 Affymetrix 流程（affy::ReadAffy + rma）\n")
  RawAffy <- ReadAffy(filenames = cel_files, celfile.path = raw_dir)
  eset <- rma(RawAffy)
  exp <- exprs(eset)
  colnames(exp) <- sapply(strsplit(colnames(exp), "\\."), `[`, 1)
  exp <- as.data.frame(exp)
  cat("探针矩阵维度:", dim(exp), "（rma 已含背景校正+quantile+log2）\n")
} else if (length(txt_files) > 0) {
  cat("走 Agilent 流程（limma::read.maimages + normexp）\n")
  RawData <- read.maimages(files = txt_files, source = "agilent",
                           path = raw_dir, green.only = TRUE)
  EList <- backgroundCorrect(RawData, method = "normexp",
                             normexp.method = "rma", offset = 50)
  mm <- range(EList$E)
  cat("最小值:", mm[1], " 最大值:", mm[2], "\n")
  exp <- EList$E
  rownames(exp) <- EList$genes$ProbeName
  colnames(exp) <- sapply(strsplit(colnames(exp), "_"), `[`, 1)
  exp <- as.data.frame(exp)
  cat("探针矩阵维度:", dim(exp), "\n")
} else {
  stop("在 ", raw_dir, " 中未找到 .CEL 或 .txt 原始文件")
}

# ============================================================
# 工具函数：子串匹配列名（不区分大小写）
# ============================================================
find_col_case_insensitive <- function(col_names, target) {
  if (is.null(target) || target == "") return(target)
  idx <- which(grepl(target, col_names, ignore.case = TRUE))
  if (length(idx) == 0) {
    cat("[列匹配警告] 未找到包含 '", target, "' 的列，可用列: ",
        paste(col_names, collapse = ", "), "\n", sep = "")
    return(target)
  }
  col_names[idx[1]]
}

# ============================================================
# 第三步：ID 转换
# ============================================================
cat("\n========== 第三步：ID转换 ==========\n")
probe2gene <- tryCatch({
  if (!isTRUE(use_idmap)) stop("已禁用 idmap() 尝试，直接使用本地GPL文件")
  cat("尝试 idmap() type='pipe'...\n")
  df <- idmap(gpl, type = "pipe"); cat("使用 idmap() type='pipe' 成功获取注释\n"); df
}, error = function(e_pipe) {
  cat("idmap() type='pipe' 调用失败:", conditionMessage(e_pipe), "\n")
  tryCatch({
    cat("尝试 idmap() type='soft'...\n")
    df <- idmap(gpl, type = "soft"); cat("使用 idmap() type='soft' 成功获取注释\n"); df
  }, error = function(e_soft) {
    cat("idmap() type='soft' 调用失败:", conditionMessage(e_soft), "\n")
    tryCatch({
      cat("尝试 idmap() type='bioc'...\n")
      df <- idmap(gpl, type = "bioc"); cat("使用 idmap() type='bioc' 成功获取注释\n"); df
    }, error = function(e_bioc) {
      cat("idmap() type='bioc' 调用失败:", conditionMessage(e_bioc), "\n")
      local_gpl_file <- {
        gpl_files <- list.files(path = output_dir,
                                pattern = paste0("^", gpl, ".*\\.txt$"),
                                full.names = TRUE)
        if (length(gpl_files) > 0) gpl_files[1] else ""
      }
      cat("切换到本地GPL文件:", local_gpl_file, "\n")
      if (!file.exists(local_gpl_file)) {
        cat("本地GPL文件不存在，尝试从NCBI下载...\n")
        gpl_url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
                          gpl, "&targ=self&form=text&view=data")
        dest_file <- file.path(output_dir, paste0(gpl, ".txt"))
        download_ok <- tryCatch({
          download.file(gpl_url, destfile = dest_file, mode = "wb", quiet = TRUE); TRUE
        }, error = function(e_dl) {
          cat("从NCBI下载GPL文件失败:", conditionMessage(e_dl), "\n"); FALSE
        })
        if (download_ok && file.exists(dest_file)) {
          local_gpl_file <- dest_file
          cat("已下载GPL文件至:", local_gpl_file, "\n")
        } else stop("本地GPL文件不存在且下载失败: ", gpl)
      }

      con <- file(local_gpl_file, "r"); skip_lines <- 0L
      repeat {
        line <- readLines(con, n = 1)
        if (length(line) == 0) break
        if (!grepl("^[#!^]", line)) break
        skip_lines <- skip_lines + 1L
      }
      close(con)
      cat("检测到注释行数:", skip_lines, "\n")

      gpl_data <- as.data.frame(fread(local_gpl_file, sep = "\t", header = TRUE,
                                      fill = TRUE, quote = "", skip = skip_lines))
      actual_probe_col <- find_col_case_insensitive(colnames(gpl_data), probe_id_col)
      actual_gene_col  <- find_col_case_insensitive(colnames(gpl_data), gene_symbol_col)
      cat("[列匹配] 探针ID列:", probe_id_col, "->", actual_probe_col, "\n")
      cat("[列匹配] 基因注释列:", gene_symbol_col, "->", actual_gene_col, "\n")
      gpl_data <- gpl_data %>% select(all_of(c(actual_probe_col, actual_gene_col)))

      gpl_data %>%
        filter(!is.na(.data[[actual_gene_col]]) & .data[[actual_gene_col]] != "---") %>%
        mutate(symbol = map_chr(str_split(.data[[actual_gene_col]], " /// "), function(x) {
          if (grepl("gene_assignment", actual_gene_col, ignore.case = TRUE)) {
            parts <- str_split(x[1], " // ")[[1]]
            if (length(parts) >= 2) str_trim(parts[2]) else NA_character_
          } else str_trim(x[1])
        })) %>%
        filter(!is.na(symbol) & symbol != tolower(symbol)) %>%
        select(all_of(actual_probe_col), symbol) %>%
        rename(ID = all_of(actual_probe_col))
    })
  })
})

probe2gene <- probe2gene %>%
  mutate(symbol = str_split(symbol, " /// ") %>% map_chr(1)) %>%
  filter(!is.na(symbol) & symbol != "") %>%
  { if ("ID" %in% colnames(.)) rename(., probe_id = ID) else . } %>%
  select(probe_id, symbol)

probe_file <- file.path(output_dir, "probe2gene.csv")
write.csv(probe2gene, file = probe_file, row.names = FALSE)
cat("注释文件已保存至", probe_file, "\n")
cat("注释文件维度:", dim(probe2gene), "\n")

# ============================================================
# 第四步：数据转换与ID映射
# ============================================================
cat("\n========== 第四步：数据转换与ID映射 ==========\n")
exp_symbol <- exp %>% as.data.frame() %>% rownames_to_column(var = "probe_id") %>%
  inner_join(probe2gene, by = "probe_id") %>% na.omit() %>%
  select(-probe_id) %>% relocate(symbol)
cat("ID映射后维度:", dim(exp_symbol), "\n")
cat("基因数量:", length(unique(exp_symbol$symbol)), "\n")

# ============================================================
# 第五步：过滤基因
# ============================================================
cat("\n========== 第五步：过滤基因 ==========\n")
exp_filtered <- exp_symbol %>%
  mutate(rowmean = rowMeans(.[, 2:ncol(.)], na.rm = TRUE)) %>%
  arrange(desc(rowmean)) %>% distinct(symbol, .keep_all = TRUE) %>%
  filter(rowmean > 0) %>% select(-rowmean) %>%
  filter(!is.na(symbol) & symbol != "") %>%
  column_to_rownames("symbol")
cat("过滤后维度:", dim(exp_filtered), "\n")

# ============================================================
# 第六步：检查是否需要log转换
# ============================================================
cat("\n========== 第六步：检查数据范围 ==========\n")
cat("转换前数据范围:", range(exp_filtered), "\n")
if (max(exp_filtered) > 20) {
  cat("数据最大值 > 20，进行log2转换\n")
  exp_filtered <- log2(exp_filtered + 1)
  cat("转换后数据范围:", range(exp_filtered), "\n")
} else {
  cat("数据无需log转换\n")
}

# ============================================================
# 第七步：保存结果
# ============================================================
cat("\n========== 第七步：保存结果 ==========\n")
exp_final <- exp_filtered %>% rownames_to_column(var = "symbol")
fwrite(exp_final, file = file.path(output_dir, "expMatrix.csv"),
       sep = ",", quote = FALSE, row.names = FALSE)
cat("表达矩阵已保存至", file.path(output_dir, "expMatrix.csv"), "\n")
cat("最终矩阵维度:", dim(exp_final), "\n")
cat("========== 分析完成 ==========\n")
