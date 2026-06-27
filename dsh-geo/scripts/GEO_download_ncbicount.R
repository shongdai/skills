#!/usr/bin/env Rscript
# =============================================================================
# GEO_download_ncbicount.R - GEO NCBI 标准化 RNA-seq counts 下载与预处理（CLI 版）
#
# 描述: 自动从 NCBI 下载官方整理的 RNA-seq counts 表达矩阵（raw_counts / FPKM /
#       TPM），完成基因 ID 到基因 symbol 的转换、基因过滤。
#       基于 gse_id 自动匹配 NCBI 整理的 RNA-seq counts 表达矩阵文件。
#       默认读取物种对应的注释文件（如 Human.GRCh38.p13.annot.tsv.gz /
#       Mouse.GRCm39.annot.tsv.gz）；若不存在则使用 clusterProfiler::bitr。
#       支持人/小鼠物种切换。
# 作者: 科研木鱼（小红书：科研木鱼）
# 版本: 2.0  —— 命令行参数化 + 物种切换
#
# 使用方法（CLI）:
#   Rscript GEO_download_ncbicount.R --gse GSE56545
#   Rscript GEO_download_ncbicount.R --gse GSE56545,GSE70089,GSE33294 --diff TRUE
#   Rscript GEO_download_ncbicount.R --gse GSE100001 --species mouse
#   Rscript GEO_download_ncbicount.R --gse GSE56545 --out ./data --proxy http://127.0.0.1:7897
#   Rscript GEO_download_ncbicount.R --help
#
# 参数:
#   --gse <ids>        必填，GSE 号（多个用逗号分隔）
#   --species <name>   human / mouse，默认 human
#   --out <dir>        输出根目录，默认 "."
#   --diff <T/F>       是否生成 group/diff 模板，默认 FALSE
#   --proxy <url>      HTTP 代理 URL，默认 http://127.0.0.1:7897；空串关闭
#   --timeout <sec>    下载超时秒数，默认 600
#   --help             显示帮助
# =============================================================================

# ---------- CLI 参数解析 ----------
parse_cli_args <- function(args) {
  defaults <- list(
    gse        = NULL,
    species    = "human",
    out        = ".",
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
      "--species"    = { defaults$species   <- val; i <- i + 2 },
      "--out"        = { defaults$out       <- val; i <- i + 2 },
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
  cat("用法: Rscript GEO_download_ncbicount.R --gse GSE_ID[,GSE_ID2] [选项]\n\n")
  cat("必填:\n")
  cat("  --gse <ids>        GSE 号（多个用逗号分隔）\n\n")
  cat("可选:\n")
  cat("  --species <name>   human / mouse（默认 human）\n")
  cat("  --out <dir>        输出根目录（默认 .）\n")
  cat("  --diff <T/F>       生成 group/diff 模板（默认 FALSE）\n")
  cat("  --proxy <url>      HTTP 代理 URL（默认 http://127.0.0.1:7897）\n")
  cat("  --timeout <sec>    下载超时（默认 600）\n")
  cat("  --help             显示帮助\n\n")
  cat("示例:\n")
  cat("  Rscript GEO_download_ncbicount.R --gse GSE56545\n")
  cat("  Rscript GEO_download_ncbicount.R --gse GSE100001 --species mouse\n")
  cat("  Rscript GEO_download_ncbicount.R --gse GSE56545,GSE70089 --diff TRUE\n")
  if (is.null(ARGS$gse) && !isTRUE(ARGS$help)) {
    cat("\n[错误] 缺少必填参数 --gse\n")
    quit(status = 1)
  }
  quit(status = 0)
}

# 解析多个 GSE
gse_id_vec <- strsplit(ARGS$gse, ",")[[1]] |> trimws() |> (\(x) x[nzchar(x)])()
species    <- ARGS$species
out_root   <- ARGS$out
diff       <- if (isTRUE(ARGS$diff)) "TRUE" else "FALSE"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# 物种配置
species_config <- list(
  human = list(
    genome_tag     = "GRCh38.p13",
    annot_filename = "Human.GRCh38.p13.annot.tsv.gz",
    org_db         = "org.Hs.eg.db"
  ),
  mouse = list(
    genome_tag     = "GRCm39",
    annot_filename = "Mouse.GRCm39.annot.tsv.gz",
    org_db         = "org.Mm.eg.db"
  )
)
if (!species %in% names(species_config)) {
  stop("[错误] 不支持的 species：", species, "（请用 human / mouse）")
}
genome_tag     <- species_config[[species]]$genome_tag
annot_filename <- species_config[[species]]$annot_filename
org_db_name    <- species_config[[species]]$org_db
cat(sprintf("[物种] species=%s -> genome_tag=%s, annot=%s, OrgDb=%s\n",
            species, genome_tag, annot_filename, org_db_name))

cat("========== CLI 参数 ==========\n")
cat("GSE 列表  :", paste(gse_id_vec, collapse = ", "), "\n")
cat("输出根目录:", normalizePath(out_root, mustWork = FALSE), "\n")
cat("species   :", species, "\n")
cat("genome_tag:", genome_tag, "\n")
cat("annot文件 :", annot_filename, "\n")
cat("org_db    :", org_db_name, "\n")
cat("diff      :", diff, "\n")
cat("proxy     :", ifelse(nzchar(ARGS$proxy), ARGS$proxy, "(未设)"), "\n")
cat("timeout   :", ARGS$timeout, "\n\n")

# ---------- 加载依赖 ----------
suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})
options(timeout = ARGS$timeout)
if (nzchar(ARGS$proxy)) {
  Sys.setenv(http_proxy = ARGS$proxy, https_proxy = ARGS$proxy)
  cat("[代理已设置]", ARGS$proxy, "\n\n")
}

# ============================================================
# 单个 GSE 完整处理
# ============================================================
process_one_gse <- function(gse) {
  cat(sprintf("\n############################################\n"))
  cat(sprintf("# 处理 GSE: %s\n", gse))
  cat(sprintf("############################################\n"))

  data_dir <- file.path(out_root, gse)
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

  clinical_csv  <- file.path(data_dir, "clinical.csv")
  exp_count_csv <- file.path(data_dir, "expMatrix_Count.csv")
  exp_fpkm_csv  <- file.path(data_dir, "expMatrix_FPKM.csv")
  exp_tpm_csv   <- file.path(data_dir, "expMatrix_TPM.csv")

  # 给定 .tsv.gz 文件名，返回对应的 expMatrix_*.csv 路径
  exp_output_for <- function(expr_file) {
    expr_type <- expr_file %>%
      gsub("^[^_]+_", "", .) %>%
      gsub(paste0("_", genome_tag, "_NCBI\\.tsv\\.gz$"), "", .) %>%
      gsub("^norm_counts_", "", .) %>%
      gsub("^raw_counts$", "Count", .)
    file.path(data_dir, paste0("expMatrix_", expr_type, ".csv"))
  }

  if (file.exists(clinical_csv) &&
      file.exists(exp_count_csv) &&
      file.exists(exp_fpkm_csv) &&
      file.exists(exp_tpm_csv)) {
    cat(sprintf("所有产物已存在，跳过整个 GSE: %s\n", gse))
    return(invisible(TRUE))
  }

  download_list <- list(
    list(
      url = sprintf(
        "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_raw_counts_%s_NCBI.tsv.gz",
        gse, gse, genome_tag
      ),
      dest = sprintf("%s_raw_counts_%s_NCBI.tsv.gz", gse, genome_tag)
    ),
    list(
      url = sprintf(
        "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_norm_counts_FPKM_%s_NCBI.tsv.gz",
        gse, gse, genome_tag
      ),
      dest = sprintf("%s_norm_counts_FPKM_%s_NCBI.tsv.gz", gse, genome_tag)
    ),
    list(
      url = sprintf(
        "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=%s&format=file&file=%s_norm_counts_TPM_%s_NCBI.tsv.gz",
        gse, gse, genome_tag
      ),
      dest = sprintf("%s_norm_counts_TPM_%s_NCBI.tsv.gz", gse, genome_tag)
    ),
    list(
      url = sprintf(
        "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts&file=%s",
        annot_filename
      ),
      dest = annot_filename
    )
  )

  # ============================================================
  # 第一步：下载 NCBI RNA-seq 数据
  # ============================================================
  cat("---------- 第一步：下载表达矩阵 ----------\n")

  for (item in download_list) {
    dest_path <- file.path(data_dir, item$dest)

    if (file.exists(dest_path)) {
      cat(sprintf("已存在本地文件，跳过下载: %s\n", dest_path))
      next
    }

    if (grepl("_(raw|norm)_counts_", item$dest)) {
      out_path <- exp_output_for(item$dest)
      if (file.exists(out_path)) {
        cat(sprintf("产物 %s 已存在，跳过下载: %s\n", basename(out_path), dest_path))
        next
      }
    }

    cat(sprintf("下载: %s\n", item$url))
    download_ok <- tryCatch({
      download.file(item$url, destfile = dest_path, method = "libcurl", mode = "wb")
      TRUE
    }, error = function(e) {
      cat(sprintf("下载失败: %s\n  错误: %s\n", item$url, conditionMessage(e)))
      FALSE
    }, warning = function(w) {
      cat(sprintf("下载警告: %s\n  信息: %s\n", item$url, conditionMessage(w)))
      file.exists(dest_path)
    })

    if (!download_ok || !file.exists(dest_path)) {
      cat(sprintf(
        "下载失败，请手动从浏览器下载: %s\n  保存到: %s\n",
        item$url, dest_path
      ))
    } else {
      cat(sprintf("下载完成: %s\n", dest_path))
    }
  }

  # ============================================================
  # 第二步：下载临床数据
  # ============================================================
  cat("\n---------- 第二步：下载临床数据 ----------\n")

  if (file.exists(clinical_csv)) {
    cat(sprintf("临床数据已存在，跳过 getGEO(): %s\n", clinical_csv))
    clinical_data <- read.csv(clinical_csv, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    gse_obj <- getGEO(gse, destdir = data_dir, getGPL = FALSE)
    clinical_data <- pData(gse_obj[[1]])
    cat("临床数据维度:", dim(clinical_data), "\n")
    write.csv(clinical_data, file = clinical_csv, row.names = FALSE)
    cat(sprintf("临床数据已保存至 %s\n", clinical_csv))
  }

  # ============================================================
  # 第三步：读取基因注释文件
  # ============================================================
  cat("\n---------- 第三步：读取基因注释文件 ----------\n")

  annot_path <- file.path(data_dir, annot_filename)
  gene_mapping <- NULL

  if (!file.exists(annot_path)) {
    cat("注释文件不存在，将使用 clusterProfiler::bitr 进行基因ID转换\n")
    suppressMessages({
      library(clusterProfiler)
      library(org_db_name, character.only = TRUE)
    })
  } else {
    gene_mapping <- annot_path %>%
      gzfile() %>%
      read.delim(sep = "\t", check.names = FALSE) %>%
      dplyr::select(GeneID, Symbol) %>%
      dplyr::mutate(GeneID = as.character(GeneID), Symbol = as.character(Symbol)) %>%
      dplyr::filter(!is.na(Symbol), Symbol != "")
    cat(sprintf("注释文件共有 %d 个基因\n", nrow(gene_mapping)))
  }

  # ============================================================
  # 第四步：处理表达矩阵文件
  # ============================================================
  cat("\n---------- 第四步：处理表达矩阵 ----------\n")

  expression_files <- c(
    sprintf("%s_raw_counts_%s_NCBI.tsv.gz", gse, genome_tag),
    sprintf("%s_norm_counts_FPKM_%s_NCBI.tsv.gz", gse, genome_tag),
    sprintf("%s_norm_counts_TPM_%s_NCBI.tsv.gz", gse, genome_tag)
  )

  for (expr_file in expression_files) {
    input_path  <- file.path(data_dir, expr_file)
    output_path <- exp_output_for(expr_file)

    if (file.exists(output_path)) {
      cat(sprintf("产物已存在，跳过处理: %s\n", output_path))
      next
    }

    if (!file.exists(input_path)) {
      cat(sprintf("文件不存在，跳过: %s\n", input_path))
      next
    }

    cat(sprintf("正在处理: %s\n", expr_file))

    expr_data <- input_path %>%
      gzfile() %>%
      read.delim(sep = "\t", check.names = FALSE) %>%
      dplyr::mutate(GeneID = as.character(GeneID))

    current_mapping <- gene_mapping
    if (is.null(current_mapping)) {
      gene_ids <- unique(expr_data$GeneID)
      bitr_result <- clusterProfiler::bitr(
        gene_ids,
        fromType = "ENTREZID",
        toType   = "SYMBOL",
        OrgDb    = get(org_db_name)
      )
      current_mapping <- bitr_result %>%
        dplyr::rename(GeneID = ENTREZID, Symbol = SYMBOL) %>%
        dplyr::mutate(GeneID = as.character(GeneID), Symbol = as.character(Symbol))
    }

    expr_data <- expr_data %>%
      dplyr::left_join(current_mapping, by = "GeneID") %>%
      dplyr::mutate(GeneID = Symbol) %>%
      dplyr::select(-Symbol) %>%
      dplyr::filter(!is.na(GeneID), GeneID != "") %>%
      dplyr::mutate(mean_expr = rowMeans(dplyr::select(., -GeneID), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(mean_expr)) %>%
      dplyr::distinct(GeneID, .keep_all = TRUE) %>%
      dplyr::filter(mean_expr > 0) %>%
      dplyr::select(-mean_expr)

    cat("表达矩阵维度:", dim(expr_data), "\n")
    cat("表达矩阵范围:", range(expr_data[, -1], na.rm = TRUE), "\n")

    write.csv(expr_data, file = output_path, row.names = FALSE)
    cat(sprintf("已保存: %s\n", output_path))
  }

  # ============================================================
  # 第五步：生成分组文件
  # ============================================================
  if (exists("diff") && diff == "TRUE") {
    cat("\n---------- 第五步：生成分组文件 ----------\n")
    group_file <- file.path(data_dir, "group.txt")
    diff_file  <- file.path(data_dir, "diff.txt")

    if (file.exists(group_file) && file.exists(diff_file)) {
      cat("分组文件与差异比较文件已存在，跳过生成\n")
    } else {
      write.table(
        data.frame(
          sample = clinical_data$geo_accession,
          group  = clinical_data$title,
          stringsAsFactors = FALSE
        ),
        file = group_file, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE
      )
      cat("分组文件已生成:", group_file, "\n")
      write.table(
        data.frame(control = character(0), treat = character(0)),
        file = diff_file, sep = "\t", quote = FALSE,
        row.names = FALSE, col.names = TRUE
      )
      cat("差异比较文件已生成:", diff_file, "\n")
    }
  }
  invisible(TRUE)
}

# ============================================================
# 主循环
# ============================================================
cat("========== GEO_download_ncbicount.R (CLI) ==========\n")
for (gse in gse_id_vec) {
  tryCatch(process_one_gse(gse), error = function(e) {
    cat("[整体跳过] ", gse, " 顶层错误: ", conditionMessage(e), "\n", sep = "")
  })
}
cat("\n全部 GSE 处理完成！\n")
