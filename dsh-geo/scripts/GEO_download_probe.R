#!/usr/bin/env Rscript
# =============================================================================
# GEO_download_probe.R - GEO 芯片 Series Matrix 下载与预处理（CLI 版）
#
# 描述: 从GEO数据库下载已归一化的Series Matrix，直接复用预处理结果，
#       完成探针ID到基因symbol的转换、基因过滤和log2转换。支持多GSE批量、
#       多平台。
# 作者: 科研木鱼（闲鱼/小红书：科研木鱼）
# 版本: 3.0  —— 命令行参数化 + bitr 二次转换
#
# 使用方法（CLI）:
#   Rscript GEO_download_probe.R --gse GSE76262
#   Rscript GEO_download_probe.R --gse GSE76262,GSE12345
#   Rscript GEO_download_probe.R --gse GSE76262 --out ./data --diff TRUE
#   Rscript GEO_download_probe.R --gse GSE76262 --use-idmap TRUE
#   Rscript GEO_download_probe.R --gse GSE8479 --use-bitr TRUE --gene-col GB_ACC --from-type REFSEQ
#   Rscript GEO_download_probe.R --help
#
# 参数:
#   --gse <ids>          必填，GSE 号（多个用逗号分隔）
#   --out <dir>          输出根目录，默认 "."
#   --use-idmap <T/F>    AnnoProbe::idmap() 在线注释，默认 FALSE
#   --probe-col <name>   本地 GPL 探针 ID 列名（模糊匹配），默认 "ID"
#   --gene-col <name>    本地 GPL 基因注释列名（模糊匹配），默认 "symbol"
#   --diff <T/F>         是否生成 group/diff 模板，默认 FALSE
#   --use-bitr <T/F>     clusterProfiler::bitr 二次转换，默认 FALSE
#   --from-type <type>   bitr 输入 ID 类型，默认 "ENTREZID"
#   --org-db <pkg>       物种注释包，默认 "org.Hs.eg.db"
#   --proxy <url>        HTTP 代理 URL，默认 http://127.0.0.1:7897；空串关闭
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
    gene_col   = "symbol",
    diff       = FALSE,
    use_bitr   = FALSE,
    from_type  = "ENTREZID",
    org_db     = "org.Hs.eg.db",
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
      "--use-bitr"   = { defaults$use_bitr  <- toupper(val) %in% c("T", "TRUE", "1", "YES"); i <- i + 2 },
      "--from-type"  = { defaults$from_type <- val; i <- i + 2 },
      "--org-db"     = { defaults$org_db    <- val; i <- i + 2 },
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
  cat("用法: Rscript GEO_download_probe.R --gse GSE_ID[,GSE_ID2] [选项]\n\n")
  cat("必填:\n")
  cat("  --gse <ids>          GSE 号（多个用逗号分隔）\n\n")
  cat("可选:\n")
  cat("  --out <dir>          输出根目录（默认 .）\n")
  cat("  --use-idmap <T/F>    AnnoProbe::idmap() 在线注释（默认 FALSE）\n")
  cat("  --probe-col <name>   本地 GPL 探针 ID 列名（默认 ID）\n")
  cat("  --gene-col <name>    本地 GPL 基因列名（默认 symbol）\n")
  cat("  --diff <T/F>         生成 group/diff 模板（默认 FALSE）\n")
  cat("  --use-bitr <T/F>     clusterProfiler::bitr 二次转换（默认 FALSE）\n")
  cat("  --from-type <type>   bitr 输入 ID 类型（默认 ENTREZID）\n")
  cat("  --org-db <pkg>       物种注释包（默认 org.Hs.eg.db）\n")
  cat("  --proxy <url>        HTTP 代理 URL（默认 http://127.0.0.1:7897）\n")
  cat("  --timeout <sec>      下载超时（默认 600）\n")
  cat("  --help               显示帮助\n\n")
  cat("示例:\n")
  cat("  Rscript GEO_download_probe.R --gse GSE76262\n")
  cat("  Rscript GEO_download_probe.R --gse GSE76262,GSE12345 --out ./data --diff TRUE\n")
  cat("  Rscript GEO_download_probe.R --gse GSE8479 --use-bitr TRUE --gene-col GB_ACC --from-type REFSEQ\n")
  if (is.null(ARGS$gse) && !isTRUE(ARGS$help)) {
    cat("\n[错误] 缺少必填参数 --gse\n")
    quit(status = 1)
  }
  quit(status = 0)
}

# 解析多个 GSE
gse_ids <- strsplit(ARGS$gse, ",")[[1]] |> trimws() |> (\(x) x[nzchar(x)])()
use_idmap       <- isTRUE(ARGS$use_idmap)
use_bitr        <- isTRUE(ARGS$use_bitr)
probe_id_col    <- ARGS$probe_col
gene_symbol_col <- ARGS$gene_col
from_type       <- ARGS$from_type
org_db          <- ARGS$org_db
diff            <- if (isTRUE(ARGS$diff)) "TRUE" else "FALSE"
out_root        <- ARGS$out
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cat("========== CLI 参数 ==========\n")
cat("GSE 列表 :", paste(gse_ids, collapse = ", "), "\n")
cat("输出根目录:", normalizePath(out_root, mustWork = FALSE), "\n")
cat("use_idmap :", use_idmap, "\n")
cat("probe-col :", probe_id_col, "\n")
cat("gene-col  :", gene_symbol_col, "\n")
cat("use_bitr  :", use_bitr, "\n")
cat("from-type :", from_type, "\n")
cat("org-db    :", org_db, "\n")
cat("diff      :", diff, "\n")
cat("proxy     :", ifelse(nzchar(ARGS$proxy), ARGS$proxy, "(未设)"), "\n")
cat("timeout   :", ARGS$timeout, "\n\n")

# ---------- 加载依赖 ----------
suppressPackageStartupMessages({
  library(limma)
  library(affy)
  library(GEOquery)
  library(AnnoProbe)
  library(devtools)
  library(data.table)
  library(tidyverse)
})
options(timeout = ARGS$timeout)
options(download.file.method.GEOquery = "libcurl")
if (nzchar(ARGS$proxy)) {
  Sys.setenv(http_proxy = ARGS$proxy, https_proxy = ARGS$proxy)
  cat("[代理已设置]", ARGS$proxy, "\n\n")
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
# 单平台处理函数
# ============================================================
process_platform <- function(gse_obj, platform_index, n_platforms, output_dir, diff_flag) {

  gpl <- annotation(gse_obj)
  cat("\n========================================\n")
  cat("处理平台 [", platform_index, "/", n_platforms, "]: ", gpl, "\n", sep = "")
  cat("========================================\n")

  suffix <- if (n_platforms == 1) "" else paste0("_", gpl)
  exp_file_check <- file.path(output_dir, paste0("expMatrix", suffix, ".csv"))
  if (file.exists(exp_file_check)) {
    cat("[跳过平台] ", gpl, " expMatrix", suffix, ".csv 已存在\n", sep = "")
    return(invisible(NULL))
  }

  exp <- exprs(gse_obj)
  cat("表达矩阵维度:", dim(exp), "\n")
  cat("表达矩阵范围:", range(exp), "\n")

  cat("\n--- 保存临床数据 ---\n")
  clinical_data <- pData(gse_obj)
  cat("临床数据维度:", dim(clinical_data), "\n")
  clinical_file <- file.path(output_dir, paste0("clinical", suffix, ".csv"))
  write.csv(clinical_data, file = clinical_file, row.names = FALSE)
  cat("临床数据已保存至", clinical_file, "\n")

  cat("\n--- ID转换 ---\n")

  load_from_local_gpl <- function() {
    local_gpl_file <- {
      gpl_files <- list.files(
        path = output_dir,
        pattern = paste0("^", gpl, ".*\\.txt$"),
        full.names = TRUE
      )
      if (length(gpl_files) > 0) gpl_files[1] else ""
    }
    cat("切换到本地GPL文件:", local_gpl_file, "\n")
    if (!file.exists(local_gpl_file)) {
      cat("本地GPL文件不存在，尝试从NCBI下载...\n")
      gpl_url <- paste0(
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
        gpl, "&targ=self&form=text&view=data"
      )
      dest_file <- file.path(output_dir, paste0(gpl, ".txt"))
      download_ok <- tryCatch({
        download.file(gpl_url, destfile = dest_file, mode = "wb", quiet = TRUE)
        TRUE
      }, error = function(e_dl) {
        cat("从NCBI下载GPL文件失败:", conditionMessage(e_dl), "\n")
        FALSE
      })
      if (download_ok && file.exists(dest_file)) {
        local_gpl_file <- dest_file
        cat("已下载GPL文件至:", local_gpl_file, "\n")
      } else {
        stop("本地GPL文件不存在且下载失败: ", gpl)
      }
    }

    con <- file(local_gpl_file, "r")
    skip_lines <- 0L
    repeat {
      line <- readLines(con, n = 1)
      if (length(line) == 0) break
      if (!grepl("^[#!^]", line)) break
      skip_lines <- skip_lines + 1L
    }
    close(con)
    cat("检测到注释行数:", skip_lines, "\n")

    gpl_data <- as.data.frame(fread(
      local_gpl_file,
      sep = "\t",
      header = TRUE,
      fill = TRUE,
      quote = "",
      skip = skip_lines
    ))

    actual_probe_col <- find_col_case_insensitive(colnames(gpl_data), probe_id_col)
    actual_gene_col  <- find_col_case_insensitive(colnames(gpl_data), gene_symbol_col)
    cat("[列匹配] 探针ID列:", probe_id_col, "->", actual_probe_col, "\n")
    cat("[列匹配] 基因注释列:", gene_symbol_col, "->", actual_gene_col, "\n")
    gpl_data <- gpl_data %>% dplyr::select(all_of(c(actual_probe_col, actual_gene_col)))

    gpl_data %>%
      filter(!is.na(.data[[actual_gene_col]]) & .data[[actual_gene_col]] != "---") %>%
      mutate(
        symbol = map_chr(str_split(.data[[actual_gene_col]], " /// "), function(x) {
          if (grepl("gene_assignment", actual_gene_col, ignore.case = TRUE)) {
            parts <- str_split(x[1], " // ")[[1]]
            if (length(parts) >= 2) str_trim(parts[2]) else NA_character_
          } else {
            str_trim(x[1])
          }
        })
      ) %>%
      dplyr::filter(!is.na(symbol) & symbol != tolower(symbol)) %>%
      dplyr::select(all_of(actual_probe_col), symbol) %>%
      rename(ID = all_of(actual_probe_col))
  }

  probe2gene <- if (isTRUE(use_idmap)) {
    tryCatch({
      cat("尝试 idmap() type='pipe'...\n")
      df <- idmap(gpl, type = "pipe")
      cat("使用 idmap() type='pipe' 成功获取注释\n")
      df
    }, error = function(e_pipe) {
      cat("idmap() type='pipe' 调用失败:", conditionMessage(e_pipe), "\n")
      tryCatch({
        cat("尝试 idmap() type='soft'...\n")
        df <- idmap(gpl, type = "soft")
        cat("使用 idmap() type='soft' 成功获取注释\n")
        df
      }, error = function(e_soft) {
        cat("idmap() type='soft' 调用失败:", conditionMessage(e_soft), "\n")
        tryCatch({
          cat("尝试 idmap() type='bioc'...\n")
          df <- idmap(gpl, type = "bioc")
          cat("使用 idmap() type='bioc' 成功获取注释\n")
          df
        }, error = function(e_bioc) {
          cat("idmap() type='bioc' 调用失败:", conditionMessage(e_bioc), "\n")
          load_from_local_gpl()
        })
      })
    })
  } else {
    cat("已禁用 idmap()，直接使用本地GPL文件\n")
    load_from_local_gpl()
  }

  probe2gene <- probe2gene %>%
    mutate(symbol = str_split(symbol, " /// ") %>% map_chr(1)) %>%
    mutate(symbol = str_split(symbol, ",") %>% map_chr(1)) %>%
    dplyr::filter(!is.na(symbol) & symbol != "" & symbol != "---") %>%
    { if ("ID" %in% colnames(.)) rename(., probe_id = ID) else . } %>%
    mutate(probe_id = as.character(probe_id)) %>%
    dplyr::select(probe_id, symbol)

  # bitr 二次转换
  if (isTRUE(use_bitr)) {
    cat("\n--- 使用 clusterProfiler::bitr 转换基因ID ---\n")
    cat(sprintf("转换: %s -> SYMBOL (%s)\n", from_type, org_db))
    suppressMessages({
      library(clusterProfiler)
      library(org_db, character.only = TRUE)
    })

    probe2gene <- probe2gene %>%
      mutate(symbol = str_remove(symbol, "\\..*$"))

    gene_ids <- unique(probe2gene$symbol)
    cat("待转换基因数:", length(gene_ids), "\n")

    bitr_result <- clusterProfiler::bitr(
      gene_ids,
      fromType = from_type,
      toType   = "SYMBOL",
      OrgDb    = get(org_db)
    )

    probe2gene <- probe2gene %>%
      dplyr::rename(!!from_type := symbol) %>%
      dplyr::inner_join(bitr_result, by = from_type) %>%
      dplyr::select(probe_id, symbol = SYMBOL) %>%
      dplyr::filter(!is.na(symbol) & symbol != "") %>%
      dplyr::distinct(probe_id, symbol, .keep_all = TRUE)

    cat("转换后 probe2gene 维度:", dim(probe2gene), "\n")
  }

  probe_file <- file.path(output_dir, paste0("probe2gene", suffix, ".csv"))
  write.csv(probe2gene, file = probe_file, row.names = FALSE)
  cat("注释文件已保存至", probe_file, "\n")
  cat("注释文件维度:", dim(probe2gene), "\n")
  head(probe2gene)

  cat("\n--- 数据转换与ID映射 ---\n")

  exp_symbol <- exp %>%
    as.data.frame() %>%
    rownames_to_column(var = "probe_id") %>%
    inner_join(probe2gene, by = "probe_id") %>%
    na.omit() %>%
    dplyr::select(-probe_id) %>%
    relocate(symbol)

  cat("ID映射后维度:", dim(exp_symbol), "\n")
  cat("基因数量:", length(unique(exp_symbol$symbol)), "\n")
  head(exp_symbol, n = 1)

  cat("\n--- 过滤基因 ---\n")

  exp_filtered <- exp_symbol %>%
    mutate(rowmean = rowMeans(.[, 2:ncol(.)], na.rm = TRUE)) %>%
    arrange(desc(rowmean)) %>%
    distinct(symbol, .keep_all = TRUE) %>%
    filter(rowmean > 0) %>%
    dplyr::select(-rowmean) %>%
    filter(!is.na(symbol) & symbol != "" & symbol != "---") %>%
    column_to_rownames("symbol")

  cat("过滤后维度:", dim(exp_filtered), "\n")
  cat("剩余基因数量:", nrow(exp_filtered), "\n")

  cat("\n--- 检查数据范围 ---\n")
  cat("转换前数据范围:", range(exp_filtered), "\n")

  if (max(exp_filtered) > 20 && min(exp_filtered) >= 0) {
    cat("数据最大值 > 20 且无负值，进行log2转换\n")
    exp_filtered <- log2(exp_filtered + 1)
    cat("转换后数据范围:", range(exp_filtered), "\n")
  } else if (max(exp_filtered) > 20 && min(exp_filtered) < 0) {
    cat("数据存在负值，跳过log2转换\n")
  } else {
    cat("数据无需log转换\n")
  }

  cat("\n--- 保存表达矩阵 ---\n")
  exp_final <- exp_filtered %>%
    rownames_to_column(var = "symbol")
  exp_file <- file.path(output_dir, paste0("expMatrix", suffix, ".csv"))
  fwrite(exp_final, file = exp_file, sep = ",", quote = FALSE, row.names = FALSE)
  cat("表达矩阵已保存至", exp_file, "\n")
  cat("最终矩阵维度:", dim(exp_final), "\n")

  if (diff_flag == "TRUE") {
    cat("\n--- 生成分组文件 ---\n")
    group_file <- file.path(output_dir, paste0("group", suffix, ".txt"))
    diff_file <- file.path(output_dir, paste0("diff", suffix, ".txt"))

    if (file.exists(group_file) && file.exists(diff_file)) {
      cat("分组文件与差异比较文件已存在，跳过生成:", group_file, ";", diff_file, "\n")
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
}

# ============================================================
# 单个 GSE 处理包装函数
# ============================================================
process_one_gse <- function(gse_id, diff_flag) {
  cat("\n############################################################\n")
  cat("# 开始处理:", gse_id, "\n")
  cat("############################################################\n")

  output_dir <- file.path(out_root, gse_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  gse <- tryCatch({
    getGEO(gse_id, destdir = output_dir, getGPL = FALSE)
  }, error = function(e) {
    cat("[跳过] ", gse_id, " 下载失败: ", conditionMessage(e), "\n", sep = "")
    return(NULL)
  })

  if (is.null(gse) || length(gse) == 0) {
    cat("[跳过] ", gse_id, " 未获取到数据\n", sep = "")
    return(invisible(FALSE))
  }

  n_platforms <- length(gse)
  cat("检测到平台数量:", n_platforms, "\n")

  for (i in seq_along(gse)) {
    cat("\n--- 处理平台 [", i, "/", n_platforms, "] ---\n", sep = "")
    result <- tryCatch({
      process_platform(gse[[i]], i, n_platforms, output_dir, diff_flag)
      TRUE
    }, error = function(e) {
      cat("[跳过平台] ", gse_id, " 平台 ", i, " 失败: ",
          conditionMessage(e), "\n", sep = "")
      FALSE
    })
    if (!result) next
  }

  cat("\n[完成] ", gse_id, " 处理结束\n", sep = "")
  invisible(TRUE)
}

# ============================================================
# 第二步：批量遍历所有 GSE
# ============================================================
cat("\n========== 第二步：批量处理所有 GSE ==========\n")

success_list <- character(0)
fail_list    <- character(0)

for (gse_id in gse_ids) {
  result <- tryCatch({
    process_one_gse(gse_id, diff)
  }, error = function(e) {
    cat("[整体跳过] ", gse_id, " 顶层错误: ",
        conditionMessage(e), "\n", sep = "")
    FALSE
  })
  if (isTRUE(result)) {
    success_list <- c(success_list, gse_id)
  } else {
    fail_list <- c(fail_list, gse_id)
  }
}

# ============================================================
# 输出下载状态统计表
# ============================================================
cat("\n========== 输出下载状态统计表 ==========\n")
status_df <- data.frame(
  GSE_ID = gse_ids,
  Status = sapply(gse_ids, function(id) {
    if (file.exists(file.path(out_root, id, "expMatrix.csv"))) "已完成" else "未完成"
  }),
  stringsAsFactors = FALSE
)
status_file <- file.path(out_root, "download_status.txt")
write.table(status_df, file = status_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat("下载状态统计表已保存至", status_file, "\n")

success_list <- status_df$GSE_ID[status_df$Status == "已完成"]
fail_list    <- status_df$GSE_ID[status_df$Status == "未完成"]

cat("\n========== 全部任务完成 ==========\n")
cat("成功: ", length(success_list), " 个 -> ",
    if (length(success_list) > 0) paste(success_list, collapse = ", ") else "无",
    "\n", sep = "")
cat("失败/跳过: ", length(fail_list), " 个 -> ",
    if (length(fail_list) > 0) paste(fail_list, collapse = ", ") else "无",
    "\n", sep = "")
