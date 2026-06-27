#!/usr/bin/env Rscript
# =============================================================================
# GEO_download_count.R - GEO 补充数据下载与预处理（CLI 版）
#
# 描述: 优先下载临床数据，再从 FTP 补充目录下载附件文件并解压。
#       仅识别并处理 *count* 文件作为表达矩阵，自动识别基因 ID 列，
#       使用 clusterProfiler::bitr 进行基因 ID 到 symbol 的转换。
#       支持人/小鼠/自定义物种切换。
# 作者: 科研木鱼（小红书：科研木鱼）
# 版本: 2.0  —— 命令行参数化 + 物种切换
#
# 使用方法（CLI）:
#   Rscript GEO_download_count.R --gse GSE152418
#   Rscript GEO_download_count.R --gse GSE266899 --species mouse --from-type ENSEMBL
#   Rscript GEO_download_count.R --gse GSE99999 --species custom --org-db org.Rn.eg.db
#   Rscript GEO_download_count.R --gse GSE99999 --keyword expression
#   Rscript GEO_download_count.R --help
#
# 参数:
#   --gse <ids>         必填，GSE 号（多个用逗号分隔）
#   --species <name>    human / mouse / custom，默认 human
#   --org-db <pkg>      custom 时必填，如 org.Rn.eg.db
#   --from-type <id>    输入基因 ID 类型，默认 ENSEMBL
#   --keyword <str>     文件名关键字，默认 count
#   --out <dir>         输出根目录，默认 "."
#   --diff <T/F>        是否生成 group/diff 模板，默认 FALSE
#   --proxy <url>       HTTP 代理 URL，默认 http://127.0.0.1:7897；空串关闭
#   --timeout <sec>     下载超时秒数，默认 600
#   --help              显示帮助
# =============================================================================

# ---------- CLI 参数解析 ----------
parse_cli_args <- function(args) {
  defaults <- list(
    gse        = NULL,
    species    = "human",
    org_db     = "",
    from_type  = "ENSEMBL",
    keyword    = "count",
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
      "--org-db"     = { defaults$org_db    <- val; i <- i + 2 },
      "--from-type"  = { defaults$from_type <- val; i <- i + 2 },
      "--keyword"    = { defaults$keyword   <- val; i <- i + 2 },
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
  cat("用法: Rscript GEO_download_count.R --gse GSE_ID[,GSE_ID2] [选项]\n\n")
  cat("必填:\n")
  cat("  --gse <ids>         GSE 号（多个用逗号分隔）\n\n")
  cat("可选:\n")
  cat("  --species <name>    human / mouse / custom（默认 human）\n")
  cat("  --org-db <pkg>      物种注释包（custom 时必填，如 org.Rn.eg.db）\n")
  cat("  --from-type <id>    输入基因 ID 类型（默认 ENSEMBL）\n")
  cat("  --keyword <str>     文件名关键字（默认 count）\n")
  cat("  --out <dir>         输出根目录（默认 .）\n")
  cat("  --diff <T/F>        生成 group/diff 模板（默认 FALSE）\n")
  cat("  --proxy <url>       HTTP 代理 URL（默认 http://127.0.0.1:7897）\n")
  cat("  --timeout <sec>     下载超时（默认 600）\n")
  cat("  --help              显示帮助\n\n")
  cat("示例:\n")
  cat("  Rscript GEO_download_count.R --gse GSE152418\n")
  cat("  Rscript GEO_download_count.R --gse GSE99999 --species mouse --from-type ENSEMBL\n")
  cat("  Rscript GEO_download_count.R --gse GSE99999 --species custom --org-db org.Rn.eg.db\n")
  if (is.null(ARGS$gse) && !isTRUE(ARGS$help)) {
    cat("\n[错误] 缺少必填参数 --gse\n")
    quit(status = 1)
  }
  quit(status = 0)
}

# 解析多个 GSE
gse_id_vec <- strsplit(ARGS$gse, ",")[[1]] |> trimws() |> (\(x) x[nzchar(x)])()
species    <- ARGS$species
from_type  <- ARGS$from_type
keyword    <- ARGS$keyword
out_root   <- ARGS$out
diff       <- if (isTRUE(ARGS$diff)) "TRUE" else "FALSE"
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

# 物种 → OrgDb 映射
org_db <- switch(
  species,
  human  = "org.Hs.eg.db",
  mouse  = "org.Mm.eg.db",
  custom = ARGS$org_db,
  stop("[错误] 未知 species: ", species, "（请用 human / mouse / custom）")
)
if (!nzchar(org_db)) {
  stop("[错误] species=custom 时必须通过 --org-db 指定 OrgDb 包名")
}
cat(sprintf("[物种] species=%s -> OrgDb=%s, fromType=%s\n", species, org_db, from_type))

cat("========== CLI 参数 ==========\n")
cat("GSE 列表  :", paste(gse_id_vec, collapse = ", "), "\n")
cat("输出根目录:", normalizePath(out_root, mustWork = FALSE), "\n")
cat("species   :", species, "\n")
cat("org_db    :", org_db, "\n")
cat("from_type :", from_type, "\n")
cat("keyword   :", keyword, "\n")
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
# 工具：解压根目录（某个 GSE 的输出目录）
# ============================================================
data_dir_for <- function(gse_id) file.path(out_root, gse_id)

# ============================================================
# 第一步：下载临床数据
# ============================================================
download_clinical <- function(gse_id) {
  data_dir <- data_dir_for(gse_id)
  if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
  cat("========== 第一步：下载临床数据 ==========\n")
  gse <- getGEO(gse_id, destdir = data_dir, getGPL = FALSE)
  clinical_data <- pData(gse[[1]])
  cat("临床数据维度:", dim(clinical_data), "\n")
  write.csv(clinical_data, file = file.path(data_dir, "clinical.csv"), row.names = FALSE)
  cat(sprintf("临床数据已保存至 %s\n", file.path(data_dir, "clinical.csv")))
  clinical_data
}

# ============================================================
# 第二步：下载 GEO 补充文件
# ============================================================
download_suppl <- function(gse_id) {
  data_dir <- data_dir_for(gse_id)
  cat("\n========== 第二步：下载GEO补充数据 ==========\n")

  gse_digits <- sub("^GSE", "", gse_id)
  n_digits   <- nchar(gse_digits)
  if (n_digits < 3) stop("GSE ID 数字部分至少 3 位")

  ftp_base <- sprintf(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE%snnn/%s/suppl/",
    substr(gse_digits, 1, n_digits - 3),
    gse_id
  )
  cat(sprintf("FTP 补充目录: %s\n", ftp_base))

  ftp_html <- tryCatch({
    paste(readLines(ftp_base, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  }, error = function(e) {
    cat(sprintf("无法访问目录: %s\n  错误: %s\n", ftp_base, conditionMessage(e)))
    ""
  })

  suppl_files <- character(0)
  if (nzchar(ftp_html)) {
    href_matches <- regmatches(ftp_html, gregexpr('href="[^"]+"', ftp_html))[[1]]
    suppl_files  <- gsub('^href="|"$', '', href_matches)
    suppl_files  <- suppl_files[!(suppl_files %in% c("", "/", "../"))]
    suppl_files  <- suppl_files[!grepl("://|^//", suppl_files)]
    suppl_files  <- suppl_files[!grepl("[?#]", suppl_files)]
    suppl_files  <- suppl_files[substr(suppl_files, 1, 1) != "/"]
    suppl_files  <- suppl_files[grepl("\\.[A-Za-z0-9]+$", suppl_files)]
  }

  if (length(suppl_files) == 0) {
    cat("未在补充目录中找到可下载文件，请检查 GSE ID 或网络连接\n")
    return(invisible(character(0)))
  }
  cat(sprintf("共发现 %d 个补充文件\n", length(suppl_files)))

  for (fname in suppl_files) {
    dest_path <- file.path(data_dir, fname)
    file_url  <- paste0(ftp_base, fname)

    if (file.exists(dest_path)) {
      cat(sprintf("已存在本地文件，跳过下载: %s\n", dest_path))
      next
    }

    cat(sprintf("下载: %s\n", file_url))
    download_ok <- tryCatch({
      download.file(file_url, destfile = dest_path, method = "libcurl", mode = "wb")
      TRUE
    }, error = function(e) {
      cat(sprintf("下载失败: %s\n  错误: %s\n", file_url, conditionMessage(e)))
      FALSE
    }, warning = function(w) {
      cat(sprintf("下载警告: %s\n  信息: %s\n", file_url, conditionMessage(w)))
      file.exists(dest_path)
    })

    if (download_ok && file.exists(dest_path)) {
      cat(sprintf("下载完成: %s\n", dest_path))
    } else {
      cat(sprintf("下载失败，请手动下载: %s\n  保存到: %s\n", file_url, dest_path))
    }
  }
  invisible(suppl_files)
}

# ============================================================
# 第三步：解压附件
# ============================================================
extract_archives <- function(gse_id) {
  data_dir <- data_dir_for(gse_id)
  cat("\n========== 第三步：解压附件 ==========\n")

  read_first_line <- function(path) {
    tryCatch({
      con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path) else file(path)
      on.exit(close(con), add = TRUE)
      line <- readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0 || is.na(line)) "(空文件)" else line[1]
    }, error = function(e) sprintf("[读取失败: %s]", conditionMessage(e)))
  }

  truncate_line <- function(line, max = 200) {
    if (is.na(line) || nchar(line) == 0) return(line)
    if (nchar(line) > max) paste0(substr(line, 1, max), "...") else line
  }

  archive_files <- list.files(
    data_dir,
    pattern    = "\\.(tar\\.gz|tgz|tar\\.bz2|tar|zip|gz)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  archive_files <- archive_files[
    !grepl("_series_matrix\\.txt\\.gz$", basename(archive_files), ignore.case = TRUE)
  ]

  if (length(archive_files) == 0) {
    cat("下载目录中未找到压缩文件\n")
    return(invisible(NULL))
  }
  cat(sprintf("共发现 %d 个压缩文件\n", length(archive_files)))

  purrr::walk(archive_files, function(f) {
    fname <- basename(f)
    cat(sprintf("\n--- %s ---\n", fname))

    if (grepl("\\.(tar\\.gz|tgz|tar\\.bz2|tar)$", fname, ignore.case = TRUE)) {
      ok <- tryCatch({ untar(f, exdir = data_dir); TRUE },
                     error = function(e) {
                       cat(sprintf("  [解压失败: %s]\n", conditionMessage(e)))
                       FALSE
                     })
      if (ok) {
        members <- tryCatch(untar(f, list = TRUE), error = function(e) character(0))
        if (length(members) == 0) {
          cat("  (空压缩包)\n")
        } else {
          cat(sprintf("  共 %d 个文件：\n", length(members)))
          for (rel in members) {
            ef <- file.path(data_dir, rel)
            first_line <- if (file.exists(ef)) read_first_line(ef) else "(未找到)"
            cat(sprintf("  - %s\n    第1行: %s\n", rel, truncate_line(first_line)))
          }
        }
      }
    } else if (grepl("\\.zip$", fname, ignore.case = TRUE)) {
      ok <- tryCatch({
        files_in_zip <- unzip(f, list = TRUE)$Name
        if (length(files_in_zip) > 0) unzip(f, files = files_in_zip, exdir = data_dir)
        TRUE
      }, error = function(e) {
        cat(sprintf("  [解压失败: %s]\n", conditionMessage(e)))
        FALSE
      })
      if (ok) {
        members <- tryCatch(unzip(f, list = TRUE)$Name, error = function(e) character(0))
        if (length(members) == 0) {
          cat("  (空压缩包)\n")
        } else {
          cat(sprintf("  共 %d 个文件：\n", length(members)))
          for (rel in members) {
            ef <- file.path(data_dir, rel)
            first_line <- if (file.exists(ef)) read_first_line(ef) else "(未找到)"
            cat(sprintf("  - %s\n    第1行: %s\n", rel, truncate_line(first_line)))
          }
        }
      }
    } else if (grepl("\\.gz$", fname, ignore.case = TRUE)) {
      original_name <- sub("\\.gz$", "", fname)
      out_path <- file.path(data_dir, original_name)
      if (file.exists(out_path)) {
        cat(sprintf("  解压目标已存在，跳过: %s\n", out_path))
      } else {
        ok <- tryCatch({
          in_con  <- gzfile(f, open = "rb")
          out_con <- file(out_path, open = "wb")
          on.exit({ close(in_con); close(out_con) }, add = TRUE)
          repeat {
            buf <- readBin(in_con, what = "raw", n = 1024 * 1024)
            if (length(buf) == 0) break
            writeBin(buf, out_con)
          }
          TRUE
        }, error = function(e) {
          cat(sprintf("  [解压失败: %s]\n", conditionMessage(e)))
          FALSE
        })
        if (ok) cat(sprintf("  已解压: %s\n", out_path))
      }
      first_line <- read_first_line(f)
      cat(sprintf("  第1行: %s\n", truncate_line(first_line)))
    }
  })
  invisible(NULL)
}

# ============================================================
# 第四步：读取 *keyword* 文件
# ============================================================
find_count_files <- function(gse_id) {
  data_dir <- data_dir_for(gse_id)
  cat("\n========== 第四步：读取 *", keyword, "* 文件 ==========\n", sep = "")

  all_data_files <- list.files(
    data_dir,
    pattern    = "\\.(tsv|csv|txt)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  count_files <- all_data_files[grepl(keyword, basename(all_data_files), ignore.case = TRUE)]
  exclude_patterns <- c("^clinical", "^expMatrix_", "^group", "^diff", "annot")
  for (pat in exclude_patterns) {
    count_files <- count_files[!grepl(pat, basename(count_files), ignore.case = TRUE)]
  }

  if (length(count_files) == 0) {
    cat(sprintf("下载目录中未找到包含 '%s' 的表达矩阵文件\n", keyword))
  } else {
    cat(sprintf("共发现 %d 个匹配文件\n", length(count_files)))
    for (f in count_files) cat("  -", basename(f), "\n")
  }
  count_files
}

# ============================================================
# 第五步：处理表达矩阵
# ============================================================
process_expr_matrices <- function(gse_id, count_files) {
  data_dir <- data_dir_for(gse_id)
  cat("\n========== 第五步：处理表达矩阵 ==========\n")

  suppressMessages({
    library(clusterProfiler)
    library(org_db, character.only = TRUE)
  })
  cat(sprintf("使用 clusterProfiler::bitr 转换基因ID (%s -> SYMBOL, %s)\n", from_type, org_db))

  for (input_path in count_files) {
    fname <- basename(input_path)
    cat(sprintf("正在处理: %s\n", fname))

    expr_data <- data.table::fread(
      input_path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      data.table = FALSE
    )
    expr_data <- expr_data %>%
      dplyr::rename(!!from_type := 1) %>%
      dplyr::mutate(!!from_type := as.character(.data[[from_type]]) %>% stringr::str_remove("\\..*$"))
    gene_ids <- expr_data[[1]] %>% unique()

    bitr_result <- clusterProfiler::bitr(
      gene_ids,
      fromType = from_type,
      toType   = "SYMBOL",
      OrgDb    = get(org_db)
    )

    expr_data <- expr_data %>%
      dplyr::left_join(bitr_result, by = from_type) %>%
      dplyr::relocate(symbol = SYMBOL, .before = everything()) %>%
      dplyr::select(-dplyr::all_of(from_type)) %>%
      dplyr::filter(!is.na(symbol), symbol != "") %>%
      dplyr::mutate(mean_expr = rowMeans(dplyr::select(., -symbol), na.rm = TRUE)) %>%
      dplyr::arrange(dplyr::desc(mean_expr)) %>%
      dplyr::distinct(symbol, .keep_all = TRUE) %>%
      dplyr::filter(mean_expr > 0) %>%
      dplyr::select(-mean_expr)

    cat("  表达矩阵维度:", dim(expr_data), "\n")
    cat("  表达矩阵范围:", range(expr_data[, -1], na.rm = TRUE), "\n")

    expr_type <- tools::file_path_sans_ext(basename(input_path))
    output_path <- file.path(data_dir, paste0("expMatrix_", expr_type, ".csv"))
    write.csv(expr_data, file = output_path, row.names = FALSE)
    cat(sprintf("已保存: %s\n", output_path))
  }
}

# ============================================================
# 第六步：生成分组文件
# ============================================================
write_group_diff <- function(gse_id, clinical_data) {
  data_dir <- data_dir_for(gse_id)
  if (exists("diff") && diff == "TRUE") {
    cat("\n========== 第六步：生成分组文件 ==========\n")
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
}

# ============================================================
# 单个 GSE 完整处理
# ============================================================
process_one_gse <- function(gse_id) {
  cat(sprintf("\n############################################\n"))
  cat(sprintf("# 处理 GSE: %s\n", gse_id))
  cat(sprintf("############################################\n"))

  clinical_data <- tryCatch(download_clinical(gse_id), error = function(e) {
    cat("[跳过] 临床数据下载失败:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(clinical_data)) return(invisible(FALSE))

  tryCatch(download_suppl(gse_id),     error = function(e) cat("[警告] 补充数据下载失败:", conditionMessage(e), "\n"))
  tryCatch(extract_archives(gse_id),   error = function(e) cat("[警告] 解压失败:", conditionMessage(e), "\n"))

  count_files <- find_count_files(gse_id)
  if (length(count_files) > 0) {
    process_expr_matrices(gse_id, count_files)
  }

  write_group_diff(gse_id, clinical_data)
  invisible(TRUE)
}

# ============================================================
# 主循环
# ============================================================
cat("========== GEO_download_count.R (CLI) ==========\n")
for (gse_id in gse_id_vec) {
  tryCatch(process_one_gse(gse_id), error = function(e) {
    cat("[整体跳过] ", gse_id, " 顶层错误: ", conditionMessage(e), "\n", sep = "")
  })
}
cat("\n全部 GSE 处理完成！\n")
