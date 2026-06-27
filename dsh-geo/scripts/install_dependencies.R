#!/usr/bin/env Rscript
# =============================================================================
# install_dependencies.R - DSH-GEO 技能依赖安装脚本（CLI 版）
#
# 描述: 一次性安装 dsh-geo 技能所需的全部 R 依赖（CRAN + Bioconductor + GitHub）。
#       支持按物种条件安装（human / mouse / all + 自定义 OrgDb）。
# 作者: 科研木鱼（闲鱼/小红书：科研木鱼）
# 版本: 2.0  —— 命令行参数化 + 按物种安装
#
# 使用方法（CLI）:
#   Rscript install_dependencies.R                          # 装人+小鼠（默认）
#   Rscript install_dependencies.R --species human          # 只装人
#   Rscript install_dependencies.R --species mouse          # 只装小鼠
#   Rscript install_dependencies.R --species all --org-db org.Rn.eg.db,org.Dr.eg.db
#   Rscript install_dependencies.R --help
#
# 参数:
#   --species <name>   human / mouse / all（默认 all = human+mouse）
#   --org-db <list>    额外 OrgDb，逗号分隔（如 org.Rn.eg.db,org.Dr.eg.db）
#   --help             显示帮助
# =============================================================================

# ---------- CLI 参数解析 ----------
parse_cli_args <- function(args) {
  defaults <- list(
    species = "all",
    org_db  = "",
    help    = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    switch(a,
      "--species" = { defaults$species <- val; i <- i + 2 },
      "--org-db"  = { defaults$org_db  <- val; i <- i + 2 },
      "--help"    = { defaults$help    <- TRUE; i <- i + 1 },
      "-h"        = { defaults$help    <- TRUE; i <- i + 1 },
      {
        cat("[警告] 未识别参数: ", a, "（已跳过）\n", sep = "")
        i <- i + 1
      }
    )
  }
  defaults
}

ARGS <- parse_cli_args(commandArgs(trailingOnly = TRUE))

if (isTRUE(ARGS$help)) {
  cat("用法: Rscript install_dependencies.R [选项]\n\n")
  cat("可选:\n")
  cat("  --species <name>   human / mouse / all（默认 all）\n")
  cat("  --org-db <list>    额外 OrgDb，逗号分隔（如 org.Rn.eg.db,org.Dr.eg.db）\n")
  cat("  --help             显示帮助\n\n")
  cat("示例:\n")
  cat("  Rscript install_dependencies.R\n")
  cat("  Rscript install_dependencies.R --species human\n")
  cat("  Rscript install_dependencies.R --species all --org-db org.Rn.eg.db,org.Dr.eg.db\n")
  quit(status = 0)
}

species <- ARGS$species
extra_org_db <- if (nzchar(ARGS$org_db)) {
  strsplit(ARGS$org_db, ",")[[1]] |> trimws() |> (\(x) x[nzchar(x)])()
} else character(0)

cat("========== install_dependencies.R (CLI) ==========\n")
cat("species :", species, "\n")
cat("extra_org_db:", if (length(extra_org_db) > 0) paste(extra_org_db, collapse = ", ") else "(无)", "\n\n")

# ============================================================
# 1. CRAN 包
# ============================================================
cran_pkgs <- c(
  "tidyverse",     # dplyr / ggplot2 / readr / ...
  "data.table",    # fread 大文件快速读取
  "archive",       # tar/zip 解压
  "devtools",      # 安装 GitHub 包
  "remotes",       # 同上
  "R.utils"        # gzip / bunzip2
)

cat("========== 1. 安装 CRAN 包 ==========\n")
cat("目标包:", paste(cran_pkgs, collapse = ", "), "\n\n")
for (pkg in cran_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  [已安装] %s\n", pkg))
  } else {
    cat(sprintf("  [安装中] %s ...\n", pkg))
    ok <- tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
      TRUE
    }, error = function(e) {
      cat(sprintf("    [错误] %s: %s\n", pkg, conditionMessage(e)))
      FALSE
    })
    cat(sprintf("  [%s] %s\n", if (ok) "完成" else "失败", pkg))
  }
}

# ============================================================
# 2. Bioconductor 包
# ============================================================
bioc_pkgs <- c(
  "limma",           # 差异分析 / 读取 agilent
  "affy",            # CEL 读取 + RMA
  "GEOquery",        # getGEO / getGEOSuppFiles
  "clusterProfiler", # bitr 基因 ID 转换
  "AnnotationDbi"    # OrgDb 依赖
)

cat("\n========== 2. 安装 Bioconductor 包 ==========\n")
cat("目标包:", paste(bioc_pkgs, collapse = ", "), "\n\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("  [安装中] BiocManager ...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/", quiet = TRUE)
}
for (pkg in bioc_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  [已安装] %s\n", pkg))
  } else {
    cat(sprintf("  [安装中] %s ...\n", pkg))
    ok <- tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
      TRUE
    }, error = function(e) {
      cat(sprintf("    [错误] %s: %s\n", pkg, conditionMessage(e)))
      FALSE
    })
    cat(sprintf("  [%s] %s\n", if (ok) "完成" else "失败", pkg))
  }
}

# ============================================================
# 3. 物种 OrgDb（按 species 决定）
# ============================================================
org_db_to_install <- c()
if (species %in% c("human", "all")) {
  org_db_to_install <- c(org_db_to_install, "org.Hs.eg.db")
}
if (species %in% c("mouse", "all")) {
  org_db_to_install <- c(org_db_to_install, "org.Mm.eg.db")
}
org_db_to_install <- c(org_db_to_install, extra_org_db)
org_db_to_install <- unique(org_db_to_install[org_db_to_install != ""])

if (length(org_db_to_install) > 0) {
  cat("\n========== 3. 安装物种 OrgDb ==========\n")
  cat("目标包:", paste(org_db_to_install, collapse = ", "), "\n\n")
  for (pkg in org_db_to_install) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(sprintf("  [已安装] %s\n", pkg))
    } else {
      cat(sprintf("  [安装中] %s ...\n", pkg))
      ok <- tryCatch({
        BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
        TRUE
      }, error = function(e) {
        cat(sprintf("    [错误] %s: %s\n", pkg, conditionMessage(e)))
        FALSE
      })
      cat(sprintf("  [%s] %s\n", if (ok) "完成" else "失败", pkg))
    }
  }
} else {
  cat("\n========== 3. 物种 OrgDb ==========\n")
  cat("  [跳过] species =", species, "，无 OrgDb 需要安装\n")
}

# ============================================================
# 4. GitHub 包
# ============================================================
cat("\n========== 4. 安装 GitHub 包 ==========\n")
gh_pkgs <- list(
  list(name = "AnnoProbe", repo = "jmzeng1314/AnnoProbe",
       desc = "GPL 在线注释 (idmap)")
)
for (pkg in gh_pkgs) {
  if (requireNamespace(pkg$name, quietly = TRUE)) {
    cat(sprintf("  [已安装] %s — %s\n", pkg$name, pkg$desc))
  } else {
    cat(sprintf("  [安装中] %s (%s) ...\n", pkg$name, pkg$repo))
    ok <- tryCatch({
      remotes::install_github(pkg$repo, upgrade = "never", quiet = TRUE)
      TRUE
    }, error = function(e) {
      cat(sprintf("    [错误] %s: %s\n", pkg$name, conditionMessage(e)))
      FALSE
    })
    cat(sprintf("  [%s] %s\n", if (ok) "完成" else "失败", pkg$name))
  }
}

# ============================================================
# 5. 安装结果汇总
# ============================================================
cat("\n========== 4. 安装结果汇总 ==========\n")
all_pkgs <- c(cran_pkgs, bioc_pkgs, org_db_to_install, sapply(gh_pkgs, `[[`, "name"))
result_table <- data.frame(
  Package = all_pkgs,
  Installed = sapply(all_pkgs, function(p) requireNamespace(p, quietly = TRUE)),
  stringsAsFactors = FALSE
)
ok_n  <- sum(result_table$Installed)
fail_n <- sum(!result_table$Installed)
cat(sprintf("共 %d 个包，已安装 %d 个，缺失 %d 个\n", length(all_pkgs), ok_n, fail_n))
if (fail_n > 0) {
  cat("\n缺失包:\n")
  print(result_table[!result_table$Installed, ])
  cat("\n请重试安装或手动安装失败包。\n")
} else {
  cat("\n所有依赖已就绪。可以使用 GEO_download_*.R 系列脚本。\n")
}
cat("========== 安装结束 ==========\n")
