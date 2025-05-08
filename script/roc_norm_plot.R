#!/usr/bin/env Rscript

# =======================================================
# Script: plot_norm_method_rocs.R
# 功能：分别针对 ORM 和 distinct 两种方法
#      画三种归一化方法（logcounts, cpm, linnorm）的 ROC
#      每种方法一个 2×2 布局的 PDF
# =======================================================

# 1. 准备
library(SingleCellExperiment)  # 读取 SCE，提取真值
library(pROC)                  # ROC 曲线函数
set.seed(123)

# 2. 通用函数：对某个方法(method)和某个数据集(type)来算 ROC 列表
process_method_rocs <- function(type, sce_path, pval_dir, method) {
  # type: de/dp/dm/db
  # method: "orm" 或 "distinct"
  sce       <- readRDS(sce_path)
  gi        <- metadata(sce)$gene_info
  truth     <- gi$category == type
  genes     <- gi$sim_gene
  norms     <- c("logcounts", "cpm", "linnorm")
  rocs      <- list()
  
  for (nm in norms) {
    # 路径里会是 de_logcounts_pval_orm.rds, de_logcounts_pval_distinct.rds 之类
    f <- file.path(pval_dir, sprintf("%s_%s_pval_%s.rds", type, nm, method))
    df <- readRDS(f)
    pv <- df[[2]]
    names(pv) <- df[[1]]
    pv <- pv[genes]
    rocs[[nm]] <- roc(
      response  = truth,
      predictor = -pv,      # 越小的 p 值越强信号，取负号
      percent   = TRUE,
      quiet     = TRUE
    )
  }
  rocs
}

# 3. 配置全局路径、数据集列表
base_dir <- "/Users/linkun/Single-Cell_Projects/Orm/data"
pval_dir <- file.path(base_dir, "P_val")
sce_dir  <- file.path(base_dir, "simdata")

datasets <- list(
  de = file.path(sce_dir, "de_sim_data.rds"),
  dp = file.path(sce_dir, "dp_sim_data.rds"),
  dm = file.path(sce_dir, "dm_sim_data.rds"),
  db = file.path(sce_dir, "db_sim_data.rds")
)

# 4. 统一颜色
cols <- c(logcounts = "red", cpm = "blue", linnorm = "darkgreen")

# 5. 分别画 ORM 和 distinct
for (method in c("orm", "distinct")) {
  out_pdf <- file.path(pval_dir, sprintf("%s_norm_methods_ROC.pdf", method))
  pdf(out_pdf, width = 10, height = 10)
  par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
  
  for (ds in names(datasets)) {
    rocs <- process_method_rocs(ds, datasets[[ds]], pval_dir, method)
    # 先画 logcounts
    plot(rocs[["logcounts"]],
         col       = cols["logcounts"],
         print.auc = TRUE,
         percent   = TRUE,
         grid      = TRUE,
         main      = sprintf("%s ROC — %s", toupper(method), toupper(ds)),
         xlab      = "False Positive Rate (%)",
         ylab      = "True Positive Rate (%)"
    )
    # 叠加 cpm / linnorm
    plot(rocs[["cpm"]],     add = TRUE, col = cols["cpm"])
    plot(rocs[["linnorm"]], add = TRUE, col = cols["linnorm"])
    # 参考线 FPR=5%
    abline(v = 5, col = "darkgray", lty = 2, lwd = 2)
    legend("bottomright",
           legend = names(cols),
           col    = cols,
           lty    = 1, cex = 0.8)
  }
  
  dev.off()
  cat("✅ 已保存", toupper(method), "ROC 图到：", out_pdf, "\n")
}
