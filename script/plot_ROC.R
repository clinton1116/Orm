# -------------------------------
# 完整示例代码：统一处理 db、de、dm、dp 四种数据集并美化 ROC 曲线（包含 scDD）
# -------------------------------

# 加载必要包
library(SingleCellExperiment)
library(distinct)
library(pROC)
library(rms)
library(dplyr)
library(scuttle)    # 提供 logNormCounts()

set.seed(123)

# 1. 清理 ORM p 值中的 gene 名称（删除 ".group" 及后缀）
clean_gene_names <- function(genes) {
  sub("\\.group.*", "", genes)
}

# 2. 根据公共基因对子对象的 metadata$gene_info 进行子集并重排
update_gene_info <- function(se_obj, common_genes) {
  gene_info <- metadata(se_obj)$gene_info
  gene_info_sub <- gene_info[gene_info$sim_gene %in% common_genes, ]
  gene_info_sub <- gene_info_sub[match(common_genes, gene_info_sub$sim_gene), ]
  metadata(se_obj)$gene_info <- gene_info_sub
  se_obj
}

# 3. 主函数：对单个数据集计算三种方法的 ROC 曲线
process_dataset <- function(type, dataset_path, orm_path, distinct_path, scdd_path) {
  cat("Processing", type, "dataset...\n")
  
  # 3.1 读入并归一化表达矩阵
  se_obj <- readRDS(dataset_path)
  se_obj <- logNormCounts(se_obj)
  
  # 3.2 读入三种方法的 p 值表
  orm_pval      <- readRDS(orm_path)
  distinct_pval <- readRDS(distinct_path)
  scdd_pval     <- readRDS(scdd_path)
  
  # 3.3 清理 ORM gene 名称
  pval_orm_genes <- clean_gene_names(orm_pval$gene)
  
  # 3.4 找到共同基因
  common_genes <- intersect(
    pval_orm_genes,
    rownames(se_obj)
  )
  cat("  Number of common genes for", type, ":", length(common_genes), "\n")
  if (length(common_genes)==0) {
    stop("No common genes found for dataset: ", type)
  }
  
  # 3.5 子集化 SingleCellExperiment
  se_sub <- se_obj[common_genes, ]
  se_sub <- update_gene_info(se_sub, common_genes)
  
  # 3.6 构造真值向量
  truth <- metadata(se_sub)$gene_info$category == type
  
  # 3.7 提取 ORM p 值子集
  p_orm_vec <- orm_pval$pvalue
  names(p_orm_vec) <- pval_orm_genes
  p_orm_common <- p_orm_vec[common_genes]
  
  # 3.8 提取 distinct p 值子集
  p_distinct_vec <- distinct_pval$p_adj.glb
  names(p_distinct_vec) <- distinct_pval$gene
  p_distinct_common <- p_distinct_vec[common_genes]
  
  # 3.9 提取 scDD p 值子集
  p_scdd_vec <- scdd_pval$p_adj.glb
  names(p_scdd_vec) <- scdd_pval$gene
  p_scdd_common <- p_scdd_vec[common_genes]
  
  # 3.10 构建结果数据框
  df <- data.frame(
    gene                  = common_genes,
    truth                 = truth,
    pvals_method_orm      = p_orm_common,
    pvals_method_distinct = p_distinct_common,
    pvals_method_scdd     = p_scdd_common,
    stringsAsFactors = FALSE
  )
  
  # 3.11 三条 ROC 曲线
  roc_orm <- roc(
    response  = df$truth,
    predictor = -df$pvals_method_orm,
    percent   = TRUE, print.auc = TRUE, grid = TRUE,
    xlab      = "False Positive Rate (%)",
    ylab      = "True Positive Rate (%)",
    main      = paste("ROC for", type, "Genes")
  )
  
  roc_distinct <- roc(
    response  = df$truth,
    predictor = -df$pvals_method_distinct,
    percent   = TRUE, print.auc = TRUE, grid = TRUE,
    add       = TRUE, col = "blue"
  )
  
  roc_scdd <- roc(
    response  = df$truth,
    predictor = -df$pvals_method_scdd,
    percent   = TRUE, print.auc = TRUE, grid = TRUE,
    add       = TRUE, col = "darkgreen"
  )
  
  # 3.12 画垂直参考线（FPR=5%）
  abline(v = 5, col = "darkgray", lwd = 2, lty = 2)
  
  # 3.13 图例
  legend(
    "bottomright",
    legend = c("ORM-based", "distinct-based", "scDD-based"),
    col    = c("red", "blue", "darkgreen"),
    lty    = 1, cex = 0.8
  )
  
  list(
    df            = df,
    roc_orm       = roc_orm,
    roc_distinct  = roc_distinct,
    roc_scdd      = roc_scdd
  )
}

# 4. 定义所有数据集路径（请酌情修改）
dataset_files <- list(
  db = list(
    dataset  = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/db_sim_data.rds",
    orm      = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_db_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_db_pval.rds",
    scdd     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/scdd_db_pval.rds"
  ),
  de = list(
    dataset  = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/de_sim_data.rds",
    orm      = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_de_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_de_pval.rds",
    scdd     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/scdd_de_pval.rds"
  ),
  dm = list(
    dataset  = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dm_sim_data.rds",
    orm      = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_dm_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_dm_pval.rds",
    scdd     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/scdd_dm_pval.rds"
  ),
  dp = list(
    dataset  = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dp_sim_data.rds",
    orm      = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_dp_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_dp_pval.rds",
    scdd     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/scdd_dp_pval.rds"
  )
)

# 5. 批量处理并保存结果
results_list <- list()
for (type in names(dataset_files)) {
  cat("----------------------------------------------------\n")
  cat("Processing dataset type:", type, "\n")
  res <- process_dataset(
    type,
    dataset_files[[type]]$dataset,
    dataset_files[[type]]$orm,
    dataset_files[[type]]$distinct,
    dataset_files[[type]]$scdd
  )
  results_list[[type]] <- res
  cat("Finished processing", type, "\n")
}

# 6. 把四张 ROC 放到 2x2 的画布里
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
for (type in names(results_list)) {
  res <- results_list[[type]]
  plot(
    res$roc_orm, col = "red",
    print.auc     = TRUE, print.auc.cex = 0.8,
    percent       = TRUE, grid = TRUE,
    main          = paste("ROC for", type, "Genes"),
    xlab          = "False Positive Rate (%)", 
    ylab          = "True Positive Rate (%)"
  )
  plot(res$roc_distinct, col = "blue",       add = TRUE)
  plot(res$roc_scdd,     col = "darkgreen",  add = TRUE)
  # FPR=5% 参考线（百分比坐标下就是 5）
  abline(v = 5, col = "darkgray", lwd = 2, lty = 2)
  legend(
    "bottomright",
    legend = c("ORM-based", "distinct-based", "scDD-based"),
    col    = c("red", "blue", "darkgreen"),
    lty    = 1, cex = 0.8
  )
}
par(mfrow = c(1, 1))

# 7. 输出三种方法的 AUC
for (type in names(results_list)) {
  roc_orm      <- results_list[[type]]$roc_orm
  roc_distinct <- results_list[[type]]$roc_distinct
  roc_scdd     <- results_list[[type]]$roc_scdd
  cat("\n", type,
      "AUC ORM:     ", auc(roc_orm),
      "; AUC distinct:", auc(roc_distinct),
      "; AUC scDD:    ", auc(roc_scdd), "\n")
}