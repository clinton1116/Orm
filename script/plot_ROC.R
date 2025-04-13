# -------------------------------
# 完整示例代码：统一处理 db、de、dm、dp 四种数据集并美化 ROC 曲线
# -------------------------------

# 加载必要包
library(SingleCellExperiment)
library(distinct)
library(pROC)
library(rms)
library(dplyr)

set.seed(123)

# 定义函数：预处理 ORM p 值的基因名称，删除 ".group" 后及其后所有内容
clean_gene_names <- function(genes) {
  sub("\\.group.*", "", genes)
}

# 定义函数：根据公共基因对 SingleCellExperiment 对象更新其 metadata$gene_info
update_gene_info <- function(se_obj, common_genes) {
  # 假设 metadata(se_obj)$gene_info 中包含 sim_gene 字段作为基因标识
  gene_info <- metadata(se_obj)$gene_info
  gene_info_sub <- gene_info[gene_info$sim_gene %in% common_genes, ]
  gene_info_sub <- gene_info_sub[match(common_genes, gene_info_sub$sim_gene), ]
  metadata(se_obj)$gene_info <- gene_info_sub
  return(se_obj)
}

# 定义主函数：对单个数据集计算 ROC 曲线及其他信息
# 输入参数：
#   type          : 数据集类别字符串，如 "db", "de", "dm", "dp"
#   dataset_path  : 数据集文件路径（RDS格式的 SingleCellExperiment 对象）
#   orm_path      : ORM 方法计算得到的 p 值文件路径
#   distinct_path : distinct 方法计算得到的 p 值文件路径
process_dataset <- function(type, dataset_path, orm_path, distinct_path) {
  cat("Processing", type, "dataset...\n")
  
  # 1. 读取数据并标准化表达矩阵
  se_obj <- readRDS(dataset_path)
  se_obj <- logNormCounts(se_obj)
  
  # 2. 读取 p 值结果
  orm_pval <- readRDS(orm_path)
  distinct_pval <- readRDS(distinct_path)
  
  # 3. 清理 ORM p 值中 gene 名称（删除“.group”及其后缀）
  pval_orm_genes <- clean_gene_names(orm_pval$gene)
  
  # 4. 取数据集行名与 p 值中的公共基因
  common_genes <- intersect(pval_orm_genes, rownames(se_obj))
  cat("  Number of common genes for", type, ":", length(common_genes), "\n")
  if (length(common_genes) == 0) {
    stop("No common genes found for dataset: ", type)
  }
  
  # 5. 子集化 SingleCellExperiment 对象
  se_sub <- se_obj[common_genes, ]
  
  # 6. 更新 metadata 中的 gene_info，使其 sim_gene 与 se_sub 的行名一致
  se_sub <- update_gene_info(se_sub, common_genes)
  
  # 7. 构建“真值”向量 truth：
  #    当 metadata(se_sub)$gene_info$category 等于当前类型时认为是真阳性
  truth <- metadata(se_sub)$gene_info$category == type
  
  # 8. 对 p 值数据进行子集化
  # 对 ORM 方法的 p 值
  pval_orm_vector <- orm_pval$pvalue
  names(pval_orm_vector) <- pval_orm_genes
  pval_orm_common <- pval_orm_vector[common_genes]
  
  # 对 distinct 方法的 p 值（假定其 gene 列已经为标准格式）
  pval_distinct_vector <- distinct_pval$p_adj.glb
  names(pval_distinct_vector) <- distinct_pval$gene
  pval_distinct_common <- pval_distinct_vector[common_genes]
  
  # 9. 构造结果数据框
  df <- data.frame(
    gene = common_genes,
    truth = truth,
    pvals_method_orm = pval_orm_common,
    pvals_method_distinct = pval_distinct_common,
    stringsAsFactors = FALSE
  )
  
  # 10. 构建 ROC 曲线
  # 由于 p 值越小越显著，而 pROC 默认认为 predictor 越大越倾向正类，因此这里取 p 值的负值作为 predictor
  roc_orm <- roc(response = df$truth, predictor = -df$pvals_method_orm, 
                 percent = TRUE, print.auc = TRUE, grid = TRUE, 
                 xlab = "False Positive Rate (%)", ylab = "True Positive Rate (%)", 
                 main = paste("ROC for", type, "Genes"))
  roc_distinct <- roc(response = df$truth, predictor = -df$pvals_method_distinct, 
                      percent = TRUE, print.auc = TRUE, grid = TRUE, add = TRUE)
  
  # 加入垂直参考线：例如在 FPR=5% (即0.05*100%)处
  abline(v = 5, col = "darkgray", lwd = 2, lty = 2)
  
  legend("bottomright", legend = c("ORM-based", "distinct-based"), 
         col = c("red", "blue"), lty = 1, cex = 0.8)
  
  # 11. 返回结果列表
  return(list(df = df, roc_orm = roc_orm, roc_distinct = roc_distinct))
}

# 定义各数据集的文件路径（请根据实际情况修改路径）
dataset_files <- list(
  db = list(
    dataset = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/db_sim_data.rds",
    orm     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_db_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_db_pval.rds"
  ),
  de = list(
    dataset = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/de_sim_data.rds",
    orm     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_de_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_de_pval.rds"
  ),
  dm = list(
    dataset = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dm_sim_data.rds",
    orm     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_dm_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_dm_pval.rds"
  ),
  dp = list(
    dataset = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dp_sim_data.rds",
    orm     = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/orm_dp_pval.rds",
    distinct = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_dp_pval.rds"
  )
)

# 处理所有数据集，并将结果存入列表
results_list <- list()
for (type in names(dataset_files)) {
  cat("----------------------------------------------------\n")
  cat("Processing dataset type:", type, "\n")
  res <- process_dataset(
    type,
    dataset_files[[type]]$dataset,
    dataset_files[[type]]$orm,
    dataset_files[[type]]$distinct
  )
  results_list[[type]] <- res
  cat("Finished processing", type, "\n")
}

# 绘图：将所有 ROC 曲线放在一个 2x2 的图中
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))
for (type in names(results_list)) {
  res <- results_list[[type]]
  # 绘制 ORM ROC 曲线，使用 pROC 的 print.auc 显示 AUC 值，单位设为百分比
  plot(res$roc_orm, col = "red", print.auc = TRUE, print.auc.cex = 0.8, 
       percent = TRUE, grid = TRUE,
       main = paste("ROC for", type, "Genes"), 
       xlab = "False Positive Rate (%)", ylab = "True Positive Rate (%)")
  plot(res$roc_distinct, col = "blue", print.auc = TRUE, print.auc.cex = 0.8, 
       percent = TRUE, add = TRUE)
  # 加入垂直参考线（例如，在 FPR=5%处）
  abline(v = 95, col = "darkgray", lwd = 2, lty = 2)
  legend("bottomright", legend = c("ORM-based", "distinct-based"),
         col = c("red", "blue"), lty = 1, cex = 0.8)
}
par(mfrow = c(1, 1))

# 输出各数据集的 AUC 值
for (type in names(results_list)) {
  auc_orm <- auc(results_list[[type]]$roc_orm)
  auc_distinct <- auc(results_list[[type]]$roc_distinct)
  cat("\n", type, "AUC ORM:", auc_orm, "; AUC distinct:", auc_distinct, "\n")
}