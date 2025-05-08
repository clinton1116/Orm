#!/usr/bin/env Rscript
# ===============================================================
# 脚本功能：
# 对 de、dp、dm、db 四种 SingleCellExperiment 数据集进行差异状态检测，
# 针对三种归一化方法（logcounts、cpm、linnorm）自动构建设计矩阵、调用 distinct_test 计算 p 值，
# 最后只保留 gene 名称和全局校正后的 p 值（p_adj.glb），并分别保存为 RDS 文件。
# ===============================================================

# 加载必要包
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(distinct)
  library(scuttle)      # logNormCounts(), calculateCPM()
  library(Linnorm)      # Linnorm()
  library(dplyr)
  library(scran)
})

# 定义计算 p 值并保存结果的函数
compute_and_save_pvals <- function(se_path, out_path,
                                   norm_method = c("logcounts","cpm","linnorm"),
                                   design_formula = "~ group_id",
                                   column_to_test = 2,
                                   min_non_zero_cells = 20,
                                   n_cores = 2) {
  norm_method <- match.arg(norm_method)
  
  # 1. 读取数据
  se_obj <- readRDS(se_path)
  
  # 2. 标准化 —— 根据 norm_method 产生相应 assay
  if (norm_method == "logcounts") {
    assay_name <- "logcounts"
  } else if (norm_method == "cpm") {
    assay_name <- "cpm"
  } else if (norm_method == "linnorm") {
    assay_name <- "linnorm"
  }
  
  # 3. 构建设计矩阵
  exp_info <- se_obj@metadata$experiment_info
  design <- model.matrix(as.formula(design_formula), data = exp_info)
  rownames(design) <- exp_info$sample_id
  
  # 4. 运行 distinct_test
  cat(sprintf("[%s] Running distinct_test on '%s' assay...\n", norm_method, assay_name))
  res <- distinct_test(
    x                       = se_obj,
    name_assays_expression  = assay_name,
    name_cluster            = "cluster_id",
    name_sample             = "sample_id",
    design                  = design,
    column_to_test          = column_to_test,
    min_non_zero_cells      = min_non_zero_cells,
    n_cores                 = n_cores
  )
  
  # 5. 保留 gene 和全局校正后的 p 值
  res <- res[, c("gene", "p_adj.glb")]
  
  # 6. 保存结果
  saveRDS(res, file = out_path)
  cat(sprintf("[%s] Results saved to: %s\n\n", norm_method, out_path))
}

# -------------------------------
# 直接写死输出根目录，不再接收命令行参数
out_base <- "/Users/linkun/Single-Cell_Projects/Orm/data/P_val"

dataset_files <- list(
  de = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/de_sim_data.rds",
  dp = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dp_sim_data.rds",
  dm = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dm_sim_data.rds",
  db = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/db_sim_data.rds"
)

norm_methods <- c("logcounts","cpm","linnorm")

for (ds in names(dataset_files)) {
  se_path <- dataset_files[[ds]]
  for (nm in norm_methods) {
    out_file <- file.path(out_base,
                          sprintf("%s_%s_pval_distinct.rds", ds, nm)
    )
    cat(sprintf("Processing dataset '%s' with norm='%s'...\n", ds, nm))
    compute_and_save_pvals(
      se_path       = se_path,
      out_path      = out_file,
      norm_method   = nm,
      design_formula= "~ group_id",
      column_to_test= 2,
      min_non_zero_cells = 20,
      n_cores       = 4
    )
  }
}