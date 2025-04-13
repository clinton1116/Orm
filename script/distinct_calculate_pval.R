# ===============================================================
# 脚本功能：
# 对 de、dp、dm、db 四种 SingleCellExperiment 数据集进行差异状态检测，
# 自动构建设计矩阵、调用 distinct_test 计算 p 值，
# 最后只保留 gene 名称和全局校正后的 p 值（p_adj.glb），并分别保存为 RDS 文件。
# ===============================================================

# 加载必要包
library(SingleCellExperiment)
library(distinct)
library(dplyr)

# 定义计算 p 值并保存结果的函数
compute_and_save_pvals <- function(file_path, out_file,
                                   design_formula = "~ group_id",
                                   column_to_test = 2,
                                   min_non_zero_cells = 20,
                                   n_cores = 2) {
  # 1. 读取数据
  se_obj <- readRDS(file_path)
  
  # 2. 标准化表达数据
  se_obj <- logNormCounts(se_obj)
  
  # 3. 提取实验信息（用于构建设计矩阵）
  exp_info <- se_obj@metadata$experiment_info
  # 注意：exp_info 应该至少包含 sample_id 和 group_id 两列
  # 构建设计矩阵，包含截距项，公式中 group_id 自动转换为哑变量
  design <- model.matrix(as.formula(design_formula), data = exp_info)
  # 设计矩阵的行名必须与 colData(se_obj)$sample_id 一致，因此这里设置为 exp_info$sample_id
  rownames(design) <- exp_info$sample_id
  
  # 4. 调用 distinct_test 进行差异状态检测
  #    注意：distinct_test 要求 x 为 SingleCellExperiment 对象，设计矩阵中行名必须与 colData(x)$sample_id 匹配
  cat("Running distinct_test...\n")
  res <- distinct_test(
    x = se_obj,
    name_assays_expression = "logcounts",
    name_cluster = "cluster_id",
    name_sample = "sample_id",
    design = design,
    column_to_test = 2,   # 例如：如果检验 group_idB 对应的效应，此处取第二列
    min_non_zero_cells = min_non_zero_cells,
    n_cores = n_cores
  )
  
  # 5. 结果处理：仅保留 gene 名称和全局校正后的 p 值（p_adj.glb）
  res <- res[, c("gene", "p_adj.glb")]
  
  # 6. 保存结果为 RDS 文件
  saveRDS(res, file = out_file)
  cat("Results saved to:", out_file, "\n")
}

# -------------------------------
# 主流程：对每种数据集依次计算 p 值并保存结果

# 定义每个数据集的文件路径（请根据实际情况修改路径）
dataset_files <- list(
  de = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/de_sim_data.rds",
  dp = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dp_sim_data.rds",
  dm = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dm_sim_data.rds",
  db = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/db_sim_data.rds"
)

# 定义输出结果的文件路径
output_files <- list(
  de = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_de_pval.rds",
  dp = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_dp_pval.rds",
  dm = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_dm_pval.rds",
  db = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/distinct_db_pval.rds"
)

# 循环处理每个数据集
for (ds in names(dataset_files)) {
  cat("Processing dataset:", ds, "\n")
  compute_and_save_pvals(file_path = dataset_files[[ds]],
                         out_file = output_files[[ds]])
  cat("Finished dataset:", ds, "\n\n")
}
