# ===============================================================
# 脚本功能：
#  对 de、dp、dm、db 四种 SCE 数据集，
#  用 scDD 计算差异分布 (DD) p 值，
#  提取全局校正后的 combined.pvalue.adj，
#  并保存 gene + p_adj.glb 两列到 RDS。
# ===============================================================

# 0. 如果没装过，请先：
#    BiocManager::install(c("scDD","scuttle","SingleCellExperiment"))
library(SingleCellExperiment)
library(scDD)         # DD 检测包
library(scuttle)      # 提供 sizeFactors() + logNormCounts()
library(dplyr)
library(BiocParallel)
# -----------------------------------------------------------------------------
# 对单个数据集跑 scDD 并保存结果
compute_and_save_pvals_scdd <- function(sce_file, out_file,
                                        condition_col = "group_id",
                                        min_nonzero   = 3,
                                        n_perms       = 0,
                                        testZeroes    = TRUE,
                                        BPPARAM       = bpparam()) {
  # 1. 读入
  sce <- readRDS(sce_file)
  
  # 2. 确保有 sizeFactors（scuttle::logNormCounts 也要它）
  if (is.null(sizeFactors(sce))) {
    # 计算
    sf <- scuttle::librarySizeFactors(sce)
    # 注入
    SingleCellExperiment::sizeFactors(sce) <- sf
  }
  
  # 3. 如果没有 normcounts assay，就手工补上
  if (!"normcounts" %in% assayNames(sce)) {
    assay(sce, "normcounts") <- sweep(counts(sce), 2, sf, "/")
  }
  
  # 4. 调用 scDD
  message("→ Running scDD on: ", sce_file)
  sce_dd <- scDD(
    SCdat        = sce,
    condition    = condition_col,
    permutations = 0,       # 0：用 KS 检验；>0 才做置换
    testZeroes   = testZeroes,
    param        = bpparam(),
    min.nonzero  = min_nonzero
  )
  
  # 5. 取结果表（存于 metadata(sce_dd)[[1]]）
  md <- results(sce_dd)
  if (!"combined.pvalue.adj" %in% colnames(md)) {
    stop("在 metadata(sce_dd)[[1]] 中找不到 combined.pvalue.adj 列，请检查 scDD 输出")
  }
  
  # 6. 拷出 gene + p_adj.glb
  out_df <- md %>%
    transmute(
      gene     = gene,
      p_adj.glb = combined.pvalue.adj
    )
  
  
  # 7. 存盘
  saveRDS(out_df, out_file)
  message("✓ Saved to ", out_file)
}

# -----------------------------------------------------------------------------
# 8. 四个数据集的输入 / 输出 路径（请按实际改）
dataset_files <- list(
  de = "/Users/linkun/Single-Cell_Projects/orm/data/simdata/de_sim_data.rds",
  dp = "/Users/linkun/Single-Cell_Projects/orm/data/simdata/dp_sim_data.rds",
  dm = "/Users/linkun/Single-Cell_Projects/orm/data/simdata/dm_sim_data.rds",
  db = "/Users/linkun/Single-Cell_Projects/orm/data/simdata/db_sim_data.rds"
)
output_files <- list(
  de = "/Users/linkun/Single-Cell_Projects/orm/data/P_val/scdd_de_pval.rds",
  dp = "/Users/linkun/Single-Cell_Projects/orm/data/P_val/scdd_dp_pval.rds",
  dm = "/Users/linkun/Single-Cell_Projects/orm/data/P_val/scdd_dm_pval.rds",
  db = "/Users/linkun/Single-Cell_Projects/orm/data/P_val/scdd_db_pval.rds"
)

# 9. 批量跑
for (ds in names(dataset_files)) {
  cat("\n====================\n处理数据集：", ds, "\n")
  compute_and_save_pvals_scdd(
    sce_file    = dataset_files[[ds]],
    out_file    = output_files[[ds]],
    condition_col = "group_id",
    min_nonzero = 20,               # 每个条件下至少 20 个非零细胞才检验
    n_perms     = 0,                # 0：快速 KS；如要置换检验可改成 e.g. 100
    testZeroes  = TRUE,
    BPPARAM     = BiocParallel::MulticoreParam(workers = 2)
  )
}