#!/usr/bin/env Rscript
#===============================================================================
# 文件名：run_seurat_DE_all.R
# 作用：对 db, de, dm, dp 四个模拟数据集并行做 Seurat 两组差异表达
# 要求：R 4.x+, Seurat, SingleCellExperiment, scuttle, dplyr, tibble
# 用法：Rscript run_seurat_DE_all.R
#===============================================================================

# —— 0. 环境准备 —— #
# 1) 安装并加载必要包
packages <- c("SingleCellExperiment","Seurat","SeuratObject","scuttle","dplyr","tibble","future")
invisible(lapply(packages, function(pkg){
  if (!requireNamespace(pkg,quietly=TRUE)) install.packages(pkg)
  library(pkg,character.only=TRUE)
}))

# 2) 并行设置（可选）
future::plan("multicore", workers = 4)

# 3) 保证结果可复现
set.seed(123)

# —— 1. 差异表达函数 —— #
compute_and_save_seurat_DE <- function(sce_path, out_file,
                                       group_col       = "group_id",
                                       assay_counts    = "counts",
                                       assay_data      = "logcounts",
                                       test.use        = "wilcox",
                                       min.pct         = 0.01,
                                       logfc.threshold = 0.1) {
  message("▶︎ 读取：", sce_path)
  sce <- readRDS(sce_path)
  sce <- logNormCounts(sce)
  
  message("  ↳ 转成 Seurat 对象并预处理")
  seu <- as.Seurat(sce, counts = assay_counts, data = assay_data)
  
  message("  ↳ 设置分组标签：", group_col)
  Idents(seu) <- seu@meta.data[[group_col]]
  grps <- levels(Idents(seu))
  if (length(grps)!=2)
    stop("group_id 必须恰好两类，用于两组比较。当前：", paste(grps, collapse=","))
  
  message("  ↳ 运行 FindMarkers: ", grps[1], " vs ", grps[2])
  markers <- FindMarkers(
    object         = seu,
    ident.1        = grps[1],
    ident.2        = grps[2],
    test.use       = test.use,
    min.pct        = min.pct,
    logfc.threshold= logfc.threshold,
    return.thresh   = 1,          # <- 返回所有基因的 p 值
    only.pos        = FALSE 
  )
  if (nrow(markers) == 0) {
    message("  ⚠️ FindMarkers 没有筛出任何基因，使用默认 p=1, log2FC=0 填充所有基因")
    all_genes <- rownames(seu)
    markers <- data.frame(
      p_val       = 1,
      p_val_adj   = 1,
      avg_log2FC  = 0,
      row.names   = all_genes,
      stringsAsFactors = FALSE
    )
  }
  
  # 整理输出
  out_df <- markers %>%
    rownames_to_column("gene") %>%
    transmute(gene, p_val, p_val_adj, avg_log2FC)
  
  message("  ↳ 保存到：", out_file)
  saveRDS(out_df, out_file)
}

# —— 2. 定义所有数据集路径 —— #
base_dir <- "/Users/linkun/Single-Cell_Projects/orm/data"
sim_dir  <- file.path(base_dir, "simdata")
pval_dir <- file.path(base_dir, "P_val")

# 四种类型及对应输入/输出
dataset_types <- c("db","de","dm","dp")
dataset_files <- setNames(
  file.path(sim_dir, paste0(dataset_types, "_sim_data.rds")),
  dataset_types
)
output_files  <- setNames(
  file.path(pval_dir, paste0("seurat_", dataset_types, "_pval.rds")),
  dataset_types
)

# —— 3. 批量并行处理 —— #
future.apply::future_lapply(dataset_types, function(ds) {
  compute_and_save_seurat_DE(
    sce_path        = dataset_files[[ds]],
    out_file        = output_files[[ds]],
    group_col       = "group_id",
    test.use        = "wilcox",
    min.pct         = 0,
    logfc.threshold = 0,
  )
})

# —— 4. 完成提示 & 会话信息 —— #
message("🎉 全部完成！四个结果保存在：\n", paste(output_files, collapse="\n"))
message("Session info:")
print(sessionInfo())