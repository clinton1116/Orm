# 加载必要的包
library(SingleCellExperiment)
library(rms)
library(dplyr)

set.seed(123)

### 1. 定义 ormPval 函数
ormPval <- function(x){
  if (length(x$fail) && x$fail) {
    return(cat("Model Did Not Converge.  No summary provided."))
  }
  intercepts <- x$non.slopes < 10
  ns <- x$non.slopes
  cik <- attr(x$coef, 'intercepts')  # 针对 fit.mult.impute 对象
  if(length(cik) && intercepts) {
    warning('intercepts=TRUE not implemented for fit.mult.impute objects')
    intercepts <- FALSE
  }
  vv <- diag(vcov(x, intercepts = if(intercepts) 'all' else 'none'))
  if(!intercepts) {
    nints <- if(!length(cik)) ns else {
      if(length(cik) == 1 && cik == 0) 0 else length(cik)
    }
    ints.to.delete <- if(ns == 0 || nints == 0) integer(0) else 1:nints
    vv <- c(rep(NA, nints), vv)
  }
  cof <- x$coef
  if(!intercepts) {
    j <- -ints.to.delete
    cof <- cof[j]
    vv <- vv[j]
  }
  obj <- list(coef = cof, se = sqrt(vv),
              aux = NULL,
              auxname = 'Penalty Scale')
  
  beta <- obj$coef
  se <- obj$se
  Z <- beta / se
  P <- 1 - pchisq(Z^2, 1)
  return(P)
}

### 2. 定义 compute_pvals() 函数
# 该函数接收一个 SingleCellExperiment 对象 se_obj 和一个可选参数 cat_filter，
# 如果 cat_filter 非 NULL，则只计算 metadata 中 gene_info$category 为该类别的基因 p 值。
compute_pvals <- function(se_obj, cat_filter = NULL) {
  # 若数据未标准化，则使用 logNormCounts 进行标准化（内部覆盖）
  se_obj <- logNormCounts(se_obj)
  
  # 提取 logcounts 表达矩阵和细胞分组信息（group_id）
  expr <- assay(se_obj, "logcounts")
  group <- colData(se_obj)$group_id
  
  # 从 metadata 中获取 gene_info，并根据 cat_filter 取要处理的基因列表
  gene_info <- metadata(se_obj)$gene_info
  if (!is.null(cat_filter)) {
    genes_to_use <- gene_info$gene[gene_info$category == cat_filter]
  } else {
    genes_to_use <- rownames(expr)
  }
  
  # 过滤掉全 NA 的基因以及表达恒定（标准差为 0）的基因
  na_genes <- apply(expr, 1, function(x) all(is.na(x)))
  constant_genes <- apply(expr, 1, function(x) {
    if (all(is.na(x))) return(FALSE)
    return(sd(x, na.rm = TRUE) == 0)
  })
  
  # 对表达矩阵进行过滤
  filtered_expr <- expr[!na_genes & !constant_genes, , drop = FALSE]
  
  # 如果指定了基因类别，则仅保留过滤后中与该类别交集的基因
  if (!is.null(cat_filter)) {
    filtered_expr <- filtered_expr[intersect(rownames(filtered_expr), genes_to_use), , drop = FALSE]
  }
  
  # 对每个基因计算 p 值（利用 sapply 遍历行名）
  pvals <- sapply(rownames(filtered_expr), function(gene) {
    gene_expr <- filtered_expr[gene, ]
    # 构造数据框，包含该基因表达和分组信息
    dat <- data.frame(expr = gene_expr, group = factor(group))
    # 拟合模型：表达 ~ group，使用 orm()（来自 rms 包）
    dd <- datadist(dat)
    assign("dd", dd, envir = .GlobalEnv)
    options(datadist = "dd")
    fit <- orm(expr ~ group, data = dat, x = TRUE, y = TRUE)
    # 计算 p 值（ormPval 返回一个向量，通常第二个值对应 group 变量的检验）
    p <- ormPval(fit)
    if (length(p) > 1) {
      return(p[2])
    } else {
      return(p)
    }
  })
  
  # 整理结果数据框：gene 名称与对应的 p 值
  results <- data.frame(gene = names(pvals), pvalue = pvals, row.names = NULL)
  return(results)
}

### 3. 主流程：处理不同差异分布数据集并保存结果
# 假设以下 RDS 文件已存在，包含对应的 SingleCellExperiment 数据集
# de: 差异表达（de）数据集
# dp: 比例差异（dp）数据集
# dm: 差异变异（dm）数据集
# db: 其他类型（db）数据集

# 请根据实际文件路径加载数据集
de <- readRDS("/Users/linkun/Single-Cell_Projects/Orm/data/simdata/de_sim_data.rds")
dp <- readRDS("/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dp_sim_data.rds")
dm <- readRDS("/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dm_sim_data.rds")
db <- readRDS("/Users/linkun/Single-Cell_Projects/Orm/data/simdata/db_sim_data.rds")

# 分别计算各数据集对应类别（de, dp, dm, db）的 p 值
res_de <- compute_pvals(de, cat_filter = NULL)
res_dp <- compute_pvals(dp, cat_filter = NULL)
res_dm <- compute_pvals(dm, cat_filter = NULL)
res_db <- compute_pvals(db, cat_filter = NULL)

# 将结果保存到各自的 RDS 文件中
saveRDS(res_de, file = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/de_pval.rds")
saveRDS(res_dp, file = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/dp_pval.rds")
saveRDS(res_dm, file = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/dm_pval.rds")
saveRDS(res_db, file = "/Users/linkun/Single-Cell_Projects/Orm/data/P_val/db_pval.rds")