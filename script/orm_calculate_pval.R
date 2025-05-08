# 加载必要的包
library(SingleCellExperiment)
library(rms)
library(dplyr)
library(BiocParallel)
library(scuttle)
set.seed(123)
# 至少两个非 NA 的值
dd <- datadist(data.frame(group = factor(c("A","B"))))
options(datadist = "dd")
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
compute_pvals <- function(se_obj, norm_method = c("logcounts", "cpm", "linnorm"),cat_filter = NULL,ncores = 4) {
  # 若数据未标准化，则使用 logNormCounts 进行标准化（内部覆盖）
  norm_method <- match.arg(norm_method)
  if(norm_method == "logcounts"){
    assay_name <- "logcounts"
  } else if(norm_method == "cpm"){
    assay_name <- "cpm"
  } else if(norm_method == "linnorm"){
    assay_name <- "linnorm"
  }
  
  # 提取 logcounts 表达矩阵和细胞分组信息（group_id）
  expr <- assay(se_obj, assay_name)
  group <- colData(se_obj)$group_id
  
  # 从 metadata 中获取 gene_info，并根据 cat_filter 取要处理的基因列表
  gene_info <- metadata(se_obj)$gene_info
  keep1 <- (rowSums(expr>0)>=20) & (apply(expr,1,var) > 0.05)
  expr <- expr[keep1,]
  # 过滤掉全 NA 的基因以及表达恒定（标准差为 0）的基因
  na_genes <- apply(expr, 1, function(x) all(is.na(x)))
  constant_genes <- apply(expr, 1, function(x) {
    if (all(is.na(x))) return(FALSE)
    return(sd(x, na.rm = TRUE) == 0)
  })
  
  # 对表达矩阵进行过滤
  filtered_expr <- expr[!na_genes & !constant_genes, , drop = FALSE]
  BPPARAM <- MulticoreParam(workers = ncores)
  # 对每个基因计算 p 值（利用 sapply 遍历行名）
  pvals_list <- bplapply(rownames(filtered_expr), function(gene) {
    gene_expr <- filtered_expr[gene, ]
    dat <- data.frame(expr = gene_expr, group = factor(group))
    
    
    fit <- orm(expr ~ group, data = dat, x=TRUE, y=TRUE)
    p <- ormPval(fit)
    if (length(p)>1) p[2] else p
  }, BPPARAM = BPPARAM)
  
  pvals <- unlist(pvals_list)
  names(pvals) <- rownames(filtered_expr)
  
  data.frame(gene = names(pvals), pvalue = pvals, row.names = NULL)
}
  
  # 整理结果数据框：gene 名称与对应的 p 值


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
datasets <- list(de = de, dp = dp, dm = dm, db = db)



norm_methods <- c("logcounts", "cpm", "linnorm")
for(ds_name in names(datasets)){
  se_obj <- datasets[[ds_name]]
  for(nm in norm_methods){
    message("Processing ", ds_name, " with ", nm, "...")
    res <- compute_pvals(se_obj, norm_method = nm, ncores = 8)
    saveRDS(res,
            file = file.path(out_dir, paste0(ds_name, "_", nm, "_pval.rds")))
  }
}

#