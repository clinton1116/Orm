# 完整脚本：基于 logcounts / cpm / Linnorm 三种归一化计算 ORM p 值并批量保存

# 1. 加载必要的包
library(SingleCellExperiment)  # assay(), SingleCellExperiment 类
library(rms)                   # orm(), datadist()
library(dplyr)                 # 管道及数据处理
library(BiocParallel)          # MulticoreParam(), bplapply()
library(scuttle)               # logNormCounts(), calculateCPM()
library(Linnorm)               # Linnorm()

set.seed(123)

# 2. 定义计算单基因 ORM p 值的辅助函数
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

# 3. 定义通用的 p 值计算函数
compute_pvals <- function(se_obj,
                          norm_method = c("logcounts", "cpm", "linnorm"),
                          ncores      = 4){
  norm_method <- match.arg(norm_method)
  # 3.1 更新 datadist，供 rms::orm 使用
  dd <<- datadist(data.frame(group = factor(colData(se_obj)$group_id)))
  options(datadist = "dd")
  
  # 3.2 归一化并选取 assay
  if(norm_method == "logcounts"){
    se_obj <- logNormCounts(se_obj, pseudo_count = 1)
    assay_name <- "logcounts"
  } else if(norm_method == "cpm"){
    assay(se_obj, "cpm") <- calculateCPM(se_obj)
    assay_name <- "cpm"
  } else if(norm_method == "linnorm"){
    cts <- as.matrix(assay(se_obj, "counts"))
    assay(se_obj, "linnorm") <- Linnorm(cts)
    assay_name <- "linnorm"
  }
  
  # 3.3 提取表达矩阵和分组信息
  expr <- assay(se_obj, assay_name)
  group <- colData(se_obj)$group_id
  
  # 3.4 基因过滤：至少 20 个非零细胞 & 方差 > 0.05
  keep <- (rowSums(expr > 0) >= 20) & (apply(expr, 1, var) > 0.05)
  expr <- expr[keep, , drop = FALSE]
  
  # 3.5 去掉全 NA 或恒定表达的基因
  na_genes  <- apply(expr, 1, function(x) all(is.na(x)))
  const_genes <- apply(expr, 1, function(x) sd(x, na.rm=TRUE) == 0)
  expr <- expr[!na_genes & !const_genes, , drop = FALSE]
  
  # 3.6 并行计算每个基因的 p 值
  BPPARAM <- MulticoreParam(workers = ncores)
  p_list <- bplapply(
    rownames(expr),
    function(gene){
      dat <- data.frame(expr = expr[gene, ], group = factor(group))
      fit <- orm(expr ~ group, data = dat, x = TRUE, y = TRUE)
      ormPval(fit)
    },
    BPPARAM = BPPARAM
  )
  pvals <- unlist(p_list)
  names(pvals) <- rownames(expr)
  
  # 3.7 返回结果数据框
  data.frame(
    gene   = names(pvals),
    pvalue = pvals,
    row.names = NULL
  )
}

# 4. 读取各模拟数据集
base_dir <- "/Users/linkun/Single-Cell_Projects/Orm/data"
de <- readRDS(file.path(base_dir, "simdata/de_sim_data.rds"))
dp <- readRDS(file.path(base_dir, "simdata/dp_sim_data.rds"))
dm <- readRDS(file.path(base_dir, "simdata/dm_sim_data.rds"))
db <- readRDS(file.path(base_dir, "simdata/db_sim_data.rds"))

datasets <- list(de = de, dp = dp, dm = dm, db = db)

# 5. 批量计算并保存所有组合的 p 值
out_dir <- file.path(base_dir, "P_val")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

message("All done!")