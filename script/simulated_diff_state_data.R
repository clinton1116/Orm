# 加载必要的包
library(muscat)
library(SingleCellExperiment)

set.seed(123)

# 读取参考数据（请根据实际文件路径修改）
sce_ref <- readRDS("/Users/linkun/Single-Cell_Projects/distinct_manuscript-main/muscat simulation and Kang real data/SNAKEMAKES/muscat-simulations/data/raw_data/ref_kang.rds")

# 设置公共参数
ng <- nrow(sce_ref)    # 使用参考数据中的所有基因数
nc <- 7200             # 模拟 7200 个细胞（两组各 3600 个）
ns <- 5                # 2 个样本（每个样本对应一个组，非配对设计）
nk <- 1                # 只模拟一个 cluster
probs_val <- list(NULL, NULL, c(0.5, 0.5))  # 每组细胞比例均等
dd_val <- TRUE         # 开启差异分布模拟
paired_val <- FALSE
p_type <- 0
lfc_val <- 1           # 对于差异基因，均值差异模拟为2倍（2^1 = 2）
rel_lfc <- NULL
phylo_tree <- NULL
phylo_pars <- c(0, 3)
force_val <- TRUE

# 模拟 DE 数据集（差异表达，de）
sim_de <- simData(sce_ref,
                  ng = ng,
                  nc = nc,
                  ns = ns,
                  nk = nk,
                  probs = probs_val,
                  dd = dd_val,
                  p_dd = c(0.8, 0, 0.2, 0, 0, 0),  # 20% de
                  paired = paired_val,
                  p_type = p_type,
                  lfc = lfc_val,
                  rel_lfc = rel_lfc,
                  phylo_tree = phylo_tree,
                  phylo_pars = phylo_pars,
                  force = force_val)
rownames(sim_de)=metadata(sim_de)$gene_info$sim_gene
saveRDS(sim_de, file = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/de_sim_data.rds")

# 模拟 DP 数据集（差异比例及/或离散性，dp）
sim_dp <- simData(sce_ref,
                  ng = ng,
                  nc = nc,
                  ns = ns,
                  nk = nk,
                  probs = probs_val,
                  dd = dd_val,
                  p_dd = c(0.8, 0, 0, 0.2, 0, 0),  # 20% dp
                  paired = paired_val,
                  p_type = p_type,
                  lfc = lfc_val,
                  rel_lfc = rel_lfc,
                  phylo_tree = phylo_tree,
                  phylo_pars = phylo_pars,
                  force = force_val)
rownames(sim_dp)=metadata(sim_dp)$gene_info$sim_gene
saveRDS(sim_dp, file = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dp_sim_data.rds")

# 模拟 DM 数据集（差异模态，dm）
sim_dm <- simData(sce_ref,
                  ng = ng,
                  nc = nc,
                  ns = ns,
                  nk = nk,
                  probs = probs_val,
                  dd = dd_val,
                  p_dd = c(0.8, 0, 0, 0, 0.2, 0),  # 20% dm
                  paired = paired_val,
                  p_type = p_type,
                  lfc = lfc_val,
                  rel_lfc = rel_lfc,
                  phylo_tree = phylo_tree,
                  phylo_pars = phylo_pars,
                  force = force_val)
rownames(sim_dm)=metadata(sim_dm)$gene_info$sim_gene
saveRDS(sim_dm, file = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/dm_sim_data.rds")

# 模拟 DB 数据集（差异变异性，db）
sim_db <- simData(sce_ref,
                  ng = ng,
                  nc = nc,
                  ns = ns,
                  nk = nk,
                  probs = probs_val,
                  dd = dd_val,
                  p_dd = c(0.8, 0, 0, 0, 0, 0.2),  # 20% db
                  paired = paired_val,
                  p_type = p_type,
                  lfc = lfc_val,
                  rel_lfc = rel_lfc,
                  phylo_tree = phylo_tree,
                  phylo_pars = phylo_pars,
                  force = force_val)
rownames(sim_db)=metadata(sim_db)$gene_info$sim_gene
saveRDS(sim_db, file = "/Users/linkun/Single-Cell_Projects/Orm/data/simdata/db_sim_data.rds")