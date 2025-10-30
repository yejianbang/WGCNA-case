# ----------------------------------------------------------------------
# 完整的 WGCNA 算例运行、生物学合理性与稳健性分析脚本
# ----------------------------------------------------------------------

# 此版本使用数字标签进行 Jaccard 相似度计算，以解决颜色标签匹配中出现的持续 0 相似度问题。

# 1. 加载必要的库
library(WGCNA)
library(cluster)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# ----------------------------------------------------------------------
# 2. 自定义函数定义 (必须放在调用代码之前)
# ----------------------------------------------------------------------

# 自定义函数 1: subsample_WGCNA (用于在子集上运行 WGCNA)
subsample_WGCNA <- function(datExpr_orig, sample_fraction = 0.8, power_val = 6) {
    # 抽取样本子集
    n_samples_orig = nrow(datExpr_orig)
    sub_samples = sample(n_samples_orig, floor(n_samples_orig * sample_fraction))
    datExpr_sub = datExpr_orig[sub_samples, ]

    # 重新运行 WGCNA
    net_sub = blockwiseModules(
        datExpr_sub,
        power = power_val,
        TOMType = "unsigned",
        minModuleSize = 30,
        reassignThreshold = 0,
        mergeCutHeight = 0.25,
        numericLabels = TRUE,
        verbose = 0,
    )

    numeric_labels <- net_sub$colors
    names(numeric_labels) <- colnames(datExpr_sub)
    return(numeric_labels)
}

# 自定义函数 2: calculate_jaccard (基于数字标签的相似度计算)
calculate_jaccard <- function(main_labels, sub_labels) {
    # main_labels: 主网络的数字标签 (命名向量)
    # sub_labels: 子网络的数字标签 (命名向量)

    # 1. 确保两个标签向量对齐 (只包含共同基因)
    common_probes <- intersect(names(main_labels), names(sub_labels))

    if (length(common_probes) == 0) return(0)

    # 提取对齐后的标签向量
    M_labels_aligned <- main_labels[common_probes]
    S_labels_aligned <- sub_labels[common_probes]

    # 2. 识别主网络中的非灰色模块 (灰色标签在 WGCNA 中通常为 0)
    module_numbers <- unique(M_labels_aligned)
    module_numbers <- module_numbers[module_numbers != 0]

    jaccard_scores <- numeric(length(module_numbers))

    # 3. 遍历主网络中的每个模块 M
    for (m_idx in 1:length(module_numbers)) {
        mod_num <- module_numbers[m_idx]

        # 提取主模块 M 的基因集 (基于数字标签)
        genes_M <- names(M_labels_aligned[M_labels_aligned == mod_num])

        max_jaccard <- 0

        # 遍历子网络中的所有非灰色模块 K
        sub_module_numbers <- unique(S_labels_aligned)
        sub_module_numbers <- sub_module_numbers[sub_module_numbers != 0]

        for (sub_mod_num in sub_module_numbers) {

            # 提取子模块 K 的基因集 (基于数字标签)
            genes_K <- names(S_labels_aligned[S_labels_aligned == sub_mod_num])

            # 计算 Jaccard 相似系数: J = |A ∩ B| / |A ∪ B|
            intersection_size <- length(intersect(genes_M, genes_K))
            union_size <- length(union(genes_M, genes_K))

            j_score <- if (union_size > 0) intersection_size / union_size else 0

            if (j_score > max_jaccard) {
                max_jaccard <- j_score
            }
        }
        jaccard_scores[m_idx] <- max_jaccard
    }

    # 返回所有模块的最大重叠度的中位数
    return(median(jaccard_scores, na.rm = TRUE))
}


# ----------------------------------------------------------------------
# A. 模拟算例数据生成
# ----------------------------------------------------------------------

cat("--- A. 模拟算例数据生成 ---\n")

# 定义模拟数据的参数
n_samples <- 500
n_genes <- 5000
n_modules <- 5

# 模拟模块结构
set.seed(42)
ME_Trait_Seed <- rnorm(n_samples, mean = 0, sd = 1)
datExpr <- matrix(rnorm(n_samples * n_genes, mean = 7, sd = 2),
                  nrow = n_samples, ncol = n_genes)

for (i in 1:n_modules) {
    start_gene <- (i - 1) * 100 + 1
    end_gene <- start_gene + 99

    ME_effect <- if (i == 1) ME_Trait_Seed else rnorm(n_samples, mean = 0, sd = 1)

    for (j in start_gene:end_gene) {
        datExpr[, j] <- datExpr[, j] + 2 * ME_effect + rnorm(n_samples, mean = 0, sd = 0.5)
    }
}

rownames(datExpr) <- paste0("Sample_", 1:n_samples)
colnames(datExpr) <- paste0("Gene_", 1:n_genes)
datExpr <- as.data.frame(datExpr)

# 模拟性状数据
TraitA <- ME_Trait_Seed + rnorm(n_samples, mean = 0, sd = 0.5)
datTraits <- data.frame(SampleName = rownames(datExpr), TraitA = TraitA)
rownames(datTraits) <- rownames(datExpr)

# 保存和加载数据
file_name <- "Simulated_WGCNA_Data.RData"
save(datExpr, datTraits, file = file_name)
load(file = file_name)

cat(paste("模拟数据加载完成。样本数:", nrow(datExpr), "，基因数:", ncol(datExpr), "\n"))
start_time <- Sys.time()
print(paste("开始运行时间:", start_time))
# ----------------------------------------------------------------------
# B. WGCNA 运行算例
# ----------------------------------------------------------------------

cat("\n--- B. WGCNA 运行算例：网络构建与模块检测 ---\n")

gsg = goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 软阈值选择 (固定为 6 用于演示)
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
power = 6

# 网络构建和模块识别
net = blockwiseModules(
    datExpr,
    power = power,
    TOMType = "unsigned",
    minModuleSize = 30,
    reassignThreshold = 0,
    mergeCutHeight = 0.25,
    numericLabels = TRUE, # 必须是 TRUE 以保证 moduleLabels 为数字
    verbose = 0,
)

moduleLabels = net$colors # 主网络数字标签 (用于稳健性分析)
moduleColors = labels2colors(net$colors) # 主网络颜色标签 (用于性状关联展示)
MEs = net$MEs

print(paste("识别到的模块总数:", length(unique(moduleColors))))

# ----------------------------------------------------------------------
# C. 结果的生物学合理性评估 (模块-性状关联)
# ----------------------------------------------------------------------

cat("\n--- C. 生物学合理性评估：模块-性状关联 ---\n")

mouseTrait <- as.data.frame(datTraits$TraitA)
names(mouseTrait) <- "TraitA"

moduleTraitCor <- cor(MEs, mouseTrait, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

MEs_no_grey <- MEs[, names(MEs) != "MEgrey"]
p_values_no_grey <- corPvalueStudent(cor(MEs_no_grey, mouseTrait, use = "p"), nSamples)

min_p_index <- which.min(p_values_no_grey)
sigME <- names(p_values_no_grey)[min_p_index]
sigModuleColor <- gsub("ME", "", sigME)

cat(paste("性状 TraitA 关联最显著的模块是:", sigModuleColor, "\n"))
cat(paste("相关系数:", round(moduleTraitCor[min_p_index], 3), "\n"))
cat(paste("P 值:", format.pval(moduleTraitPvalue[min_p_index], digits = 2), "\n"))

cat("\n==> 生物学合理性结论: \n")
cat("如果关键模块与性状显著相关 (P < 0.05 且相关系数较高)，则表明结果具有初步合理性。\n")


# ----------------------------------------------------------------------
# D. 结果稳健性分析：子采样重现性 (使用 Jaccard 相似系数)
# ----------------------------------------------------------------------

cat("\n--- D. 结果稳健性分析：子采样重现性 (使用 Jaccard 相似系数) ---\n")

# 2. 执行多次子采样
N_SUBSAMPLES = 5
subsample_results = list()
all_probes <- colnames(datExpr) # 所有基因名称

for (i in 1:N_SUBSAMPLES) {
    # 修复 set.seed 错误：使用循环变量 i 设置种子
    set.seed(i * 100 + 42)

    cat(paste("运行第", i, "次子采样 (80% 样本)..."))

    # subsample_results 包含数字标签
    subsample_results[[i]] = subsample_WGCNA(datExpr, power_val = power)
    cat(" 完成。\n")
}

# 3. 评估稳健性
cat("\n评估模块重叠度 (稳健性指标):\n")

jaccard_medians <- numeric(N_SUBSAMPLES)

for (i in 1:N_SUBSAMPLES) {

    sub_labels <- subsample_results[[i]]

    # 修复 unused arguments 错误：使用位置参数调用修正后的函数
    median_jaccard <- calculate_jaccard(
        moduleLabels,      # main_labels (主网络数字标签)
        sub_labels         # sub_labels (子网络数字标签)
    )
    jaccard_medians[i] <- median_jaccard

    cat(paste0(
        "子采样 ", i, " vs 主网络：模块中位 Jaccard 相似度 = ",
        round(median_jaccard, 3), "\n"
    ))
}

final_robustness_score <- mean(jaccard_medians)

cat(paste0("\n==> 最终稳健性得分 (中位 Jaccard 相似度平均值): ", round(final_robustness_score, 3), "\n"))
cat("稳健性结论: 如果得分高于 0.7-0.8，则模块重现性好，结果稳健。\n")

# ----------------------------------------------------------------------
# E. 完成
# ----------------------------------------------------------------------

cat("\n完整的 WGCNA 算例运行和分析完成。\n")
end_time <- Sys.time()
print(paste("结束运行时间:", end_time))
run_time <- end_time - start_time
run_time_seconds <- as.numeric(run_time, units = "secs")
print(paste("总运行时长:", round(run_time_seconds, 2), "秒"))
