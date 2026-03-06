################################################################################
# This code is written by Zhiwen Wang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com
#
# Group BSRN stations and perform quality control for each group
# 将BSRN站点进行分组，再对分组后的站点分别进行质量控制
################################################################################

# Clear workspace and load packages
# 清理工作台加载包
rm(list = ls(all = TRUE))
library(dplyr)
libs <- c("foreach", "isotree", "data.table")
invisible(lapply(libs, require, character.only = TRUE))

options(digits = 5)

################################################################################
# Global variables and paths
# 全局变量与路径
################################################################################
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data")

################################################################################
# Extract data for each cluster
# 提取每个分类下的数据
################################################################################
setwd(dir_data)
cluster <- read.csv("cluster.csv", header = TRUE, sep = ",")
# cluster_names <- c("PSA", "SSC", "TMF", "MSU", "MWT", "HMP", "HIC")
cluster_names <- unique(cluster$Group)

setwd(file.path(dir_data, "BSRNtxt"))
groups <- split(cluster$stn, cluster$Group) # Group stations / 将站点按照聚类分组

# Iterate through each group, read corresponding txt files, and merge
# 遍历每个分组，读取对应 txt 文件并合并
group_data_list <- lapply(names(groups), function(g) {
  stns <- groups[[g]]
  # Read all station files for this group
  # 读取该组所有站点文件
  data_list <- lapply(stns, function(s) {
    file_name <- paste0(s, ".txt")
    if (file.exists(file_name)) {
      # Use fread to speed up the loading of txt files
      # 使用fread加速txt文件的读取
      df <- fread(file_name, header = TRUE, sep = "\t", data.table = FALSE)

      df <- df %>%
        mutate(Station = s) %>%
        dplyr::select("Time", "Gh", "Dh", "Bn", "Z", "ETR", "kt", "kd", "kb", "Station") %>%
        filter(
          Z < 90,
          !is.na(kt), !is.na(kb), !is.na(kd),
          !is.na(Gh), !is.na(Dh), !is.na(Bn),
          !is.nan(kt), !is.nan(kb), !is.nan(kd)
        ) %>%
        return(df)
    } else {
      return(NULL)
    }
  })
  group_df <- bind_rows(data_list)
  return(group_df)
})

# Data sampling via 1-degree Zenith Bins (Top 5% for Gh, Bn, Dh)
# 按照天顶角1度进行分组采样（分别保留Gh, Bn, Dh最高的前5%）
sampled_group_list <- lapply(group_data_list, function(df) {
  # Assign each point to an integer degree bin based on Zenith (Z)
  # 天顶角Z取整
  df_binned <- df %>% mutate(Z_bin = floor(Z))

  # Extract top 5% for each irradiance component
  # 分别提取Gh, Bn, Dh的前5%
  top_gh <- df_binned %>%
    group_by(Z_bin) %>%
    slice_max(order_by = Gh, prop = 0.05, with_ties = FALSE) %>%
    ungroup()
  top_bn <- df_binned %>%
    group_by(Z_bin) %>%
    slice_max(order_by = Bn, prop = 0.05, with_ties = FALSE) %>%
    ungroup()
  top_dh <- df_binned %>%
    group_by(Z_bin) %>%
    slice_max(order_by = Dh, prop = 0.05, with_ties = FALSE) %>%
    ungroup()

  # Add a further 10% random sampling per bin
  # 额外从每个bin中随机抽取10%
  set.seed(123)
  random_samp <- df_binned %>%
    group_by(Z_bin) %>%
    slice_sample(prop = 0.1) %>%
    ungroup()

  # Combine uniquely / 合并并去重
  sampled_df <- bind_rows(top_gh, top_bn, top_dh, random_samp) %>%
    distinct() %>%
    dplyr::select(-Z_bin) # Cleanup temporary columns / 删除临时列

  return(sampled_df)
})

# F1-score function
# F1-score函数
calculate_f1_score <- function(cal, data, st, weight_z_gt_50 = FALSE) {
  cal <- cal %>%
    mutate(
      T = ifelse(data <= st, 1, 0), # T: Predicted positive samples / 预测为正的样本
      F = ifelse(data > st, 1, 0), # F: Predicted negative samples / 预测为负的样本
      W = ifelse(rep(weight_z_gt_50, n()) & .data$Z > 50, 5, 1) # Weighting for Z > 50
    )
  total_weight <- sum(cal$W)

  # Calculate TP, FP, TN, FN with optional weighting
  # 计算 TP, FP, TN, FN
  metrics <- cal %>%
    summarise(
      TP = sum(.data$W * (.data$P == 1 & .data$T == 1)) / total_weight * 100,
      FP = sum(.data$W * (.data$P == 0 & .data$T == 1)) / total_weight * 100,
      TN = sum(.data$W * (.data$N == 1 & .data$F == 1)) / total_weight * 100,
      FN = sum(.data$W * (.data$N == 0 & .data$F == 1)) / total_weight * 100
    )

  evaluation <- metrics %>%
    mutate(
      precision = .data$TP / (.data$TP + .data$FP), # Precision / 精确率
      recall = .data$TP / (.data$TP + .data$FN) # Recall / 召回率
    )

  f1_score <- 2 * (evaluation$precision * evaluation$recall) / (evaluation$precision + evaluation$recall)
  return(f1_score)
}

# Result storage
# 结果存储
results_list <- list()

# Iterate through all clusters
# 遍历所有聚类
for (group_index in seq_along(names(groups))) {
  cat(sprintf("\nProcessing Group %d: %s...\n", group_index, cluster_names[group_index]))

  data.cal <- sampled_group_list[[group_index]] %>%
    as_tibble()
  if (nrow(data.cal) == 0) {
    cat("Warning: No data for group, skipping...\n")
    next
  }

  ################################################################################
  # Physically Possible Limit (PPL) Filter
  # 物理可能极限(PPL)过滤
  ################################################################################
  data.cal <- data.cal %>%
    mutate(
      cosZ = pmax(0, cos(Z * pi / 180))
    ) %>%
    filter(
      Gh >= -4 & Gh <= ETR * 1.5 * (cosZ^1.2) + 100,
      Dh >= -4 & Dh <= ETR * 0.95 * (cosZ^1.2) + 50,
      Bn >= -4 & Bn <= ETR
    ) %>%
    dplyr::select(-cosZ)

  if (nrow(data.cal) < 10) {
    cat("Warning: Insufficient data after PPL filter, skipping...\n")
    next
  }

  ################################################################################
  # Isolate anomalous values with Isolation Forest
  # 孤立森林分离异常值
  ################################################################################
  # Isolation Forest algorithm
  # 孤立森林算法
  iso_features <- data.cal %>% dplyr::select("kt", "kd", "kb") # Extract parameters / 提取总直散相关参数
  iso <- isotree::isolation.forest(iso_features, ndim = 3, ntrees = 100, nthreads = 3)
  pred <- predict(iso, iso_features) # Calculate anomaly score / 计算孤立森林异常得分

  # Append score
  # 评分
  cal <- data.cal %>% mutate(score = pred)

  # Calculate KDE of scores
  # 得分计算KDE
  kde_result <- density(cal$score)
  log_den <- log10(kde_result$y + 1e-10)

  # Find threshold
  # 找到孤立森林的阈值
  threshold.index <- max(which(log_den > 0))
  threshold <- kde_result$x[threshold.index + 1]

  # Extract P/N values
  # 孤立森林PN值提取
  cal <- cal %>%
    mutate(marker = ifelse(score > threshold, "N", "P")) %>%
    mutate(
      P = ifelse(marker == "P", 1, 0),
      N = ifelse(marker == "N", 1, 0)
    )

  ################################################################################
  # BNI Optimization
  # BNI 最优化
  ################################################################################
  max_Bn <- cal %>%
    mutate(Z_rounded = round(Z, 1)) %>%
    group_by(Z_rounded) %>%
    slice_max(order_by = Bn, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(max_Bn = Bn) %>%
    filter(!is.na(max_Bn))

  fitness_function_bni_init <- function(params) {
    a <- params[1]
    b <- params[2]
    c <- params[3]
    fitted_values <- a * max_Bn$ETR * (cos(max_Bn$Z_rounded / 180 * pi)^b) + c

    # Assign higher penalty weight to Zenith angles > 50
    # 给天顶角(Z > 50)赋予更大的误差惩罚权重
    weights <- ifelse(max_Bn$Z_rounded > 50, 5, 1)

    # Calculate weighted Mean Squared Error
    mse <- sum(weights * (max_Bn$max_Bn - fitted_values)^2, na.rm = TRUE) / sum(weights)
    return(mse)
  }

  initial_params <- c(a = 1, b = 1, c = 20)
  initial <- optim(
    par = initial_params, fn = fitness_function_bni_init, method = "L-BFGS-B",
    lower = c(0.1, 0.1, 0), upper = c(2, 0.5, 10)
  )

  fitness_function_bni <- function(params) {
    a <- params[1]
    b <- params[2]
    c <- params[3]
    bni <- cal %>%
      mutate(bn.lo = a * ETR * (cos(Z / 180 * pi)^b) + c) %>%
      dplyr::select(Z, Bn, bn.lo, P, N)
    f1_score <- calculate_f1_score(bni, bni$Bn, bni$bn.lo, weight_z_gt_50 = TRUE)
    return(-f1_score)
  }

  # Safe boundaries preventing lb > ub collisions
  lb <- c(max(0.1, initial$par[1] - 0.05), max(0.05, initial$par[2] - 0.05), max(0, initial$par[3] - 10))
  ub <- c(min(2, initial$par[1] + 0.05), min(0.5, initial$par[2] + 0.05), min(10, initial$par[3] + 10))


  result_bni <- optim(
    par = initial$par, fn = fitness_function_bni, method = "L-BFGS-B",
    lower = lb, upper = ub
  )

  best_params.bni <- result_bni$par
  best_f1_score.bni <- -result_bni$value

  ################################################################################
  # DHI Optimization
  # DHI 最优化
  ################################################################################
  max_Dh <- cal %>%
    filter(P == 1) %>%
    filter(Z > (min(cal$Z) + 10)) %>%
    mutate(Z_rounded = round(Z)) %>%
    group_by(Z_rounded) %>%
    slice_max(order_by = Dh, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(max_Dh = Dh) %>%
    filter(!is.na(max_Dh))

  fitness_function_dhi_init <- function(params) {
    a <- params[1]
    b <- params[2]
    c <- params[3]
    fitted_values <- a * max_Dh$ETR * (cos(max_Dh$Z_rounded / 180 * pi)^b) + c
    mse <- mean((max_Dh$max_Dh - fitted_values)^2, na.rm = TRUE)
    return(mse)
  }

  initial <- optim(
    par = initial_params, fn = fitness_function_dhi_init, method = "L-BFGS-B",
    lower = c(0.1, 0.1, 0), upper = c(2, 2, 30)
  )

  fitness_function_dhi <- function(params) {
    a <- params[1]
    b <- params[2]
    c <- params[3]
    dhi <- cal %>%
      mutate(dh.lo = a * ETR * (cos(Z / 180 * pi)^b) + c) %>%
      dplyr::select(Z, Dh, dh.lo, P, N)
    f1_score <- calculate_f1_score(dhi, dhi$Dh, dhi$dh.lo)
    return(-f1_score)
  }

  lb <- c(max(0.1, initial$par[1] - 1), max(0.1, initial$par[2] - 1), max(0, initial$par[3] - 10))
  ub <- c(min(2, initial$par[1] + 1), min(2, initial$par[2] + 1), min(30, initial$par[3] + 10))

  result_dhi <- optim(
    par = initial$par, fn = fitness_function_dhi, method = "L-BFGS-B",
    lower = lb, upper = ub
  )

  best_params.dhi <- result_dhi$par
  best_f1_score.dhi <- -result_dhi$value

  ################################################################################
  # GHI Optimization
  # GHI 最优化
  ################################################################################
  max_Gh <- cal %>%
    filter(P == 1) %>%
    filter(Z > (min(cal$Z) + 5)) %>%
    mutate(Z_rounded = round(Z, 1)) %>%
    group_by(Z_rounded) %>%
    slice_max(order_by = Gh, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    rename(max_Gh = Gh) %>%
    filter(!is.na(max_Gh))

  fitness_function_ghi_init <- function(params) {
    a <- params[1]
    b <- params[2]
    c <- params[3]
    fitted_values <- a * max_Gh$ETR * (cos(max_Gh$Z_rounded / 180 * pi)^b) + c
    mse <- mean((max_Gh$max_Gh - fitted_values)^2, na.rm = TRUE)
    return(mse)
  }

  initial <- optim(
    par = initial_params, fn = fitness_function_ghi_init, method = "L-BFGS-B",
    lower = c(0.1, 0.1, 0), upper = c(2, 2, 50)
  )

  fitness_function_ghi <- function(params) {
    a <- params[1]
    b <- params[2]
    c <- params[3]
    ghi <- cal %>%
      mutate(gh.lo = a * ETR * (cos(Z / 180 * pi)^b) + c) %>%
      dplyr::select(Z, Gh, gh.lo, P, N)
    f1_score <- calculate_f1_score(ghi, ghi$Gh, ghi$gh.lo)
    return(-f1_score)
  }

  lb <- c(max(0.1, initial$par[1] - 0.05), max(0.1, initial$par[2] - 0.05), max(0, initial$par[3] + best_params.dhi[3] - 10))
  ub <- c(min(2, initial$par[1] + 0.05), min(2, initial$par[2] + 0.05), min(50, initial$par[3] + best_params.dhi[3] + 10))

  result_ghi <- optim(
    par = initial$par, fn = fitness_function_ghi, method = "L-BFGS-B",
    lower = lb, upper = ub
  )

  best_params.ghi <- result_ghi$par
  best_f1_score.ghi <- -result_ghi$value

  # Compile row for this cluster
  # 编译当前聚类结果行
  cluster_row <- data.frame(
    Cluster_Name = cluster_names[group_index],
    Stations = paste(groups[[group_index]], collapse = ", "),
    gh_a = best_params.ghi[1],
    gh_b = best_params.ghi[2],
    gh_c = best_params.ghi[3],
    F1_gh = best_f1_score.ghi,
    bn_a = best_params.bni[1],
    bn_b = best_params.bni[2],
    bn_c = best_params.bni[3],
    F1_bn = best_f1_score.bni,
    dh_a = best_params.dhi[1],
    dh_b = best_params.dhi[2],
    dh_c = best_params.dhi[3],
    F1_dh = best_f1_score.dhi
  )

  results_list[[group_index]] <- cluster_row
}

# Write results to CSV
# 将结果写入CSV
results_df <- bind_rows(results_list)
output_csv <- file.path(dir_data, "cluster_parameter.csv")
write.csv(results_df, output_csv, row.names = FALSE)

cat("\nAll operations successfully completed. Results saved to:", output_csv, "\n")
