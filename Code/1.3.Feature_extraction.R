################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Extract distribution features from BSRN files (e.g. skew-normal mixtures)
# 从BSRN处理后的txt文件中提取分布特征（如偏态正态混合分布参数等）
################################################################################

# Clear workspace and load necessary packages
# 清空工作区并加载所需的库
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "mixsmsn", "bbmle", "sn", "geoTS", "foreach", "doSNOW")
invisible(lapply(libs, lapply, require, character.only = TRUE))

################################################################################
# Global input and functions
# 全局输入与函数
################################################################################

# Log-likelihood function for power and two-component skew-normal mixture
# 幂律分布与两个偏态正态混合分布的对数似然函数
power_sn2_loglik <- function(logit_p_power, logit_p_sn1,
                             log_a_power,
                             xi_sn1, log_omega_sn1, alpha_sn1,
                             xi_sn2, log_omega_sn2, alpha_sn2,
                             x) {
  # Transform parameters
  # 转换参数
  p_power <- plogis(logit_p_power)
  p_sn1 <- plogis(logit_p_sn1)
  a_power <- exp(log_a_power)
  omega_sn1 <- exp(log_omega_sn1)
  omega_sn2 <- exp(log_omega_sn2)

  # Ensure probabilities sum to 1
  # 确保概率之和为1
  p_sn2 <- 1 - p_power - p_sn1
  if (p_sn2 < 0.001 || p_sn2 > 0.999) {
    return(1e10) # Penalize invalid probabilities
    # 惩罚不合理的概率
  }

  # Calculate component densities
  # 计算组分密度函数
  power_dens <- a_power * x^(a_power - 1)
  sn1_dens <- sn::dsn(x, xi = xi_sn1, omega = omega_sn1, alpha = alpha_sn1)
  sn2_dens <- sn::dsn(x, xi = xi_sn2, omega = omega_sn2, alpha = alpha_sn2)

  # Mixture density
  # 混合密度
  mixture_dens <- p_power * power_dens + p_sn1 * sn1_dens + p_sn2 * sn2_dens

  # Handle very small values
  # 处理极小值防止log报错
  mixture_dens <- pmax(mixture_dens, .Machine$double.eps)

  # Return negative log-likelihood
  # 返回负对数似然
  -sum(log(mixture_dens))
}

# Improved initial parameter estimation using K-means
# 使用 K-means 的改进初始参数估计
get_smart_initial_params <- function(data) {
  # Estimate mixing proportions using k-means clustering
  # 使用k-means聚类估计混合比例
  km <- kmeans(data, centers = 3, nstart = 25)
  props <- table(km$cluster) / length(data)

  # Sort clusters by their means
  # 按均值排序聚类簇
  cluster_means <- tapply(data, km$cluster, mean)
  sorted_clusters <- order(cluster_means)

  # Assign components based on cluster characteristics:
  # Power distribution (typically has many small values)
  # 幂律分布通常包含许多较小的值
  power_cluster <- sorted_clusters[1]
  power_data <- data[km$cluster == power_cluster]

  # Skew-normal clusters
  # 偏态正态簇
  sn1_cluster <- sorted_clusters[2]
  sn2_cluster <- sorted_clusters[3]
  sn1_data <- data[km$cluster == sn1_cluster]
  sn2_data <- data[km$cluster == sn2_cluster]

  # Estimate power distribution parameter
  # 估计幂律分布参数
  a_power_est <- if (mean(power_data) > 0) {
    1 / (1 - mean(power_data)) # Rough estimate
    # 粗略估计
  } else {
    0.8
  }

  # Estimate skew-normal parameters using method of moments
  # 矩方法估计偏态正态参数
  estimate_sn_params <- function(x) {
    m <- mean(x)
    s <- sd(x)
    skew <- mean((x - m)^3) / s^3
    alpha <- sign(skew) * min(10, abs(skew) * 2) # Regularized estimate
    # 正则化估计
    list(xi = m, omega = s, alpha = alpha)
  }

  sn1_params <- estimate_sn_params(sn1_data)
  sn2_params <- estimate_sn_params(sn2_data)

  return(list(
    logit_p_power = qlogis(props[power_cluster]),
    logit_p_sn1 = qlogis(props[sn1_cluster]),
    log_a_power = log(pmax(0.1, pmin(5, a_power_est))),
    xi_sn1 = sn1_params$xi,
    log_omega_sn1 = log(pmax(0.05, sn1_params$omega)),
    alpha_sn1 = sn1_params$alpha,
    xi_sn2 = sn2_params$xi,
    log_omega_sn2 = log(pmax(0.05, sn2_params$omega)),
    alpha_sn2 = sn2_params$alpha
  ))
}

################################################################################
# Path Settings and Metadata Loading
# 路径设置和元数据读取
################################################################################

# Base directory
# 基础目录
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation/Data"

# Read in metadata for BSRN sites
# 读取BSRN站点地理信息
bsrn_info_file <- file.path(dir0, "BSRN information.csv")
bsrn <- read.csv(bsrn_info_file, header = TRUE)
bsrn$stn <- tolower(bsrn$stn)

# Set working directory to the processed txt files
# 设置工作区为处理后的txt数据目录
bsrn_txt_dir <- file.path(dir0, "BSRNtxt")
setwd(bsrn_txt_dir)
files <- list.files(pattern = "\\.txt$")

features <- tibble(
  stn = bsrn$stn,
  lat = bsrn$Latitude,
  lon = bsrn$Longitude,
  f1.Z.r = NA, # range for zenith angle
  f2.kt.h = NA, # the area of the high-kt component PDF
  f3.kt.l = NA, # the area of the low-kt component PDF
  f4.kb.h = NA, # the area of the high-kb component PDF
  f5.kb.l = NA, # the area of the low-kb component PDF
  f6.kd.h = NA, # the area of the high-kd component PDF
  f7.kd.l = NA, # the area of the low-kd component PDF
  f8.G.c0 = NA, # intercept for harmonic fit for Gh
  f9.G.r = NA, # range for harmonic fit for Gh
  f10.B.c0 = NA, # intercept for harmonic fit for Bh
  f11.B.r = NA, # range for harmonic fit for Bh
  f12.D.c2 = NA, # intercept for harmonic fit for Dh
  f13.D.r = NA # range for harmonic fit for Dh
)

# Extract station names accurately from available processed texts
# 从可用的处理后的文本中准确提取站点名称
stns <- substr(files, 1, 3)

# Setup parallel processing
# 设置并行计算
cat("\nStarting parallel feature extraction...\n")
n_cores <- parallel::detectCores() - 4
cl <- snow::makeCluster(n_cores)
doSNOW::registerDoSNOW(cl)

pb <- txtProgressBar(max = length(stns), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run feature extraction in parallel
# 并行提取特征
extracted_features_list <- foreach(stn = stns, .packages = libs, .options.snow = opts, .verbose = FALSE) %dopar% {
  # Try-catch to handle errors without stopping the loop
  # 尝试-捕获错误，以防止错误中断循环
  result <- tryCatch(
    {
      # Read in txt data efficiently, filter early, and convert Time
      # 高效读取txt数据，提前过滤数据并转换时间格式
      data <- read.table(paste0(stn, ".txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
        as_tibble() %>%
        filter(
          Z < 85,
          Gh > 0, Dh > 0, Bn > 0,
          Gh < 1361, Dh < 1361, Bn < 1361,
          kt > 0, kb > 0, kd > 0,
          kt < 1.2, kb < 0.9, kd < 0.7,
          !is.na(kt), !is.na(kb), !is.na(kd),
          !is.na(Gh), !is.na(Dh), !is.na(Bn),
          !is.nan(kt), !is.nan(kb), !is.nan(kd)
        ) %>%
        mutate(Time = lubridate::ymd_hms(Time))

      # Performance optimization: Subsample to 1 year of daytime equivalent (approx 262,800 points)
      # 性能优化：随机下采样至大约1年的白天等效数据量（约262,800个数据点）用于分布特征提取
      max_points <- 365 * 12 * 60
      if (nrow(data) > max_points) {
        data_sub <- data %>% sample_n(max_points)
      } else {
        data_sub <- data
      }

      # Zenith angle range
      # 天顶角范围
      f1_Z_r <- min(data$Z)

      # Density-based features
      # 密度分布特征

      # Use three-component skew-normal mixture for kt
      # 对kt使用三分量偏态正态混合分布
      kt <- data_sub$kt
      mixsn3 <- mixsmsn::smsn.mix(kt, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE)
      delta <- mixsn3$shape / (sqrt(1 + (mixsn3$shape)^2))
      mean_val <- mixsn3$mu + sqrt(mixsn3$sigma2) * delta * sqrt(2 / pi)
      f2_kt_h <- mixsn3$pii[order(mean_val)[3]]
      f3_kt_l <- mixsn3$pii[order(mean_val)[1]]

      # Use power-sn-sn mixture for kb
      # 对kb使用幂律分布与双生偏态正态混合分布
      kb <- data_sub$kb
      smart_params <- get_smart_initial_params(kb)
      mixpsn2 <- bbmle::mle2(
        power_sn2_loglik,
        start = smart_params,
        data = list(x = kb),
        method = "L-BFGS-B",
        lower = list(
          logit_p_power = -20, logit_p_sn1 = -20, log_a_power = log(0.001),
          xi_sn1 = min(kb) - 0.1, log_omega_sn1 = log(0.001), alpha_sn1 = -30,
          xi_sn2 = min(kb) - 0.1, log_omega_sn2 = log(0.001), alpha_sn2 = -30
        ),
        upper = list(
          logit_p_power = 20, logit_p_sn1 = 20, log_a_power = log(100),
          xi_sn1 = max(kb) + 0.1, log_omega_sn1 = log(10), alpha_sn1 = 30,
          xi_sn2 = max(kb) + 0.1, log_omega_sn2 = log(10), alpha_sn2 = 30
        ),
        control = list(maxit = 1000, trace = 0, factr = 1e9)
      )

      cf <- coef(mixpsn2)
      p_pwr <- plogis(as.numeric(cf["logit_p_power"]))
      p_sn1 <- plogis(as.numeric(cf["logit_p_sn1"]))

      mixpsn2_vec <- c(
        p_power = p_pwr,
        p_sn1 = p_sn1,
        p_sn2 = 1 - p_pwr - p_sn1,
        a_power = exp(as.numeric(cf["log_a_power"])),
        xi_sn1 = as.numeric(cf["xi_sn1"]),
        omega_sn1 = exp(as.numeric(cf["log_omega_sn1"])),
        alpha_sn1 = as.numeric(cf["alpha_sn1"]),
        xi_sn2 = as.numeric(cf["xi_sn2"]),
        omega_sn2 = exp(as.numeric(cf["log_omega_sn2"])),
        alpha_sn2 = as.numeric(cf["alpha_sn2"])
      )
      delta1 <- mixpsn2_vec["alpha_sn1"] / (sqrt(1 + (mixpsn2_vec["alpha_sn1"])^2))
      mean1 <- mixpsn2_vec["xi_sn1"] + mixpsn2_vec["omega_sn1"] * delta1 * sqrt(2 / pi)
      delta2 <- mixpsn2_vec["alpha_sn2"] / (sqrt(1 + (mixpsn2_vec["alpha_sn2"])^2))
      mean2 <- mixpsn2_vec["xi_sn2"] + mixpsn2_vec["omega_sn2"] * delta2 * sqrt(2 / pi)
      f4_kb_h <- mixpsn2_vec[order(c(mean1, mean2))[2] + 1]
      f5_kb_l <- mixpsn2_vec["p_power"]

      # Use three-component skew-normal mixture for kd
      # 对kd使用三分量偏态正态混合分布
      kd <- data_sub$kd
      mixsn3_d <- mixsmsn::smsn.mix(kd, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE)
      delta_d <- mixsn3_d$shape / (sqrt(1 + (mixsn3_d$shape)^2))
      mean_val_d <- mixsn3_d$mu + sqrt(mixsn3_d$sigma2) * delta_d * sqrt(2 / pi)
      f6_kd_h <- mixsn3_d$pii[order(mean_val_d)[3]]
      f7_kd_l <- mixsn3_d$pii[order(mean_val_d)[1]]

      # Harmonic-based features
      # 谐波分布特征
      daily <- data %>%
        mutate(Bh = Bn * cos(insol::radians(Z)), Time = lubridate::floor_date(Time, "day")) %>%
        group_by(Time) %>%
        summarise(Gh = mean(Gh), Bh = mean(Bh), Dh = mean(Dh), .groups = "drop")

      harmR_Gh <- geoTS::haRmonics(y = daily$Gh, method = "harmR", numFreq = 25, delta = 0.1)
      f8_G_c0 <- harmR_Gh$amplitude[1]
      f9_G_r <- max(harmR_Gh$fitted) - min(harmR_Gh$fitted)

      harmR_Bh <- geoTS::haRmonics(y = daily$Bh, method = "harmR", numFreq = 25, delta = 0.1)
      f10_B_c0 <- harmR_Bh$amplitude[1]
      f11_B_r <- max(harmR_Bh$fitted) - min(harmR_Bh$fitted)

      harmR_Dh <- geoTS::haRmonics(y = daily$Dh, method = "harmR", numFreq = 25, delta = 0.1)
      f12_D_c2 <- harmR_Dh$amplitude[1]
      f13_D_r <- max(harmR_Dh$fitted) - min(harmR_Dh$fitted)

      # Return data frame row
      # 返回包含特征的单行数据框
      tibble(
        stn = stn,
        f1.Z.r = f1_Z_r,
        f2.kt.h = as.numeric(f2_kt_h),
        f3.kt.l = as.numeric(f3_kt_l),
        f4.kb.h = as.numeric(f4_kb_h),
        f5.kb.l = as.numeric(f5_kb_l),
        f6.kd.h = as.numeric(f6_kd_h),
        f7.kd.l = as.numeric(f7_kd_l),
        f8.G.c0 = as.numeric(f8_G_c0),
        f9.G.r = as.numeric(f9_G_r),
        f10.B.c0 = as.numeric(f10_B_c0),
        f11.B.r = as.numeric(f11_B_r),
        f12.D.c2 = as.numeric(f12_D_c2),
        f13.D.r = as.numeric(f13_D_r)
      )
    },
    error = function(e) {
      cat(sprintf("\nError processing station %s: %s\n", stn, e$message))
      return(NULL)
    }
  )

  if (!is.null(result)) {
    return(result)
  } else {
    return(NULL)
  }
}

close(pb)
snow::stopCluster(cl)

# Combine the list of results into a single tibble
# 将提取的特征列表合并为单个数据框
extracted_features <- dplyr::bind_rows(extracted_features_list)

# Merge extracted features back to final structure
# 将提取的特征合并回最终结构
features <- features %>%
  dplyr::select(stn, lat, lon) %>%
  left_join(extracted_features, by = "stn")

# Output the extracted features to CSV
# 导出提取的所有特征至CSV文件
output_file <- file.path(dir0, "extraction_features.csv")
write.csv(features, output_file, row.names = FALSE)
cat(sprintf("\nExtraction complete! Data saved to: %s\n", output_file))
