################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Compare different distribution fitting results for BSRN data for one station
# 在一个站点比较BSRN数据的不同分布拟合结果
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "fitdistrplus", "mixsmsn", "mixtools", "sn", "bbmle", "geoTS", "ggplot2", "ggthemes", "scales", "ggpubr")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# global input
# 全局输入
################################################################################
plot.size <- 8
line.size <- 0.2
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14
# dir0 <- "/Users/dyang/Dropbox/Working papers/QC"
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
################################################################################

# Set working directory to processed txt files
# 设置工作区为处理后的txt数据目录
setwd(file.path(dir0, "Data/BSRNtxt"))

# Load and process data for BOS station
# 读取并处理BOS站点数据
data <- read.table("bos.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  mutate(Time = ymd_hms(Time)) %>%
  filter(
    Z < 85,
    Gh > 0, Dh > 0, Bn > 0,
    Gh < 1361, Dh < 1361, Bn < 1361,
    !is.na(kt), !is.na(kb), !is.na(kd),
    !is.na(Gh), !is.na(Dh), !is.na(Bn),
    !is.nan(kt), !is.nan(kb), !is.nan(kd)
  )

# Performance optimization: Subsample to 1 year of daytime equivalent (approx 262,800 points)
# 性能优化：随机下采样至大约1年的白天等效数据量（约262,800个数据点）
max_points <- 365 * 12 * 60
if (nrow(data) > max_points) {
  data <- data %>% sample_n(max_points)
}

################################################################################
# plot kt fits
# 拟合并绘制kt的分布
################################################################################

kt <- pull(data, "kt")
# two-component Gaussian mixture
mixn2 <- normalmixEM(kt, mu = c(0.2, 0.7), sigma = c(0.1, 0.1), k = 2)
# three-component Gaussian mixture
mixn3 <- normalmixEM(kt, mu = c(0.2, 0.4, 0.7), sigma = c(0.1, 0.2, 0.1), k = 3)
# three-component skew-normal mixture
mixsn3 <- smsn.mix(kt, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE)

# construct data for plotting
data.plot.kt <- NULL
threshold <- seq(0, 1, 0.005)
data.plot.kt <- data.plot.kt %>%
  bind_rows(., tibble(x = threshold, y = mixn2$lambda[1] * dnorm(threshold, mean = mixn2$mu[1], sd = mixn2$sigma[1]) + mixn2$lambda[2] * dnorm(threshold, mean = mixn2$mu[2], sd = mixn2$sigma[2]), group = "norm-norm mixture", quantity = "kt")) %>%
  bind_rows(., tibble(x = threshold, y = mixn3$lambda[1] * dnorm(threshold, mean = mixn3$mu[1], sd = mixn3$sigma[1]) + mixn3$lambda[2] * dnorm(threshold, mean = mixn3$mu[2], sd = mixn3$sigma[2]) + mixn3$lambda[3] * dnorm(threshold, mean = mixn3$mu[3], sd = mixn3$sigma[3]), group = "norm-norm-norm mixture", quantity = "kt")) %>%
  bind_rows(., tibble(x = seq(0, 1, 0.005), y = mixsmsn:::d.mixedSN(seq(0, 1, 0.005), mixsn3$pii, mixsn3$mu, mixsn3$sigma2, mixsn3$shape), group = "sn-sn-sn mixture", quantity = "kt"))

data.plot.kt$group <- factor(data.plot.kt$group, levels = c("norm-norm mixture", "norm-norm-norm mixture", "sn-sn-sn mixture"))

p1 <- ggplot() +
  geom_histogram(aes(x = kt, y = ..density..), data = data, bins = 60, fill = "grey80", color = "gray30", linewidth = line.size) +
  geom_line(aes(x = x, y = y, color = group), data = data.plot.kt, linewidth = line.size * 2, alpha = 0.8) +
  scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  scale_x_continuous(name = expression(paste("Clearness index, ", italic(k)[italic(t)], " [dimensionless]")), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste("Density [dimensionless]"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.2, 0, 0.2), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "inside", legend.position.inside = c(0.3, 0.7), legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "vertical")

p1

################################################################################
# plot kd fits
# 拟合并绘制kd的分布
################################################################################
kd <- pull(data, "kd")

# two-component Gaussian mixture
# 双峰高斯混合
mixn2.d <- normalmixEM(kd, mu = c(0.15, 0.45), sigma = c(0.1, 0.1), k = 2)
# three-component Gaussian mixture
# 三峰高斯混合
mixn3.d <- normalmixEM(kd, mu = c(0.15, 0.3, 0.45), sigma = c(0.1, 0.2, 0.1), k = 3)
# three-component skew-normal mixture
mixsn3.d <- smsn.mix(kd, nu = 3, g = 3, get.init = TRUE, criteria = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE)

# construct data for plotting
data.plot.kd <- NULL
threshold <- seq(0, 1, 0.005)
data.plot.kd <- data.plot.kd %>%
  bind_rows(., tibble(x = threshold, y = mixn2.d$lambda[1] * dnorm(threshold, mean = mixn2.d$mu[1], sd = mixn2.d$sigma[1]) + mixn2.d$lambda[2] * dnorm(threshold, mean = mixn2.d$mu[2], sd = mixn2.d$sigma[2]), group = "norm-norm mixture", quantity = "kd")) %>%
  bind_rows(., tibble(x = threshold, y = mixn3.d$lambda[1] * dnorm(threshold, mean = mixn3.d$mu[1], sd = mixn3.d$sigma[1]) + mixn3.d$lambda[2] * dnorm(threshold, mean = mixn3.d$mu[2], sd = mixn3.d$sigma[2]) + mixn3.d$lambda[3] * dnorm(threshold, mean = mixn3.d$mu[3], sd = mixn3.d$sigma[3]), group = "norm-norm-norm mixture", quantity = "kd")) %>%
  bind_rows(., tibble(x = seq(0, 1, 0.005), y = mixsmsn:::d.mixedSN(seq(0, 1, 0.005), mixsn3.d$pii, mixsn3.d$mu, mixsn3.d$sigma2, mixsn3.d$shape), group = "sn-sn-sn mixture", quantity = "kd"))

data.plot.kd$group <- factor(data.plot.kd$group, levels = c("norm-norm mixture", "norm-norm-norm mixture", "sn-sn-sn mixture"))

p2 <- ggplot() +
  geom_histogram(aes(x = kd, y = ..density..), data = data, bins = 60, fill = "grey80", color = "gray30", linewidth = line.size) +
  geom_line(aes(x = x, y = y, color = group), data = data.plot.kd, linewidth = line.size * 2, alpha = 0.8) +
  scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  scale_x_continuous(name = expression(paste("Diffuse transmittance, ", italic(k)[italic(d)], " [dimensionless]")), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste("Density [dimensionless]"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.2, 0, 0.2), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "inside", legend.position.inside = c(0.7, 0.7), legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "vertical")

p2

################################################################################
# plot kb fits
# 拟合并绘制kb的分布
################################################################################
kb <- pull(data, "kb")

# power and sn mixture
# 幂律分布与偏态正态混合分布的对数似然函数
power_sn_loglik <- function(logit_p_power, log_a_power,
                            xi_sn, log_omega_sn, alpha_sn,
                            x = kb) {
  # Transform parameters
  # 转换参数
  p_power <- plogis(logit_p_power)
  a_power <- exp(log_a_power)
  omega_sn <- exp(log_omega_sn)

  # Calculate component densities
  # 计算组分密度函数
  # Power distribution: f(x) = a * x^(a-1)
  # 幂律分布：f(x) = a * x^(a-1)
  power_dens <- a_power * x^(a_power - 1)

  # Skew-normal distribution (truncated to avoid extreme values)
  # 偏态正态分布（截断以避免极值）
  sn_dens <- sn::dsn(x, xi = xi_sn, omega = omega_sn, alpha = alpha_sn)

  # Mixture density
  # 混合密度
  mixture_dens <- p_power * power_dens + (1 - p_power) * sn_dens

  # Handle very small values and numerical issues
  # 处理极小值和数值问题
  mixture_dens <- pmax(mixture_dens, .Machine$double.eps)

  # Return negative log-likelihood
  -sum(log(mixture_dens))
}
initial_params <- list(
  logit_p_power = qlogis(0.4),
  log_a_power = log(0.8),
  xi_sn = 0.3,
  log_omega_sn = log(0.15),
  alpha_sn = 3
)
mixpsn <- mle2(
  power_sn_loglik,
  start = initial_params,
  method = "L-BFGS-B",
  lower = list(
    logit_p_power = -10, log_a_power = log(0.01),
    xi_sn = -1, log_omega_sn = log(0.01), alpha_sn = -10
  ),
  upper = list(
    logit_p_power = 10, log_a_power = log(10),
    xi_sn = 2, log_omega_sn = log(2), alpha_sn = 20
  ),
  control = list(maxit = 1000)
)
mixpsn <- c(
  p_power = plogis(as.numeric(coef(mixpsn)["logit_p_power"])),
  a_power = exp(as.numeric(coef(mixpsn)["log_a_power"])),
  xi_sn = as.numeric(coef(mixpsn)["xi_sn"]),
  omega_sn = exp(as.numeric(coef(mixpsn)["log_omega_sn"])),
  alpha_sn = as.numeric(coef(mixpsn)["alpha_sn"])
)

# power and beta mixture
# 幂律分布与贝塔混合分布的对数似然函数
power_beta_loglik <- function(logit_p_power, log_a_power,
                              log_alpha_beta, log_beta_beta,
                              x = kb) {
  # Transform parameters
  # 转换参数
  p_power <- plogis(logit_p_power)
  a_power <- exp(log_a_power)
  alpha_beta <- exp(log_alpha_beta)
  beta_beta <- exp(log_beta_beta)

  # Calculate component densities
  # 计算组分密度函数
  # Power distribution: f(x) = a * x^(a-1) for x in [0,1]
  # 幂律分布：x在[0,1]之间时，f(x) = a * x^(a-1)
  power_dens <- a_power * x^(a_power - 1)

  # Beta distribution
  # 贝塔分布
  beta_dens <- dbeta(x, alpha_beta, beta_beta)

  # Mixture density
  # 混合密度
  mixture_dens <- p_power * power_dens + (1 - p_power) * beta_dens

  # Handle very small values and numerical issues
  # 处理极小值和数值问题
  mixture_dens <- pmax(mixture_dens, .Machine$double.eps)

  # Return negative log-likelihood
  -sum(log(mixture_dens))
}
initial_params <- list(
  logit_p_power = qlogis(0.4),
  log_a_power = log(0.8),
  log_alpha_beta = log(2.5),
  log_beta_beta = log(2.0)
)
mixpbe <- mle2(
  power_beta_loglik,
  start = initial_params,
  method = "L-BFGS-B",
  lower = list(
    logit_p_power = -10, log_a_power = log(0.01),
    log_alpha_beta = log(0.01), log_beta_beta = log(0.01)
  ),
  upper = list(
    logit_p_power = 10, log_a_power = log(10),
    log_alpha_beta = log(20), log_beta_beta = log(20)
  ),
  control = list(maxit = 1000)
)
mixpbe <- c(
  p_power = plogis(as.numeric(coef(mixpbe)["logit_p_power"])),
  a_power = exp(as.numeric(coef(mixpbe)["log_a_power"])),
  alpha_beta = exp(as.numeric(coef(mixpbe)["log_alpha_beta"])),
  beta_beta = exp(as.numeric(coef(mixpbe)["log_beta_beta"]))
)

# Power and two-component sn mixture
# 幂律分布与两个偏态正态混合分布的对数似然函数
power_sn2_loglik <- function(logit_p_power, logit_p_sn1,
                             log_a_power,
                             xi_sn1, log_omega_sn1, alpha_sn1,
                             xi_sn2, log_omega_sn2, alpha_sn2,
                             x = kb) {
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
smart_params <- get_smart_initial_params(kb)
mixpsn2 <- mle2(
  power_sn2_loglik,
  start = smart_params,
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
mixpsn2 <- c(
  p_power = plogis(as.numeric(coef(mixpsn2)["logit_p_power"])),
  p_sn1 = plogis(as.numeric(coef(mixpsn2)["logit_p_sn1"])),
  p_sn2 = 1 - plogis(as.numeric(coef(mixpsn2)["logit_p_power"])) - plogis(as.numeric(coef(mixpsn2)["logit_p_sn1"])),
  a_power = exp(as.numeric(coef(mixpsn2)["log_a_power"])),
  xi_sn1 = as.numeric(coef(mixpsn2)["xi_sn1"]),
  omega_sn1 = exp(as.numeric(coef(mixpsn2)["log_omega_sn1"])),
  alpha_sn1 = as.numeric(coef(mixpsn2)["alpha_sn1"]),
  xi_sn2 = as.numeric(coef(mixpsn2)["xi_sn2"]),
  omega_sn2 = exp(as.numeric(coef(mixpsn2)["log_omega_sn2"])),
  alpha_sn2 = as.numeric(coef(mixpsn2)["alpha_sn2"])
)


# Construct data for plotting
# 构建用于绘图的数据
data.plot.kb <- NULL
threshold <- seq(0, 1, 0.005)
data.plot.kb <- data.plot.kb %>%
  bind_rows(., tibble(x = threshold, y = mixpbe["p_power"] * (mixpbe["a_power"] * threshold^(mixpbe["a_power"] - 1)) + (1 - mixpbe["p_power"]) * dbeta(threshold, shape1 = mixpbe["alpha_beta"], shape2 = mixpbe["beta_beta"]), group = "power-beta mixture", quantity = "kb")) %>%
  bind_rows(., tibble(x = threshold, y = mixpsn["p_power"] * (mixpsn["a_power"] * threshold^(mixpsn["a_power"] - 1)) + (1 - mixpsn["p_power"]) * sn::dsn(threshold, xi = mixpsn["xi_sn"], omega = mixpsn["omega_sn"], alpha = mixpsn["alpha_sn"]), group = "power-sn mixture", quantity = "kb")) %>%
  bind_rows(., tibble(x = threshold, y = mixpsn2["p_power"] * (mixpsn2["a_power"] * threshold^(mixpsn2["a_power"] - 1)) + mixpsn2["p_sn1"] * sn::dsn(threshold, xi = mixpsn2["xi_sn1"], omega = mixpsn2["omega_sn1"], alpha = mixpsn2["alpha_sn1"]) + mixpsn2["p_sn2"] * sn::dsn(threshold, xi = mixpsn2["xi_sn2"], omega = mixpsn2["omega_sn2"], alpha = mixpsn2["alpha_sn2"]), group = "power-sn-sn mixture", quantity = "kb"))


data.plot.kb$group <- factor(data.plot.kb$group, levels = c("power-beta mixture", "power-sn mixture", "power-sn-sn mixture"))

p3 <- ggplot() +
  geom_histogram(aes(x = kb, y = ..density..), data = data, bins = 60, fill = "grey80", color = "gray30", linewidth = line.size) +
  geom_line(aes(x = x, y = y, color = group), data = data.plot.kb, linewidth = line.size * 2, alpha = 0.8) +
  # scale_linetype_manual(values=c("dashed", "dashed", "solid")) +
  scale_color_manual(values = colorblind_pal()(8)[2:4]) +
  # facet_wrap(~quantity, labeller = label_parsed, nrow = 1) +
  scale_x_continuous(name = expression(paste("Beam transmittance, ", italic(k)[italic(b)], " [dimensionless]")), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(name = expression(paste("Density [dimensionless]")), limits = c(0, 5)) +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.2, 0, 0.2), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")), panel.spacing = unit(0, "lines"), legend.position = "inside", legend.position.inside = c(0.35, 0.7), legend.title = element_blank(), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines"), legend.direction = "vertical", legend.spacing.x = unit(-0.2, "lines"))

p3

p <- ggpubr::ggarrange(p1, p3, p2, ncol = 1, align = "h", labels = c("(a)", "(b)", "(c)"), heights = c(1, 1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

# Save plot to PDF
# 将图像保存至PDF文件
setwd(file.path(dir0, "tex"))
ggsave(filename = "CompareFit.pdf", plot = p, width = 80, height = 140, unit = "mm")
