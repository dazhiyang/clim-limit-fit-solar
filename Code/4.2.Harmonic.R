################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Daily irradiance visualization and harmonic analysis/fitting example
# 逐日太阳辐射可视化与谐波分析/拟合示例
################################################################################

# Clear workspace and load necessary packages
# 清空工作区并加载必须的R包
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "ggplot2", "fitdistrplus", "mixsmsn", "geoTS", "ggthemes", "scales")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# Global Input and Path Settings
# 全局输入与路径设置
################################################################################
plot.size <- 8
line.size <- 0.2
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14

# Base directory for the project
# 项目基础目录
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
# Data and output directories
# 数据与输出目录
dir_data <- file.path(dir0, "Data/BSRNtxt")
dir_tex <- file.path(dir0, "tex")

################################################################################
# Load Data and Data Processing
# 数据加载与处理
################################################################################

# Get all station files from the processed BSRNtxt folder
# 从处理后的 BSRNtxt 文件夹获取所有站点文件
setwd(dir_data)
files <- list.files(pattern = "\\.txt$")

# Select a station for plotting (BOS station)
# 选择一个站点进行绘图 (BOS 站点)
stn_target <- "BOS"
file_target <- paste0(stn_target, ".txt")

if (!file.exists(file_target)) {
  stop(sprintf("Station file %s not found in %s", file_target, dir_data))
}

cat(sprintf("\nProcessing station: %s for harmonic visualization...\n", stn_target))
data <- tibble::as_tibble(read.table(file_target, header = TRUE, sep = "\t")) %>%
  mutate(Time = lubridate::ymd_hms(Time)) %>%
  filter(Z < 85) %>%
  filter(!is.na(kt) & !is.na(kb) & !is.na(kd)) %>%
  filter(!is.na(Gh) & !is.na(Dh) & !is.na(Bn)) %>%
  filter(Gh > 0 & Dh > 0 & Bn > 0)

# Aggregate data to daily means
# 将数据聚合为逐日平均值
daily <- data %>%
  mutate(Bh = Bn * cos(insol::radians(Z))) %>%
  mutate(Time = lubridate::round_date(Time, "day")) %>%
  group_by(Time) %>%
  summarise(across(everything(), mean, .names = "{.col}"), .groups = "drop")

################################################################################
# Harmonic Analysis and Fitting
# 谐波分析与拟合
################################################################################

variable <- c("Gh", "Bh", "Dh")
data.plot <- NULL
harm.feature <- numeric(length(variable))

for (i in seq_along(variable)) {
  irradiance <- dplyr::pull(daily, variable[i])
  # Store observation data
  # 存储观测数据
  data.plot <- data.plot %>%
    dplyr::bind_rows(tibble::tibble(Time = daily$Time, y = irradiance, quantity = variable[i], group = "Observation"))

  # Perform harmonic fitting using geoTS package
  # 使用 geoTS 包进行谐波拟合
  harmR <- geoTS::haRmonics(y = irradiance, method = "harmR", numFreq = 25, delta = 0.1)

  # Store fitted data
  # 存储拟合数据
  data.plot <- data.plot %>%
    dplyr::bind_rows(tibble::tibble(Time = daily$Time, y = harmR$fitted, quantity = variable[i], group = "Fitted"))

  # Extract the first amplitude as a feature
  # 提取第一振幅作为特征
  harm.feature[i] <- harmR$amplitude[1]
}

# Print a summary of extracted features
# 打印提取特征的简要说明
stn_id <- stn_target
print(tibble::tibble(stn = stn_id, G.harm = harm.feature[1], B.harm = harm.feature[2], D.harm = harm.feature[3]))

# Prepare labels for the plot (mathematical expressions)
# 准备绘图标签（数学公式）
data.plot$quantity <- factor(data.plot$quantity,
  levels = c("Gh", "Bh", "Dh"),
  labels = c("Daily~italic(G)[italic(h)]", "Daily~italic(B)[italic(h)]", "Daily~italic(D)[italic(h)]")
)

################################################################################
# Visualization
# 数据可视化
################################################################################

p1 <- ggplot2::ggplot() +
  ggplot2::geom_line(ggplot2::aes(x = Time, y = y),
    data = data.plot %>% dplyr::filter(group == "Observation"),
    linewidth = line.size / 2
  ) +
  ggplot2::geom_line(ggplot2::aes(x = Time, y = y),
    data = data.plot %>% dplyr::filter(group == "Fitted"),
    linewidth = line.size * 2,
    color = ggthemes::colorblind_pal()(8)[2]
  ) +
  ggplot2::facet_wrap(~quantity, labeller = ggplot2::label_parsed, nrow = 1) +
  ggplot2::scale_x_datetime(
    name = "Date",
    breaks = scales::date_breaks("1 year"),
    labels = scales::date_format("%Y"),
    expand = c(0.07, 0)
  ) +
  ggplot2::scale_y_continuous(
    name = expression(paste("Irradiance [W ", m^-2, "]")),
    expand = c(0, 0)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    plot.margin = unit(c(0.1, 0.2, 0, 0.0), "lines"),
    text = ggplot2::element_text(family = "Times", size = plot.size),
    axis.text = ggplot2::element_text(family = "Times", size = plot.size),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(colour = "gray80", linewidth = line.size / 2),
    panel.background = ggplot2::element_rect(fill = "transparent"),
    plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
    strip.text.x = ggplot2::element_text(family = "Times", size = plot.size, margin = unit(c(0.1, 0, 0.1, 0), "lines")),
    strip.text.y = ggplot2::element_text(family = "Times", size = plot.size, margin = unit(c(0, 0.1, 0, 0.1), "lines")),
    panel.spacing = unit(0, "lines"),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(family = "Times", size = plot.size),
    legend.background = ggplot2::element_rect(fill = "transparent"),
    legend.key = ggplot2::element_rect(fill = "transparent"),
    legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines")
  )

# Preview plot
# p1

# Save plot to PDF
# 保存图形为 PDF
if (!dir.exists(dir_tex)) dir.create(dir_tex, recursive = TRUE)
output_pdf <- file.path(dir_tex, "HAexample.pdf")
ggplot2::ggsave(filename = output_pdf, plot = p1, width = 160, height = 50, unit = "mm")

cat(sprintf("\nHarmonic analysis plot successfully saved to: %s\n", output_pdf))
