################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Illustrate standard BSRN Extremely Rare Limits (ERL) for different components
# 展示不同辐射分量的 BSRN 标准“极罕见极限” (ERL)
################################################################################

# Clear workspace and load libraries
# 清空工作区并加载库
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "SolarData", "scattermore", "ggthemes", "insol")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# Global Input and Functions
# 全局输入与函数定义
################################################################################

# Function to calculate solar geometry (Zenith angle and ETR)
# 计算太阳几何参数（天顶角和地外辐射）的函数
calZen <- function(Tm, lat, lon, alt = 0) {
  jd <- insol::JD(Tm) # Julian day
  # 儒略日

  sunv <- insol::sunvector(jd, lat, lon, timezone = 0) # Solar vector
  # 太阳向量

  azi <- round(insol::sunpos(sunv)[, 1], 3) # Azimuth
  # 方位角
  zen <- round(insol::sunpos(sunv)[, 2], 3) # Zenith angle
  # 天顶角

  doy <- insol::daydoy(Tm) # Day of year
  # 积日
  da <- (2 * pi / 365) * (doy - 1) # Day angle
  # 日角

  # Eccentricity correction factor
  # 轨道偏心率修正系数
  re <- 1.000110 + 0.034221 * cos(da) + 0.001280 * sin(da) + 0.00719 * cos(2 * da) + 0.000077 * sin(2 * da)
  # Extraterrestrial direct normal irradiance
  # 地外法向直接辐射
  E0n <- round(1361.1 * re, 3)

  solpos <- list(zen, azi, E0n)
  names(solpos) <- c("zenith", "azimuth", "E0n")
  return(solpos)
}

# Path Settings
# 路径设置
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data")
dir_tex <- file.path(dir0, "tex")

# Station details for plotting (e.g., QIQ)
# 绘图使用的站点信息
stn <- "QIQ"
lat <- 47.7957
lon <- 124.4852
alt <- 170

# Plotting parameters
# 绘图参数
plot.size <- 8
line.size <- 0.2
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14

################################################################################
# Load BSRN Data
# 加载 BSRN 数据
################################################################################

# Set working directory to the station's raw data folder
# 设置工作目录至站点的原始数据文件夹
setwd(file.path(dir_data, "BSRN", tolower(stn)))
files <- list.files(pattern = "\\.gz$")
# Sort files chronologically by month
# 按照月份顺序排序文件
files <- files[order(paste0(substr(files, 6, 7), substr(files, 4, 5)))]

cat(sprintf("\nLoading data for station %s...\n", stn))
data_list <- list()
pb <- txtProgressBar(max = length(files), style = 3)

for (i in seq_along(files)) {
  # Read the BSRN format data (SolarData package)
  # 使用 SolarData 包读取 BSRN 格式数据
  tmp <- BSRN.read(files[i], use.qc = FALSE, use.agg = FALSE, directory = getwd())

  # Select time and radiation variables
  # 筛选时间与辐射变量
  tmp <- tmp %>%
    rename(GHI = dw_solar, DNI = direct_n, DIF = diffuse) %>%
    dplyr::select(Time, GHI, DNI, DIF)

  data_list[[i]] <- tmp
  setTxtProgressBar(pb, i)
}
close(pb)

data_combined <- bind_rows(data_list) %>%
  filter(lubridate::year(Time) == 2024)

if (nrow(data_combined) == 0) {
  stop("No data found for year 2024! Please check the station files.")
}

################################################################################
# Solar Positioning and QC Limits
# 太阳位置计算与 QC 极限应用
################################################################################

# Calculate solar geometry (center of 1-min interval)
# 计算太阳位置（取 1 分钟间隔的中心点）
solpos <- calZen(data_combined$Time - 30, lat = lat, lon = lon, alt = alt)

data_combined <- data_combined %>%
  mutate(Z = solpos$zenith, ETR = solpos$E0n) %>%
  filter(Z <= 90)

# Compute standard BSRN Extremely Rare Limits (ERL)
# 计算标准的 BSRN “极罕见极限” (ERL)
data_combined <- data_combined %>%
  mutate(
    Glim = 1.2 * ETR * (cos(Z * pi / 180))^1.2 + 50,
    Dlim = 0.75 * ETR * (cos(Z * pi / 180))^1.2 + 30,
    Blim = 0.95 * ETR * (cos(Z * pi / 180))^0.2 + 10
  )

################################################################################
# Data Preparation for Visualization
# 可视化数据准备
################################################################################

# Reshape data for faceted plotting
# 为分面绘图转换数据格式
data.plot <- rbind(
  data.frame(x = data_combined$Z, y = data_combined$GHI, ext.lim = data_combined$Glim, group = "GHI"),
  data.frame(x = data_combined$Z, y = data_combined$DNI, ext.lim = data_combined$Blim, group = "BNI"),
  data.frame(x = data_combined$Z, y = data_combined$DIF, ext.lim = data_combined$Dlim, group = "DHI")
)

# Set factor levels and labels for mathematical expressions
# 设置因子水平和用于数学公式展示的标签
data.plot$group <- factor(data.plot$group,
  levels = c("GHI", "BNI", "DHI"),
  labels = c("italic(G)[italic(h)]", "italic(B)[italic(n)]", "italic(D)[italic(h)]")
)

################################################################################
# Generate Plot
# 生成图形
################################################################################

p <- ggplot(data.plot) +
  geom_scattermore(aes(x = x, y = y), pointsize = 0, alpha = 0.5) +
  geom_scattermore(aes(x = x, y = ext.lim), color = colorblind_pal()(8)[2], pointsize = 0) +
  facet_wrap(~group, labeller = label_parsed, nrow = 1) +
  scale_x_continuous(name = expression(paste("Zenith angle [", degree, "]"))) +
  scale_y_continuous(name = expression(paste("Irradiance [W ", m^-2, "]"))) +
  theme_bw() +
  theme(
    plot.margin = unit(c(0.1, 0.2, 0, 0.0), "lines"),
    text = element_text(family = "Times", size = plot.size),
    axis.text = element_text(family = "Times", size = plot.size),
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Times", size = plot.size),
    panel.spacing = unit(0.5, "lines")
  )

# Preview plot (if interactive)
# p

# Export plot to PDF
# 导出 PDF 图形
if (!dir.exists(dir_tex)) dir.create(dir_tex, recursive = TRUE)
output_pdf <- file.path(dir_tex, "ERLexample.pdf")

ggsave(
  filename = output_pdf, plot = p,
  width = 160, height = 60, units = "mm"
)

cat(sprintf("\nSuccessfully saved the ERL example plot to: %s\n", output_pdf))
