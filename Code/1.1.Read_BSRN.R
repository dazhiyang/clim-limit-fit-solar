################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Read multi-month BSRN data, calculate Zenith angle and other parameters, and generate corresponding txt files
# 读取BSRN多月数据、计算天顶角等太阳几何参数、生成对应的txt文件
################################################################################

# Clear workspace and load necessary packages
# 清理工作台并加载必须的R包
rm(list = ls(all = TRUE))
libs <- c("SolarData", "insol", "dplyr")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# Global Functions
# 全局函数定义
################################################################################

# Function to calculate solar geometry parameters (Zenith angles, extraterrestrial irradiance)
# 计算太阳几何参数的函数（天顶角、地外太阳辐射等）
calZen <- function(Tm, lat, lon, tz = 0, alt = 0) {
  require(insol)
  jd <- JD(Tm) # Julian day
  # 儒略日
  sunv <- sunvector(jd, lat, lon, tz) # Solar vector
  # 太阳向量

  azi <- round(sunpos(sunv)[, 1], 3) # Azimuth of the sun
  # 太阳方位角
  zen <- round(sunpos(sunv)[, 2], 3) # Zenith angle
  # 太阳天顶角

  doy <- daydoy(Tm) # Day of year
  # 积日
  da <- (2 * pi / 365) * (doy - 1) # Day angle
  # 日角

  # Eccentricity correction factor of the Earth's orbit
  # 地球轨道偏心率修正系数
  re <- 1.000110 + 0.034221 * cos(da) + 0.001280 * sin(da) + 0.00719 * cos(2 * da) + 0.000077 * sin(2 * da)

  E0n <- round(1361.1 * re, 3) # Extraterrestrial direct normal irradiance
  # 地外法向直接辐射
  E0 <- round(1361.1 * re * cos(radians(zen))) # Horizontal extraterrestrial irradiance
  # 地外水平面辐射
  E0 <- ifelse(zen >= 90, 0, E0) # Set negative radiance (night) to 0
  # 将夜晚辐射设置为0

  out <- list(zen, azi, E0, E0n)
  names(out) <- c("zenith", "azimuth", "E0", "E0n")
  return(out)
}

# Conversion between radians and degrees
# 弧度与角度转换
r2d <- function(x) {
  x * 180 / pi
}
d2r <- function(x) {
  x * pi / 180
}

################################################################################
# Global Input and Path Settings
# 全局输入和路径设置
################################################################################
options(digits = 3)

# Directories for data loading and saving
# 数据读取和保存的文件夹路径
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation/Data"
dir.out <- file.path(dir0, "BSRNtxt") # Processed data output folder
# 处理后数据保存的文件夹

################################################################################
# Main Data Processing Loop
# 主数据处理循环
################################################################################

# Read BSRN station metadata
# 读取BSRN站点地理信息
bsrn <- read.csv(file.path(dir0, "BSRN information.csv"), header = TRUE)
bsrn$stn <- tolower(bsrn$stn)

# Get all station names from folders
# 从文件夹中获取所有站点的名称列表
bsrn_dir <- file.path(dir0, "BSRN")
stn <- list.dirs(bsrn_dir, full.names = FALSE, recursive = FALSE)
for (i in seq_along(stn)) {
  # List all monthly files (.gz) for the current station
  # 获取当前站点的所有月度数据文件 (.gz)
  dir.tmp <- list.files(file.path(bsrn_dir, stn[i]), pattern = "\\.gz$", full.names = FALSE, recursive = TRUE)

  if (length(dir.tmp) == 0) {
    message(paste("No data found for station:", stn[i], "- Skipping / 未找到数据，跳过"))
    next
  }

  data.stn <- NULL
  cat(sprintf("\nProcessing station: %s (%d of %d) / 正在处理站点: %s\n", stn[i], i, length(stn), stn[i]))

  # Setup progress bar
  # 设置进度条
  pb <- txtProgressBar(max = length(dir.tmp), style = 3)

  for (j in seq_along(dir.tmp)) {
    # Read the data using the SolarData package
    # 使用SolarData库读取BSRN数据格式
    tmp <- BSRN.read(dir.tmp[j], directory = file.path(bsrn_dir, stn[i]), use.qc = FALSE, use.agg = FALSE)

    # Select only the shortwave variables and time
    # 筛选出时间与短波辐射变量
    tmp <- tmp %>%
      rename(Gh = dw_solar, Dh = diffuse, Bn = direct_n) %>%
      dplyr::select(one_of("Time", "Gh", "Dh", "Bn"))

    # Append monthly data into one big tibble
    # 将当前月度数据拼接到整个站点的结果中
    data.stn <- bind_rows(data.stn, tmp)

    setTxtProgressBar(pb, j)
  }
  close(pb)

  # Sort timeframe chronologically
  # 按时间排序
  data.stn <- data.stn %>% arrange(Time)

  # Find matching metadata row for current station
  # 获取当前站点的经纬度和海拔
  stn_meta <- bsrn %>% filter(stn == !!stn[i])
  if (nrow(stn_meta) == 0) {
    warning(paste("Metadata not found for station:", stn[i], "- Cannot calculate zenith. Skipping / 缺失地理信息，跳过"))
    next
  }

  # Perform solar positioning calculations
  # 进行太阳位置和辐射参量计算
  # Note: Time - 30s is used to center the 1-minute interval
  # 时间减去30秒以居中表示1分钟数据
  solpos <- calZen(
    Tm = data.stn$Time - 30,
    lat = stn_meta$Latitude[1],
    lon = stn_meta$Longitude[1],
    tz = 0,
    alt = stn_meta$Elevation[1]
  )

  # Replace infinite values with NA
  # 将无穷大值替换为NA
  replace_inf <- function(x) ifelse(is.finite(x), x, NA)

  # Compute clearness indices and format data
  # 计算晴空指数(clearness indices)并格式化
  data <- data.stn %>%
    mutate(
      Z = solpos$zenith,
      ETR = solpos$E0n,
      kt = round(Gh / solpos$E0, 3),
      kd = round(Dh / solpos$E0, 3),
      kb = round(Bn / ETR, 3)
    ) %>%
    mutate_at(vars(kt, kd, kb), replace_inf) %>%
    mutate(Time = as.character(format(Time)))

  # Export combined station data to txt
  # 导出站点的汇编txt文件
  output_file <- file.path(dir.out, paste0(stn[i], ".txt"))
  write.table(data, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE)
}

cat("\nAll operations successfully completed. / 所有站点的数据合并处理已完成。\n")
