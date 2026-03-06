################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Diagnostic Multiplot for Isolated Climate Regimes
# 汇集各气候分离区域模型观测图进行多图组合诊断可视化
################################################################################

# Clear workspace and load libraries
# 清空工作区并加载库
# nolint start
rm(list = ls(all = TRUE))
libs <- c("dplyr", "ggplot2", "lubridate", "data.table", "gridExtra", "insol", "ggthemes", "tidyr", "scattermore", "magrittr")
invisible(lapply(libs, require, character.only = TRUE))

################################################################################
# Global Input and Path Settings
# 全局输入与路径设置
################################################################################
# Base directory for the project
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data/BSRNtxt")
dir_tex <- file.path(dir0, "tex")
setwd(dir0) # Force root directory for relative paths / 强制根目录路径

plot.size <- 8
line.size <- 0.2
point.size <- 0.5
legend.size <- 0.4
text.size <- plot.size * 5 / 14

cat("\nPlot objects loaded. Assembling multiplot...\n")

################################################################################
# Regime 1 (Equatorial Monsoon Bifurcation)
# 区域1（赤道季风分岔）
# Target Station: dar (Darwin, Australia) from Cluster 1
# 目标站点：达尔文（澳大利亚）
# Characteristic: High theoretical clear-sky $G_h$ max, severe seasonal variance.
# 特征：理论晴空 $G_h$ 最大值高，季节性方差大。
################################################################################
cat("Generating Subplot 1 (Regime 1: Darwin, AUS)...\n")

# Load Actual Observations
# 加载实际观测数据
dar_obs <- data.table::fread("Data/BSRNtxt/dar.txt") %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time),
        Year = lubridate::year(Time_Parsed),
        DOY = lubridate::yday(Time_Parsed)
    ) %>%
    filter(Year == 2014, Gh > 0) %>%
    group_by(DOY) %>%
    summarise(Max_Observed_Gh = max(Gh, na.rm = TRUE), .groups = "drop")

# Load McClear Clear-Sky Simulations
# 加载 McClear 晴空模拟数据
dar_mcclear <- data.table::fread("Data/McClear/dar/dar_2014.csv") %>%
    rename(
        Clear_Gh = ghi_clear
    ) %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time),
        DOY = lubridate::yday(Time_Parsed)
    ) %>%
    group_by(DOY) %>%
    summarise(Max_Clear_Gh = max(Clear_Gh, na.rm = TRUE), .groups = "drop")

# Merge Datasets for 2014 Daily Maximums
# 合并数据获取 2014 年每日最大值
dar_merged <- dar_obs %>%
    inner_join(dar_mcclear, by = "DOY")

# Create Subplot 1
# 创建子图 1
p1 <- ggplot(dar_merged, aes(x = DOY)) +
    # Observed Variance Envelope (Shaded + Points)
    # 观测方差包络（阴影 + 散点）
    geom_ribbon(aes(ymin = 0, ymax = Max_Observed_Gh), fill = "gray80", alpha = 0.7) +
    geom_point(aes(y = Max_Observed_Gh), color = colorblind_pal()(8)[3], size = point.size, alpha = 0.7, stroke = 0) +
    # Theoretical McClear Clear-Sky Ceiling
    # 理论 McClear 晴空上限
    geom_line(aes(y = Max_Clear_Gh), color = colorblind_pal()(8)[2], linewidth = line.size) +
    labs(
        x = "Day of year",
        y = expression(paste("Maximum daily ", italic(G)[italic(h)], " [W ", m^-2, "]")),
        color = NULL
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.2, 0.2, 0, 0.0), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        strip.text.x = element_text(family = "Times", size = plot.size, margin = margin(0.1, 0, 0.1, 0, "lines")),
        strip.text.y = element_text(family = "Times", size = plot.size, margin = margin(0, 0.1, 0, 0.1, "lines")),
        panel.spacing = unit(0, "lines"),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size = plot.size),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.box.margin = unit(c(-0.7, 0, 0, 0.0), "lines")
    )

# 4. Assemble and save the multiplot
# 使用 gridExtra::grid.arrange() 组合并输出图形
# We will use gridExtra::grid.arrange() to define the layout and output.
cat("Saving diagnostic multiplot...\n")
ggsave("tex/Regime1.pdf", plot = p1, width = 85, height = 55, units = "mm", dpi = 300)

################################################################################
# Regime 2 (Polar and Sub-Polar, PSP)
# 区域2（极地和副极地）
# Target Station: nya (Ny-Ålesund, Svalbard, Norway) from Cluster 2
# 目标站点：新奥尔松（挪威斯瓦尔巴群岛）
# Characteristic: Extreme seasonal cycles, Polar Night (zero irradiance), low max G_h
# 特征：极端的季节周期，极夜（零辐射），低理论G_h上限
################################################################################
cat("\nGenerating Subplot 2 (Regime 2: Ny-Ålesund, NOR)...\n")

# Load High-Resolution McClear Clear-Sky Simulations for nya
# 加载 McClear 高分辨率晴空模拟数据
nya_mcclear <- data.table::fread("Data/McClear/nya/nya_2021.csv") %>%
    rename(
        Clear_Gh = ghi_clear
    ) %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time),
        DOY = lubridate::yday(Time_Parsed),
        Hour = lubridate::hour(Time_Parsed),
        # Calculate precise Solar Elevation mathematically (h = 90 - Z)
        # 计算精确的太阳高度角 (h = 90 - Z)
        JD = insol::JD(Time_Parsed),
        SunV = insol::sunvector(JD, latitude = 78.9227, longitude = 11.9273, timezone = 0),
        Zenith = insol::sunpos(SunV)[, 2],
        Elevation = 90 - Zenith
    ) %>%
    # Calculate hourly average Clear-Sky Gh and Elevation to create smooth heatmap pixels
    # 计算每小时平均晴空 Gh 和高度角以创建平滑的热图像素
    group_by(DOY, Hour) %>%
    summarise(
        Mean_Clear_Gh = mean(Clear_Gh, na.rm = TRUE),
        Mean_Elevation = mean(Elevation, na.rm = TRUE),
        .groups = "drop"
    )

# Create Subplot 2: Diurnal Heatmap (The Hourglass)
# 创建子图 2: 日变化热图（沙漏图）
p2 <- ggplot(nya_mcclear, aes(x = DOY, y = Hour)) +
    geom_raster(aes(fill = Mean_Clear_Gh), interpolate = TRUE) +
    # Add Solar Elevation Contours
    # 添加太阳高度角等值线
    geom_contour(aes(z = Mean_Elevation), color = "white", alpha = 0.5, binwidth = 10, linewidth = line.size) +
    # Highlight the horizon (Elevation = 0), marking the boundary of Polar Day/Night
    # 突出显示地平线（高度角 = 0），标记极昼极夜的边界
    geom_contour(aes(z = Mean_Elevation), breaks = 0, color = "white", linewidth = line.size * 2, linetype = "dashed") +
    # Use default Viridis for continuous color mapping
    # 使用 Viridis 提供连续色彩映射
    scale_fill_viridis_c(
        option = "viridis",
        name = expression(paste(italic(G)[italic(h)], " [W ", m^-2, "]"))
    ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 24, 6)) +
    labs(
        x = "Day of year",
        y = "Hour of day (UTC)"
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.1, 0.2, 0, 0.1), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        panel.grid = element_blank(), # Remove grid lines over the raster
        legend.position = "bottom",
        legend.key.height = unit(0.3, "cm"),
        legend.key.width = unit(1.0, "cm"),
        legend.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "lines"),
        legend.box.margin = margin(t = -0.5, b = 0, l = 0, r = 0, unit = "lines"),
        legend.title = element_text(family = "Times", size = plot.size),
        legend.text = element_text(family = "Times", size = plot.size)
    )

cat("Saving diagnostic plot for Regime 2...\n")
ggsave("tex/Regime2.pdf", plot = p2, width = 85, height = 55, units = "mm", dpi = 300)

################################################################################
# Regime 3 (Arid and Semi-Arid Deserts)
# 区域3（干旱和半干旱沙漠）
# Target Station: dra (Desert Rock, USA) from Cluster 3
# 目标站点：沙漠岩（美国）
# Characteristic: Highly stable, dry air masses, predominantly clear-sky
# 特征：高度稳定、干燥的气团，主要为晴空
################################################################################
cat("\nGenerating Subplot 3 (Regime 3: Desert Rock, USA)...\n")

# Load 1-minute Actual Observations for dra (Year 2021)
# 加载 dra 的 1 分钟实际观测数据 (2021年)
dra_obs <- data.table::fread("Data/BSRNtxt/dra.txt") %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time),
        Year = lubridate::year(Time_Parsed),
        Month = lubridate::month(Time_Parsed, label = TRUE, abbr = TRUE)
    ) %>%
    filter(Year == 2021, Gh >= 0) %>%
    select(Time_Parsed, Month, Gh)

# Load McClear Clear-Sky Simulations for dra
# 加载 McClear 晴空模拟数据
dra_mcclear <- data.table::fread("Data/McClear/dra/dra_2021.csv") %>%
    rename(
        Clear_Gh = ghi_clear
    ) %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time),
        # Calculate precise Solar Zenith parameter mathematically
        # 数学计算精确太阳天顶角
        JD = insol::JD(Time_Parsed),
        SunV = insol::sunvector(JD, latitude = 36.626, longitude = -116.018, timezone = 0),
        Zenith = insol::sunpos(SunV)[, 2]
    ) %>%
    select(Time_Parsed, Clear_Gh, Zenith)

# Merge and Calculate Clear-Sky Index (Kc)
# 合并并计算晴空指数 (Kc)
dra_merged <- dra_obs %>%
    inner_join(dra_mcclear, by = "Time_Parsed") %>%
    # Filter highly oblique angles to prevent massive Kc artifact stretching (Z < 85)
    # 过滤掉高倾角以防止 Kc 产生极大伪影 (Z < 85)
    filter(Zenith < 85, Clear_Gh > 10) %>%
    mutate(
        Kc = Gh / Clear_Gh
    ) %>%
    # Filter extreme outliers just for density visual clarity (0 to 1.5 strictly)
    # 仅为了分布图清晰，限制极大离群值范围 (0 - 1.5)
    filter(Kc >= 0, Kc <= 1.5)

# Create Subplot 3: Violin Density Plot of Kc
# 创建子图 3: Kc 月度小提琴图
p3 <- ggplot(dra_merged, aes(x = Month, y = Kc)) +
    # Use violin plots to prove stability (massive horizontal spike at 1.0)
    # 使用小提琴图证明稳定性（在 1.0 处有巨大的峰值）
    geom_violin(fill = colorblind_pal()(8)[2], alpha = 0.7, scale = "width", linewidth = line.size / 2) +
    # Add horizontal reference line exactly at Kc = 1.0 (Theoretical perfect clear-sky)
    # 在 Kc = 1.0 处添加水平参考线（理论完美晴空）
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "gray30", linewidth = line.size) +
    scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, 0.5)) +
    labs(
        x = "Month",
        y = expression(paste("Clear-sky index, ", italic(kappa)))
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.2, 0.2, 0, 0.0), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)
    )

cat("Saving diagnostic plot for Regime 3...\n")
ggsave("tex/Regime3.pdf", plot = p3, width = 85, height = 55, units = "mm", dpi = 300)

################################################################################
# Regime 4 (Marine and Oceanic)
# 区域4（海洋区域）
# Target Stations: ler (Lerwick, UK - Marine) vs fpe (Fort Peck, USA - Cont)
# 目标站点：勒威克 (英国 - 海洋) 对比 福特佩克 (美国 - 大陆)
# Characteristic: Persistent marine stratocumulus producing higher diffuse fraction
# 特征：持续的海洋层积云产生更高的散射比例
################################################################################
cat("\nGenerating Subplot 4 (Regime 4: Marine vs Continental Diffuse)... \n")

# Load 1-minute Actual Observations for ler (Marine, All Years)
# 加载 ler 的所有 1 分钟实际观测数据
ler_obs <- data.table::fread("Data/BSRNtxt/ler.txt") %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time)
    ) %>%
    # Filter mathematically valid daytime observations
    # 过滤有明确物理意义的日间数据
    filter(Gh > 10, Dh >= 0, Dh <= Gh) %>%
    mutate(
        Station = "Lerwick (Marine)",
        Regime = "Regime 4",
        kd = Dh / Gh
    ) %>%
    select(Time_Parsed, Station, Regime, kd)

# Load 1-minute Actual Observations for fpe (Continental, All Years)
# 加载 fpe 的所有 1 分钟实际观测数据
fpe_obs <- data.table::fread("Data/BSRNtxt/fpe.txt") %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time)
    ) %>%
    # Filter mathematically valid daytime observations
    # 过滤有明确物理意义的日间数据
    filter(Gh > 10, Dh >= 0, Dh <= Gh) %>%
    mutate(
        Station = "Fort Peck (Continental)",
        Regime = "Regime 5",
        kd = Dh / Gh
    ) %>%
    select(Time_Parsed, Station, Regime, kd)

# Bind both datasets linearly
# 线性合并两个数据集
diffuse_merged <- bind_rows(ler_obs, fpe_obs)

# Create Subplot 4: Overlapping Diffuse Fraction Densities
# 创建子图 4: 叠加的散射比例核密度图
p4 <- ggplot(diffuse_merged, aes(x = kd, fill = Station, color = Station)) +
    # Use density distributions to prove bulk variance shifts
    # 使用核密度分布证明整体方差偏移
    geom_density(alpha = 0.5, linewidth = line.size) +
    scale_fill_manual(values = c("Lerwick (Marine)" = colorblind_pal()(8)[2], "Fort Peck (Continental)" = colorblind_pal()(8)[3])) +
    scale_color_manual(values = c("Lerwick (Marine)" = colorblind_pal()(8)[2], "Fort Peck (Continental)" = colorblind_pal()(8)[3])) +
    scale_x_continuous(limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2)) +
    labs(
        x = expression(paste("Diffuse fraction, ", italic(k))),
        y = "Density"
    ) +
    guides(
        fill = guide_legend(ncol = 1),
        color = guide_legend(ncol = 1)
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.1, 0.2, 0, 0.1), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.5, 0.5),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size = plot.size),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "lines"),
        legend.key.size = unit(0.3, "cm")
    )

cat("Saving diagnostic plot for Regime 4...\n")
ggsave("tex/Regime4.pdf", plot = p4, width = 85, height = 55, units = "mm", dpi = 300)

################################################################################
# Regime 5 (Continental Temperate Mid-Latitudes)
# 区域5（大陆性温带中纬度）
# Target Station: pay (Payerne, Switzerland) from Cluster 5
# 目标站点：帕耶讷 (瑞士)
# Characteristic: Distinct seasonal cycle superimposed with severe day-to-day variability
# 特征：明显的季节周期叠加剧烈的日际波动 (受天气系统交替影响)
################################################################################
cat("\nGenerating Subplot 5 (Regime 5: Temperate High-Variance Dual Scatter)... \n")

# Load 1-minute Actual Observations for pay (Continental, Year 2021)
# 加载 pay 的 1 分钟实际观测数据 (以 2021 年作为展示年份)
pay_obs <- data.table::fread("Data/BSRNtxt/pay.txt") %>%
    mutate(
        Time_Parsed = lubridate::as_datetime(Time),
        Year = lubridate::year(Time_Parsed),
        DOY = lubridate::yday(Time_Parsed)
    ) %>%
    # Filter for the full target year and physically valid data
    # 过滤完整的目标年份并剔除物理上不合理的数据
    filter(Year == 2021, Gh >= 0, Bn >= 0, Dh >= 0)

# Aggregate 1-minute irradiance into Daily Mean (W/m^2) mathematically equivalent to daily integral energy
# 将 1 分钟辐照度聚合为日均值 (W/m^2)，这在数学上等同于日累积总能量
pay_daily <- pay_obs %>%
    group_by(DOY) %>%
    # Taking the mean across 24 hours mathematically scales linearly to total megajoules
    # 取 24 小时的平均值在数学上等同于总焦耳数
    summarise(
        Daily_Bn = mean(Bn, na.rm = TRUE),
        Daily_Dh = mean(Dh, na.rm = TRUE)
    )

# Reshape from Wide to Long format to explicitly overlay Direct and Diffuse on one graph
# 将宽表重塑为长表格式，以便在同一张图上显式叠加直射和散射
pay_long <- tidyr::pivot_longer(pay_daily,
    cols = c(Daily_Bn, Daily_Dh),
    names_to = "Component",
    values_to = "Irradiance"
) %>%
    mutate(
        # Keep underlying structural labels intact since we will mathematically format them in the legend scale
        Component = case_when(
            Component == "Daily_Bn" ~ "Daily_Bn",
            Component == "Daily_Dh" ~ "Daily_Dh"
        )
    )

# Create Subplot 5: 365-Day Dual Scatter Plot mapping massive daily variance
# 创建子图 5: 365 天双散点图，映射巨大的日均方差
p5 <- ggplot(pay_long, aes(x = DOY, y = Irradiance, color = Component, fill = Component)) +
    # Use thin connecting lines to definitively trace the "day-to-day wild zig-zagging"
    # 使用极细的连线明确追踪"剧烈的日际锯齿波动"
    geom_line(alpha = 0.3, linewidth = line.size) +
    # Map explicit daily occurrence points proving the distinct variance
    # 映射特定的每日发生点证明明显的分散度
    geom_point(alpha = 0.8, size = point.size, shape = 21, color = "white", stroke = 0.1) +
    scale_fill_manual(
        values = c("Daily_Bn" = colorblind_pal()(8)[6], "Daily_Dh" = colorblind_pal()(8)[2]),
        labels = c("Daily_Bn" = expression(paste("BNI (", italic(B)[italic(n)], ")")), "Daily_Dh" = expression(paste("DHI (", italic(D)[italic(h)], ")")))
    ) +
    scale_color_manual(
        values = c("Daily_Bn" = colorblind_pal()(8)[6], "Daily_Dh" = colorblind_pal()(8)[2]),
        labels = c("Daily_Bn" = expression(paste("BNI (", italic(B)[italic(n)], ")")), "Daily_Dh" = expression(paste("DHI (", italic(D)[italic(h)], ")")))
    ) +
    scale_x_continuous(limits = c(0, 366), breaks = seq(0, 360, 60)) +
    labs(
        x = "Day of year",
        y = expression(paste("Daily mean irradiance (W ", m^{
            -2
        }, ")"))
    ) +
    guides(
        fill = guide_legend(ncol = 2),
        color = guide_legend(ncol = 2)
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.2, 0.2, 0, 0.0), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size = plot.size),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        legend.margin = margin(t = -0.5, b = 0, l = 0, r = 0, unit = "lines"),
        legend.key.size = unit(0.5, "cm")
    )

cat("Saving diagnostic plot for Regime 5...\n")
ggsave("tex/Regime5.pdf", plot = p5, width = 85, height = 55, units = "mm", dpi = 300)

################################################################################
# Regime 6 (High Aerosol/Polluted)
# 区域6（高气溶胶/污染区域）
# Target Stations: gan (Gandhinagar) and gur (Gurgaon) from Cluster 6
# 目标站点：甘地讷格尔 (gan) 和 古尔冈 (gur)
# Characteristic: Extreme aerosol optical depths (AOD) leading to Bn attenuation and high Dh
# 特征：极高的 AOD 导致直射 (Bn) 显著衰减，散射 (Dh) 比例升高
################################################################################
cat("\nGenerating Subplot 6 (Regime 6: Aerosol Attenuation Scatter)... \n")

# Define processing function for Regime 6 contrast: Indian stations vs Clean Baselines
# 定义区域6对比处理函数：印度站点 vs 清洁基准站点
process_regime6_contrast <- function(stn, lat, lon, year) {
    # Load 1-minute observations
    obs_file <- paste0("Data/BSRNtxt/", stn, ".txt")
    obs <- data.table::fread(obs_file) %>%
        mutate(
            Time_Parsed = lubridate::as_datetime(Time),
            Year = lubridate::year(Time_Parsed),
            DOY = lubridate::yday(Time_Parsed)
        ) %>%
        # Filter for selected year and valid daytime data
        filter(Year == year, Gh > 10, Bn >= 0) %>%
        # Calculate solar geometry
        mutate(
            jd = insol::JD(Time_Parsed),
            sv = insol::sunvector(jd, lat, lon, 0), # UTC
            sp = insol::sunpos(sv),
            Zenith = sp[, 2],
            # Calculate ERL limit for Bn (BSRN Extremely Rare Limit)
            # 计算 Bn 的 BSRN 极罕见极限 (ERL)
            I0n = 1367 * (1 + 0.033 * cos(2 * pi * DOY / 365.25)),
            Bn_ERL = 0.95 * I0n * (cos(Zenith * pi / 180)^0.2) + 10
        ) %>%
        filter(Zenith < 85 & Zenith > 0) %>%
        mutate(
            Station = case_when(
                stn == "gan" ~ "GAN (India, 23.1°N)",
                stn == "gur" ~ "GUR (India, 28.4°N)",
                stn == "iza" ~ "IZA (Spain, 28.3°N)",
                stn == "yus" ~ "YUS (Taiwan, 23.5°N)"
            ),
            Type = ifelse(stn %in% c("gan", "gur"), "Regime 6 (Polluted)", "Baseline (Clean)")
        ) %>%
        select(Time_Parsed, Station, Type, Bn, Zenith, Bn_ERL)
}

# Process GAN (2015), GUR (2015), IZA (2021), YUS (2019)
# 处理 GAN (23°N), GUR (28°N), IZA (28°N, 清洁), YUS (23°N, 清洁)
reg6_data <- bind_rows(
    process_regime6_contrast("gan", 23.1101, 72.6276, 2015),
    process_regime6_contrast("gur", 28.4249, 77.156, 2018),
    process_regime6_contrast("iza", 28.3093, -16.4993, 2021),
    process_regime6_contrast("yus", 23.4876, 120.9595, 2019)
)

# Create Subplot 6: Latitudinal Contrast of Bn vs Zenith
# 创建子图 6: 相似纬度下的直射辐照度对比
p6 <- ggplot(reg6_data, aes(x = Zenith, y = Bn)) +
    # Use scattermore for high-density 1-minute scatter points
    scattermore::geom_scattermore(color = "black", alpha = 0.5, pointsize = 0) +
    # Add the ERL limit curve as points (Extremely Rare Limit) using project standard color
    # 添加 ERL 极限曲线作为散点 (极罕见极限)，使用项目标准颜色 (Orange from colorblind_pal)
    scattermore::geom_scattermore(aes(y = Bn_ERL), color = colorblind_pal()(8)[2], pointsize = 0) +
    facet_wrap(~Station, ncol = 2) +
    scale_x_continuous(limits = c(0, 90), breaks = seq(0, 90, 30)) +
    scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 500)) +
    labs(
        x = expression(paste("Solar zenith angle, ", italic(Z), " (", degree, ")")),
        y = expression(paste("Observed BNI, ", italic(B)[italic(n)], " (W ", m^{
            -2
        }, ")"))
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.2, 0.2, 0, 0.2), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        strip.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)
    )

cat("Saving diagnostic plot for Regime 6 (with BSRN ERL limit)...\n")
ggsave("tex/Regime6.pdf", plot = p6, width = 85, height = 75, units = "mm", dpi = 300)

################################################################################
# Regime 7 (Low-latitude anomalous)
# 区域7（低纬度异常区域）
# Target Stations: nau, ohy, sms, tir from Cluster 7
# 目标站点：瑙鲁 (nau), 豪拉基 (ohy), 圣玛丽亚 (sms), 蒂鲁伯蒂 (tir)
# Characteristic: Statistical similarity across geographically distant low-latitude sites
# 特征：地理位置遥远但统计特性高度相似的低纬度站点
################################################################################
cat("\nGenerating Subplot 7 (Regime 7: Statistical Distribution Overlay)... \n")

# Function to load all available data for Regime 7 stations
process_regime7_stat <- function(stn) {
    obs_file <- paste0("Data/BSRNtxt/", stn, ".txt")
    data.table::fread(obs_file) %>%
        mutate(Station = toupper(stn)) %>%
        select(Station, kt, kd, kb)
}

# Load all data for the four stations
reg7_data <- bind_rows(
    process_regime7_stat("nau"),
    process_regime7_stat("ohy"),
    process_regime7_stat("sms"),
    process_regime7_stat("tir")
)

# Prepare data for distribution plotting (Pivot to long format)
reg7_dist_long <- reg7_data %>%
    tidyr::pivot_longer(cols = c(kt, kd, kb), names_to = "Parameter", values_to = "Value") %>%
    mutate(Parameter = factor(Parameter, levels = c("kt", "kd", "kb"), labels = c("italic(k)[italic(t)]", "italic(k)[italic(d)]", "italic(k)[italic(b)]")))

# Create Subplot 7: Statistical Distribution Overlay
p7 <- ggplot(reg7_dist_long, aes(x = Value, color = Station, fill = Station)) +
    geom_density(alpha = 0.2, linewidth = line.size) +
    facet_wrap(~Parameter, scales = "free", labeller = label_parsed, ncol = 1) +
    scale_x_continuous(limits = c(0, 1.2), breaks = seq(0, 1, 0.5)) +
    scale_color_manual(values = colorblind_pal()(8)[c(2:4, 7)]) +
    scale_fill_manual(values = colorblind_pal()(8)[c(2:4, 7)]) +
    labs(
        x = "Parameter value",
        y = "Density"
    ) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        strip.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size = plot.size - 1),
        legend.key.size = unit(0.3, "cm"),
        legend.margin = margin(t = -0.5, b = 0, l = 0, r = 0, unit = "lines"),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA)
    )

cat("Saving diagnostic plot for Regime 7 Statistical Overlay...\n")
ggsave("tex/Regime7.pdf", plot = p7, width = 85, height = 90, units = "mm", dpi = 300)

# nolint end
