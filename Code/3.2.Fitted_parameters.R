################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Generating LaTeX Table for Fitted Cluster Parameters
# 将聚类拟合参数转换为 LaTeX 表格
################################################################################

# Load libraries
# 加载库
library(dplyr)
library(tidyr)
library(xtable)

# Set directories
# 设置目录
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
setwd(dir0)

# Load the cluster parameters
# 加载聚类参数
cat("Reading cluster_parameter.csv...\n")
df <- read.csv("Data/cluster_parameter.csv")

# 1. Exclude the 'Stations' column
# 排除 'Stations' 列
df_clean <- df %>% select(-Stations)

# 2. Format numerical values
# 将数值格式化
# Cluster_Name should be an integer, others two-digit decimals
df_formatted <- df_clean %>%
    mutate(Cluster_Name = as.character(as.integer(Cluster_Name))) %>%
    mutate(across(-Cluster_Name, ~ formatC(.x, format = "f", digits = 2)))

# Define cleaner column names for LaTeX
# 定义更整洁的 LaTeX 列名
colnames(df_formatted) <- c(
    "Cluster",
    "$a_G$", "$b_G$", "$c_G$", "$F_{1,G}$",
    "$a_B$", "$b_B$", "$c_B$", "$F_{1,B}$",
    "$a_D$", "$b_D$", "$c_D$", "$F_{1,D}$"
)

# 3. Convert to LaTeX table
# 转换为 LaTeX 表格
cat("\nGenerating LaTeX table content:\n\n")

# Using xtable for a clean LaTeX output
latex_table <- xtable(df_formatted,
    align = c("l", "c", rep("r", ncol(df_formatted) - 1)),
    caption = "Optimized parameters for the proposed QC limits across seven climatic regimes.",
    label = "tab:fitted_parameters"
)

# Print the table code
# We use sanitize.text.function = identity to keep the $math$ mode intact
print(latex_table,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    booktabs = TRUE,
    floating = TRUE,
    table.placement = "htbp",
    caption.placement = "top"
)

# Optionally save to a .tex file
cat("\nSaving table to tex/fitted_parameters.tex...\n")
sink("tex/fitted_parameters.tex")
print(latex_table,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    booktabs = TRUE,
    floating = TRUE,
    table.placement = "htbp",
    caption.placement = "top"
)
sink()

cat("Finished.\n")
