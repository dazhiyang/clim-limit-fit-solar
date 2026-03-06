################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# Harbin Institute of Technology
# emails: stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Generate a LaTeX table of BSRN Station Metadata
# 生成BSRN站点地理信息的LaTeX表格
################################################################################

# Clear workspace and load necessary packages
# 清理工作台并加载必须的R包
rm(list = ls(all = TRUE))
libs <- c("dplyr", "xtable")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# Global Input and Path Settings
# 全局输入与路径设置
################################################################################

# Base directory for data
# 基础数据获取目录
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation/Data"

################################################################################
# Data Processing
# 数据处理
################################################################################

# Read in metadata for BSRN sites
# 读取BSRN站点地理信息
bsrn_info_file <- file.path(dir0, "BSRN information.csv")
bsrn <- read.csv(bsrn_info_file, header = TRUE, stringsAsFactors = FALSE)

# Clean and select relevant columns for the table
# 清洗并选择用于表格的相关列
if (ncol(bsrn) > 9) {
    bsrn <- bsrn[, 1:9]
}

# Delete the first column "Station" and the third column "Location"
# 保留：Code, Lat, Lon, Elev, Months, Start, End
bsrn <- bsrn %>% select(stn, Latitude, Longitude, Elevation, month, start, end)

# Ensure stn is character string
bsrn$stn <- as.character(bsrn$stn)


# Formats negative signs to use LaTeX math mode (e.g., -23.80 -> $-$23.80)
# 将负数转换为数学模式LaTeX格式，保留两位小数
format_negative_math <- function(x) {
    ifelse(x < 0, paste0("$-$", sprintf("%.2f", abs(x))), sprintf("%.2f", x))
}

bsrn$Latitude <- format_negative_math(as.numeric(bsrn$Latitude))
bsrn$Longitude <- format_negative_math(as.numeric(bsrn$Longitude))

# Rename columns for the LaTeX table
# 重命名LaTeX表格的列名
colnames(bsrn) <- c("Code", "Lat.", "Lon.", "Elev.", "\\# Mon.", "Start", "End")

# Split the dataframe into two side-by-side tables
# 将表切割成左右两列以缩短表格长度
n_rows <- nrow(bsrn)
half_rows <- ceiling(n_rows / 2)

df1 <- bsrn[1:half_rows, ]
df2 <- bsrn[(half_rows + 1):n_rows, ]

# Pad the second half with NAs if the number of rows is odd
# 如果行数为奇数，则使用NA填充
if (nrow(df2) < half_rows) {
    df2[nrow(df2) + 1, ] <- NA
}

# Bind them together sideway
# 左右拼接成一个新表格
bsrn_wide <- cbind(df1, df2)

# Create the xtable object (It now has 14 columns, plus 1 for rownames = 15 align items)
align_str <- c("l", rep(c("l", "r", "r", "r", "c", "c", "c"), 2))
latex_table <- xtable(bsrn_wide,
    caption = "Metadata for BSRN stations used in this study.",
    label = "tab:metadata",
    align = align_str
)

# Define the output directory for tex files
# 定义tex文件的输出路径
tex_dir <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation/tex"
if (!dir.exists(tex_dir)) {
    dir.create(tex_dir)
}

# Export the LaTeX table
# 导出LaTeX表格
output_file <- file.path(tex_dir, "BSRN_metadata_table.tex")

print(latex_table,
    type = "latex",
    file = output_file,
    include.rownames = FALSE,
    booktabs = TRUE,
    sanitize.text.function = identity, # prevents escaping the math mode $ signs and \#
    caption.placement = "top"
)

cat("\nLaTeX table successfully generated at:", output_file, "\n")
