################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Generating LaTeX Table for Component-wise Rejection Rates (Old ERL vs New Limits)
# 生成旧 ERL 限制与新限制的组件分项（G, B, D）拒绝率 LaTeX 表格
################################################################################

# Load libraries
library(dplyr)
library(xtable)

# Set directories
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
setwd(dir0)

# Load the rejection summary
cat("Reading rejection_summary.csv...\n")
df <- read.csv("Data/rejection_summary.csv")

# Prepare data for the rejection rate table
# Format: "Abs rej numbers" on top and "(percentage%)" below using \makecell
cat("Formatting table data with multi-line cells...\n")
df_formatted <- df %>%
    mutate(
        Points = format(Total_Points, big.mark = ","),
        G_Old = paste0("\\makecell{", format(G_Old_Rej, big.mark = ","), " \\\\ (", formatC(G_Old_Rate, format = "f", digits = 2), "\\%)}"),
        G_New = paste0("\\makecell{", format(G_New_Rej, big.mark = ","), " \\\\ (", formatC(G_New_Rate, format = "f", digits = 2), "\\%)}"),
        B_Old = paste0("\\makecell{", format(B_Old_Rej, big.mark = ","), " \\\\ (", formatC(B_Old_Rate, format = "f", digits = 2), "\\%)}"),
        B_New = paste0("\\makecell{", format(B_New_Rej, big.mark = ","), " \\\\ (", formatC(B_New_Rate, format = "f", digits = 2), "\\%)}"),
        D_Old = paste0("\\makecell{", format(D_Old_Rej, big.mark = ","), " \\\\ (", formatC(D_Old_Rate, format = "f", digits = 2), "\\%)}"),
        D_New = paste0("\\makecell{", format(D_New_Rej, big.mark = ","), " \\\\ (", formatC(D_New_Rate, format = "f", digits = 2), "\\%)}")
    ) %>%
    select(Cluster, Points, G_Old, G_New, B_Old, B_New, D_Old, D_New)

# Define cleaner column names with LaTeX math mode
colnames(df_formatted) <- c(
    "Cluster", "Total points",
    "Global (ERL)", "Global (New)",
    "Beam (ERL)", "Beam (New)",
    "Diffuse (ERL)", "Diffuse (New)"
)

# Convert to LaTeX table
cat("\nGenerating LaTeX table content:\n\n")

latex_table <- xtable(df_formatted,
    align = c("l", "l", "r", "r", "r", "r", "r", "r", "r"),
    caption = "Component-wise Rejection Counts and Rates: Absolute number (top) and Percentage (bottom).",
    label = "tab:component_rejection_rates"
)

# Print the table code
print(latex_table,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    booktabs = TRUE,
    floating = TRUE,
    table.placement = "htbp",
    caption.placement = "top",
    hline.after = c(-1, 0, nrow(df_formatted))
)

# Save to a .tex file
cat("\nSaving table to tex/rejection_rates_component.tex...\n")
if (!dir.exists("tex")) dir.create("tex")
sink("tex/rejection_rates_component.tex")
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
