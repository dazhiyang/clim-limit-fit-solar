################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# False Negative Analysis (Regime 6)
# Aim: Demonstrate that data accepted by ERL is physically invalid (direct component too high)
#      by showing poor closure relationship or being outside physical limits,
#      while being caught by the New limit.
################################################################################

rm(list = ls(all = TRUE))
libs <- c("dplyr", "ggplot2", "lubridate", "data.table", "gridExtra", "insol", "ggthemes", "tidyr")
invisible(lapply(libs, require, character.only = TRUE))

# Global Input and Path Settings
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data")
dir_tex <- file.path(dir0, "tex")
setwd(dir_data)

plot.size <- 8
line.size <- 0.2
point.size <- 0.5
legend.size <- 0.4
text.size <- plot.size * 5 / 14

# 1. Load the False Negative Discrepancy Cases for Regime 6
cat("Loading discrepancy cases...\n")
regime6_cases <- data.table::fread("Regime6_New_only.txt") %>%
    mutate(Time = as_datetime(Time))

# 2. Find the absolute best false negative day
cat("Scanning for the optimal false negative case...\n")
fn_summary <- regime6_cases %>%
    mutate(
        Date = as.Date(Time),
        G_calc = Dh + Bn * cos(SZA * pi / 180),
        Closure_Error = abs(Gh - G_calc) / Gh * 100
    ) %>%
    filter(SZA <= 85, Gh > 10) %>%
    group_by(Station, Date) %>%
    summarise(
        FN_Count = sum(Flag_New == TRUE & Flag_ERL == FALSE, na.rm = TRUE),
        Mean_Error_FN = mean(Closure_Error[Flag_New == TRUE & Flag_ERL == FALSE], na.rm = TRUE),
        Max_Error_FN = max(Closure_Error[Flag_New == TRUE & Flag_ERL == FALSE], na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # At least one point caught by New and missed by ERL
    filter(FN_Count >= 1) %>%
    # Rank by the single worst closure error missed by ERL
    arrange(desc(Max_Error_FN))

cat("\nTop 20 False Negative Candidates (Ranked by Max Error of Missed Points):\n")
print(head(fn_summary, 20))

# Automatically select the best candidate (prefer station GUR for visual stability if top)
target_stn <- fn_summary$Station[1]
target_date <- fn_summary$Date[1]

cat(sprintf("\n=> Automatically selected optimal day: Station %s on %s\n\n", toupper(target_stn), target_date))

target_data <- regime6_cases %>%
    filter(Station == target_stn, as.Date(Time) == target_date)

if (nrow(target_data) == 0) {
    stop("Target data not found. Please check station and date.")
}

# 3. Calculate Closure and Components
cat("Calculating closure relationships...\n")
plot_data <- target_data %>%
    mutate(
        Hour = hour(Time) + minute(Time) / 60,
        G_calc = Dh + Bn * cos(SZA * pi / 180),
        Closure_Error = abs(Gh - G_calc) / Gh * 100
    ) %>%
    filter(SZA <= 85)

# 4. Generate Visualizations
cat("Generating dual-panel visualization...\n")

# Top Panel: Closure Error Time Series
p_top <- ggplot(plot_data, aes(x = Hour)) +
    geom_line(aes(y = Closure_Error, color = "Relative closure error"), linewidth = line.size, na.rm = TRUE) +
    geom_point(aes(y = Closure_Error, fill = "Relative closure error"), shape = 21, color = "white", size = point.size * 1.25, stroke = 0, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(
        values = c("Relative closure error" = "black"),
        labels = c("Relative closure error" = "Relative closure error")
    ) +
    scale_fill_manual(
        values = c("Relative closure error" = colorblind_pal()(8)[2]),
        labels = c("Relative closure error" = "Relative closure error")
    ) +
    labs(
        x = NULL,
        y = "Closure error [%]"
    ) +
    scale_x_continuous(limits = c(0, 24), expand = c(0, 0), breaks = seq(0, 24, 6)) +
    guides(
        color = "none",
        fill = guide_legend(ncol = 1, order = 1, override.aes = list(size = 3))
    ) +
    theme_bw(base_size = plot.size, base_family = "Times") +
    theme(
        plot.margin = unit(c(0.2, 0.4, 0, 0.0), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = c(0.25, 0.6),
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size = plot.size),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(0.1, 0, 0, 0, "lines"),
        legend.spacing.y = unit(0.05, "lines"),
        legend.box.margin = margin(-0.7, 0, 0, 0, "lines")
    )

# Bottom Panel: Gh vs Limits (Showcasing flags based on global irradiance)
# Create long format for limits
limits_long <- plot_data %>%
    select(Hour, Gh, Glim_ERL, Glim_New, Flag_New) %>%
    tidyr::pivot_longer(cols = c(Glim_ERL, Glim_New), names_to = "Limit_Type", values_to = "Limit_Value") %>%
    mutate(Limit_Type = ifelse(Limit_Type == "Glim_ERL", "ERL", "New limit")) %>%
    filter(!is.na(Gh), !is.na(Limit_Value), !is.na(Flag_New))

p_bot <- ggplot(limits_long, aes(x = Hour)) +
    geom_point(aes(y = Gh, fill = Flag_New), shape = 21, color = "white", size = point.size * 1.25, stroke = 0, na.rm = TRUE) +
    geom_line(aes(y = Limit_Value, linetype = Limit_Type), color = "gray30", linewidth = line.size, na.rm = TRUE) +
    scale_fill_manual(
        values = c("FALSE" = "black", "TRUE" = colorblind_pal()(8)[3]),
        labels = c("FALSE" = "Passed", "TRUE" = "Flagged by New"),
        na.translate = FALSE
    ) +
    labs(
        x = "Hour of day (UTC)",
        y = expression(paste(italic(G)[italic(h)], " [W ", m^-2, "]")),
        fill = "New test status",
        linetype = "Limit parameter"
    ) +
    scale_x_continuous(limits = c(0, 24), expand = c(0, 0), breaks = seq(0, 24, 6)) +
    guides(
        fill = guide_legend(position = "inside", ncol = 1, order = 1, override.aes = list(size = 3)),
        linetype = guide_legend(nrow = 1, order = 2)
    ) +
    theme_bw(base_size = plot.size, base_family = "Times") +
    theme(
        plot.margin = unit(c(0.2, 0.4, 0, 0.0), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "gray80", linewidth = line.size / 2),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.position = "bottom",
        legend.justification.inside = c(0.15, 0.7),
        legend.box = "vertical",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times", size = plot.size),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        legend.key.size = unit(0.5, "cm"),
        legend.margin = margin(0.1, 0, 0, 0, "lines"),
        legend.spacing.y = unit(0.05, "lines"),
        legend.box.margin = margin(-0.5, 0, 0, 0, "lines")
    )

# Combine plots with categorical alignment
cat("Aligning and saving plot...\n")
library(gtable)
library(grid)

g_top <- ggplotGrob(p_top)
g_bot <- ggplotGrob(p_bot)

# Align widths
max_width <- unit.pmax(g_top$widths, g_bot$widths)
g_top$widths <- max_width
g_bot$widths <- max_width

# Arrange with specific heights
combined_plot <- gridExtra::arrangeGrob(g_top, g_bot, ncol = 1, heights = c(1, 1.3))

# Save
setwd(dir0)
ggsave("tex/FalseNegative_Regime6.pdf", plot = combined_plot, width = 85, height = 90, units = "mm", dpi = 300)

# Print Summary
summary_stats <- plot_data %>%
    summarise(
        Mean_Closure_Error_Pct = mean(Closure_Error, na.rm = TRUE),
        Max_Closure_Error_Pct = max(Closure_Error, na.rm = TRUE),
        Total_Points = n(),
        Flagged_Points = sum(Flag_New, na.rm = TRUE)
    )

cat("\nAnalysis Summary:\n")
print(summary_stats)
cat("Done.\n")
