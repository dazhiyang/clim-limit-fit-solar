################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# False Positive Analysis (Regime 2)
# Aim: Demonstrate that data rejected by ERL is physically valid by showing
#      excellent closure relationship (G = D + B*cos(Z)).
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

# 1. Load the False Positive Discrepancy Cases for Regime 2
cat("Loading discrepancy cases...\n")
regime2_cases <- data.table::fread("Regime2_ERL_only.txt") %>%
    mutate(Time = as_datetime(Time))

# 2. Find the absolute best false positive day: Highest FP Count with Lowest Closure Error
cat("Scanning for the optimal false positive case...\n")
fp_summary <- regime2_cases %>%
    # Calculate closure relationship for all points first to find the best day
    mutate(
        Date = as.Date(Time),
        G_calc = Dh + Bn * cos(SZA * pi / 180),
        Closure_Error = abs(Gh - G_calc) / Gh * 100
    ) %>%
    # Filter for standard daylight conditions to ensure fair evaluation
    filter(SZA <= 85, Gh > 10) %>%
    group_by(Station, Date) %>%
    summarise(
        FP_Count = sum(Flag_ERL == TRUE & Flag_New == FALSE, na.rm = TRUE),
        Mean_Closure_Error_Pct = mean(Closure_Error, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    # We want a day with excellent closure (< 1.5% average error) but many false positives
    filter(Mean_Closure_Error_Pct < 1.5) %>%
    arrange(desc(FP_Count), Mean_Closure_Error_Pct)

cat("\nTop 10 Optimal False Positive Candidates (Excellent Closure + High Rejection):\n")
print(head(fp_summary, 10))

# Automatically select the absolute best candidate
target_stn <- fp_summary$Station[1]
target_date <- fp_summary$Date[1]

cat(sprintf("\n=> Automatically selected optimal day: Station %s on %s\n\n", toupper(target_stn), target_date))

target_data <- regime2_cases %>%
    filter(Station == target_stn, as.Date(Time) == target_date)

if (nrow(target_data) == 0) {
    stop("Target data not found. Please check station and date.")
}

# 3. Calculate Closure and Components
# G_obs = Gh
# G_calc = Dh + Bn * cos(SZA * pi / 180)
# Closure Error (%) = (G_obs - G_calc) / G_obs * 100
cat("Calculating closure relationships...\n")
plot_data <- target_data %>%
    mutate(
        Hour = hour(Time) + minute(Time) / 60,
        G_calc = Dh + Bn * cos(SZA * pi / 180),
        Closure_Error = abs(Gh - G_calc) / Gh * 100
    ) %>%
    filter(SZA <= 85) # Filter extreme angles for cleaner visualization

# 4. Generate Visualizations
cat("Generating dual-panel visualization...\n")

# Top Panel: Closure Relationship (G_obs vs G_calc)
p_top <- ggplot(plot_data, aes(x = Hour)) +
    geom_line(aes(y = Gh, color = "Observed Gh"), linewidth = line.size, na.rm = TRUE) +
    geom_point(aes(y = G_calc, fill = "Calculated (Dh + Bn*cosZ)"), shape = 21, color = "white", size = point.size * 1.25, stroke = 0, alpha = 0.8, na.rm = TRUE) +
    scale_color_manual(
        values = c("Observed Gh" = "black"),
        labels = c("Observed Gh" = expression(italic(G)[italic(h)]))
    ) +
    scale_fill_manual(
        values = c("Calculated (Dh + Bn*cosZ)" = colorblind_pal()(8)[2]),
        labels = c("Calculated (Dh + Bn*cosZ)" = expression(italic(D)[italic(h)] + italic(B)[italic(n)] * cos * italic(Z)))
    ) +
    labs(
        x = NULL,
        y = expression(paste(italic(G)[italic(h)], " [W ", m^-2, "]"))
    ) +
    scale_x_continuous(limits = c(0, 24), expand = c(0, 0), breaks = seq(0, 24, 6)) +
    guides(
        color = guide_legend(ncol = 1, order = 1),
        fill = guide_legend(ncol = 1, order = 2, override.aes = list(size = 3))
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

# Bottom Panel: Dh vs Limits
# Create long format for limits
limits_long <- plot_data %>%
    select(Hour, Dh, Dlim_ERL, Dlim_New, Flag_ERL) %>%
    tidyr::pivot_longer(cols = c(Dlim_ERL, Dlim_New), names_to = "Limit_Type", values_to = "Limit_Value") %>%
    mutate(Limit_Type = ifelse(Limit_Type == "Dlim_ERL", "ERL", "New limit")) %>%
    filter(!is.na(Dh), !is.na(Limit_Value), !is.na(Flag_ERL))

p_bot <- ggplot(limits_long, aes(x = Hour)) +
    geom_point(aes(y = Dh, fill = Flag_ERL), shape = 21, color = "white", size = point.size * 1.25, stroke = 0, na.rm = TRUE) +
    geom_line(aes(y = Limit_Value, linetype = Limit_Type), color = "gray30", linewidth = line.size, na.rm = TRUE) +
    scale_fill_manual(
        values = c("FALSE" = "black", "TRUE" = colorblind_pal()(8)[3]),
        labels = c("FALSE" = "Passed", "TRUE" = "Flagged by ERL"),
        na.translate = FALSE
    ) +
    labs(
        x = "Hour of day (UTC)",
        y = expression(paste(italic(D)[italic(h)], " [W ", m^-2, "]")),
        fill = "ERL test status",
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

# Combine plots
combined_plot <- gridExtra::grid.arrange(p_top, p_bot, ncol = 1, heights = c(1, 1.3))

# Save
cat("Saving plot...\n")
setwd(dir0)
ggsave("tex/FalsePositive_Regime2.pdf", plot = combined_plot, width = 85, height = 90, units = "mm", dpi = 300)

# Print Summary
summary_stats <- plot_data %>%
    summarise(
        Mean_Closure_Error_Pct = mean(Closure_Error, na.rm = TRUE),
        Max_Closure_Error_Pct = max(Closure_Error, na.rm = TRUE),
        Total_Points = n(),
        Flagged_Points = sum(Flag_ERL, na.rm = TRUE)
    )

cat("\nAnalysis Summary:\n")
print(summary_stats)
cat("Done.\n")
