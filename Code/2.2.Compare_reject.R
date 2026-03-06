################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
#
# Compare new cluster-based limits vs old ERL limits using Component Sum
################################################################################
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "xtable", "foreach", "data.table")
invisible(lapply(libs, require, character.only = TRUE))

options(digits = 5)

################################################################################
# Global Variables & Functions
################################################################################
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data")
dir_tex <- file.path(dir0, "tex")

# Function to calculate ETR correction factor
get_RE <- function(Tm) {
  doy <- lubridate::yday(Tm)
  da <- (2 * pi / 365) * (doy - 1)
  re <- 1.000110 + 0.034221 * cos(da) + 0.001280 * sin(da) + 0.00719 * cos(2 * da) + 0.000077 * sin(2 * da)
  return(re)
}

################################################################################
# Load Data & Perform Cluster-based Analysis
################################################################################
setwd(dir_data)
cluster_meta <- read.csv("cluster.csv", header = TRUE, sep = ",")
params <- read.csv("cluster_parameter.csv", header = TRUE, sep = ",")

cluster_names <- c("TSS", "PSP", "ASD", "MCM", "TML", "HAA", "LLA")
groups <- split(cluster_meta$stn, cluster_meta$Group)

setwd(file.path(dir_data, "BSRNtxt"))

summary_list <- list()

# Loop through all seven clusters
for (group_index in 1:7) {
  cat(sprintf("\nProcessing Group %d: %s...\n", group_index, cluster_names[group_index]))

  stns <- groups[[group_index]]
  group_data_list <- list()

  for (s in stns) {
    file_name <- paste0(s, ".txt")
    if (file.exists(file_name)) {
      df <- data.table::fread(file_name, header = TRUE, sep = "\t")
      # Sample 1/10th
      set.seed(123)
      n <- nrow(df)
      sample_size <- ceiling(n / 10)
      df <- df[sample(1:n, sample_size), ]

      df <- df %>%
        mutate(Time = lubridate::ymd_hms(Time)) %>%
        mutate(RE = get_RE(Time)) %>%
        mutate(ETR = 1361.1 * RE) %>%
        rename(SZA = Z) %>%
        mutate(Station = s) %>%
        dplyr::select(Station, Time, Gh, Dh, Bn, SZA, ETR) %>%
        filter(SZA <= 90) %>%
        filter(Gh > 0) # Basic filter to avoid div by zero in relative error

      group_data_list[[s]] <- df
    }
  }

  if (length(group_data_list) == 0) {
    cat("Warning: No data for group, skipping...\n")
    next
  }

  all_data <- bind_rows(group_data_list)
  rm(group_data_list) # Free memory

  # Parameters for this cluster
  gh_a <- params$gh_a[group_index]
  gh_b <- params$gh_b[group_index]
  gh_c <- params$gh_c[group_index]
  dh_a <- params$dh_a[group_index]
  dh_b <- params$dh_b[group_index]
  dh_c <- params$dh_c[group_index]
  bn_a <- params$bn_a[group_index]
  bn_b <- params$bn_b[group_index]
  bn_c <- params$bn_c[group_index]

  # Apply Limits
  all_data <- all_data %>%
    mutate(
      Glim_ERL = 1.2 * ETR * (cos(SZA * pi / 180))^1.2 + 50,
      Dlim_ERL = 0.75 * ETR * (cos(SZA * pi / 180))^1.2 + 30,
      Blim_ERL = 0.95 * ETR * (cos(SZA * pi / 180))^0.2 + 10,
      Glim_New = gh_a * ETR * (cos(SZA * pi / 180))^gh_b + gh_c,
      Dlim_New = dh_a * ETR * (cos(SZA * pi / 180))^dh_b + dh_c,
      Blim_New = bn_a * ETR * (cos(SZA * pi / 180))^bn_b + bn_c
    ) %>%
    mutate(
      G_Flag_ERL = (Gh > Glim_ERL),
      D_Flag_ERL = (Dh > Dlim_ERL),
      B_Flag_ERL = (Bn > Blim_ERL),
      Flag_ERL = G_Flag_ERL | D_Flag_ERL | B_Flag_ERL,
      G_Flag_New = (Gh > Glim_New),
      D_Flag_New = (Dh > Dlim_New),
      B_Flag_New = (Bn > Blim_New),
      Flag_New = G_Flag_New | D_Flag_New | B_Flag_New
    )

  # Calculate Component Sum Error & BSRN Standard Closure Test
  all_data <- all_data %>%
    mutate(
      Sum_SW = Bn * cos(SZA * pi / 180) + Dh,
      Rel_Error = abs(Gh - Sum_SW) / Gh,

      # BSRN Standard Closure Test (Long & Dutton 2002)
      # 8% for SZA < 75, 15% for 75 < SZA < 93, only if Sum_SW > 50
      BSRN_Closure_Limit = ifelse(SZA < 75, 0.08, 0.15),
      Fail_BSRN_Closure = (Sum_SW > 50) & (Rel_Error > BSRN_Closure_Limit)
    )

  # Categorize Points
  all_data <- all_data %>%
    mutate(
      Category = case_when(
        !Flag_ERL & !Flag_New ~ "Clean",
        !Flag_ERL & Flag_New ~ "Caught by New (Missed by ERL)",
        Flag_ERL & !Flag_New ~ "Caught by ERL (Missed by New)",
        Flag_ERL & Flag_New ~ "Caught by Both"
      )
    )

  # Summary Statistics for this cluster: Head-to-head comparison
  # We want to know: Of the points the NEW limits flag, how many fail the standard BSRN Closure Test?
  new_flags <- all_data %>% filter(Flag_New)
  old_flags <- all_data %>% filter(Flag_ERL)

  summary_stats <- data.frame(
    Cluster = cluster_names[group_index],
    Total_Points = nrow(all_data),

    # Absolute Counts
    G_Old_Rej = sum(all_data$G_Flag_ERL, na.rm = TRUE),
    B_Old_Rej = sum(all_data$B_Flag_ERL, na.rm = TRUE),
    D_Old_Rej = sum(all_data$D_Flag_ERL, na.rm = TRUE),
    Ext_Old_Rej = nrow(old_flags),
    G_New_Rej = sum(all_data$G_Flag_New, na.rm = TRUE),
    B_New_Rej = sum(all_data$B_Flag_New, na.rm = TRUE),
    D_New_Rej = sum(all_data$D_Flag_New, na.rm = TRUE),
    Ext_New_Rej = nrow(new_flags),

    # Rejection Rates (%)
    G_Old_Rate = (sum(all_data$G_Flag_ERL, na.rm = TRUE) / nrow(all_data)) * 100,
    B_Old_Rate = (sum(all_data$B_Flag_ERL, na.rm = TRUE) / nrow(all_data)) * 100,
    D_Old_Rate = (sum(all_data$D_Flag_ERL, na.rm = TRUE) / nrow(all_data)) * 100,
    Ext_Old_Rate = (nrow(old_flags) / nrow(all_data)) * 100,
    G_New_Rate = (sum(all_data$G_Flag_New, na.rm = TRUE) / nrow(all_data)) * 100,
    B_New_Rate = (sum(all_data$B_Flag_New, na.rm = TRUE) / nrow(all_data)) * 100,
    D_New_Rate = (sum(all_data$D_Flag_New, na.rm = TRUE) / nrow(all_data)) * 100,
    Ext_New_Rate = (nrow(new_flags) / nrow(all_data)) * 100,

    # Total points caught by each limit that fail the official BSRN Closure test
    Old_Fail_BSRN_Closure = sum(old_flags$Fail_BSRN_Closure, na.rm = TRUE),
    New_Fail_BSRN_Closure = sum(new_flags$Fail_BSRN_Closure, na.rm = TRUE),

    # Percentage of rejections that are "Hard Errors" according to BSRN
    Percent_Old_Fail = (sum(old_flags$Fail_BSRN_Closure, na.rm = TRUE) / nrow(old_flags)) * 100,
    Percent_New_Fail = (sum(new_flags$Fail_BSRN_Closure, na.rm = TRUE) / nrow(new_flags)) * 100
  )

  summary_list[[group_index]] <- summary_stats

  rm(all_data) # Free memory before next cluster loop
}

# Save Rejection Summary
# 保存拒绝率汇总
final_summary <- bind_rows(summary_list)
cat("\nFinal Summary Table:\n")
print(as.data.frame(final_summary))

cat("\nSaving rejection summary to Data/rejection_summary.csv...\n")
write.csv(final_summary, file = file.path(dir_data, "rejection_summary.csv"), row.names = FALSE)
cat("Finished.\n")
