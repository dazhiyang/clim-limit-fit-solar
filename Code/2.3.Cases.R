################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
#
# Aim: Identify specific cases where ERL and New Clustered Limits disagree.
# 1. Regime 2 (SSC/Polar): Points rejected by ERL but passed by New Limits.
# 2. Regime 6 (HAA/Polluted): Points rejected by New Limits but passed by ERL.
################################################################################

# nolint start
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "data.table")
invisible(lapply(libs, require, character.only = TRUE))

# Global Variables
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data")
setwd(dir_data)

# Load Meta and Parameters
cluster_meta <- read.csv("cluster.csv", header = TRUE, sep = ",")
params <- read.csv("cluster_parameter.csv", header = TRUE, sep = ",")
cluster_names <- c("TSS", "PSP", "ASD", "MCM", "TML", "HAA", "LLA")
groups <- split(cluster_meta$stn, cluster_meta$Group)

# Function to extract discrepancy cases for a specific regime
extract_cases <- function(group_idx, case_type) {
    cat(sprintf(
        "\nProcessing Regime %d (Cluster %s) for Case: %s\n",
        group_idx, cluster_names[group_idx], case_type
    ))

    stns <- groups[[group_idx]]
    case_data_list <- list()

    # Parameters for this cluster
    p <- params[group_idx, ]

    for (s in stns) {
        file_path <- file.path("BSRNtxt", paste0(s, ".txt"))
        if (file.exists(file_path)) {
            cat(sprintf("  Reading %s...\n", s))
            df <- data.table::fread(file_path, header = TRUE, sep = "\t", colClasses = list(character = "Time"))

            # Note: No subsampling here, as we want to find specific cases.

            df <- df %>%
                mutate(Time = as_datetime(Time)) %>%
                rename(SZA = Z) %>%
                filter(SZA <= 90) %>%
                filter(Gh > 0)

            # Apply both limits and round radiation to integers
            df <- df %>%
                mutate(
                    Gh = round(Gh),
                    Dh = round(Dh),
                    Bn = round(Bn),
                    Glim_ERL = round(1.2 * ETR * (cos(SZA * pi / 180))^1.2 + 50),
                    Dlim_ERL = round(0.75 * ETR * (cos(SZA * pi / 180))^1.2 + 30),
                    Blim_ERL = round(0.95 * ETR * (cos(SZA * pi / 180))^0.2 + 10),
                    Glim_New = round(p$gh_a * ETR * (cos(SZA * pi / 180))^p$gh_b + p$gh_c),
                    Dlim_New = round(p$dh_a * ETR * (cos(SZA * pi / 180))^p$dh_b + p$dh_c, 2), # Note: Keep 2 digits for New limits if they are very small? User said integers. I will use round() for all.
                    Blim_New = round(p$bn_a * ETR * (cos(SZA * pi / 180))^p$bn_b + p$bn_c)
                ) %>%
                mutate(
                    Dlim_New = round(p$dh_a * ETR * (cos(SZA * pi / 180))^p$dh_b + p$dh_c)
                ) %>%
                mutate(
                    Flag_ERL = (Gh > Glim_ERL) | (Dh > Dlim_ERL) | (Bn > Blim_ERL),
                    Flag_New = (Gh > Glim_New) | (Dh > Dlim_New) | (Bn > Blim_New)
                )

            # Identify discrepancy dates
            if (case_type == "ERL_only") {
                # Rejected by ERL, Passed by New
                discrepancy_indices <- which(df$Flag_ERL == TRUE & df$Flag_New == FALSE)
            } else if (case_type == "New_only") {
                # Rejected by New, Passed by ERL
                discrepancy_indices <- which(df$Flag_New == TRUE & df$Flag_ERL == FALSE)
            }

            if (length(discrepancy_indices) > 0) {
                # Get unique dates with discrepancies
                discrepancy_dates <- unique(as.Date(df$Time[discrepancy_indices]))

                # Filter full day data for those dates
                full_day_data <- df %>%
                    filter(as.Date(Time) %in% discrepancy_dates) %>%
                    mutate(Station = s)

                case_data_list[[s]] <- full_day_data
            }
        }
    }

    return(bind_rows(case_data_list))
}

# 1. Aim 1: Regime 2 (SSC/Polar) - ERL Rejected but New Passed
regime2_cases <- extract_cases(2, "ERL_only")

# 2. Aim 2: Regime 6 (HAA/Polluted) - New Rejected but ERL Passed
regime6_cases <- extract_cases(6, "New_only")

# Save Results
cat("\nSaving discrepancy cases to files...\n")
write.table(regime2_cases, file = "Regime2_ERL_only.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(regime6_cases, file = "Regime6_New_only.txt", sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nSummary:")
cat(sprintf("\nRegime 2 (ERL rejected, New passed): %d points", nrow(regime2_cases)))
cat(sprintf("\nRegime 6 (New rejected, ERL passed): %d points", nrow(regime6_cases)))
cat("\n\nFinished.\n")

# nolint end
