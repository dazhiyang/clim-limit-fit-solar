################################################################################
# This code is written by Zhiwen Wang and Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: Stevenwangzw@gmail.com; yangdazhi.nus@gmail.com
#
# Generates a 7x3 scatterplot grid visualizing QoS Zenith angle limits
# across all dynamically identified clusters. Depicts both the legacy ERL
# limits and the locally optimized parametric ERL curves.
################################################################################

# Clear workspace and load libraries
rm(list = ls(all = TRUE))
libs <- c("dplyr", "ggplot2", "lubridate", "scattermore", "ggthemes", "readr", "data.table")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# Global Input and Functions
################################################################################

# Path Settings
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
dir_data <- file.path(dir0, "Data")
dir_tex <- file.path(dir0, "tex")

# Plotting typography settings based on previous modules
plot.size <- 8

################################################################################
# 1. Load Cluster Structure & Limits
################################################################################
cat("Loading cluster definitions and custom optimizations...\n")

cluster_file <- file.path(dir_data, "cluster.csv")
if (!file.exists(cluster_file)) {
    stop("cluster.csv not found!")
}
clusters <- read.csv(cluster_file, stringsAsFactors = FALSE)

params_file <- file.path(dir_data, "cluster_parameter.csv")
if (!file.exists(params_file)) {
    stop("cluster_parameter.csv not found!")
}
params <- read.csv(params_file, stringsAsFactors = FALSE)

num_clusters <- nrow(params)

################################################################################
# 2. Extract Subsampled Scatter Data
################################################################################
cat("Collecting station point data for visualization...\n")

# We need to extract the exact identical dataset used for fitting in 2.1
data_list <- list()

groups <- split(clusters$stn, clusters$Group)

set.seed(123) # Reproducibility exactly mirroring 2.1's sampling strategy

for (i in seq_len(num_clusters)) {
    c_name <- params$Cluster_Name[i]
    stns <- groups[[c_name]]

    cat(sprintf("Collecting and identically sampling all actual %d stations for Cluster %s...\n", length(stns), c_name))

    cluster_df_list <- lapply(stns, function(s) {
        file_path <- file.path(dir_data, "BSRNtxt", paste0(s, ".txt"))
        if (file.exists(file_path)) {
            df <- fread(file_path, header = TRUE, sep = "\t", data.table = FALSE) %>%
                mutate(Station = s) %>%
                dplyr::select("Time", "Gh", "Dh", "Bn", "Z", "ETR", "kt", "kd", "kb", "Station") %>%
                filter(
                    Z < 90,
                    !is.na(kt), !is.na(kb), !is.na(kd),
                    !is.na(Gh), !is.na(Dh), !is.na(Bn),
                    !is.nan(kt), !is.nan(kb), !is.nan(kd)
                )
            return(df)
        } else {
            return(NULL)
        }
    })

    c_df <- bind_rows(cluster_df_list)

    if (nrow(c_df) == 0) {
        cat(sprintf("Warning: No viable points for Cluster %s. Skipping.\n", c_name))
        next
    }

    # Data sampling via 1-degree Zenith Bins (Top 5% for Gh, Bn, Dh) exactly mirroring 2.1
    c_df_binned <- c_df %>% mutate(Z_bin = floor(Z))

    top_gh <- c_df_binned %>%
        group_by(Z_bin) %>%
        slice_max(order_by = Gh, prop = 0.05, with_ties = FALSE) %>%
        ungroup()
    top_bn <- c_df_binned %>%
        group_by(Z_bin) %>%
        slice_max(order_by = Bn, prop = 0.05, with_ties = FALSE) %>%
        ungroup()
    top_dh <- c_df_binned %>%
        group_by(Z_bin) %>%
        slice_max(order_by = Dh, prop = 0.05, with_ties = FALSE) %>%
        ungroup()

    # Combine uniquely
    # 合并并去重
    random_samp <- c_df_binned %>%
        group_by(Z_bin) %>%
        slice_sample(prop = 0.1) %>%
        ungroup()

    # Further random sampling for visualization clarity
    # 为了清晰可视化，进一步随机采样
    c_df <- bind_rows(top_gh, top_bn, top_dh, random_samp) %>%
        distinct() %>%
        dplyr::select(-Z_bin) %>%
        slice_sample(prop = 0.2)

    # Apply identical Physical Limit boundaries before extracting columns
    # 物理可能极限(PPL)过滤
    c_df <- c_df %>%
        mutate(cosZ = pmax(0, cos(Z * pi / 180))) %>%
        filter(
            Gh >= -4 & Gh <= ETR * 1.5 * (cosZ^1.2) + 100,
            Dh >= -4 & Dh <= ETR * 0.95 * (cosZ^1.2) + 50,
            Bn >= -4 & Bn <= ETR
        ) %>%
        dplyr::select(-cosZ)

    # Finally generate the structural components required for plotting
    df <- c_df %>%
        mutate(
            Cluster = factor(c_name, levels = params$Cluster_Name),

            # Old Empirical Limits (Legacy BSRN)
            Gh_lim_old = 1.2 * ETR * (cos(Z * pi / 180))^1.2 + 50,
            Bn_lim_old = 0.95 * ETR * (cos(Z * pi / 180))^0.2 + 10,
            Dh_lim_old = 0.75 * ETR * (cos(Z * pi / 180))^1.2 + 30,

            # New Fitted Parametric Limits using ETR logic
            Gh_lim_new = params$gh_a[i] * ETR * (cos(Z * pi / 180))^params$gh_b[i] + params$gh_c[i],
            Bn_lim_new = params$bn_a[i] * ETR * (cos(Z * pi / 180))^params$bn_b[i] + params$bn_c[i],
            Dh_lim_new = params$dh_a[i] * ETR * (cos(Z * pi / 180))^params$dh_b[i] + params$dh_c[i]
        )

    data_list[[i]] <- df
}

master_data <- bind_rows(data_list)

if (nrow(master_data) == 0) {
    stop("Critical: Failed to load any station data for plotting.")
}

################################################################################
# 3. Shape Data for 7x3 Grid
################################################################################
cat("Formatting visual panels into tidy dimensions...\n")

# We restructure Gh, Bn, Dh into one variable mapped against its respective limits
data.plot <- rbind(
    data.frame(
        Z = master_data$Z,
        Irradiance = master_data$Gh,
        Lim_Old = master_data$Gh_lim_old,
        Lim_New = master_data$Gh_lim_new,
        Cluster = master_data$Cluster,
        Component = "GHI"
    ),
    data.frame(
        Z = master_data$Z,
        Irradiance = master_data$Bn,
        Lim_Old = master_data$Bn_lim_old,
        Lim_New = master_data$Bn_lim_new,
        Cluster = master_data$Cluster,
        Component = "BNI"
    ),
    data.frame(
        Z = master_data$Z,
        Irradiance = master_data$Dh,
        Lim_Old = master_data$Dh_lim_old,
        Lim_New = master_data$Dh_lim_new,
        Cluster = master_data$Cluster,
        Component = "DHI"
    )
)

data.plot$Component <- factor(
    data.plot$Component,
    levels = c("GHI", "BNI", "DHI"),
    labels = c("italic(G)[italic(h)]", "italic(B)[italic(n)]", "italic(D)[italic(h)]")
)

################################################################################
# 4. Generate Plot Grid
################################################################################
cat("Rendering scatter plots and mathematical curves mapping 7x3...\n")

# Color palette safely handling dynamic clusters
if (num_clusters <= 8) {
    cluster_colors <- ggthemes::colorblind_pal()(8)[c(2:8)]
} else {
    cluster_colors <- viridis::viridis(num_clusters)
}

# The user explicitly asked for "two limits" and "7 rows by 3 columns"
p <- ggplot(data.plot) +
    # The actual scatter points
    geom_scattermore(aes(x = Z, y = Irradiance, color = Cluster), pointsize = 0, alpha = 0.5) +
    # Legacy limit explicitly drawn in a contrasting style (e.g. gray line)
    geom_scattermore(aes(x = Z, y = Lim_Old), color = "gray50", pointsize = 0) +
    # The new clustered limit drawn visibly over it
    geom_scattermore(aes(x = Z, y = Lim_New), color = "black", pointsize = 0) +
    # Grid construction
    facet_grid(Cluster ~ Component, labeller = label_parsed, scales = "free_y") +
    scale_x_continuous(name = expression(paste("Zenith angle [", degree, "]"))) +
    scale_y_continuous(name = expression(paste("Irradiance [W ", m^-2, "]"))) +
    scale_color_manual(values = cluster_colors) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.1, 0.2, 0, 0.0), "lines"),
        text = element_text(family = "Times", size = plot.size),
        axis.text = element_text(family = "Times", size = plot.size - 1),
        panel.grid.minor = element_blank(),
        strip.text = element_text(family = "Times", size = plot.size),
        panel.spacing = unit(0.5, "lines"),
        legend.position = "none" # Subplots perfectly denote clusters by row automatically
    )

################################################################################
# 5. Output Graphics
################################################################################
if (!dir.exists(dir_tex)) dir.create(dir_tex, recursive = TRUE)
output_pdf <- file.path(dir_tex, "ScatterVis.pdf")

# We scale the image taller vertically since it is 7 rows vs 1 row in 4.1
ggsave(
    filename = output_pdf, plot = p,
    width = 160, height = 170, units = "mm"
)

cat(sprintf("\nSuccessfully saved the grouped QoS scatter grid to: %s\n", output_pdf))
