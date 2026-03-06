################################################################################
# 4.5 World Map Visualization
# Depict the cluster distribution globally leveraging the aesthetic style
# from Fig.KGC_BSRN_laylout.R
################################################################################
rm(list = ls(all = TRUE))

libs <- c(
    "dplyr", "ggplot2", "ggthemes", "RColorBrewer", "tidyverse",
    "ggrepel", "sf", "sp", "rnaturalearth", "raster",
    "rnaturalearthdata"
)
invisible(lapply(libs, function(x) {
    if (!require(x, character.only = TRUE, quietly = TRUE)) {
        warning(paste("Package", x, "not installed, please install it."))
    }
}))

################################################################################
# Global variables and paths
################################################################################
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
data_dir <- file.path(dir0, "Data")
out_dir <- file.path(dir0, "tex")
if (!dir.exists(out_dir)) dir.create(out_dir)

plot.size <- 8
line.size <- 0.2
point.size <- 0.05
legend.size <- 0.4
text.size <- plot.size * 5 / 14

################################################################################
# 1. Load Data
################################################################################
cat("Loading cluster data...\n")
cluster_file <- file.path(data_dir, "cluster.csv")
if (!file.exists(cluster_file)) {
    stop("cluster.csv not found! Please run 1.4.PCA_hclust.R first.")
}
bsrn <- readr::read_csv(cluster_file)
bsrn$Group <- as.factor(bsrn$Group)

# Determine number of clusters properly
num_clusters <- length(unique(bsrn$Group))

# Dynamically construct matching cluster palette combining colorblind + viridis natively if needed
if (num_clusters <= 8) {
    cluster_colors <- ggthemes::colorblind_pal()(8)[c(2:8)]
} else {
    cluster_colors <- viridis::viridis(num_clusters)
}

# Optional Coastline Geometry (must use countries to get closed polygons for fill)
coastlines_sim2_df <- ne_countries(scale = "medium", returnclass = "sf")

################################################################################
# 2. Render Geographic World Map
################################################################################
cat("Rendering geographic world map applying laylout aesthetics...\n")

p <- ggplot() +
    # Render clean background map boundaries using classic SF with shaded land area
    geom_sf(data = coastlines_sim2_df, fill = "gray80", color = "gray60", size = line.size) +
    # Map actual points styled like laylout.R
    geom_point(data = bsrn, aes(lon, lat, color = Group), size = point.size * 50) +
    # Set typography labels utilizing repulsion directly
    ggrepel::geom_text_repel(
        data = bsrn,
        aes(lon, lat, label = stn),
        family = "Times",
        max.overlaps = 20,
        size = text.size,
        segment.size = 0.2,
        point.padding = 0,
        force = 0.5
    ) +
    labs(color = "Regime") +
    scale_color_manual(values = cluster_colors) +
    coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) +
    theme_bw() +
    theme(
        plot.margin = unit(c(0.1, 0.1, 0, 0), "lines"),
        panel.spacing = unit(0.02, "lines"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        text = element_text(family = "Times", size = plot.size),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(family = "Times", size = plot.size),
        legend.key.height = unit(0.2, "lines"),
        legend.key.width = unit(0.4, "lines"),
        legend.box.background = element_rect(fill = "transparent", color = "transparent"),
        legend.background = element_rect(fill = "transparent", colour = "transparent"),
        legend.position = "bottom",
        legend.title = element_text(family = "Times", size = plot.size),
        legend.box.margin = unit(c(-1.0, 0, 0, 0), "lines"),
        plot.background = element_rect(fill = "transparent", colour = "transparent")
    ) +
    guides(
        color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 3))
    )

################################################################################
# 3. Save the Plot
################################################################################
cat("Saving Map Output to tex directory...\n")

pdf_file <- file.path(out_dir, "WorldMap.pdf")
ggsave(filename = pdf_file, plot = p, scale = 1, width = 160, height = 85, unit = "mm", dpi = 300)

cat(sprintf("World Map generation complete! Saved to %s\n", pdf_file))
