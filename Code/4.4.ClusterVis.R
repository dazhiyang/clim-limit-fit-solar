################################################################################
# 4.4 Cluster Visualization
# This script processes the clustering distance matrix into a paper-grade Heatmap
################################################################################
rm(list = ls(all = TRUE))

libs <- c("dplyr", "tibble", "ggplot2", "ComplexHeatmap", "viridis", "circlize", "ggthemes")
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
out_dir <- file.path(dir0, "tex") # Updated to output exactly to tex folder
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################################################
# 1. Load Data
################################################################################
cat("Loading PCA and clustering data...\n")
feature_file <- file.path(data_dir, "extraction_features.csv")
data <- as_tibble(read.csv(feature_file, stringsAsFactors = FALSE))
data <- na.omit(data)
vars <- data[, 4:16]

# Recalculate PCA and distance matrix just like 1.4b
pca_fit <- prcomp(vars, center = TRUE, scale. = TRUE)
pca_scores <- as.data.frame(pca_fit$x)

percentage_of_var <- round(pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100, 3)
threshold <- 90
total <- 0
i <- 1
while (total < threshold && i <= length(percentage_of_var)) {
    total <- total + percentage_of_var[i]
    i <- i + 1
}
n <- i - 1
pcadata.new <- pca_scores[, 1:n]

dist_matrix_obj <- dist(pcadata.new, method = "euclidean")
dist_matrix <- as.matrix(dist_matrix_obj)
rownames(dist_matrix) <- data$stn
colnames(dist_matrix) <- data$stn

################################################################################
# 2. Get the clustering that heatmap uses
################################################################################
cat("Processing Heatmap order based on optimal tree cuts...\n")

# Load cluster assignments directly from cluster.csv instead of recalculating
cat("Loading cluster assignments from cluster.csv...\n")
cluster_file <- file.path(data_dir, "cluster.csv")
if (!file.exists(cluster_file)) {
    stop("cluster.csv not found! Please run 1.4.PCA_hclust.R first.")
}
cluster_data <- read.csv(cluster_file, stringsAsFactors = FALSE)

# Ensure the ordering matches the dataset exactly
cluster_data <- cluster_data[match(data$stn, cluster_data$stn), ]
cluster_assignments_raw <- setNames(cluster_data$Group, cluster_data$stn)
num_clusters <- length(unique(cluster_assignments_raw))

cat(sprintf("Loaded %d optimal clusters from cluster.csv\n", num_clusters))

# Generate the hierarchical clustering dendrogram tree
hc_fit <- hclust(dist_matrix_obj, method = "ward.D2")

cat("Extracting visual Heatmap slice sequence...\n")
# We use ComplexHeatmap's engine to tell us exactly how it visually sorts the slices!
# We do a silent draw to extract the true visual slice sequence from left-to-right.
pdf(NULL)
ht_temp <- Heatmap(dist_matrix, cluster_columns = as.dendrogram(hc_fit), column_split = num_clusters)
ht_temp <- draw(ht_temp)
co <- column_order(ht_temp)
dev.off()

# Extract the TRUE cluster ID (from cluster.csv) for each visual block left-to-right
visual_group_order <- sapply(co, function(idx) {
    # Get the station name from the first item in this visual block
    stn_name <- colnames(dist_matrix)[idx[1]]
    as.character(cluster_assignments_raw[stn_name])
})

################################################################################
# 3. Plot the heatmap using ComplexHeatmap
################################################################################
cat("Rendering paper-grade Heatmap...\n")

# Create quantile-based color function
create_quantile_colors <- function(matrix, n_colors = 256) {
    vals <- as.vector(matrix)
    quantiles <- quantile(vals, probs = seq(0, 1, length.out = n_colors), na.rm = TRUE)
    colors <- viridis::viridis(n_colors)
    circlize::colorRamp2(quantiles, colors)
}

col_fun <- create_quantile_colors(dist_matrix)

# Safely dynamically allocate colors matching 4.5.WorldMap.R safely
if (num_clusters <= 8) {
    base_colors <- ggthemes::colorblind_pal()(8)[c(2:8)]
} else {
    base_colors <- viridis::viridis(num_clusters)
}

# The dendrogram blocks are assembled left-to-right, so we map our colors
# and labels to correspond with the sequence the groups appear in that tree.
group_labels <- paste0(visual_group_order)
group_colors <- base_colors[as.numeric(visual_group_order)]

plot.size <- 7 # Scaled down for 85mm canvas

# Render PDF for the Paper Version (85mm x 85mm)
pdf(file.path(out_dir, "HierarchicalResult.pdf"), width = 85 / 25.4, height = 80 / 25.4)
ht_paper <- Heatmap(dist_matrix,
    name = "WD",
    col = col_fun,
    # Change dendrogram line size
    row_dend_width = unit(1.3, "cm"), # Width of row dendrogram
    column_dend_height = unit(1.3, "cm"), # Height of column dendrogram
    row_dend_gp = gpar(lwd = 0.3), # Row dendrogram line width
    column_dend_gp = gpar(lwd = 0.3), # Column dendrogram line width
    # Disable row and column names
    show_row_names = FALSE,
    show_column_names = FALSE,
    # Top annotation (columns)
    top_annotation = HeatmapAnnotation(
        foo = anno_block(
            gp = gpar(fill = group_colors),
            labels = group_labels,
            height = unit(0.3, "cm"),
            labels_gp = gpar(col = "white", fontsize = plot.size - 1, fontfamily = "Times")
        )
    ),
    # Left annotation (rows) - SYMMETRIC to top
    left_annotation = rowAnnotation(
        foo = anno_block(
            gp = gpar(fill = group_colors),
            labels = group_labels,
            width = unit(0.3, "cm"),
            labels_gp = gpar(col = "white", fontsize = plot.size - 1, fontfamily = "Times")
        )
    ),
    cluster_rows = as.dendrogram(hc_fit),
    cluster_columns = as.dendrogram(hc_fit),
    column_split = num_clusters,
    row_split = num_clusters,
    column_title = NULL,
    row_title = NULL,
    # Same font settings for both axes
    column_names_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times"),
    row_names_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times"),
    # Legend settings
    heatmap_legend_param = list(
        title_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times", fontface = "plain"),
        labels_gp = gpar(fontsize = plot.size - 1, fontfamily = "Times"),
        legend_height = unit(4, "cm")
    )
)
draw(ht_paper)
dev.off()

cat("Paper Heatmap saved to tex/HierarchicalResult.pdf\n")
cat("Heatmap generation complete!\n")
