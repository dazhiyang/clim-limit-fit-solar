################################################################################
# This code implements PCA followed by Ward's Hierarchical Clustering
# It includes Silhouette Analysis and Latitudinal Validation
################################################################################
# Clear workspace and load necessary packages
rm(list = ls(all = TRUE))
libs <- c("dplyr", "tibble", "ggplot2", "cluster", "factoextra")
invisible(lapply(libs, function(x) {
    if (!require(x, character.only = TRUE, quietly = TRUE)) {
        warning(paste("Package", x, "not installed, please install it."))
    }
}))

options(digits = 5)

################################################################################
# Global variables and paths
################################################################################
dir0 <- "/Volumes/Macintosh Research/Data/BSRN_QC_isolation"
data_dir <- file.path(dir0, "Data")
out_dir <- file.path(dir0, "tex")
if (!dir.exists(out_dir)) dir.create(out_dir)

################################################################################
# 1. Data Processing and PCA
################################################################################
cat("Reading extracted features...\n")
feature_file <- file.path(data_dir, "extraction_features.csv")
data <- as_tibble(read.csv(feature_file, stringsAsFactors = FALSE))
data <- na.omit(data) # Remove rows with NA to prevent distance errors
vars <- data[, 4:16]

cat("Performing PCA...\n")
pca_fit <- prcomp(vars, center = TRUE, scale. = TRUE)
pca_scores <- pca_fit$x

# Determine number of PCs explaining >95% variance
percentage_of_var <- round(pca_fit$sdev^2 / sum(pca_fit$sdev^2) * 100, 3)
threshold <- 90
total <- 0
i <- 1
while (total < threshold && i <= length(percentage_of_var)) {
    total <- total + percentage_of_var[i]
    i <- i + 1
}
n <- i - 1
cat(sprintf("Using first %d PCs explaining %.2f%% of total variance.\n", n, total))

pcadata.new <- as.data.frame(pca_scores)[, 1:n]

################################################################################
# 2. Hierarchical Clustering (Ward's Method)
################################################################################
cat("Performing Ward's Hierarchical Clustering...\n")
# Calculate Euclidean distance
dist_matrix <- dist(pcadata.new, method = "euclidean")
# Ward's method (ward.D2 is mathematically the original Ward's criterion)
hc_fit <- hclust(dist_matrix, method = "ward.D2")

# Determine optimal k using Average Silhouette Width
k_values <- 4:10
sil_widths <- numeric(length(k_values))

for (idx in seq_along(k_values)) {
    k <- k_values[idx]
    cluster_assignments <- cutree(hc_fit, k = k)
    sil <- silhouette(cluster_assignments, dist_matrix)
    sil_widths[idx] <- mean(sil[, "sil_width"])
}

# Select k that maximizes average silhouette width
best_k <- k_values[which.max(sil_widths)]
cat(sprintf("Optimal number of clusters based on Average Silhouette Width: %d\n", best_k))

# Cut tree at optimal k
final_clusters <- cutree(hc_fit, k = best_k)

################################################################################
# 3. Save Results
################################################################################
cluster_df <- tibble(
    stn = data$stn,
    lat = data$lat,
    lon = data$lon,
    Group = as.factor(final_clusters)
)
output_file <- file.path(data_dir, "cluster.csv")
write.csv(cluster_df, output_file, row.names = FALSE)
cat(sprintf("\nClustering complete! Results saved to %s\n", output_file))
