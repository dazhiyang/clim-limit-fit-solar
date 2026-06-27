################################################################################
# This code implements PCA followed by Ward's Hierarchical Clustering
# It includes Silhouette Analysis and Latitudinal Validation
################################################################################
# Clear workspace and load necessary packages
rm(list = ls(all = TRUE))
libs <- c("dplyr", "tibble", "ggplot2", "cluster", "factoextra", "xtable")
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
revision_dir <- file.path(dir0, "Revision1")
if (!dir.exists(out_dir)) dir.create(out_dir)
if (!dir.exists(revision_dir)) dir.create(revision_dir)

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
# 3. Within- and Inter-Cluster Distances
################################################################################
cat("Computing within- and inter-cluster distances...\n")

dist_mat <- as.matrix(dist_matrix)
rownames(dist_mat) <- data$stn
colnames(dist_mat) <- data$stn

cluster_ids <- sort(unique(final_clusters))
k <- length(cluster_ids)

mean_within <- function(indices) {
    if (length(indices) < 2) {
        return(NA_real_)
    }
    sub <- dist_mat[indices, indices, drop = FALSE]
    mean(sub[upper.tri(sub)])
}

mean_between <- function(idx_i, idx_j) {
    mean(as.vector(dist_mat[idx_i, idx_j, drop = FALSE]))
}

within_df <- tibble(
    Cluster = as.integer(cluster_ids),
    N_stations = vapply(cluster_ids, function(g) sum(final_clusters == g), integer(1)),
    Mean_within = vapply(cluster_ids, function(g) {
        mean_within(which(final_clusters == g))
    }, numeric(1))
)

dist_matrix_summary <- matrix(NA_real_, nrow = k, ncol = k)
rownames(dist_matrix_summary) <- paste0("Cluster ", cluster_ids)
colnames(dist_matrix_summary) <- paste0("Cluster ", cluster_ids)

for (i in seq_len(k)) {
    idx_i <- which(final_clusters == cluster_ids[i])
    for (j in seq_len(k)) {
        idx_j <- which(final_clusters == cluster_ids[j])
        if (i == j) {
            dist_matrix_summary[i, j] <- mean_within(idx_i)
        } else {
            dist_matrix_summary[i, j] <- mean_between(idx_i, idx_j)
        }
    }
}

within_pairs <- unlist(lapply(cluster_ids, function(g) {
    idx <- which(final_clusters == g)
    if (length(idx) < 2) {
        return(numeric(0))
    }
    sub <- dist_mat[idx, idx, drop = FALSE]
    sub[upper.tri(sub)]
}))

between_pairs <- unlist(lapply(combn(seq_len(k), 2, simplify = FALSE), function(pair) {
    idx_i <- which(final_clusters == cluster_ids[pair[1]])
    idx_j <- which(final_clusters == cluster_ids[pair[2]])
    as.vector(dist_mat[idx_i, idx_j, drop = FALSE])
}))

summary_df <- tibble(
    Metric = c("Mean within-cluster distance", "Mean inter-cluster distance"),
    Value = c(mean(within_pairs), mean(between_pairs))
)

write.csv(within_df, file.path(data_dir, "cluster_within_distances.csv"), row.names = FALSE)
write.csv(dist_matrix_summary, file.path(data_dir, "cluster_distance_matrix.csv"))
write.csv(summary_df, file.path(data_dir, "cluster_distance_summary.csv"), row.names = FALSE)

cat(sprintf(
    "  Mean within-cluster distance: %.3f\n  Mean inter-cluster distance: %.3f\n",
    summary_df$Value[1], summary_df$Value[2]
))

################################################################################
# 4. Save Results
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

################################################################################
# 5. LaTeX Table for Revision (7 x 7 distance matrix)
################################################################################
cat("Generating LaTeX table...\n")

dist_formatted <- matrix(
    formatC(as.vector(dist_matrix_summary), format = "f", digits = 3),
    nrow = k,
    ncol = k,
    byrow = TRUE,
    dimnames = list(
        as.character(cluster_ids),
        as.character(cluster_ids)
    )
)
for (i in seq_len(k)) {
    min_j <- which.min(dist_matrix_summary[i, ])
    dist_formatted[i, min_j] <- paste0("\\textbf{", dist_formatted[i, min_j], "}")
}
dist_formatted <- as.data.frame(dist_formatted, stringsAsFactors = FALSE)
dist_formatted <- cbind(Cluster = as.character(cluster_ids), dist_formatted)

matrix_table <- xtable(
    dist_formatted,
    align = c("l", rep("r", k + 1)),
    caption = paste0(
        "Mean pairwise Euclidean distances in PCA space ",
        "(first ", n, " PCs explaining ", formatC(total, format = "f", digits = 2),
        "\\% variance). Diagonal: within-cluster; off-diagonal: inter-cluster."
    ),
    label = "tab:cluster_distance_matrix"
)

latex_file <- file.path(revision_dir, "cluster_distances.tex")
sink(latex_file)
print(
    matrix_table,
    include.rownames = FALSE,
    sanitize.text.function = identity,
    booktabs = TRUE,
    floating = TRUE,
    table.placement = "htbp",
    caption.placement = "top"
)
sink()

cat(sprintf("LaTeX table saved to %s\n", latex_file))
cat("Distance CSVs saved to Data/cluster_within_distances.csv, ",
    "Data/cluster_distance_matrix.csv, Data/cluster_distance_summary.csv\n",
    sep = ""
)
