setwd("C:/Users/thorp/Downloads")
rm(list=ls())

library(cluster)

tbprofiler <- read.delim("tbprofilerall.txt")
names(tbprofiler)[names(tbprofiler) =="sample"] <- "Sample"

matrix_data <- read.delim("matrix.tsv", row.names=1)

distance <- as.dist(matrix_data)

cluster_analysis <- function(threshold, distance) {
  clusters <- agnes(distance, diss = TRUE, method = "average")
  clusters <- as.data.frame(cutree(as.hclust(clusters), h = threshold))
  colnames(clusters) <- "cluster_id"
  clusterIDs <- unique(subset(clusters, duplicated(clusters$cluster_id))$cluster_id)
  clusters <- subset(clusters, cluster_id %in% clusterIDs)
  cluster_sizes <- table(clusters$cluster_id)
  total_clusters <- length(unique(clusters$cluster_id))
  median_size <- median(cluster_sizes)
  range_size <- range(cluster_sizes)
  
  return(list(total_clusters=total_clusters, median_size=median_size, range_size=range_size))
}
thresholds <- c(1, 5, 10, 20)

results <- lapply(thresholds, function(t) cluster_analysis(t, distance))

for (i in 1:length(thresholds)) {
  print(paste("For threshold", thresholds[i], "SNPs:"))
  print(paste("  Total number of clusters:", results[[i]]$total_clusters))
  print(paste("  Median cluster size:", results[[i]]$median_size))
  print(paste("  Range of cluster sizes:", results[[i]]$range_size[1], "to", results[[i]]$range_size[2]))
}

study_sample_names <- list()

for (t in thresholds) {
  current_clusters <- as.data.frame(cutree(as.hclust(agnes(distance, diss = TRUE, method = "average")), h = t))
  colnames(current_clusters) <- "cluster_id"
  current_clusters$Sample <- rownames(current_clusters)
  merged_data_current <- merge(tbprofiler, current_clusters, by="Sample", all=TRUE)

  study_samples <- subset(merged_data_current, 
                          !is.na(cluster_id) & 
                            cluster_id != 1 & 
                            grepl("^MTB", Sample) & 
                            cluster_sizes[as.numeric(cluster_id)] >= 2)
  
  study_sample_names[[as.character(t)]] <- study_samples$Sample
}

for (t in thresholds) {
  print(paste("For threshold", t, "SNPs:"))
  print(paste("  Study sample names in clusters (2 or more):", toString(study_sample_names[[as.character(t)]])))
}
