setwd("C:/Users/thorp/Downloads")
rm(list=ls())

library(cluster)
##read in metadata
tbprofiler <- read.delim("tbprofilerall.txt")
##reformat sample
names(tbprofiler)[names(tbprofiler) =="sample"] <- "Sample"
##read in matrix data
matrix_data <- read.delim("matrix.tsv", row.names=1)
##get distances from  matrix
distance <- as.dist(matrix_data)
##define a function 'cluster_analysis' for clustering analysis
cluster_analysis <- function(threshold, distance) {
  ##perform agglomerative hierarchical clustering
  clusters <- agnes(distance, diss = TRUE, method = "average")
  ##convert the dendrogram to clusters based on the specified height threshold
  clusters <- as.data.frame(cutree(as.hclust(clusters), h = threshold))
  colnames(clusters) <- "cluster_id"
  ##extract unique cluster IDs
  clusterIDs <- unique(subset(clusters, duplicated(clusters$cluster_id))$cluster_id)
  clusters <- subset(clusters, cluster_id %in% clusterIDs)
  ##calculate cluster sizes and other statistics
  cluster_sizes <- table(clusters$cluster_id)
  total_clusters <- length(unique(clusters$cluster_id))
  median_size <- median(cluster_sizes)
  range_size <- range(cluster_sizes)
  
  ##return a list of cluster statistics
  return(list(total_clusters=total_clusters, median_size=median_size, range_size=range_size))
}
##define thresholds for clustering
thresholds <- c(1, 5, 10, 20)
##apply the cluster analysis function
results <- lapply(thresholds, function(t) cluster_analysis(t, distance))
##for each threshold print out the clustering results
for (i in 1:length(thresholds)) {
  print(paste("For threshold", thresholds[i], "SNPs:"))
  print(paste("  Total number of clusters:", results[[i]]$total_clusters))
  print(paste("  Median cluster size:", results[[i]]$median_size))
  print(paste("  Range of cluster sizes:", results[[i]]$range_size[1], "to", results[[i]]$range_size[2]))
}
##initialise a list to store sample names for each study
study_sample_names <- list()

##loop over the thresholds to identify study samples
for (t in thresholds) {
  ## Perform clustering at the current threshold
  current_clusters <- as.data.frame(cutree(as.hclust(agnes(distance, diss = TRUE, method = "average")), h = t))
  colnames(current_clusters) <- "cluster_id"
  current_clusters$Sample <- rownames(current_clusters)
  ## Merge cluster data with tbprofiler metadata
  merged_data_current <- merge(tbprofiler, current_clusters, by="Sample", all=TRUE)
  
  ## Calculate cluster sizes here
  cluster_sizes <- table(current_clusters$cluster_id)
  
  ## Subset the data to identify study samples based on cluster size and IDs
  study_samples <- subset(merged_data_current, 
                          !is.na(cluster_id) & 
                            cluster_id != 1 & 
                            grepl("^MTB", Sample) & 
                            cluster_sizes[as.numeric(cluster_id)] >= 3)
  
  ## Store the study sample names for each threshold
  study_sample_names[[as.character(t)]] <- study_samples$Sample
}

## Loop over the thresholds to print study sample names
for (t in thresholds) {
  print(paste("For threshold", t, "SNPs:"))
  print(paste("  Study sample names in clusters (3 or more):", toString(study_sample_names[[as.character(t)]])))
}
